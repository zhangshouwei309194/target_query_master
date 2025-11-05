import pandas as pd
import json
import os
import time
from collections import OrderedDict

class GTExParser:
    """
    面向对象的GTEx (Genotype-Tissue Expression) 数据解析器
    封装了从GTEx数据文件提取基因表达信息、解析和转换的功能
    """
    
    def __init__(self, ensembl_id=None, 
                 sample_expr_file=None, 
                 tissue_expr_file=None, 
                 metadata_file=None):
        """
        初始化GTEx解析器
        
        Args:
            ensembl_id (str, optional): Ensembl基因ID
            sample_expr_file (str): 样本水平表达数据文件路径
            tissue_expr_file (str): 组织水平表达数据文件路径  
            metadata_file (str): 样本元数据文件路径
        """
        self.ensembl_id = ensembl_id
        self.sample_expr_file = sample_expr_file
        self.tissue_expr_file = tissue_expr_file
        self.metadata_file = metadata_file
        self.raw_sample_data = None
        self.raw_tissue_data = None
        self.parsed_data = None
        self.processing_time = None
        
        # 验证文件存在性
        if sample_expr_file and not os.path.exists(sample_expr_file):
            raise FileNotFoundError(f"样本表达文件不存在: {sample_expr_file}")
        if tissue_expr_file and not os.path.exists(tissue_expr_file):
            raise FileNotFoundError(f"组织表达文件不存在: {tissue_expr_file}")
        if metadata_file and not os.path.exists(metadata_file):
            raise FileNotFoundError(f"元数据文件不存在: {metadata_file}")
    
    def extract_gene_data_memory_safe(self, file_path, target_gene):
        """
        内存安全的基因数据提取方法
        
        Args:
            file_path: 数据文件路径
            target_gene: 目标基因ID
            use_cluster: 是否使用分布式集群
            
        Returns:
            dict: 包含基因数据和元信息的字典
        """
        start_time = time.time()
        client = None
        
        try:
            # 逐行扫描文件避免内存问题
            target_lines = []
            line_count = 0
            found_count = 0
            
            with open(file_path, 'r') as f:
                for _ in range(2):  # 跳过GCT文件头
                    next(f)
                    line_count += 1
                
                # 读取列名
                columns = next(f).strip().split('\t')
                
                # 逐行扫描目标基因
                for line in f:
                    line_count += 1
                    if line_count % 1000000 == 0:
                        print(f"📈 已扫描 {line_count} 行，找到 {found_count} 个匹配")
                    
                    parts = line.split('\t', 1)
                    if not parts:
                        continue
                    
                    gene_id_full = parts[0]
                    gene_id = gene_id_full.split('.')[0] if '.' in gene_id_full else gene_id_full
                    
                    if gene_id == target_gene:
                        target_lines.append(line.strip())
                        found_count += 1
            
            if not target_lines:
                return {
                    'data': pd.DataFrame(),
                    'gene_found': False,
                    'message': f'基因 {target_gene} 未找到'
                }
            
            # 构建结果DataFrame
            data_rows = []
            for line in target_lines:
                values = line.split('\t')
                if len(values) > len(columns):
                    values = values[:len(columns)]
                elif len(values) < len(columns):
                    values.extend([''] * (len(columns) - len(values)))
                data_rows.append(values)
            
            result_df = pd.DataFrame(data_rows, columns=columns)
            
            return {
                'data': result_df,
                'gene_found': True,
                'lines_scanned': line_count,
                'message': f'成功提取基因 {target_gene} 的 {len(result_df)} 行数据'
            }
            
        except Exception as e:
            return {
                'data': pd.DataFrame(),
                'gene_found': False,
                'error': str(e)
            }
    
    def process_sample_expression(self, gene_info_df, sample_info_df):
        """
        处理样本水平表达数据
        
        Args:
            gene_info_df: 基因表达数据框
            sample_info_df: 样本信息数据框
            
        Returns:
            pd.DataFrame: 处理后的表达数据
        """
        if gene_info_df.empty:
            return pd.DataFrame()
        
        sample_columns = gene_info_df.columns[2:-1]  # 跳过Name, Description和GeneID列
        
        expression_long_list = []
        for sample in sample_columns:
            expression_long_list.append({
                'Name': gene_info_df['Name'].iloc[0],
                'Description': gene_info_df['Description'].iloc[0],
                'Sample': sample,
                'Expression': gene_info_df[sample].iloc[0]
            })
        
        expression_long = pd.DataFrame(expression_long_list)
        
        # 合并样本信息
        sample_info_subset = sample_info_df[['SAMPID', 'SMTS', 'SMTSD']].copy()
        
        expression_long['Sample_Clean'] = expression_long['Sample'].str.strip()
        sample_info_subset['SAMPID_Clean'] = sample_info_subset['SAMPID'].str.strip()
        
        merged_df = pd.merge(
            expression_long, 
            sample_info_subset, 
            left_on='Sample_Clean', 
            right_on='SAMPID_Clean', 
            how='left'
        )
        
        merged_df.drop(['Sample_Clean', 'SAMPID_Clean'], axis=1, inplace=True)
        final_columns = ['Name', 'Description', 'Sample', 'Expression', 'SAMPID', 'SMTS', 'SMTSD']
        
        return merged_df[final_columns]
    
    def process_tissue_expression(self, median_gene_info_df, target_gene):
        """
        处理组织特异性表达数据
        
        Args:
            median_gene_info_df: 组织表达数据框
            target_gene: 目标基因ID
            
        Returns:
            dict: 包含组织表达和特异性信息字典
        """
        median_gene_info_df['GeneID'] = median_gene_info_df['Name'].str.split('.').str[0]
        target_gene_info_df = median_gene_info_df[median_gene_info_df['GeneID'] == target_gene]
        target_gene_info_df = target_gene_info_df.reset_index(drop=True)
        
        if target_gene_info_df.empty:
            return {'Tissue_expression': pd.DataFrame(), 'Tissue_specificity': {}}
        
        tissue_specificity_columns = target_gene_info_df.columns[2:5]
        tissue_columns = target_gene_info_df.columns[5:-1] 
        
        expression_long_list = []
        for tissue in tissue_columns:
            expression_long_list.append({
                'Name': target_gene_info_df['Name'].iloc[0],
                'Description': target_gene_info_df['Description'].iloc[0],
                'Tissue': tissue,
                'Expression': target_gene_info_df[tissue].iloc[0]
            })
        
        expression_long_df = pd.DataFrame(expression_long_list)

        # 组织特异性信息
        tissue_specificity_dict = {}
        for tissue_specificity in tissue_specificity_columns:
            tissue_specificity_dict[tissue_specificity] = target_gene_info_df[tissue_specificity].iloc[0]

        return {
            'Tissue_expression': expression_long_df,
            'Tissue_specificity': tissue_specificity_dict
        }
    
    def _categorize_expression(self, expression_value):
        """根据表达值分类表达水平"""
        if expression_value == 0:
            return "Not detected"
        elif expression_value < 1:
            return "Low"
        elif expression_value < 10:
            return "Medium"
        else:
            return "High"
    
    def _calculate_expression_stats(self, sample_df, tissue_df):
        """计算表达统计信息"""
        sample_expressions = sample_df['Expression'].astype(float)
        tissue_expressions = tissue_df['Expression'].astype(float)
        
        stats = {
            "sample_level": {
                "mean_expression": float(sample_expressions.mean()),
                "median_expression": float(sample_expressions.median()),
                "expression_range": {
                    "min": float(sample_expressions.min()),
                    "max": float(sample_expressions.max())
                },
                "detected_samples": int((sample_expressions > 0).sum()),
                "total_samples": len(sample_expressions)
            },
            "tissue_level": {
                "mean_expression": float(tissue_expressions.mean()),
                "median_expression": float(tissue_expressions.median()),
                "expression_range": {
                    "min": float(tissue_expressions.min()),
                    "max": float(tissue_expressions.max())
                },
                "highly_expressed_tissues": int((tissue_expressions >= 10).sum())
            }
        }
        
        return stats
    
    def create_structured_json(self, sample_expression_df, tissue_expression_df, tissue_specificity_dict):
        """
        创建结构化JSON输出[9](@ref)
        
        Args:
            sample_expression_df: 样本表达数据
            tissue_expression_df: 组织表达数据  
            tissue_specificity_dict: 组织特异性字典
            
        Returns:
            str: JSON格式字符串
        """
        if sample_expression_df.empty or tissue_expression_df.empty:
            return json.dumps({"error": "无有效基因表达数据"}, indent=2)
        
        # 基础基因信息
        gene_info = {
            "gene_id": sample_expression_df['Name'].iloc[0],
            "gene_name": sample_expression_df['Description'].iloc[0],
        }
        
        # 样本级别表达数据
        sample_expression_data = []
        for _, row in sample_expression_df.iterrows():
            sample_data = {
                "sample_id": row['Sample'],
                "expression_value": float(row['Expression']),
                "tissue_type": row['SMTS'],
                "tissue_subtype": row['SMTSD'],
                "sample_accession": row['SAMPID']
            }
            sample_expression_data.append(sample_data)
        
        # 组织级别表达数据
        tissue_expression_data = []
        for _, row in tissue_expression_df.iterrows():
            tissue_data = {
                "tissue_name": row['Tissue'],
                "expression_value": float(row['Expression']),
                "expression_category": self._categorize_expression(float(row['Expression']))
            }
            tissue_expression_data.append(tissue_data)
        
        # 组织特异性信息
        specificity_info = {
            "tau_score": float(tissue_specificity_dict.get('Tau', 0)),
            "specificity_type": tissue_specificity_dict.get('Specificity_Type', 'Unknown'),
            "specific_tissues": tissue_specificity_dict.get('Specific_Tissues', '')
        }
        
        # 表达统计信息
        expression_stats = self._calculate_expression_stats(sample_expression_df, tissue_expression_df)
        
        # 构建完整数据结构
        structured_data = OrderedDict([
            ("gene_information", gene_info),
            ("sample_level_expression", {
                "total_samples": len(sample_expression_data),
                "samples": sample_expression_data
            }),
            ("tissue_level_expression", {
                "tissues_analyzed": len(tissue_expression_data),
                "tissue_data": tissue_expression_data
            }),
            ("tissue_specificity_analysis", specificity_info),
            ("expression_statistics", expression_stats),
            ("data_metadata", {
                "data_source": "GTEx Analysis",
                "processing_date": pd.Timestamp.now().strftime("%Y-%m-%d"),
                "units": "TPM (Transcripts Per Million)",
                "ensembl_id": self.ensembl_id
            })
        ])
        
        return json.dumps(structured_data, indent=2, ensure_ascii=False)
    
    def parse(self, ensembl_id=None, strategy='auto'):
        """
        主解析方法：执行完整的GTEx数据解析流程
        
        Args:
            ensembl_id (str, optional): 基因ID，如果提供则覆盖初始化参数
            strategy (str): 处理策略 ('auto', 'safe', 'dask')
            
        Returns:
            str: 格式化后的JSON字符串
        """
        start_time = time.time()
        
        if ensembl_id:
            self.ensembl_id = ensembl_id
        
        if not self.ensembl_id:
            raise ValueError("必须提供Ensembl基因ID")
        
        if not all([self.sample_expr_file, self.tissue_expr_file, self.metadata_file]):
            raise ValueError("必须提供所有必要的文件路径")
        
        try:
            print(f"🔍 开始解析基因 {self.ensembl_id} 的GTEx数据...")
            
            # 1. 提取样本水平表达数据
            sample_expr_result = self.extract_gene_data_memory_safe(
                self.sample_expr_file, self.ensembl_id
            )
            
            if not sample_expr_result['gene_found']:
                return json.dumps({
                    "error": f"未找到基因 {self.ensembl_id} 的表达数据",
                    "message": sample_expr_result.get('message', '')
                }, indent=2)
            
            sample_expr_df = sample_expr_result['data']

            
            # 2. 加载样本元数据
            sample_metadata_df = pd.read_csv(self.metadata_file, sep='\t', low_memory=False)
            
            # 3. 处理样本表达数据
            processed_sample_expr = self.process_sample_expression(sample_expr_df, sample_metadata_df)
            
            # 4. 处理组织表达数据
            tissue_expr_df = pd.read_csv(self.tissue_expr_file, sep='\t')
            tissue_result = self.process_tissue_expression(tissue_expr_df, self.ensembl_id)
            # 5. 生成结构化JSON
            json_output = self.create_structured_json(
                processed_sample_expr,
                tissue_result['Tissue_expression'],
                tissue_result['Tissue_specificity']
            )
            # 记录处理时间
            self.processing_time = time.time() - start_time
            print(f"✅ 基因 {self.ensembl_id} 数据解析完成，耗时: {self.processing_time:.2f}秒")
            
            self.parsed_data = json.loads(json_output)
            return json_output
            
        except Exception as e:
            error_msg = f"解析基因 {self.ensembl_id} 数据时发生错误: {str(e)}"
            print(f"❌ {error_msg}")
            return json.dumps({"error": error_msg}, indent=2)
    
    def get_processing_time(self):
        """获取处理时间"""
        return self.processing_time
    
    def get_parsed_data(self):
        """获取解析后的数据"""
        return self.parsed_data
    
    def save_to_file(self, filename, ensembl_id=None):
        """
        将解析结果保存到文件
        
        Args:
            filename (str): 输出文件名
            ensembl_id (str, optional): 基因ID
        """
        if not self.parsed_data:
            if ensembl_id:
                self.parse(ensembl_id)
            else:
                self.parse()
        
        json_output = json.dumps(self.parsed_data, indent=2, ensure_ascii=False)
        with open(filename, 'w', encoding='utf-8') as f:
            f.write(json_output)
        print(f"💾 数据已保存到: {filename}")
    
    def __str__(self):
        """字符串表示"""
        if self.parsed_data:
            return f"GTExParser(基因: {self.ensembl_id}, 状态: 已解析, 耗时: {self.processing_time:.2f}秒)"
        elif self.ensembl_id:
            return f"GTExParser(基因: {self.ensembl_id}, 状态: 未解析)"
        else:
            return "GTExParser(状态: 未初始化)"