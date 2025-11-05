import pandas as pd
import json
import os
import time
from collections import OrderedDict

class MET500Parser:
    """
    面向对象的MET500 (MET500 cohort) 数据解析器
    封装了从MET500数据文件提取基因表达信息、解析和转换的功能
    """
    
    def __init__(self, ensembl_id=None,
                 genename=None,
                 sample_expr_file=None, 
                 metadata_file=None):
        """
        初始化MET500解析器
        
        Args:
            ensembl_id (str, optional): Ensembl基因ID
            sample_expr_file (str): 样本水平表达数据文件路径  
            metadata_file (str): 样本元数据文件路径
        """
        self.ensembl_id = ensembl_id
        self.genename = genename
        self.sample_expr_file = sample_expr_file
        self.metadata_file = metadata_file
        self.raw_sample_data = None
        self.parsed_data = None
        self.processing_time = None
        
        # 验证文件存在性
        if sample_expr_file and not os.path.exists(sample_expr_file):
            raise FileNotFoundError(f"样本表达文件不存在: {sample_expr_file}")
        if metadata_file and not os.path.exists(metadata_file):
            raise FileNotFoundError(f"元数据文件不存在: {metadata_file}")
    
    def extract_gene_data_memory_safe(self, file_path, target_gene):
        """
        内存安全的基因数据提取方法
        
        Args:
            file_path: 数据文件路径
            target_gene: 目标基因ID
            
        Returns:
            dict: 包含基因数据和元信息的字典
        """
        start_time = time.time()
        
        try:
            # 逐行扫描文件避免内存问题
            target_lines = []
            line_count = 0
            found_count = 0
            
            with open(file_path, 'r') as f:
                # 读取列名
                columns = next(f).strip().split('\t')
                columns[0] = 'GeneID'
                
                # 逐行扫描目标基因
                for line in f:
                    line_count += 1
                    if line_count % 10000 == 0:
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
                    print (f"{columns[0]}列名的个数和值不匹配，值的个数多于样本数匹配")
                    values = values[:len(columns)]
                elif len(values) < len(columns):
                    print (f"{columns[0]}列名的个数和值不匹配，值的个数少于样本数匹配")
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
        
        sample_columns = gene_info_df.columns[1:-1]  # 跳过GeneID列
        
        expression_long_list = []
        for sample in sample_columns:
            expression_long_list.append({
                'GeneID': gene_info_df['GeneID'].iloc[0],
                'Sample_id': sample,
                'Expression': gene_info_df[sample].iloc[0]
            })
        
        expression_long = pd.DataFrame(expression_long_list)
     
        # 合并样本信息
        sample_info_subset = sample_info_df[['Sample_id', 'sample_type', 'tissue','cohort','biopsy_tissue','tc']].copy()
        
        merged_df = pd.merge(
            expression_long, 
            sample_info_subset, 
            on='Sample_id', 
            how='left'
        )
        
        final_columns = ['GeneID', 'Sample_id', 'sample_type', 'tissue','cohort','biopsy_tissue','tc','Expression']
        merged_df = merged_df[final_columns]
        merged_df = merged_df.rename(columns={'tc':'tumor_content'})
        return merged_df
    
    def _calculate_expression_stats(self, sample_df):
        """计算表达统计信息"""
        sample_expressions = sample_df['Expression'].astype(float)
        
        stats = {
            "mean_expression": float(sample_expressions.mean()),
            "median_expression": float(sample_expressions.median()),
            "expression_range": {
                "min": float(sample_expressions.min()),
                "max": float(sample_expressions.max())
            },
            "detected_samples": int((sample_expressions > 0).sum()),
            "total_samples": len(sample_expressions)
        }
        return stats
    
    def create_structured_json(self, sample_expression_df,target_gene_name):
        """
        创建结构化JSON输出
        
        Args:
            sample_expression_df: 样本表达数据
            
        Returns:
            str: JSON格式字符串
        """
        if sample_expression_df.empty:
            return json.dumps({"error": "无有效基因表达数据"}, indent=2)
        
        # 基础基因信息
        gene_info = {
            "gene_id": sample_expression_df['GeneID'].iloc[0],
            "gene_name": target_gene_name,
        }
        
        # 样本级别表达数据
        sample_expression_data = []
        for _, row in sample_expression_df.iterrows():
            sample_data = {
                "sample_id": row['Sample_id'],
                "expression_value": float(row['Expression']),
                'sample_type':row['sample_type'],
                'tissue':row['tissue'],
                'cohort':row['cohort'],
                'biopsy_tissue':row['biopsy_tissue'],
                'tumor_content':row['tumor_content']
            }
            sample_expression_data.append(sample_data)
        
        # 表达统计信息
        expression_stats = self._calculate_expression_stats(sample_expression_df)
        
        # 构建完整数据结构
        structured_data = OrderedDict([
            ("gene_information", gene_info),
            ("sample_level_expression", {
                "total_samples": len(sample_expression_data),
                "samples": sample_expression_data
            }),
            ("expression_statistics", expression_stats),
            ("data_metadata", {
                "data_source": "MET500 Analysis",
                "processing_date": pd.Timestamp.now().strftime("%Y-%m-%d"),
                "units": "FPKM (Fragments Per Kilobase of exon model per Million mapped fragments)",
                "ensembl_id": self.ensembl_id,
                "gene_name":self.genename
            })
        ])
        
        return json.dumps(structured_data, indent=2, ensure_ascii=False)
    
    def parse(self, ensembl_id=None):
        """
        主解析方法：执行完整的GTEx数据解析流程
        
        Args:
            ensembl_id (str, optional): 基因ID，如果提供则覆盖初始化参数
            
        Returns:
            str: 格式化后的JSON字符串
        """
        start_time = time.time()
        
        if ensembl_id:
            self.ensembl_id = ensembl_id
        
        if not self.ensembl_id:
            raise ValueError("必须提供Ensembl基因ID")
        
        if not all([self.sample_expr_file, self.metadata_file]):
            raise ValueError("必须提供所有必要的文件路径")        
        
        try:
            print(f"🔍 开始解析基因 {self.ensembl_id} 的MET500数据...")
            
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

            # 4. 生成结构化JSON
            json_output = self.create_structured_json(
                processed_sample_expr,
                self.genename
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
            return f"MET500Parser(基因: {self.ensembl_id}, 状态: 已解析, 耗时: {self.processing_time:.2f}秒)"
        elif self.ensembl_id:
            return f"MET500Parser(基因: {self.ensembl_id}, 状态: 未解析)"
        else:
            return "MET500Parser(状态: 未初始化)"