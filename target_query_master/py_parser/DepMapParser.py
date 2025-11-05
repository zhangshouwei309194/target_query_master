import pandas as pd
import numpy as np
import json
import os
import time
from collections import OrderedDict
import dask.dataframe as dd
from typing import Dict, Union

class DepMapParser:
    """
    面向对象的DepMap数据解析器
    封装了从DepMap数据文件提取基因表达信息、解析和转换的功能
    """
    
    def __init__(self,
                 genename=None,
                 sample_expr_file=None, 
                 metadata_file=None):
        """
        初始化TCGA解析器
        
        Args:
            genename (str, optional): gene symbol
            sample_expr_file (str): 样本水平表达数据文件路径
            metadata_file (str): 样本元数据所有癌症细胞系样本类型文件路径
        """
        self.genename = genename
        self.sample_expr_file = sample_expr_file
        self.metadata_file = metadata_file
        self.parsed_data = None
        self.processing_time = None
        
        # 验证文件存在性
        if sample_expr_file and not os.path.exists(sample_expr_file):
            raise FileNotFoundError(f"样本表达文件不存在: {sample_expr_file}")
        if metadata_file and not os.path.exists(metadata_file):
            raise FileNotFoundError(f"所有癌症样本类型元数据文件不存在: {metadata_file}")
    
    def extract_gene_expression_dask(self,file_path: str, target_gene: str) -> Dict[str, Union[dd.DataFrame, bool, str]]:
        """
        使用Dask高效读取基因表达数据并提取目标基因的表达信息
        
        参数:
            file_path: CSV文件路径
            target_gene: 目标基因名称(如"TNMD")
        
        返回:
            字典包含:
            - 'data': Dask DataFrame包含四列: Symbol, NCBI_Gene, ModelID, Expression_log
            - 'gene_found': 布尔值，表示是否成功找到基因
            - 'error': 错误信息字符串，如果没有错误则为空
        """
        
        # 初始化返回字典
        result = {
            'data': pd.DataFrame(),  # 空DataFrame
            'gene_found': False,
            'error': ''
        }
        
        try:
            # 1. 使用Pandas安全读取列名（避免Dask样本大小问题）
            try:
                print("使用Pandas读取列名...")
                df_sample = pd.read_csv(file_path, nrows=0)
                all_columns = df_sample.columns.tolist()
                print(f"成功读取 {len(all_columns)} 列")
            except Exception as pd_error:
                result['error'] = f"无法读取文件列名: {pd_error}"
                return result
            
            # 2. 改进列名解析逻辑，适应实际文件格式
            gene_column_info = {}  # 存储基因列信息: {基因名: (列索引, NCBI_ID)}
            for i, col_name in enumerate(all_columns):
                # 跳过第一列(ModelID列)
                if i == 0:
                    continue
                    
                # 改进的正则表达式，适应实际列名格式
                # 从图片看，列名格式为: "TSPAN6 (7105)"（注意有空格）
                try:
                    gene_name = col_name.split(' (')[0]
                    ncbi_id = col_name.split(' (')[1].split(')')[0]
                    gene_column_info[gene_name] = (i, ncbi_id)
                except Exception:
                    # 如果解析失败，跳过此列
                    continue
            
            print(f"成功解析 {len(gene_column_info)} 个基因列")
            
            # 3. 查找目标基因（不区分大小写）
            target_gene_upper = target_gene.upper()
            found_gene = None
            for gene_name in gene_column_info:
                if gene_name.upper() == target_gene_upper:
                    found_gene = gene_name
                    break
            
            if found_gene is None:
                available_genes = list(gene_column_info.keys())[:10]
                result['error'] = f"未找到基因 '{target_gene}'。可用的基因包括: {available_genes}"
                return result
            
            target_col_index, target_ncbi_id = gene_column_info[found_gene]
            print(f"定位到目标基因: {found_gene} (NCBI: {target_ncbi_id})")
            
            # 4. 确定需要读取的列
            columns_to_read = [all_columns[0], all_columns[target_col_index]]
            print(f"将读取列: {columns_to_read}")
            
            # 5. 使用Dask读取指定列，增加样本大小解决长列名问题
            dtypes = {
                all_columns[0]: 'object',    # ModelID列为字符串
                all_columns[target_col_index]: 'float64'  # 表达值为浮点数
            }
            
            try:
                # 方法1: 使用Dask并增加样本大小
                gene_data = dd.read_csv(
                    file_path,
                    usecols=columns_to_read,
                    dtype=dtypes,
                    sample=50 * 1024 * 1024,  # 增加样本大小到50MB
                    blocksize='64MB',
                    assume_missing=True
                )
                print("使用Dask成功读取数据")
                
            except Exception as e:
                print(f"Dask读取失败: {e}")
                # 方法2: 使用Pandas读取，然后转换为Dask
                print("使用Pandas读取数据...")
                try:
                    pandas_data = pd.read_csv(file_path, usecols=columns_to_read)
                    gene_data = dd.from_pandas(pandas_data, npartitions=4)
                    print("使用Pandas成功读取并转换为Dask DataFrame")
                except Exception as pd_error:
                    result['error'] = f"Pandas读取也失败: {pd_error}"
                    return result
            
            # 6. 重命名列并添加基因信息
            gene_data = gene_data.rename(columns={
                all_columns[0]: 'ModelID',
                all_columns[target_col_index]: 'Expression_log'
            })
            
            # 添加基因符号和NCBI基因ID列
            gene_data['Symbol'] = found_gene
            gene_data['NCBI_Gene'] = target_ncbi_id
            
            # 重新排列列顺序
            final_columns = ['Symbol', 'NCBI_Gene', 'ModelID', 'Expression_log']
            gene_data = gene_data[final_columns]
            gene_data = gene_data.compute()
            gene_data['Expression'] = gene_data['Expression_log'].apply(lambda x:(2**float(x)-1))
            
            # 更新返回结果
            result['data'] = gene_data
            result['gene_found'] = True
            
            return result
            
        except Exception as e:
            result['error'] = f"发生未知错误: {str(e)}"
            return result
    
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

        focus_col = ['ModelID', 'PatientID', 'CellLineName', 'StrippedCellLineName', 
                     'DepmapModelType', 'OncotreeLineage', 'OncotreePrimaryDisease', 
                     'OncotreeSubtype', 'OncotreeCode', 'PatientSubtypeFeatures',
                     'PrimaryOrMetastasis','SampleCollectionSite','ModelType',
                     'GrowthPattern','ModelSubtypeFeatures']
        sample_info_subset = sample_info_df[focus_col]
        # 合并样本信息

        gene_info_df['ModelID'] = gene_info_df['ModelID'].str.strip()
        
        sample_info_subset['ModelID'] = sample_info_df['ModelID'].str.strip()
        
        merged_df = pd.merge(
            gene_info_df, 
            sample_info_subset, 
            on='ModelID',
            how='left'
        )
        
        return merged_df
    
    def _calculate_expression_stats(self, sample_df):
        """计算表达统计信息"""
        sample_expressions = sample_df['Expression'].astype(float)
        
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
            }
        }
        return stats
    
    def create_structured_json(self, sample_expression_df, sample_metadata_df):
        """
        创建结构化JSON输出
        
        Args:
            sample_expression_df: 细胞系表达数据
            sample_metadata_df: 细胞系样本信息  
            
        Returns:
            str: JSON格式字符串
        """
        if sample_expression_df.empty or sample_metadata_df.empty:
            return json.dumps({"error": "无有效基因表达数据"}, indent=2)
        
        # 基础基因信息
        gene_info = {
            "gene_id": sample_expression_df['NCBI_Gene'].iloc[0],
            "gene_name": sample_expression_df['Symbol'].iloc[0],
        }
        # 样本级别表达数据 
        sample_expression_data = []
        for _, row in sample_expression_df.iterrows():
            sample_data = {
                "ModelID": row['ModelID'],
                "PatientID":  row['PatientID'],
                "expression_value": float(row['Expression']),
                "CellLineName": row['CellLineName'],
                "StrippedCellLineName": row['StrippedCellLineName'],
                "DepmapModelType": row['DepmapModelType'],
                "OncotreeLineage":row['OncotreeLineage'],
                "OncotreePrimaryDisease":row['OncotreePrimaryDisease'],
                "OncotreeSubtype":row['OncotreeSubtype'],
                "OncotreeCode":row['OncotreeCode'],
                "PatientSubtypeFeatures":row['PatientSubtypeFeatures'],
                "PrimaryOrMetastasis":row['PrimaryOrMetastasis'],
                "SampleCollectionSite":row['SampleCollectionSite'],
                "ModelType":row['ModelType'],
                "GrowthPattern":row['GrowthPattern'],
                "ModelSubtypeFeatures":row['ModelSubtypeFeatures']
            }
            sample_expression_data.append(sample_data)

        # 表达统计信息
        expression_stats = self._calculate_expression_stats(sample_expression_df)

        # DepMap细胞系表型信息
        samples_info_list = []
        columns_to_extract = sample_metadata_df.columns.tolist()
        # 遍历每一行（每个样本）
        for index, row in sample_metadata_df.iterrows():
            # 使用OrderedDict确保键的顺序
            sample_dict = OrderedDict()
            # 按照columns_to_extract的顺序提取数据
            for col in columns_to_extract:
                value = row[col]
                if pd.isna(value):
                    sample_dict[col] = None
                else:
                    sample_dict[col] = value
            samples_info_list.append(sample_dict)

        # 构建完整数据结构
        structured_data = OrderedDict([
            ("gene_information", gene_info),
            ("sample_level_expression", {
                "total_samples": len(sample_expression_data),
                "samples": sample_expression_data
            }),
            ("expression_statistics", expression_stats),
            ("cellline_metadata",{'samples':samples_info_list,'sample_count':len(samples_info_list)}),
            ("data_metadata", {
                "data_source": "DepMap Analysis",
                "processing_date": pd.Timestamp.now().strftime("%Y-%m-%d"),
                "units": "TPM (Transcripts Per Million)",
                "gene_id": sample_expression_df['NCBI_Gene'].iloc[0],
                "gene_name":self.genename
            })
        ])
        
        return json.dumps(structured_data, indent=2, ensure_ascii=False)
    
    def parse(self, genename=None):
        """
        主解析方法：执行完整的DepMap数据解析流程
        
        Args:
            genename (str, optional): 基因名字，如果提供则覆盖初始化参数
            
        Returns:
            str: 格式化后的JSON字符串
        """
        start_time = time.time()
        
        if genename:
            self.genename = genename
        
        if not self.genename:
            raise ValueError("必须提供基因名称（Symbol,如NECTIN4）")
      
        if not all([self.sample_expr_file, self.metadata_file]):
            raise ValueError("必须提供所有必要的文件路径")
        
        try:
            print(f"🔍 开始解析基因 {self.genename} 的DepMap数据...")
            
            # 1. 提取样本水平表达数据
            sample_expr_result = self.extract_gene_expression_dask(
                self.sample_expr_file, self.genename
            )
            
            if not sample_expr_result['gene_found']:
                return json.dumps({
                    "error": f"未找到基因 {self.genename} 的表达数据",
                    "message": sample_expr_result.get('message', '')
                }, indent=2)
            sample_expr_df = sample_expr_result['data']
            
            # 2. 加载样本元数据
            sample_metadata_df = pd.read_csv(self.metadata_file, low_memory=False)
            
            # 3. 处理样本表达数据
            processed_sample_expr = self.process_sample_expression(sample_expr_df, sample_metadata_df)
            
            # 4. 生成结构化JSON
            json_output = self.create_structured_json(
                processed_sample_expr,
                sample_metadata_df
            )
            # 记录处理时间
            self.processing_time = time.time() - start_time
            print(f"✅ 基因 {self.genename} 数据解析完成，耗时: {self.processing_time:.2f}秒")
            
            self.parsed_data = json.loads(json_output)
            return json_output
            
        except Exception as e:
            error_msg = f"解析基因 {self.genename} 数据时发生错误: {str(e)}"
            print(f"❌ {error_msg}")
            return json.dumps({"error": error_msg}, indent=2)
    
    def get_processing_time(self):
        """获取处理时间"""
        return self.processing_time
    
    def get_parsed_data(self):
        """获取解析后的数据"""
        return self.parsed_data
    
    def save_to_file(self, filename, genename=None):
        """
        将解析结果保存到文件
        
        Args:
            filename (str): 输出文件名
            genename (str, optional): 基因名称
        """
        if not self.parsed_data:
            if genename:
                self.parse(genename)
            else:
                self.parse()
        
        json_output = json.dumps(self.parsed_data, indent=2, ensure_ascii=False)
        with open(filename, 'w', encoding='utf-8') as f:
            f.write(json_output)
        print(f"💾 数据已保存到: {filename}")
    
    def __str__(self):
        """字符串表示"""
        if self.parsed_data:
            return f"DepMapParser(基因: {self.genename}, 状态: 已解析, 耗时: {self.processing_time:.2f}秒)"
        elif self.genename:
            return f"DepMapParser(基因: {self.genename}, 状态: 未解析)"
        else:
            return "DepMapParser(状态: 未初始化)"