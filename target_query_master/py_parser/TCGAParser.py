import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
import json
import os
import time
from collections import OrderedDict

class TCGAParser:
    """
    面向对象的TCGA数据解析器
    封装了从TCGA数据文件提取基因表达信息、解析和转换的功能
    """
    
    def __init__(self, ensembl_id=None,
                 genename=None,
                 sample_expr_file=None, 
                 metadata_file=None,
                 metadata_all_json_file=None):
        """
        初始化TCGA解析器
        
        Args:
            ensembl_id (str, optional): Ensembl基因ID
            genename: gene symbol
            sample_expr_file (str): 样本水平表达数据文件路径
            metadata_file (str): 样本元数据所有癌症样本类型文件路径
            metadata_all_json_file (str): 所有TCGA样本的临床信息，以json的格式存储
        """
        self.ensembl_id = ensembl_id
        self.genename = genename
        self.sample_expr_file = sample_expr_file
        self.metadata_file = metadata_file
        self.metadata_all_json_file = metadata_all_json_file
        self.parsed_data = None
        self.processing_time = None
        
        # 验证文件存在性
        if sample_expr_file and not os.path.exists(sample_expr_file):
            raise FileNotFoundError(f"样本表达文件不存在: {sample_expr_file}")
        if metadata_file and not os.path.exists(metadata_file):
            raise FileNotFoundError(f"所有癌症样本类型元数据文件不存在: {metadata_file}")
        if metadata_all_json_file and not os.path.exists(metadata_all_json_file):
            raise FileNotFoundError(f"所有癌症临床信息json格式文件不存在: {metadata_all_json_file}")
    
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
                # 读取列名
                columns = next(f).strip().split('\t')
                
                # 逐行扫描目标基因
                for line in f:
                    line_count += 1
                    if line_count % 100000 == 0:
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
                    continue
                elif len(values) < len(columns):
                    print (f"{columns[0]}列名的个数和值不匹配，值的个数少于样本数匹配")
                    continue
                else:
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
                'Gene_id': gene_info_df['Gene'].iloc[0],
                'Gene_name': self.genename,
                'sampleID': sample,
                'Expression': gene_info_df[sample].iloc[0]
            })
        
        expression_long = pd.DataFrame(expression_long_list)
  
        # 合并样本信息
        sample_info_subset = sample_info_df[['sampleID', 'patient', 'primary_disease', 'cancer_type', 'sample_type']].copy()

        expression_long['sampleID'] = expression_long['sampleID'].str.strip()
        sample_info_subset['sampleID'] = sample_info_subset['sampleID'].str.strip()
        
        merged_df = pd.merge(
            expression_long, 
            sample_info_subset, 
            on='sampleID',
            how='left'
        )
        
        return merged_df
    
    def process_sampletype_expression(self, sample_expr_df,sam_type_list):
        """
        处理TCGA不同癌症样本类型（关注样本类型:Primary Tumor|Solid Tissue Normal）,并计算不同癌症样本类型的中位值
        
        Args:
            sample_expr_df: TCGA所有癌症样本相应基因的表达数据
            sam_type_list: ['Primary Tumor','Solid Tissue Normal']
            
        Returns:
            dataframe: 包含不同癌症类型不同样本类型中位表达信息数据框
        """
        sample_expr_subset_df = sample_expr_df[sample_expr_df['sample_type'].isin(sam_type_list)].copy(deep=True)
        sample_expr_subset_df = sample_expr_subset_df.reset_index(drop=True)
        sample_expr_subset_gdf = sample_expr_subset_df.groupby(['cancer_type','sample_type'])
        cancertype_samtype_median_expr = []
        for name,group in sample_expr_subset_gdf:
            cancer = name[0]
            samtype = name[1]
            geneid = group['Gene_id'].tolist()[0]
            genename = group['Gene_name'].tolist()[0]
            group['Expression'] = pd.to_numeric(group['Expression'], errors='coerce')
            median_expression = np.median(group['Expression'])
            primary_disease = group['primary_disease'].tolist()[0]
            cancertype_samtype_median_expr.append([geneid,genename,cancer,primary_disease,samtype,median_expression])
            
        median_expression_long_df = pd.DataFrame(cancertype_samtype_median_expr,columns=['Gene_id','Gene_name','cancer_type','primary_disease','sample_type','Expression'])
        return median_expression_long_df

    def cancer_specificty(self,cancer_level_expression_df,tpm_threshold=1.0):
        """
        处理TCGA不同癌症（关注样本类型:Primary Tumor）组织特异性指数Tau (过滤TPM < 1的数据点)，并筛选特异性表达癌症
        
        Args:
            cancer_level_expression_df: TCGA所有癌症类型，基于所有癌症类型筛选特异性表达癌症
            
        Returns:
            dict: 返回该基因的组织特异性指数Tau及特异性表达组织，及Top1,Top2,Top3基因表达，相应的组织及75分位性表达值，中位表达值
        """        
        # 过滤掉TPM < threshold的值（设为0）
        cancer_level_expression_filtered_df = cancer_level_expression_df.copy()
        cancer_level_expression_filtered_df = cancer_level_expression_df[cancer_level_expression_df['sample_type']=='Primary Tumor']
        cancer_level_expression_filtered_df['Expression'] =  cancer_level_expression_filtered_df['Expression'].apply(lambda x:float(x))
        expr_vector = cancer_level_expression_df['Expression'].tolist()
        filtered_expr = expr_vector.copy()
        filtered_expr = [ value if value >= tpm_threshold else 0 for value in filtered_expr ]
        filtered_expr = np.array(filtered_expr)
        n_total_tissues = len(filtered_expr)
        # 找到最大表达值（忽略0值）
        max_expr = np.max(filtered_expr)
        
        # 如果最大表达值为0，Tau设为0
        if max_expr == 0:
            tau = 0.0
        else:
            # 计算相对表达量（所有组织都参与计算，过滤的为0）
            relative_expr = filtered_expr / max_expr
            
            # 计算组织特异性分数（1 - 相对表达量）
            specificity_scores = 1 - relative_expr
            
            # 计算Tau指数：分母使用总组织数-1
            if n_total_tissues <= 1:
                tau = 0.0
            else:
                tau = np.sum(specificity_scores) / (n_total_tissues - 1)

        ## 基于tau值，筛选特异性表达组织
        # 按表达值从高到低排序
        cancer_level_expression_filtered_df = cancer_level_expression_filtered_df.sort_values(by=['Expression'],ascending=[False])
        cancer_level_expression_filtered_df = cancer_level_expression_filtered_df.reset_index(drop=True)
        # 获取前三个最高表达的组织和表达值
        cancer_level_expression_filtered_df['primary_disease_str'] = cancer_level_expression_filtered_df.apply(lambda x:x['primary_disease']+"("+x['cancer_type']+")",axis=1)
        top_tissues = cancer_level_expression_filtered_df['primary_disease_str'].tolist()
        top_values = cancer_level_expression_filtered_df['Expression'].tolist()
        top1, top2, top3 = top_values[0], top_values[1], top_values[2]
        tissue1, tissue2, tissue3 = top_tissues[0], top_tissues[1], top_tissues[2]
        tissue_75percentile = np.percentile(top_values, 75) 
        tissue_median = np.median(top_values)
        #Tau阈值，默认0.8
        tau_threshold = 0.8
        if tau < tau_threshold:
            specificity_types = 'No'
            specific_tissues = ''
        else:
            # 防止除0错误，设置最小值
            epsilon = 1e-10
            top2_adj = top2 if top2 > 0 else epsilon
            top3_adj = top3 if top3 > 0 else epsilon
            
            # 计算log2比值
            log2_ratio1 = np.log2(top1 / top2_adj) if top1 > 0 else -np.inf
            log2_ratio2 = np.log2(top2_adj / top3_adj) if top2_adj > 0 else -np.inf
            # 判断组织特异性类型
            if log2_ratio1 > 1:
                specificity_types='Single tissue specificity'
                specific_tissues=tissue1
            elif log2_ratio2 > 0.5:
                specificity_types='Double tissue specificity'
                specific_tissues=f"{tissue1}, {tissue2}"
            else:
                specificity_types='No'
                specific_tissues=''
        top1, top2, top3 = top_values[0], top_values[1], top_values[2]
        tissue1, tissue2, tissue3 = top_tissues[0], top_tissues[1], top_tissues[2]
        
        cancer_expression_specificity_dict =  OrderedDict()
        cancer_expression_specificity_dict['tau_score'] = tau
        cancer_expression_specificity_dict['specificity_type'] = specificity_types
        cancer_expression_specificity_dict['specific_tissues'] = specific_tissues
        cancer_expression_specificity_dict['cancer.top1.expression'] = top1
        cancer_expression_specificity_dict['cancer.top1'] = tissue1
        cancer_expression_specificity_dict['cancer.top2.expression'] = top2
        cancer_expression_specificity_dict['cancer.top2'] = tissue2
        cancer_expression_specificity_dict['cancer.top3.expression'] = top3
        cancer_expression_specificity_dict['cancer.top3'] = tissue3
        cancer_expression_specificity_dict['cancer.75percentile.expression'] = tissue_75percentile
        cancer_expression_specificity_dict['cancer.median.expression'] = tissue_median
        return cancer_expression_specificity_dict
    
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

    def wilcoxon_ranksum_test(self,feature_value_df,group_df,group_col,sample_col,g1,g2):
        """在指定的2组间进行检验"""
        g1_samlist = group_df[group_df[group_col]==g1][sample_col].tolist()
        g2_samlist = group_df[group_df[group_col]==g2][sample_col].tolist()
        # 初始化结果存储
        results = []
    
        # 遍历每个基因进行检验
        for gene in feature_value_df.columns:
            # 提取两组数据
            feature_value_samlist = feature_value_df.index.tolist()
            g1_samlist = list(set(feature_value_samlist)&set(g1_samlist))
            g2_samlist = list(set(feature_value_samlist)&set(g2_samlist))
            group1 = feature_value_df.loc[g1_samlist, gene]
            group2 = feature_value_df.loc[g2_samlist, gene]
            # 执行Mann-Whitney U检验
            try:
                stat, pval = mannwhitneyu(group1, group2, alternative='two-sided')
            except ValueError:
                # 处理无法计算的情况（如数据全相同）
                stat, pval = np.nan, np.nan
    
            num_group1 = len(g1_samlist)
            num_group2 = len(g2_samlist)
            ## 计算两组间的平均值
            group1_mean = np.mean(group1)
            group2_mean = np.mean(group2)
            feature_exp_name = "Expression_type"
            if group1_mean > group2_mean:
                feature_exp = 'High'
            elif group2_mean > group1_mean :
                feature_exp = 'Low'
            else:
                feature_exp = 'Equal'
            if group2_mean == 0:
                if group1_mean !=0:
                    log2fc = np.inf
                else:
                    log2fc = 0
            else:
                if group1_mean == 0:
                    log2fc = -np.inf
                else:
                    log2fc = np.log2(group1_mean/group2_mean)
            results.append({
                'Gene': gene,
                'Statistic': stat,
                'Mean_TPM:'+g1:group1_mean,
                'Mean_TPM:'+g2:group2_mean,
                'Number_of_samples'+g1:num_group1,
                'Number_of_samples'+g2:num_group2,
                'Log2FC':log2fc,
                feature_exp_name:feature_exp,
                'PValue': pval
            })
        
        # 转换为DataFrame
        result_df = pd.DataFrame(results)
        result_df['Significant'] = result_df.apply(lambda x:"Up" if x['PValue'] < 0.05 and x[feature_exp_name]=="High" else ( "Down" if x['PValue'] < 0.05 and x[feature_exp_name]=="Low" else "No" ),axis=1)
        return result_df  

    def run_clin_wilcoxon_ranksum_test(self,clinical_df,cancertype_col,group_col,sample_col,primary_group,feature_value_df):
        cancertype_list = set(clinical_df[cancertype_col].to_list())
        cancertype_wilcoxon_test_dflist = []
        for cancertype in cancertype_list:
            cancer_clinical_df = clinical_df[clinical_df[cancertype_col]==cancertype]
            primary_disease = clinical_df['primary_disease'].tolist()[0]
            if len(set(cancer_clinical_df[group_col].tolist())) < 2: continue
            group_list = list(set(cancer_clinical_df[group_col].tolist()))
            group_list.remove(primary_group)
            normal_group = group_list[0]
            wilcoxon_test_df = self.wilcoxon_ranksum_test(feature_value_df,cancer_clinical_df,group_col,sample_col,primary_group,normal_group)
            wilcoxon_test_df['cancer_type'] = cancertype
            wilcoxon_test_df['primary_disease'] = primary_disease
            cancertype_wilcoxon_test_dflist.append(wilcoxon_test_df)
        cancertype_wilcoxon_test_df = pd.concat(cancertype_wilcoxon_test_dflist)
        cancertype_wilcoxon_test_df = cancertype_wilcoxon_test_df.reset_index(drop=True)
        return cancertype_wilcoxon_test_df

    def run_clin_wilcoxon_ranksum_test_do(self,sample_expression_df,sam_type_list):
        """
        针对Primary Tumor vs Solid Tissue Normal比较组，进行差异分析

        Args:
            sample_expression_df: 样本表达数据
            sam_type_list: 样本类型 ["Primary Tumor","Solid Tissue Normal"]
        Returns:
            dataframe: pandas dataframe
        """
        sample_expression_subset_df = sample_expression_df[sample_expression_df['sample_type'].isin(sam_type_list)].copy(deep=True)
        sample_expression_subset_df['Expression'] = sample_expression_subset_df['Expression'].apply(lambda x:float(x))
        clinical_df = sample_expression_subset_df[['sampleID','cancer_type','primary_disease','sample_type']]
        feature_value_df = sample_expression_subset_df[['sampleID','Gene_id','Gene_name','Expression']]
        feature_value_wide_df = feature_value_df.pivot(index='sampleID', columns='Gene_id', values='Expression')
        primary_group = "Primary Tumor"
        allcancertype_wilcoxon_test_df = self.run_clin_wilcoxon_ranksum_test(clinical_df,"cancer_type","sample_type","sampleID",primary_group,feature_value_wide_df)
        return allcancertype_wilcoxon_test_df
    
    def create_structured_json(self, sample_expression_df, cancer_samtype_expr_df,allcancertype_wilcoxon_test_df,allcancertype_clin_info,cancer_expression_specificity_dict):
        """
        创建结构化JSON输出
        
        Args:
            sample_expression_df: 样本表达数据
            cancer_samtype_expr_df: 不同癌症不同样本类型中位表达数据
            allcancertype_wilcoxon_test_df: 不同癌症不同样本间的wilcoxon差异分析结果 
            allcancertype_clin_info: TCGA所有癌症样本结构化表型信息
            cancer_expression_specificity_dict：癌症特异性表达信息
        Returns:
            str: JSON格式字符串
        """
        if sample_expression_df.empty or cancer_samtype_expr_df.empty:
            return json.dumps({"error": "无有效基因表达数据"}, indent=2)
        
        # 基础基因信息
        gene_info = {
            "gene_id": sample_expression_df['Gene_id'].iloc[0],
            "gene_name": sample_expression_df['Gene_name'].iloc[0],
        }
        
        # 样本级别表达数据
        sample_expression_data = []
        for _, row in sample_expression_df.iterrows():
            sample_data = {
                "sampleID": row['sampleID'],
                "patient":  row['patient'],
                "expression_value": float(row['Expression']),
                "primary_disease": row['primary_disease'],
                "cancer_type": row['cancer_type'],
                "sample_type": row['sample_type']
            }
            sample_expression_data.append(sample_data)

        
        # 基于癌症类型不同样本类型表达中位值
        cancer_samtype_expression_data = []
        for _, row in cancer_samtype_expr_df.iterrows():
            cancer_samtype_data = {
                "cancer_type": row['cancer_type'],
                "primary_disease": row['primary_disease'],
                "sample_type":row['sample_type'],
                "expression_value":float(row['Expression']),
                "expression_category": self._categorize_expression(float(row['Expression']))
            }
            cancer_samtype_expression_data.append(cancer_samtype_data)
        
        # 表达统计信息
        expression_stats = self._calculate_expression_stats(sample_expression_df, cancer_samtype_expr_df)

        # 癌症特异性表达分析信息
        specificity_info = cancer_expression_specificity_dict
        
        # 不同癌症特定样本类型间的wilcoxon差异分析
        cancer_samtype_wilcoxon_data = []
        for _, row in allcancertype_wilcoxon_test_df.iterrows():
            cancertype_wilcoxon_data = {
                "cancer_type": row['cancer_type'],
                "primary_disease":  row['primary_disease'],
                "Mean_TPM:Primary Tumor": float(row['Mean_TPM:Primary Tumor']),
                "Mean_TPM:Primary Tumor": float(row['Mean_TPM:Solid Tissue Normal']),
                "Number_of_samples in Primary Tumor": float(row['Number_of_samplesPrimary Tumor']),
                "Number_of_samples in Solid Tissue Normal": float(row['Number_of_samplesSolid Tissue Normal']),
                "Log2FC":float(row['Log2FC']),
                "Expression_type":row['Expression_type'],
                "PValue":float(row['PValue']),
                "Significant":row['Significant']
            }
            cancer_samtype_wilcoxon_data.append(cancertype_wilcoxon_data)
             
        
        # 构建完整数据结构
        structured_data = OrderedDict([
            ("gene_information", gene_info),
            ("sample_level_expression", {
                "total_samples": len(sample_expression_data),
                "samples": sample_expression_data
            }),
            ("cancer_sample_type_level_expression", {
                "cancer_sample_type_analyzed": len(cancer_samtype_expression_data),
                "cancer_sample_type_data": cancer_samtype_expression_data,
                "cancer_specificity_analysis": specificity_info
            }),
            ("expression_statistics", expression_stats),
            ("cancer_sample_type_deg_data",cancer_samtype_wilcoxon_data),
            ("cancer_clin_data",allcancertype_clin_info),
            ("data_metadata", {
                "data_source": "TCGA Analysis",
                "processing_date": pd.Timestamp.now().strftime("%Y-%m-%d"),
                "units": "TPM (Transcripts Per Million)",
                "ensembl_id": self.ensembl_id,
                "gene_name":self.genename
            })
        ])
        
        return json.dumps(structured_data, indent=2, ensure_ascii=False)
    
    def parse(self, ensembl_id=None, strategy='auto'):
        """
        主解析方法：执行完整的TCGA数据解析流程
        
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
      
        if not all([self.sample_expr_file, self.metadata_file, self.metadata_all_json_file]):
            raise ValueError("必须提供所有必要的文件路径")
        
        try:
            print(f"🔍 开始解析基因 {self.ensembl_id} 的TCGA数据...")
            
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
            # 4. 针对不同癌症类型，针对Primary Tumor vs Solid Tissue Normal间进行wilcoxon差异分析
            sam_type_list = ['Primary Tumor','Solid Tissue Normal']
            processed_allcancertype_wilcoxon_test_df = self.run_clin_wilcoxon_ranksum_test_do(processed_sample_expr,sam_type_list)
            # 5. 按不同癌症样本类型计算相应基因表达的中位值
            processed_cancer_samtype_expr = self.process_sampletype_expression(processed_sample_expr,sam_type_list)
            # 6.基于癌症Primary Tumor计算特异性高表达癌症
            cancer_expression_specificity_dict = self.cancer_specificty(processed_cancer_samtype_expr)
            # 7. 导入TCGA json表型数据
            with open(self.metadata_all_json_file, 'r', encoding='utf-8') as file:
                processed_allcancertype_clin_info = json.load(file)
            # 8. 生成结构化JSON
            json_output = self.create_structured_json(
                processed_sample_expr,
                processed_cancer_samtype_expr,
                processed_allcancertype_wilcoxon_test_df,
                processed_allcancertype_clin_info,
                cancer_expression_specificity_dict
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
            return f"TCGAParser(基因: {self.ensembl_id}, 状态: 已解析, 耗时: {self.processing_time:.2f}秒)"
        elif self.ensembl_id:
            return f"TCGAParser(基因: {self.ensembl_id}, 状态: 未解析)"
        else:
            return "TCGAParser(状态: 未初始化)"