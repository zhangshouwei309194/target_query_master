import pandas as pd
import numpy as np
import json
import os
import time
from collections import OrderedDict

class GeneInfoAggregate:
    """
    面向对象的GeneInfoAggregate 基因信息注释程序
    封装了HGNC，NCBI，Uniprot和OncoKB的基因相关注释信息    
    """

    def __init__(self,
                 hgnc_gene_file=None, 
                 ncbi_gene_file=None,
                 unipro_gene_file=None,
                 oncokb_gene_file=None):
        """
        初始化GeneInfoAggregate解析器
        
        Args:
            hgnc_gene_file (str): HGNC 基因信息（human gene nomenclature,https://www.genenames.org/）
            ncbi_gene_file (str): NCBI基因信息（https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info） 
            unipro_gene_file (str): UniProt（https://www.uniprot.org/uniprotkb?query=*&facets=model_organism%3A9606）
            oncokb_gene_file (str): OncoKB Oncogene/TSG（https://www.oncokb.org/cancer-genes）
        """
        self.hgnc_gene_file = hgnc_gene_file
        self.ncbi_gene_file = ncbi_gene_file
        self.unipro_gene_file = unipro_gene_file
        self.oncokb_gene_file = oncokb_gene_file

            # 初始化属性
        self.parsed_data = None
        self.processing_time = None
        
        # 验证文件存在性
        if hgnc_gene_file and not os.path.exists(hgnc_gene_file):
            raise FileNotFoundError(f"HGNC基因信息文件不存在: {sample_expr_file}")
        if ncbi_gene_file and not os.path.exists(ncbi_gene_file):
            raise FileNotFoundError(f"NCBI基因信息文件不存在: {tissue_expr_file}")
        if unipro_gene_file and not os.path.exists(unipro_gene_file):
            raise FileNotFoundError(f"UniProt信息文件不存在: {metadata_file}")
        if oncokb_gene_file and not os.path.exists(oncokb_gene_file):
            raise FileNotFoundError(f"OncoKB信息文件不存在: {oncokb_gene_file}")        

    ## 提取NCBI中dbXrefs信息标签
    def _extract_tag_value(self, input_string, tag_name):
        """
        从格式为"key1:value1|key2:value2|..."的字符串中提取指定标签对应的值
        
        参数:
            input_string (str): 包含多个键值对的字符串
            tag_name (str): 要提取的标签名（如"HGNC"、"Ensembl"等）
            
        返回:
            str 或 None: 成功则返回标签对应的值，未找到则返回None
        """
        # 构建动态正则表达式模式
        pattern = rf'{tag_name}:([^|]+)'
        
        match = re.search(pattern, input_string)
        return match.group(1) if match else np.nan

    def _merge_synonyms(self,row):
        """
        合并Previous symbols, Alias symbols和NCBI.Synonyms三列为Synonyms列
        如果三列都为空，则标记为缺失值NaN
        
        参数:
            row: pandas DataFrame的一行数据
        
        返回:
            合并后的字符串或NaN
        """
        # 获取三列的值
        previous = row['Previous symbols']
        alias = row['Alias symbols']
        ncbi = row['NCBI.Synonyms']
        
        # 收集所有非空值
        synonyms_parts = []
        
        # 检查每列是否有有效值（非NaN且非空字符串）
        if pd.notna(previous) and str(previous).strip() != '':
            prvious_list = previous.split(', ')
            synonyms_parts += prvious_list
        
        if pd.notna(alias) and str(alias).strip() != '':
            alias_list = alias.split(', ')
            synonyms_parts += alias_list
        
        if pd.notna(ncbi) and str(ncbi).strip() != '':
            if ncbi != "-":
                ncbi_list = ncbi.split('|')
                synonyms_parts += ncbi_list
        
        if synonyms_parts:
            return ','.join(synonyms_parts)
        else:
            return np.nan
    
    def integrate_ncbi_to_hgnc(self, hgnc_df, ncbi_df):
        """
        将NCBI数据整合到HGNC数据中
        
        参数:
            hgnc_df: HGNC基因数据信息
            ncbi_df: NCBI基因数据信息
        
        返回:
            整合后的HGNC数据框
        """
        # 创建结果数据框的副本
        result_df = hgnc_df.copy()
        hgnc_gene_info_df = hgnc_df[hgnc_df['Status']=='Approved']
        ncbi_df['Ensembl'] = ncbi_df['dbXrefs'].apply(lambda x:self._extract_tag_value(x,"Ensembl"))
        ncbi_df['HGNC ID'] = ncbi_df['dbXrefs'].apply(lambda x:self._extract_tag_value(x,"HGNC"))
        ncbi_gene_info_simplify_df = ncbi_df[['GeneID','Synonyms','type_of_gene','Ensembl','HGNC ID']]
        ncbi_gene_info_simplify_df = ncbi_gene_info_simplify_df.rename(columns={'Synonyms':'NCBI.Synonyms'})
    
        hgnc_gene_info_addncbi_df = pd.merge(hgnc_gene_info_df,ncbi_gene_info_simplify_df,on='HGNC ID',how='left')
        hgnc_gene_info_addncbi_df['Protein_coding'] = hgnc_gene_info_addncbi_df['type_of_gene'].apply(lambda x:"YES" if x=="protein-coding" else "NO" )

        ## 合并gene alias信息
        hgnc_gene_info_addncbi_df['Synonyms'] = hgnc_gene_info_addncbi_df.apply(self._merge_synonyms, axis=1)
    
        col_order = ['HGNC ID',
                     'Approved Symbol',
                     'Approved name',
                     'Status',
                     'Synonyms',
                     'Ensembl',
                     'GeneID',
                     'Chromosome',
                     'RefSeq IDs',
                     'Gene group name',
                     'OMIM ID(supplied by OMIM)',
                     'type_of_gene',
                     'Protein_coding']
    
        hgnc_gene_info_addncbi_df = hgnc_gene_info_addncbi_df[col_order]
        return hgnc_gene_info_addncbi_df

    def integrate_uniprot_to_hgnc(self, hgnc_df, uniprot_df):
        """
        优化版本 - 将UniProt数据整合到HGNC数据中
        
        参数:
            hgnc_df: HGNC数据框
            uniprot_df: UniProt数据框
        
        返回:
            整合后的HGNC数据框和匹配统计信息
        """
        # 创建结果数据框的副本
        result_df = hgnc_df.copy()
        
        # 添加新的列，初始值为NaN
        result_df['Uniprot Entry'] = np.nan
        result_df['Uniprot Reviewed'] = np.nan
        result_df['Uniprot Protein names'] = np.nan
        result_df['Uniprot Protein length'] = np.nan
        result_df['Match Type'] = 'No Match'
        
        # 预处理UniProt数据 - 创建基因名到记录的映射
        gene_to_uniprot = {}
        
        for idx, row in uniprot_df.iterrows():
            gene_names = str(row['Gene Names']).strip()
            if gene_names and gene_names != 'nan':
                for gene in re.split(r'[;\s]+', gene_names):
                    gene = gene.strip()
                    if gene and gene not in gene_to_uniprot:
                        gene_to_uniprot[gene] = {
                            'Entry': row['Entry'],
                            'Reviewed': row['Reviewed'],
                            'Protein names': row['Protein names'],
                            'Length': row['Length']
                        }
        
        # 统计信息
        match_stats = {
            'Approved Symbol Matches': 0,
            'Synonyms Matches': 0,
            'No Matches': 0
        }
        
        # 第一步：基于Approved Symbol进行匹配
        approved_symbols = result_df['Approved Symbol'].astype(str).str.strip()
        
        # 使用向量化操作查找匹配 - 修复：正确处理NaN值
        approved_matches = approved_symbols.map(gene_to_uniprot)
        
        # 处理匹配结果 - 修复：使用pd.notna()检查而非is not None
        for idx, match_data in approved_matches.items():
            if pd.notna(match_data):  # 修复：使用pd.notna()而不是is not None
                result_df.loc[idx, 'Uniprot Entry'] = match_data['Entry']
                result_df.loc[idx, 'Uniprot Reviewed'] = match_data['Reviewed']
                result_df.loc[idx, 'Uniprot Protein names'] = match_data['Protein names']
                result_df.loc[idx, 'Uniprot Protein length'] = match_data['Length']
                result_df.loc[idx, 'Match Type'] = 'Approved Symbol'
                match_stats['Approved Symbol Matches'] += 1
        
        # 第二步：基于Synonyms进行匹配
        unmatched_mask = result_df['Match Type'] == 'No Match'
        unmatched_indices = result_df[unmatched_mask].index
        
        for idx in unmatched_indices:
            synonyms = result_df.loc[idx, 'Synonyms']  # 修复：使用loc而非at保持一致性
            
            if pd.isna(synonyms) or str(synonyms).strip() == '':
                match_stats['No Matches'] += 1
                continue
            
            synonym_list = [s.strip() for s in str(synonyms).split(',')]
            matched = False
            
            for synonym in synonym_list:
                if synonym and synonym in gene_to_uniprot:
                    match_data = gene_to_uniprot[synonym]
                    result_df.loc[idx, 'Uniprot Entry'] = match_data['Entry']
                    result_df.loc[idx, 'Uniprot Reviewed'] = match_data['Reviewed']
                    result_df.loc[idx, 'Uniprot Protein names'] = match_data['Protein names']
                    result_df.loc[idx, 'Uniprot Protein length'] = match_data['Length']
                    result_df.loc[idx, 'Match Type'] = 'Synonyms'
                    match_stats['Synonyms Matches'] += 1
                    matched = True
                    break
            
            if not matched:
                match_stats['No Matches'] += 1
        del result_df['Match Type']
        return result_df, match_stats    

    def integrate_oncokb_to_hgnc(self, hgnc_df, oncokb_df):
        """
        优化版本 - 将UniProt数据整合到HGNC数据中
        
        参数:
            hgnc_df: HGNC数据框
            oncokb_df: Oncokb基因属性数据框(Oncogene/TSG)
        
        返回:
            整合后的HGNC数据框和匹配统计信息
        """
        # 创建结果数据框的副本
        result_df = hgnc_df.copy() 
        
        # 添加新的列，初始值为NaN
        result_df['OncoKB Gene Type'] = np.nan
        result_df['OncoKB GRCh38 Isoform'] = np.nan
        result_df['OncoKB GRCh38 RefSeq'] = np.nan
        result_df['Match Type'] = 'No Match'
        
        # 预处理Oncokb数据 - 创建基因名到记录的映射
        gene_to_oncokb = {}
        oncokb_df['Gene Names'] = oncokb_df.apply(lambda x:','.join([x['Hugo Symbol'],x['Gene Aliases']]) if not pd.isnull(x['Gene Aliases']) else x['Hugo Symbol'],axis=1)
        
        for idx, row in oncokb_df.iterrows():
            gene_names = str(row['Gene Names']).strip()
            if gene_names and gene_names != 'nan':
                for gene in re.split(r'[,\s]+', gene_names):
                    gene = gene.strip()
                    if gene and gene not in gene_to_oncokb:
                        gene_to_oncokb[gene] = {
                            'Gene Type': row['Gene Type'],
                            'GRCh38 Isoform': row['GRCh38 Isoform'],
                            'GRCh38 RefSeq': row['GRCh38 RefSeq']
                        }
        
        # 统计信息
        match_stats = {
            'Approved Symbol Matches': 0,
            'Synonyms Matches': 0,
            'No Matches': 0
        }
        
        # 第一步：基于Approved Symbol进行匹配
        approved_symbols = result_df['Approved Symbol'].astype(str).str.strip()
        
        # 使用向量化操作查找匹配 - 修复：正确处理NaN值
        approved_matches = approved_symbols.map(gene_to_oncokb)
        
        # 处理匹配结果 - 修复：使用pd.notna()检查而非is not None
        for idx, match_data in approved_matches.items():
            if pd.notna(match_data):  # 修复：使用pd.notna()而不是is not None
                result_df.loc[idx, 'OncoKB Gene Type'] = match_data['Gene Type']
                result_df.loc[idx, 'OncoKB GRCh38 Isoform'] = match_data['GRCh38 Isoform']
                result_df.loc[idx, 'OncoKB GRCh38 RefSeq'] = match_data['GRCh38 RefSeq']
                result_df.loc[idx, 'Match Type'] = 'Approved Symbol'
                match_stats['Approved Symbol Matches'] += 1
        
        # 第二步：基于Synonyms进行匹配
        unmatched_mask = result_df['Match Type'] == 'No Match'
        unmatched_indices = result_df[unmatched_mask].index
        
        for idx in unmatched_indices:
            synonyms = result_df.loc[idx, 'Synonyms']  # 修复：使用loc而非at保持一致性
            
            if pd.isna(synonyms) or str(synonyms).strip() == '':
                match_stats['No Matches'] += 1
                continue
            
            synonym_list = [s.strip() for s in str(synonyms).split(',')]
            matched = False
            
            for synonym in synonym_list:
                if synonym and synonym in gene_to_oncokb:
                    match_data = gene_to_oncokb[synonym]
                    result_df.loc[idx, 'OncoKB Gene Type'] = match_data['Gene Type']
                    result_df.loc[idx, 'OncoKB GRCh38 Isoform'] = match_data['GRCh38 Isoform']
                    result_df.loc[idx, 'OncoKB GRCh38 RefSeq'] = match_data['GRCh38 RefSeq']
                    result_df.loc[idx, 'Match Type'] = 'Synonyms'
                    match_stats['Synonyms Matches'] += 1
                    matched = True
                    break
            
            if not matched:
                match_stats['No Matches'] += 1
        del result_df['Match Type']
        return result_df, match_stats
    
    
    def parse(self):
        """
        主解析方法：执行完整的GeneInfoAggregate数据解析流程
            
        Returns:
            dataframe: 注释后的基因信息数据框文件
        """
        start_time = time.time()
        
        if not all([self.hgnc_gene_file, self.ncbi_gene_file, self.unipro_gene_file,self.oncokb_gene_file]):
            raise ValueError("必须提供所有必要的文件路径")
        
        try:
            print(f"🔍 开始基因信息注释")        
            
            # 1. 加载基因注释信息数据
            hgnc_gene_info_df = pd.read_table(self.hgnc_gene_file,header=0,sep='\t',low_memory=False)
            ncbi_gene_info_df = pd.read_table(self.ncbi_gene_file,header=0,sep='\t',low_memory=False)
            uniprot_geneinfo_df = pd.read_table(self.unipro_gene_file,header=0,sep='\t')          
            oncokb_gene_info_df = pd.read_table(self.oncokb_gene_file,header=0,sep='\t',low_memory=False)
            
            # 2. 注释NCBI的管线信息
            hgnc_gene_info_anno_df = self.integrate_ncbi_to_hgnc(hgnc_gene_info_df, ncbi_gene_info_df)
            
            # 3. 注释UniProt信息
            hgnc_gene_info_anno_df,uniprot_map_status = self.integrate_uniprot_to_hgnc(hgnc_gene_info_anno_df, uniprot_geneinfo_df)
            ##
            print ("UniProt注释")
            for key, value in uniprot_map_status.items():
                print(f"{key}: {value}")

            # 4. 注释OncoKB信息
            hgnc_gene_info_anno_df,oncokb_map_status = self.integrate_oncokb_to_hgnc(hgnc_gene_info_anno_df, oncokb_gene_info_df)
            ##
            print ("OncoKB注释")
            for key, value in oncokb_map_status.items():
                print(f"{key}: {value}")

            hgnc_gene_info_anno_df['GeneID'] = hgnc_gene_info_anno_df['GeneID'].astype('Int64')
            hgnc_gene_info_anno_df['Uniprot Protein length'] = hgnc_gene_info_anno_df['Uniprot Protein length'].astype('Int64')

            self.processing_time = time.time() - start_time
            print(f"✅ 基因信息注释完成，耗时: {self.processing_time:.2f}秒")            
            
            self.parsed_data = hgnc_gene_info_anno_df
            return hgnc_gene_info_anno_df
        except Exception as e:
            error_msg = f"注释基因信息时发生错误: {str(e)}"
            print(f"❌ {error_msg}")
            return json.dumps({"error": error_msg}, indent=2)
            
    def get_processing_time(self):
        """获取处理时间"""
        return self.processing_time
    
    def get_parsed_data(self):
        """获取解析后的数据"""
        return self.parsed_data
    
    def save_to_file(self, filename):
        """
        将解析结果保存到数据表中
        
        Args:
            filename (str): 输出文件名
        """
        if self.parsed_data is None:
            self.parse()

        anno_table = self.parsed_data
        anno_table.to_csv(filename,header=True,index=False,sep='\t')
        print(f"💾 数据已保存到: {filename}")      
    
    def __str__(self):
        """字符串表示"""
        if self.parsed_data:
            return f"GeneInfoAggregate(状态: 已解析, 耗时: {self.processing_time:.2f}秒)"
        else:
            return "GeneInfoAggregate(状态: 未初始化)"
