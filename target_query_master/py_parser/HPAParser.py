import pandas as pd
import json
import requests
from collections import OrderedDict

class HPAParser:
    """
    面向对象的HPA (Human Protein Atlas) 数据解析器
    封装了从HPA API获取数据、解析和转换的功能
    """
    
    def __init__(self, ensembl_id=None):
        """
        初始化解析器
        
        Args:
            ensembl_id (str, optional): Ensembl基因ID，如果提供则自动获取数据
        """
        self.ensembl_id = ensembl_id
        self.raw_data = None
        self.parsed_data = None
        
        if ensembl_id:
            self.fetch_data(ensembl_id)
    
    def fetch_data(self, ensembl_id):
        """
        从HPA API获取数据
        
        Args:
            ensembl_id (str): Ensembl基因ID
            
        Returns:
            bool: 获取是否成功
        """
        self.ensembl_id = ensembl_id
        # HPA API URL - 使用实际的HPA数据端点
        url = f"https://www.proteinatlas.org/{ensembl_id}.tsv"
        
        try:
            # 使用pandas读取TSV数据
            df = pd.read_csv(url, sep='\t')
            
            # 检查DataFrame是否为空
            if df.empty:
                print(f"警告: 未找到Ensembl ID {ensembl_id} 的数据")
                return False
            
            # 转换为JSON格式并提取第一条记录
            json_data_str = df.to_json(orient='records', indent=2, force_ascii=False)
            json_data = json.loads(json_data_str)
            self.raw_data = json_data[0] if json_data else {}
            
            return True
            
        except pd.errors.EmptyDataError:
            print(f"错误: 接收到空数据文件")
            return False
        except pd.errors.ParserError:
            print(f"错误: 解析TSV数据失败")
            return False
        except Exception as e:
            print(f"意外错误: {str(e)}")
            return False
    
    def _get_structure_dict(self):
        """
        获取结构化信息定义
        
        Returns:
            dict: 结构化信息字典
        """
        return {
            'INFORMATION': ['Gene', 'Gene synonym', 'Ensembl', 'Chromosome', 'Position', 
                           'Gene description', 'Protein class', 'Evidence', 'Disease involvement'],
            'PROTEIN EXPRESSION AND LOCALIZATION': ['Subcellular location', 'Subcellular main location', 
                                                   'Subcellular additional location', 'Secretome location', 
                                                   'Secretome function'],
            'TISSUE RNA EXPRESSION': ['RNA tissue specificity', 'Tissue expression cluster', 
                                     'RNA tissue distribution', 'RNA tissue specificity score', 
                                     'RNA tissue specific nTPM', 'RNA brain regional specificity', 
                                     'Brain expression cluster', 'RNA brain regional distribution', 
                                     'RNA brain regional specificity score', 'RNA brain regional specific nTPM'],
            'CELL TYPE RNA EXPRESSION': ['RNA single cell type specificity', 'Single cell expression cluster', 
                                        'RNA single cell type distribution', 'RNA single cell type specificity score', 
                                        'RNA single cell type specific nTPM', 'RNA tissue cell type enrichment'], 
            'PROTEINS IN BLOOD': ['Blood expression cluster', 'RNA blood cell specificity', 
                                 'RNA blood cell distribution', 'RNA blood cell specificity score', 
                                 'RNA blood cell specific nTPM', 'RNA blood lineage specificity', 
                                 'RNA blood lineage distribution', 'RNA blood lineage specificity score', 
                                 'RNA blood lineage specific nTPM', 'Blood concentration - Conc. blood IM [pg/L]', 
                                 'Blood concentration - Conc. blood MS [pg/L]'],
            'CANCER & CELL LINES': ['Cancer prognostics', 'RNA cancer specificity', 'RNA cancer distribution', 
                                   'RNA cancer specificity score', 'RNA cancer specific FPKM', 
                                   'Cell line expression cluster', 'RNA cell line specificity', 
                                   'RNA cell line distribution', 'RNA cell line specificity score', 
                                   'RNA cell line specific nTPM'],
            'PROTEIN FUNCTION': ['Biological process', 'Molecular function']
        }
    
    def _process_cancer_prognostics(self):
        """
        特殊处理CANCER & CELL LINES类别中的Cancer prognostics字段
        
        Returns:
            str: 处理后的Cancer prognostics信息
        """
        cancer_prognostics_list = []
        
        # 遍历所有键，寻找包含"Cancer prognostics"的键
        for key, value in self.raw_data.items():
            if "Cancer prognostics" in key and "unprognostic" not in str(value):
                # 提取癌症类型信息（去掉"Cancer prognostics - "前缀）
                cancer_type = key.replace("Cancer prognostics - ", "")
                # 创建格式化的字符串
                formatted_entry = f"{cancer_type}: {value}"
                cancer_prognostics_list.append(formatted_entry)
        
        # 如果有多个结果，用逗号连接
        if cancer_prognostics_list:
            return ", ".join(cancer_prognostics_list)
        else:
            return "No significant prognostic data"
    
    def _categorize_data(self):
        """
        根据结构化信息对HPA数据进行归类
        
        Returns:
            OrderedDict: 归类后的数据
        """
        structure_dict = self._get_structure_dict()
        categorized_data = OrderedDict()
        
        # 处理每个类别
        for category, fields in structure_dict.items():
            category_data = OrderedDict()
            
            # 特殊处理CANCER & CELL LINES类别中的Cancer prognostics
            if category == "CANCER & CELL LINES" and "Cancer prognostics" in fields:
                category_data["Cancer prognostics"] = self._process_cancer_prognostics()
                
                # 处理该类别中的其他字段
                for field in fields:
                    if field != "Cancer prognostics" and field in self.raw_data:
                        category_data[field] = self.raw_data[field]
            
            else:
                # 处理其他类别
                for field in fields:
                    if field in self.raw_data:
                        category_data[field] = self.raw_data[field]
            
            categorized_data[category] = category_data
        categorized_data['THE HUMAN PROTEIN ATLAS'] = f'https://www.proteinatlas.org/{self.ensembl_id}-{self.raw_data['Gene']}'
        
        return categorized_data
    
    def _convert_to_json(self, ensure_ascii=False):
        """
        将OrderedDict转换为JSON字符串
        
        Args:
            ensure_ascii (bool): 是否确保ASCII编码
            
        Returns:
            str: JSON格式字符串
        """
        try:
            json_string = json.dumps(
                self.parsed_data, 
                indent=2, 
                ensure_ascii=ensure_ascii,
                sort_keys=False
            )
            return json_string
        except TypeError as e:
            print(f"序列化错误: {e}")
            
            def fallback_serializer(obj):
                if hasattr(obj, '__dict__'):
                    return obj.__dict__
                return str(obj)
            
            return json.dumps(
                self.parsed_data, 
                indent=2, 
                ensure_ascii=ensure_ascii,
                default=fallback_serializer
            )
    
    def parse(self):
        """
        主解析方法：执行完整的数据解析流程
        
        Returns:
            str: 格式化后的JSON字符串
        """
        if not self.raw_data:
            if not self.ensembl_id:
                raise ValueError("未提供Ensembl ID且无原始数据")
            if not self.fetch_data(self.ensembl_id):
                raise Exception("数据获取失败")
        
        # 解析数据
        self.parsed_data = self._categorize_data()
        
        return self._convert_to_json()
    
    def get_raw_data(self):
        """获取原始数据"""
        return self.raw_data
    
    def get_parsed_data(self):
        """获取解析后的数据"""
        return self.parsed_data
    
    def save_to_file(self, filename, ensure_ascii=False):
        """
        将解析结果保存到文件
        
        Args:
            filename (str): 文件名
            ensure_ascii (bool): 是否确保ASCII编码
        """
        if not self.parsed_data:
            self.parse()
        
        json_output = self._convert_to_json(ensure_ascii)
        with open(filename, 'w', encoding='utf-8') as f:
            f.write(json_output)
    
    def __str__(self):
        """字符串表示"""
        if self.parsed_data:
            return f"HPAParser(ID: {self.ensembl_id}, 已解析)"
        elif self.raw_data:
            return f"HPAParser(ID: {self.ensembl_id}, 原始数据已加载)"
        else:
            return f"HPAParser(ID: {self.ensembl_id}, 未初始化)"