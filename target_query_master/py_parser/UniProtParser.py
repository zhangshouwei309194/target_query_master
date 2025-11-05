import requests
import json
from collections import OrderedDict

class UniProtParser:
    """
    面向对象的UniProt数据解析器
    封装了从UniProt API获取数据、解析和转换的功能
    """
    
    def __init__(self, uniprot_id=None):
        """
        初始化解析器
        
        Args:
            uniprot_id (str, optional): UniProt ID，如果提供则自动获取数据
        """
        self.uniprot_id = uniprot_id
        self.raw_data = None
        self.parsed_data = None
        
        if uniprot_id:
            self.fetch_data(uniprot_id)
    
    def fetch_data(self, uniprot_id):
        """
        从UniProt API获取数据
        
        Args:
            uniprot_id (str): UniProt ID
            
        Returns:
            bool: 获取是否成功
        """
        self.uniprot_id = uniprot_id
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}"
        
        try:
            response = requests.get(url)
            response.raise_for_status()
            self.raw_data = response.json()
            return True
        except requests.exceptions.RequestException as e:
            print(f"请求错误: {e}")
            return False
        except json.JSONDecodeError as e:
            print(f"JSON解析错误: {e}")
            return False
    
    @staticmethod
    def _update_dict_with_concatenation(a_dict, new_data):
        """
        静态方法：更新字典，支持值拼接
        
        Args:
            a_dict (dict): 目标字典
            new_data: 新数据
            
        Returns:
            dict: 更新后的字典
        """
        if isinstance(new_data, dict):
            items = new_data.items()
        else:
            items = new_data
        
        for key, value in items:
            if key in a_dict:
                a_dict[key] += value
            else:
                a_dict[key] = value
        return a_dict
    
    def _parse_go_domain(self):
        """
        解析GO功能域信息（分子功能、生物过程、细胞组分、结构域）
        
        Returns:
            OrderedDict: GO功能域字典
        """
        go_domain_dict = OrderedDict()
        key_content = ['Molecular function', 'Biological process', 
                      'Cellular component', 'Domain']
        
        # 初始化空列表
        for key in key_content:
            go_domain_dict[key] = []
        
        # 解析关键词信息
        uniprot_info_keywords = self.raw_data['keywords']      
        for keywords_dict in uniprot_info_keywords:
            if keywords_dict['category'] in key_content:
                key_info_dict = {
                    keywords_dict['category']: [keywords_dict['name']]
                }
                go_domain_dict = self._update_dict_with_concatenation(
                    go_domain_dict, key_info_dict
                )
        return go_domain_dict
    
    def _parse_comment_section(self, key_content_dict):
        """
        解析comment信息
        
        Args:
            key_content_dict (dict): 键值映射字典
            
        Returns:
            OrderedDict: 解析后的评论信息
        """
        info_dict = OrderedDict()
        uniprot_info = self.raw_data['comments']
        
        for subpart_uniprot_dict in uniprot_info:
            if subpart_uniprot_dict['commentType'] in key_content_dict:
                parts_info_dict = {}
                search_dict = key_content_dict[subpart_uniprot_dict['commentType']]
                
                for key in search_dict:
                    if key != "note":
                        content_key = key
                        value_key = search_dict[key]
                        if content_key == 'reaction':
                            parts_info_dict[subpart_uniprot_dict['commentType']] = [
                                subpart_uniprot_dict[content_key][value_key]
                            ]
                        else:
                            parts_info_dict[subpart_uniprot_dict['commentType']] = [
                                subpart_uniprot_dict[content_key][0][value_key]
                            ]
                    else:
                        content1_key = key
                        for subkey in search_dict[key]:
                            content2_key = subkey
                            value_key = search_dict[key][content2_key]
                            if content1_key in subpart_uniprot_dict:
                                parts_info_dict[subpart_uniprot_dict['commentType']] = [
                                    subpart_uniprot_dict[content1_key][content2_key][0][value_key]
                                ]
                
                info_dict = self._update_dict_with_concatenation(info_dict, parts_info_dict)
        
        return info_dict
    
    def _parse_subcellular_location(self):
        """
        解析亚细胞定位信息
        
        Returns:
            OrderedDict: 亚细胞定位信息
        """
        info_dict = OrderedDict()
        info_dict.setdefault('Total', {})
        info_dict['Total'].setdefault('note', [])
        info_dict['Total'].setdefault('location', [])
        info_dict.setdefault('Isoform', {})
        
        uniprot_info = self.raw_data['comments']
        
        for subpart_uniprot_dict in uniprot_info:
            if subpart_uniprot_dict['commentType'] == "SUBCELLULAR LOCATION":
                if 'molecule' not in subpart_uniprot_dict:
                    # 总定位信息
                    if 'note' in subpart_uniprot_dict:
                        info_dict['Total']['note'].append(
                            subpart_uniprot_dict['note']['texts'][0]['value']
                        )
                    if 'subcellularLocations' in subpart_uniprot_dict:
                        for location in subpart_uniprot_dict['subcellularLocations']:
                            info_dict['Total']['location'].append(
                                location['location']['value']
                            )
                    else:
                        info_dict['Total']['location'].append('')
                elif 'molecule' in subpart_uniprot_dict:
                    # 异构体特定定位信息
                    isoform_id = subpart_uniprot_dict['molecule']
                    info_dict['Isoform'].setdefault(isoform_id, {})
                    info_dict['Isoform'][isoform_id].setdefault('note', [])
                    info_dict['Isoform'][isoform_id].setdefault('location', [])
                    
                    if 'note' in subpart_uniprot_dict:
                        info_dict['Isoform'][isoform_id]['note'].append(
                            subpart_uniprot_dict['note']['texts'][0]['value']
                        )
                    if 'subcellularLocations' in subpart_uniprot_dict:
                        for location in subpart_uniprot_dict['subcellularLocations']:
                            info_dict['Isoform'][isoform_id]['location'].append(
                                location['location']['value']
                            )
                    else:
                        info_dict['Isoform'][isoform_id]['location'].append('')
        
        return info_dict
    
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
            if not self.uniprot_id:
                raise ValueError("未提供UniProt ID且无原始数据")
            if not self.fetch_data(self.uniprot_id):
                raise Exception("数据获取失败")
        
        # 1. 解析基础信息
        basis_info_dict = OrderedDict()
        basis_info_dict['Protein'] = self.raw_data['proteinDescription']['recommendedName']['fullName']['value']
        basis_info_dict['primaryAccession'] = self.raw_data['primaryAccession'] 
        basis_info_dict['Amino acids'] = self.raw_data['sequence']['length']
        basis_info_dict['Molweight'] = self.raw_data['sequence']['molWeight']
        basis_info_dict['Protein existence'] = self.raw_data['proteinExistence']
        basis_info_dict['Status'] = self.raw_data['entryType']
        basis_info_dict['Annotation score'] = self.raw_data['annotationScore']
        
        # 2. 解析GO功能域
        go_domain_dict = self._parse_go_domain()
        
        # 3. 解析function部分
        function_key_content_dict = {
            "FUNCTION": {"texts": "value"},
            "CAUTION": {"texts": "value"},
            "CATALYTIC ACTIVITY": {"reaction": "name"},
            "COFACTOR": {"cofactors": "name", "note": {"texts": "value"}},
            "ACTIVITY REGULATION": {"texts": "value"},
            'PATHWAY': {"texts": "value"}
        }
        function_comment_dict = self._parse_comment_section(function_key_content_dict)
        
        # 4. 解析亚细胞定位
        subcellular_dict = self._parse_subcellular_location()
        
        # 5. 解析表达信息
        expression_key_content_dict = {
            "TISSUE SPECIFICITY": {"texts": "value"},
            "INDUCTION": {"texts": "value"}
        }
        expression_comment_dict = self._parse_comment_section(expression_key_content_dict)
        
        # 6. 获取序列信息
        sequence_dict = {'Sequence': self.raw_data['sequence']['value']}
        
        # 7. 整合所有数据
        self.parsed_data = OrderedDict()
        self.parsed_data['Basis Info'] = basis_info_dict
        
        # 合并功能信息和GO信息
        for go_category in ['Molecular function', 'Biological process', 'Cellular component']:
            function_comment_dict[go_category] = go_domain_dict[go_category]
        
        self.parsed_data['Function'] = function_comment_dict
        self.parsed_data['Subcellulary Location'] = subcellular_dict
        self.parsed_data['Expression'] = expression_comment_dict
        self.parsed_data['Domain'] = go_domain_dict['Domain']
        self.parsed_data['Sequence'] = sequence_dict['Sequence']
        self.parsed_data['UniProt'] = 'https://www.uniprot.org/uniprotkb/' + self.uniprot_id
        
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
            return f"UniProtParser(ID: {self.uniprot_id}, 已解析)"
        elif self.raw_data:
            return f"UniProtParser(ID: {self.uniprot_id}, 原始数据已加载)"
        else:
            return f"UniProtParser(ID: {self.uniprot_id}, 未初始化)"