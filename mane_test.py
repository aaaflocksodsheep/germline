import pandas as pd
import pysam
import copy
import csv

class ApplyMane:
    
    '''
    Date : 2023-10-06
    Description : vep annotation된 vcf 파일의 CSQ 필드를 mane transcript 고려하여 작성
    '''
    
    def __init__(self):
        self.mane = pd.read_csv('/data/script/lims/Hereditary_Disease/Gene_Transcript.csv')
        self.ghr = pd.read_csv('/data/ref/hg19/hereditary_disease/vep/GHR/ghr_summary.csv')
        self.omim = pd.read_csv('/data/ref/hg19/hereditary_disease/vep/OMIM/omim.csv')

    def mane_dict(self):
        # mane dataframe 중복제거, Transcript 컬럼이 NaN인 경우, '.'인 경우 제거
        gene_transcript = self.mane.drop_duplicates()
        gene_transcript = gene_transcript.loc[(gene_transcript['Transcript'].notna()) & (gene_transcript['Transcript'] != '.')]

        # mane transcript가 1개인 gene(=>gene_counts_1)과 1개 이상인 gene(=>gene_counts_multi) 분리
        gene_counts = gene_transcript.Gene_Name.value_counts()
        gene_mane_1 = gene_counts[gene_counts == 1].index.to_list()
        gene_mane_multi = gene_counts[gene_counts > 1].index.to_list()
        #print(f'mane 1개 : {len(gene_mane_1)}, mane 1개 이상 : {len(gene_mane_multi)}')

        # mane transcript가 1개인 gene(=>gene_counts_1)과 1개 이상인 gene(=>gene_counts_multi)에 대한 dictionary 생성
        gene_mane_1_dict = {}
        gene_mane_multi_dict = {}
        
        for index, row in gene_transcript.iterrows():
            # gene_transcript에서 Gene_Name이 gene_mane_1에 해당하는 것만 추출
            if row['Gene_Name'] in gene_mane_1:
                gene_mane_1_dict[row['Gene_Name']] = row['Transcript']
    
            # gene_transcript에서 Gene_Name이 gene_counts_multi에 해당하는 것만 추출
            elif row['Gene_Name'] in gene_mane_multi:
                if not row['Gene_Name'] in gene_mane_multi_dict:
                    gene_mane_multi_dict[row['Gene_Name']] = []
                gene_mane_multi_dict[row['Gene_Name']].append(row['Transcript'])
        
        self.mane_dict_1 = gene_mane_1_dict
        self.mane_dict_multi = gene_mane_multi_dict

    def parse_vep(self, record_info, new_df):
        '''
        Define a function to parse VEP annotations from the INFO field of a VCF record
        '''
        
        info_field_vep = record_info.get('CSQ', None)
        #print(">> info_field_vep list : ", info_field_vep)
    
        if info_field_vep:
            for csq_entry in info_field_vep:
                csq_fields = csq_entry.split('|')
                impact = csq_fields[2]
                gene_name = csq_fields[3]
                transcript_id = csq_fields[6]
                    
                if impact != 'MODIFIER':
                    ## 1) mane가 1개인 gene dictionary에서 key와 value가 일치하는 경우
                    if transcript_id == self.mane_dict_1.get(gene_name):
                        new_df['IMPACT_MANE'].append(impact)
                        if len(new_df['CSQ_MANE']) < 5:
                            new_df['CSQ_MANE'].append(csq_entry)                       
                
                    ## 2) gene_counts_max_dict의 key와 value가 일치하는 경우
                    elif self.mane_dict_multi.get(gene_name) is not None and transcript_id in self.mane_dict_multi.get(gene_name):
                        new_df['IMPACT_MANE'].append(impact)
                        if len(new_df['CSQ_MANE']) < 5:
                            new_df['CSQ_MANE'].append(csq_entry)
                        
                    ## 3) mane transcript가 없는 경우                    
                    else:
                        new_df['IMPACT_NOT_MANE'].append(impact)
                        if len(new_df['CSQ_NOT_MANE']) < 5:
                            new_df['CSQ_NOT_MANE'].append(csq_entry)
                                   
            # CSQ 하나만 선택
            if new_df['CSQ_MANE']:
                new_df['CSQ_selected'] = new_df['CSQ_MANE'][0]
            elif new_df['CSQ_NOT_MANE']:
                new_df['CSQ_selected'] = new_df['CSQ_NOT_MANE'][0]
            else:
                return
                
            # new_df에 column명은 CSQ_description, 값은 CSQ_selected를 "|"로 분리한 값을 넣어줌
            if new_df['CSQ_selected']:
                for i in range(len(self.CSQ_description)):
                    if self.CSQ_description[i] == 'EXON' or self.CSQ_description[i] == 'INTRON':
                        new_df[self.CSQ_description[i]] = "'"+new_df['CSQ_selected'].split('|')[i]
                    else:
                        new_df[self.CSQ_description[i]] = new_df['CSQ_selected'].split('|')[i]
                
            del new_df['CSQ_selected']

            # CSQ 리스트를 string으로 변환
            new_df['CSQ_MANE'] = 'CSQ='.join([''] + new_df['CSQ_MANE'])
            new_df['IMPACT_MANE'] = ';'.join(new_df['IMPACT_MANE'])
            new_df['CSQ_NOT_MANE'] = 'CSQ='.join([''] + new_df['CSQ_NOT_MANE'])
            new_df['IMPACT_NOT_MANE'] = ';'.join(new_df['IMPACT_NOT_MANE'])
                
            # new_df을 return
            return new_df.values()

    def apply_parse_vep(self, vcf_file, output_file):
        '''
        Define a function to apply the parse_vep function to a VCF file
        '''
        self.mane_dict()
        
        # Open the VCF file
        vcf_file = pysam.VariantFile(vcf_file, 'r')
        
        # vcf_file의 header에 있는 CSQ 정보를 가져옴 -> (output file에 header로 사용, parse_vep 함수에서 사용)
        CSQ_description = vcf_file.header.info.get('CSQ').description
        self.CSQ_description = CSQ_description.split('Format: ')[1].split('|')
                
        # Write the header to the output file
        with open(output_file, mode='w', newline='') as f:
            wf = csv.writer(f, quoting=csv.QUOTE_NONE)
            header = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'GENOTYPE', 'DEPTH', 'VAF', 'INFO', 'IMPACT_MANE', 'IMPACT_NOT_MANE', 'CSQ_MANE', 'CSQ_NOT_MANE']
            header += self.CSQ_description
            wf.writerow((map(str, header)))
        
        # Iterate through the records in the VCF file
        for record in vcf_file:

            ## 1) record.samples에서 [ GENOTYPE, DEPTH, VAF ] 작성을 위한 정보 추출
            for sample in record.samples:
                sample_name = sample
                sample_values = record.samples[sample_name].values()
                sample_GT = '/'.join(map(str, record.samples[sample_name]['GT']))
                # ./. 이면 unknown, 0/1 이면 HET, 1/1 이면 HOM 를 dictionary 사용해서 적용
                GT_dict = {"./.": "unknown", "0/1": "HET", "1/1": "HOM"}
                sample_GT = GT_dict.get(sample_GT)
                sample_AD = record.samples[sample_name]['AD'][1]               
                sample_DP = record.samples[sample_name]['DP']

            ## 2) record.info를 복사하여 CSQ를 제외한 정보 추출 -> [ INFO ] 작성을 위한 정보 추출
            re_info = dict(record.info)
            if 'CSQ' in re_info:
                del re_info['CSQ']
            re_info = ';'.join([f'{k}={round(v[0],3)}' if isinstance(v, tuple) else f'{k}={round(v,3)}' for k, v in re_info.items()])

            ## 3) record에서 [ CHROM, POS, ID, REF, ALT, QUAL, FILTER, GENOTYPE, DEPTH, VAF, INFO(CSQ 제외)] 정보 추출하여 작성
            new_df = {'CHROM': record.chrom, 
                      'POS': record.pos, 
                      'ID' : record.id,
                      'REF' : record.ref,
                      'ALT' : record.alts[0],
                      'QUAL' : round(record.qual, 2),
                      'FILTER' : record.filter.keys()[0],
                      'GENOTYPE' : sample_GT,
                      'DEPTH' : sample_DP,
                      'VAF' : 0 if sample_DP == 0 else sample_AD/sample_DP*100,
                      'INFO' : re_info,
                      'IMPACT_MANE': [], 'IMPACT_NOT_MANE': [], 'CSQ_MANE' : [], 'CSQ_NOT_MANE': []}
            
            ## 4) record.info에서 [ CSQ ] 정보 추출하여 파싱 (parse_vep 함수 사용)
            with open(output_file, mode='a', newline='') as f:
                wf = csv.writer(f, quoting=csv.QUOTE_NONE, escapechar='\\')
                # 해당 line의 INFO 필드에서 CSQ 정보 파싱
                vep_annotations = self.parse_vep(record.info, new_df)
                # vep_annotations가 None이 아닌 경우에만 write
                if vep_annotations: 
                    wf.writerow(map(str,vep_annotations))
            
        # Close the VCF file
        vcf_file.close()

    def ghr_omim_annotation(self, input_file):
        vep_df = pd.read_csv(input_file)
        # GHR annotation
        annotation_df = pd.merge(vep_df, self.ghr, left_on = "SYMBOL", right_on = "gene", how="left")
        # Omim annotation
        omim_df = self.omim.replace('\n', ' ', regex=True)
        annotation_df = pd.merge(annotation_df, omim_df, left_on = "SYMBOL", right_on = "Gene", how="left")
        annotation_df.to_csv(input_file, sep='\t', index=False)

if __name__ == '__main__':
    
    #### input, output file ####
    vcf_file = '/data/script/namhee21/test_1004/NA12878-1.bwa.RG.markdup.deldecoy.hc.merged.LATV.chr.vep_test.vcf' 
    vep_output_file = '/data/script/namhee21/test_1004/output.csv'
    
    ApplyMane = ApplyMane()
    ApplyMane.apply_parse_vep(vcf_file, vep_output_file)
    ApplyMane.ghr_omim_annotation(vep_output_file)
    
    print(f'mane 1개 dictionary : {len(ApplyMane.mane_dict_1)}, mane 1개 이상 dictionary : {len(ApplyMane.mane_dict_multi)}')
    print(f'gene : A1BG, mane transcript : {ApplyMane.mane_dict_1.get("A1BG")}')
    print(f'gene : BMP1, mane transcript : {ApplyMane.mane_dict_multi.get("BMP1")}')
