import pandas as pd
import re
import xml.etree.ElementTree as ET
import tqdm
import openpyxl
from openpyxl.utils.dataframe import dataframe_to_rows
from openpyxl import load_workbook
from openpyxl.utils import get_column_letter
import sys, os
from ftplib import FTP

variation_dict = {}
# xml parsing
def xml_parsing(ftp):
    target_directory = '/pub/clinvar/xml/clinvar_variation/'
    ftp.cwd(target_directory)

    # 파일 목록 출력
    last_updated_date_file_name = sorted([i for i in ftp.nlst() if i.endswith('.xml.gz')])[-1]
    filedate = last_updated_date_file_name.split('.')[0].split('_')[-1]
    # wget xml file    
    if os.path.isfile(f'/data/ref/hg19/hereditary_disease/vep/clinvar/ClinVarVariationRelease_{filedate}.xml') :
        print(f'ClinVarVariationRelease_{filedate}.xml is already exist')
    else :
        print(f'ClinVarVariationRelease_{filedate}.xml is not exist')
        os.system(f'wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/clinvar_variation/ClinVarVariationRelease_{filedate}.xml.gz -o /data/ref/hg19/hereditary_disease/vep/clinvar/ClinVarVariationRelease_{filedate}.xml.gz')
        os.system(f'gunzip -dc /data/ref/hg19/hereditary_disease/vep/clinvar/ClinVarVariationRelease_{filedate}.xml.gz > /data/ref/hg19/hereditary_disease/vep/clinvar/ClinVarVariationRelease_{filedate}.xml')
    # xml 파싱
    xml_file_name = f'/data/ref/hg19/hereditary_disease/vep/clinvar/ClinVarVariationRelease_{filedate}.xml'
    xml_file = open(xml_file_name, 'r')
    XRefList = []
    variation_id = ''
    # class 선언
    class VariantClass :
        def __init__(self, root):            
            variation_id = (root.attrib['VariationID'])
            self.variation_id = root.attrib['VariationID']
            self.variation_type = root.attrib['VariationType']
            self.variation_name = root.attrib['VariationName']
            self.number_of_submissions = root.attrib['NumberOfSubmissions']
            self.number_of_submitters = root.attrib['NumberOfSubmitters']
            self.datelastupdated = root.attrib['DateLastUpdated']
            self.gene = ''
            
            if root.attrib['RecordType'] == 'included' :
                record = root.find('IncludedRecord')
            elif root.attrib['RecordType'] == 'interpreted' :
                record = root.find('InterpretedRecord')
            
            if record.find('SimpleAllele') :
                vartype = record.find('SimpleAllele')
            elif record.find('Haplotype') :
                vartype = record.find('Haplotype').find('SimpleAllele')
            elif record.find('Genotype') :
                if record.find('Genotype').find('Haplotype') :
                    vartype = record.find('Genotype').find('Haplotype').find('SimpleAllele')
                else :
                    vartype = record.find('Genotype').find('SimpleAllele')

            self.allele_id = vartype.attrib['AlleleID']

            if vartype.find('GeneList'):
                self.gene_list = []
                for gene in vartype.find('GeneList').findall('Gene'):
                    self.gene_list.append(gene.attrib['Symbol'])
                    if gene.attrib['Symbol'] in self.variation_name :
                        self.gene = gene.attrib['Symbol']

                self.gene_list = '/'.join(sorted(self.gene_list))

            if vartype.find('Location') != None:
                if vartype.find('Location').find('SequenceLocation') != None:
                    self.chromosome = vartype.find('Location').find('SequenceLocation').attrib['Chr'] if 'Chr' in vartype.find('Location').find('SequenceLocation').attrib else None
                    self.strand = vartype.find('Location').find('SequenceLocation').attrib['Strand'] if 'Strand' in vartype.find('Location').find('SequenceLocation').attrib else None
                if vartype.find('Location').find('CytogeneticLocation') != None:
                    self.cytogenetic_location = vartype.find('Location').find('CytogeneticLocation').text
                
            if vartype.find('HGVSlist') != None :
                for HGVS in vartype.find('HGVSlist').findall('HGVS'):
                    if HGVS.attrib['Type'] == "genomic, top-level" :
                        if 'Assembly' in HGVS.attrib :
                            if HGVS.attrib['Assembly'] == 'GRCh37' and HGVS.find('NucleotideExpression') != None:
                                self.grch37_hgvs_g = HGVS.find('NucleotideExpression').find('Expression').text 
                            elif HGVS.attrib['Assembly'] == 'GRCh38' and HGVS.find('NucleotideExpression') != None: 
                                self.grch38_hgvs_g = HGVS.find('NucleotideExpression').find('Expression').text

                    elif HGVS.attrib['Type'] == "coding" :
                        if HGVS.find('NucleotideExpression') != None and 'MANESelect' in HGVS.find('NucleotideExpression').attrib :
                            self.hgvs_c = HGVS.find('NucleotideExpression').find('Expression').text
                            if HGVS.find('ProteinExpression') != None :
                                self.hgvs_p = HGVS.find('ProteinExpression').find('Expression').text

            if vartype.find('XRefList') != None :
                for XRef in vartype.find('XRefList').findall('XRef'): 
                    if XRef.attrib['DB'] == 'dbSNP' :
                        self.dbsnp_id = XRef.attrib['ID']
                    elif XRef.attrib['DB'] == 'dbVar' :
                        self.dbvar_id = XRef.attrib['ID']
                    elif XRef.attrib['DB'] == 'ClinGen' :
                        self.clingen = XRef.attrib['ID']
                    elif XRef.attrib['DB'] == 'OMIM' :
                        self.omim_id = XRef.attrib['ID']           

            if vartype.find('AlleleFrequencyList') != None:
                self.af_list = []
                self.max_af = 0
                for AlleleFrequency in vartype.find('AlleleFrequencyList').findall('AlleleFrequency'):
                    if self.max_af < float(AlleleFrequency.attrib['Value']) :
                        self.max_af = float(AlleleFrequency.attrib['Value'])
                    self.af_list.append(AlleleFrequency.attrib['Source'] + ':' + AlleleFrequency.attrib['Value']) 
                self.af_list = '/'.join(self.af_list)
                        
            self.review_status = record.find('ReviewStatus').text        
            self.interpretation = []
            self.trait_list = []
            for Interpretation in record.find('Interpretations').findall('Interpretation'):
                if not Interpretation.find('Description').text in self.interpretation :
                    self.interpretation.append(Interpretation.find('Description').text)            
                if Interpretation.find('ConditionList') != None :
                    for traitset in Interpretation.find('ConditionList').findall('TraitSet') :
                        for trait in traitset.findall('Trait') :                            
                            for traitname in trait.findall('Name') :
                                self.trait_list.append(traitname.find('ElementValue').text)                            

            self.interpretation = '/'.join(self.interpretation)
            self.trait_list = '/'.join(self.trait_list)            

            self.clinical_assertion_localkey = {}                 
            if record.find('ClinicalAssertionList') != None :
                for ClinicalAssertion in record.find('ClinicalAssertionList').findall('ClinicalAssertion'):
                    if ClinicalAssertion.find('ClinVarSubmissionID').attrib['localKey'] != None :
                        localkey = ClinicalAssertion.find('ClinVarSubmissionID').attrib['localKey']
                        # localkey value is interpretation, sumbitername, dateupdated, inheritance
                        self.clinical_assertion_localkey[localkey] = ['Unknown','Unknown','Unknown','Unknown']
                        if ClinicalAssertion.find('Interpretation').find('Description') != None :
                            self.clinical_assertion_localkey[localkey][0] = ClinicalAssertion.find('Interpretation').find('Description').text
                        if ClinicalAssertion.find('ClinVarAccession').attrib['SubmitterName'] != None :
                            self.clinical_assertion_localkey[localkey][1] = ClinicalAssertion.find('ClinVarAccession').attrib['SubmitterName']
                        if 'DateUpdated' in ClinicalAssertion.find('ClinVarAccession').attrib and ClinicalAssertion.find('ClinVarAccession').attrib['DateUpdated'] is not None:
                            self.clinical_assertion_localkey[localkey][2] = ClinicalAssertion.find('ClinVarAccession').attrib['DateUpdated']
                        for attributeset in  ClinicalAssertion.findall('AttributeSet') :
                            if attributeset.find('Attribute').attrib['Type'] == 'ModeOfInheritance':
                                self.clinical_assertion_localkey[localkey][3] = attributeset.find('Attribute').text

            self.clinical_assertion_localkey = '/'.join(['('+ key + '_' + '_'.join(value) + ')' for key, value in self.clinical_assertion_localkey.items()])            

    VariationArchive = False
    pbar = tqdm.tqdm(total=2198285)
    count = 0
    while True :
        line = xml_file.readline().strip()
        if not line :
            line = xml_file.readline().strip()
            if not line:
                break
            
        if line.startswith('<VariationArchive') :
            VariationArchive = True
            variant_contents = line.strip()

        elif line.startswith('</VariationArchive') :
            VariationArchive = False
            variant_contents += line.strip()
            root = ET.fromstring(variant_contents)
            vclass = VariantClass(root)
            variation_dict[vclass.allele_id] = vclass
            pbar.update(1)
            if count == 10000 :
                pass
            else :
                count += 1

        elif VariationArchive:
            variant_contents += line.strip()

    xml_file.close()

def update_vcf(ftp):
    target_directory = '/pub/clinvar/vcf_GRCh37/weekly/'
    ftp.cwd(target_directory)

    # 파일 목록 출력
    last_updated_date_file_name = sorted([i for i in ftp.nlst() if i.endswith('.vcf.gz') and not 'papu' in i])[-1]
    filedate = last_updated_date_file_name.split('.')[0].split('_')[-1]

    if os.path.isfile(f'/data/ref/hg19/hereditary_disease/vep/clinvar/clinvar_{filedate}.vcf') :
        print(f'clinvar_{filedate}.vcf is already exist')
    else :
        print(f'clinvar_{filedate}.vcf is not exist')
        os.system(f'wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/weekly/clinvar_{filedate}.vcf.gz -o /data/ref/hg19/hereditary_disease/vep/clinvar/clinvar_{filedate}.vcf.gz')
        os.system(f'gunzip /data/ref/hg19/hereditary_disease/vep/clinvar/clinvar_{filedate}.vcf.gz')

    outfile = open(f'/data/ref/hg19/hereditary_disease/vep/clinvar/clinvar_{filedate}_xml_update.vcf', 'w')
    with open(f'/data/ref/hg19/hereditary_disease/vep/clinvar/clinvar_{filedate}.vcf') as f:   
        for line in f:
            if line.startswith('#') :
                outfile.write(line)
                continue
            else :
                line_sp = line.strip().split('\t')
                info = line_sp[7]
                allele_id = [i.split('=')[-1] for i in info.split(';') if i.startswith('ALLELEID=')][0]
                add_info = []
                if allele_id in variation_dict :
                    for key, value in vars(variation_dict[allele_id]).items():
                        if key in ['af_list'] : continue
                        if type(value) == str : value = value.replace('=','%3D')
                        add_info.append(f'{key}={value}')
                    add_info = [''] + add_info
                    line_sp[7] = line_sp[7] + ';'.join(add_info)
                outfile.write('\t'.join(line_sp)+'\n')
    outfile.close()
    os.system(f'bgzip /data/ref/hg19/hereditary_disease/vep/clinvar/clinvar_{filedate}.vcf')

if __name__ == "__main__":
    if len(sys.argv) != 1:
        print('Usage: python clinvar_parsing.py')
        sys.exit()
    elif len(sys.argv) == 1:
        # FTP 서버에 연결
        ftp = FTP('ftp.ncbi.nlm.nih.gov')  # FTP 서버 주소 입력
        # 로그인
        ftp.login()  # FTP 서버 로그인 정보 입력

        xml_parsing(ftp)
        update_vcf(ftp)

        # FTP 연결 종료
        ftp.quit()