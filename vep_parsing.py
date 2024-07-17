
import os
import sys
import re
import json
import argparse


if len(sys.argv) != 3:
    print("usage: python vep_parsing.py run_id sample_id")
    sys.exit(1)
sys.argv[1] = run_id
sys.argv[2] = sample_id

sample_dict = {}

sample_dict["sample_preprocess_dir"] = f'/data/analysis/project/hd/{run_id}/{sample_id}/1.preprocessing'
sample_dict["vep_vcf"] = f'{sample_dict["sample_preprocess_dir"]}/{sample_id}.bwa.RG.markdup.deldecoy.hc.merged.LATV.chr.vep.vcf'

def parsing_vep():
    #vep 에서 나온 결과 중에 CSQ 기준으로 지정된 transcript ID에 해당하는 정보만 가져옴.
    summary_df = pd.read_csv("/data/script/lims/Hereditary_Disease/Gene_Transcript.csv")
    
    summary_dict = {}
    for index, row in summary_df.iterrows():
        gene_name = row['Gene_Name']
        transcript_name = row['Transcript']
        gene_id = row['HGNC_ID']
        summary_dict[gene_name] = [transcript_name, gene_id]
    #vep 결과를 parsing하여 필요한 정보만 추출
    vep_parsing_output_file = open(f'{sample_dict["vep_vcf"]}.txt', 'w')
    with open(sample_dict["vep_vcf"]) as f:
        for line in f:
            if line.startswith('##VEP'):
                vep_version = line.strip()
            elif line.startswith('##INFO=<ID=CSQ'):
                csq_header = [i for i in line.strip().split(',') if i.startswith('Description=')][0]            
                csq_header = csq_header.split('"')[1].split(':')[-1].strip().split('|')
                csq_df = pd.DataFrame(columns=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'GENOTPYPE', 'DEPTH', 'VAF', 'INFO'] +csq_header+ ['FORMAT', 'SAMPLE'])
                vep_parsing_output_file.write('\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'GENOTPYPE', 'DEPTH', 'VAF', 'INFO']+csq_header+['FORMAT', 'SAMPLE'])+'\n')
            elif line.startswith('#') :
                continue
            else :
                chrom, pos, id, ref, alt, qual, filter, info, format, sample = line.strip().split('\t')
                # genotype, depth, vaf는 sample 정보를 기반으로 따로 분리하여 표시
                genotype = sample.split(':')[0]
                if genotype == '0/1' : genotype = 'HET'
                elif genotype == '1/1' : genotype = 'HOM'
                elif genotype == '0/0' : genotype = 'REF'
                else : genotype = 'OTH'
                allele_depth = sample.split(':')[1].split(',')[1]
                depth = sample.split(':')[2]
                try :
                    vaf = str(round(int(allele_depth)/int(depth)*100, 2))
                except:
                    vaf = '0'
                new_info = []
                for i in info.split(';') :
                    # vep는 '=' 기준으로 split 하기 때문에 '='이 들어간 값은 임시로 다른 문자로 대체하여 split 후 다시 원래 문자로 대체
                    i = i.replace('%3D', '=').replace('%2C', ',').replace('%26', '&').replace('%7C', '|')
                    if i.startswith('CSQ') :
                        csq_list = i[4:].split(',')
                        csq_info = ''              
                        for csq in csq_list :
                            csq_sp = csq.split('|')
                            gene_id = csq_sp[3]
                            transcript_id = csq_sp[6]
                            # exon에서 /가 들어가면 excel에서 오류가 발생하므로 /를 _로 대체
                            csq_sp[8] = csq_sp[8].replace('/','_')
                            csq_sp[9] = csq_sp[9].replace('/','_')
                            if csq_info == '' :
                                csq_info = csq_sp
                            if gene_id in summary_dict :
                                if transcript_id in summary_dict[gene_id] :
                                    csq_info = csq_sp
                    else :
                        new_info.append(i)
                # IMPACT가 MODIFIER라면 제외
                if 'MODIFIER' in csq_info :
                    continue
                new_line = '\t'.join([chrom, pos, id, ref, alt, qual, filter, genotype, depth, vaf, ';'.join(new_info)] + csq_info + [format, sample])
                vep_parsing_output_file.write(new_line+'\n')
    vep_parsing_output_file.close()


def main()
    parsing_vep()

if __name__ == '__main__':
    main()