import sys
import os
import subprocess
import argparse
import logging
import concurrent.futures
import glob
import re
from pathlib import Path
import time
import resource
import shutil
import json
import pandas as pd
from datetime import datetime

now = datetime.now()
formatted_time = now.strftime("%y%m%d_%H%M%S")

hg19 = "/data/ref/hg19/hereditary_disease/hs37d5.fa.gz"
hg19_fai = "/data/ref/hg19/hereditary_disease/hs37d5.fa.gz.fai"

dbsnp = "/data/ref/hg19/dbsnp_138.hg19.vcf"
contig_mapping ="/data/ref/hg19/hereditary_disease/contig_mapping.txt"

sample_dict = {}

parser = argparse.ArgumentParser(description='In-house Pipeline')
parser.add_argument('run_id', help='run_id (example: 220723_M70544_0051_000000000-G9C6V)')
parser.add_argument('sample_id', help='sample_id (example: 190966)')
parser.add_argument('-o', dest='output_dir', default=f'/data/analysis/project/hd/', help='output directory (default: %(default)s)')
parser.add_argument('-f1', dest='raw_fastq1', default=f'', help='output directory (default: %(default)s)')
parser.add_argument('-f2', dest='raw_fastq2', default=f'', help='output directory (default: %(default)s)')
parser.add_argument('-t', dest='threads', type=int, default=12, help='threads (default: %(default)s)')
parser.add_argument('-b', dest='bwa', default='quay.io/biocontainers/bwa:0.7.17--hed695b0_7', help='docker image (default: %(default)s)')
parser.add_argument('-s', dest='samtools', default='quay.io/biocontainers/samtools:1.16.1--h6899075_1', help='docker image (default: %(default)s)')
parser.add_argument('-g', dest='gatk', default='broadinstitute/gatk:4.3.0.0', help='docker image (default: %(default)s)')
parser.add_argument('-c', dest='picard', default='broadinstitute/picard:2.27.5', help='docker image (default: %(default)s)')
parser.add_argument('-f', dest='fastqc', default='biocontainers/fastqc:v0.11.9_cv8', help='docker image (default: %(default)s)')
parser.add_argument('-tabix', dest='tabix', default='biocontainers/tabix:v1.9-11-deb_cv1', help='docker image (default: %(default)s)')
parser.add_argument('-bcftools', dest='bcftools', default='biocontainers/bcftools:v1.9-1-deb_cv1', help='docker image (default: %(default)s)')
parser.add_argument('-v', dest='vep', default='ensemblorg/ensembl-vep:release_108.0', help='docker image (default: %(default)s)')
parser.add_argument('-j', dest='jvarkit', default='jvarkit:0f51971fd08b70247f44', help='docker image (default: %(default)s)')
parser.add_argument('--rbed', dest='region_bed', default='/data/ref/hg19/hereditary_disease/TSOE_MANE_exon.bed', help='docker image (default: %(default)s)')
parser.add_argument('--gbed', dest='gene_bed', default='/data/ref/hg19/hereditary_disease/TSOE_MANE_transcript.bed', help='docker image (default: %(default)s)')
parser.add_argument('--interval', dest='interval', default='/data/ref/hg19/hereditary_disease/TSOne_Expanded.change_chr.intervallist', help='docker image (default: %(default)s)')
args = parser.parse_args()

sample_dict["sample_dir"] = f'{args.output_dir}/{args.run_id}/{args.sample_id}'
os.makedirs(sample_dict["sample_dir"], exist_ok=True)

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter('[%(levelname)5s - %(asctime)s] %(message)s',"%Y-%m-%d %H:%M:%S")

stream_handler = logging.StreamHandler(sys.stdout)
stream_handler.setFormatter(formatter)
stream_handler.setLevel(logging.INFO)
logger.addHandler(stream_handler)

file_handler = logging.FileHandler(f'{sample_dict["sample_dir"]}/run_{formatted_time}.log','w')
file_handler.setFormatter(formatter)
file_handler.setLevel(logging.DEBUG)
logger.addHandler(file_handler)

start_time = time.time()

def init_env_check():
    logging.info(f'{"-"*20}{"germline variant call In-house Pipeline":^50}{"-"*20}')

    logging.info(f'{"cmd":<15}: {" ".join(sys.argv)}')
    logging.info(f'{"log file":<15}: {sample_dict["sample_dir"]}/preprocess.log')
    #fastq 파일이 존재하는지 확인
    if args.raw_fastq1 == '' and args.raw_fastq2 == '' :
        sample_dict["raw_fastq1"] = glob.glob(f'/data/raw_data/{args.run_id}/{args.sample_id}*R1*.fastq.gz')[0]
        sample_dict["raw_fastq2"] = glob.glob(f'/data/raw_data/{args.run_id}/{args.sample_id}*R2*.fastq.gz')[0]                
    else :
        sample_dict["raw_fastq1"] = args.raw_fastq1
        sample_dict["raw_fastq2"] = args.raw_fastq2

    if sample_dict["raw_fastq1"] == '' or sample_dict["raw_fastq2"] == '' :
        print('file not exist')
        exit()
    #저장 기간에 따라 폴더를 분리.
    #fastq = 2년, prerocess = 용량이 크므로 1년, result = 영구보관
    sample_dict["sample_fastq_dir"] = f'{args.output_dir}/{args.run_id}/{args.sample_id}/0.fastq'
    sample_dict["sample_preprocess_dir"] = f'{args.output_dir}/{args.run_id}/{args.sample_id}/1.preprocessing'
    sample_dict["sample_result_dir"] = f'{args.output_dir}/{args.run_id}/{args.sample_id}/2.result'
    os.makedirs(sample_dict["sample_fastq_dir"], exist_ok=True)
    os.makedirs(sample_dict["sample_preprocess_dir"], exist_ok=True)
    os.makedirs(sample_dict["sample_result_dir"], exist_ok=True)
    shutil.copy(sample_dict["raw_fastq1"], f'{sample_dict["sample_dir"]}/0.fastq/')
    shutil.copy(sample_dict["raw_fastq2"], f'{sample_dict["sample_dir"]}/0.fastq/')

def finish_summary():
    #spend time as a min
    end_time = time.time()
    spend_time = (end_time - start_time)/60    

    logging.info(f'{"-"*20}{"Summary":^50}{"-"*20}')
    logging.info(f'{"sample_id":<15}: {args.sample_id}')
    logging.info(f'{"run_id":<15}: {args.run_id}')
    logging.info(f'{"fastq1":<15}: {sample_dict["raw_fastq1"]}')
    logging.info(f'{"fastq2":<15}: {sample_dict["raw_fastq2"]}')
    logging.info(f'{"sample_fastq_dir":<15}: {sample_dict["sample_fastq_dir"]}')
    logging.info(f'{"sample_preprocess_dir":<15}: {sample_dict["sample_preprocess_dir"]}')
    logging.info(f'{"sample_result_dir":<15}: {sample_dict["sample_result_dir"]}')
    logging.info(f'{"spend time":<15}: {spend_time:.2f} min')
    logging.info(f'{"pipeline ver":<15}: 1.0.0')
    logging.info(f'{"-"*20}{"Pipeline Finished":^50}{"-"*20}')

def log_subprocess_output(step, cmd):
    logging.info(f'{"processing":<15}: {step:<30}')
    logging.debug(f'$ {cmd}')
    #subprocess 이후 진행 상황을 file_handler에 저장
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    while True:
        line = process.stdout.readline()
        if not line :
            break
        logging.debug('> %s', line.decode('utf-8').strip())
    
    process.wait()

def run_bwa():
    #bwa 맵핑.
    sample_dict["bwa_sam"] = f'{sample_dict["sample_preprocess_dir"]}/{args.sample_id}.bwa.sam'
    cmd = f'docker run -u {os.getuid()}:{os.getgid()} --rm -v {sample_dict["sample_dir"]}:{sample_dict["sample_dir"]} -v /data/ref/hg19:/data/ref/hg19 -v /data/raw_data/{args.run_id}:/data/raw_data/{args.run_id} {args.bwa} bwa mem \
        -M \
        -t {args.threads}\
        {hg19} \
        {sample_dict["raw_fastq1"]}  \
        {sample_dict["raw_fastq2"]} \
        -o {sample_dict["bwa_sam"]} '    
    log_subprocess_output(f'bwa mem ({args.sample_id})', cmd)

def run_samtools():
    #bwa 결과를 sam에서 bam으로 변환
    sample_dict["bwa_bam"] = f'{sample_dict["sample_preprocess_dir"]}/{args.sample_id}.bwa.bam'
    cmd = f'docker run -u {os.getuid()}:{os.getgid()} --rm -v {sample_dict["sample_dir"]}:{sample_dict["sample_dir"]} {args.samtools} samtools view \
        -@ {args.threads} \
        -h \
        -bS \
        {sample_dict["bwa_sam"]} \
        > {sample_dict["bwa_bam"]} '    
    log_subprocess_output(f'samtools view ({args.sample_id})', cmd)

def run_gatk_preprocess():
    #bwa 결과에 read group 추가
    sample_dict["bwa_bam"] = f'{sample_dict["sample_preprocess_dir"]}/{args.sample_id}.bwa.bam'
    sample_dict["gatk_RG_bam"] = f'{sample_dict["sample_preprocess_dir"]}/{args.sample_id}.bwa.RG.bam'
    cmd = f'docker run -u {os.getuid()}:{os.getgid()} --rm -v {sample_dict["sample_dir"]}:{sample_dict["sample_dir"]} {args.gatk} gatk AddOrReplaceReadGroups \
        -I {sample_dict["bwa_bam"]} \
        -O {sample_dict["gatk_RG_bam"]} \
        --SORT_ORDER coordinate \
        --CREATE_INDEX true \
        --RGLB TSOE \
        --RGPL ILLUMINA \
        --RGPU {args.run_id} \
        --RGSM {args.sample_id}'
    log_subprocess_output(f'gatk AddOrReplaceReadGroups ({args.sample_id})', cmd)
    sample_dict["gatk_RG_tx_bamstats05"] = f'{sample_dict["sample_preprocess_dir"]}/{args.sample_id}.bwa.RG.bam.tx_bamstats05.txt'
    cmd = f'docker run -u {os.getuid()}:{os.getgid()} --rm -d -v /data/ref/hg19:/data/ref/hg19 -v {sample_dict["sample_dir"]}:{sample_dict["sample_dir"]} {args.jvarkit} java -jar /opt/jvarkit/dist/jvarkit.jar bamstats05 \
        -o {sample_dict["gatk_RG_tx_bamstats05"]} \
        -B {args.gene_bed} \
        -m 10 \
        {sample_dict["gatk_RG_bam"]} --merge'
    log_subprocess_output(f'jvarkit bamstats05 tx ({args.sample_id})', cmd) 
    # bamstats05로 엑손 별 depth, coverage 등을 계산
    sample_dict["gatk_RG_exon_bamstats05"] = f'{sample_dict["sample_preprocess_dir"]}/{args.sample_id}.bwa.RG.bam.exon_bamstats05.txt'
    cmd = f'docker run -u {os.getuid()}:{os.getgid()} --rm -v /data/ref/hg19:/data/ref/hg19 -v {sample_dict["sample_dir"]}:{sample_dict["sample_dir"]} {args.jvarkit} java -jar /opt/jvarkit/dist/jvarkit.jar bamstats05 \
        -o {sample_dict["gatk_RG_exon_bamstats05"]} \
        -B {args.region_bed} \
        -m 10 \
        {sample_dict["gatk_RG_bam"]}'
    log_subprocess_output(f'jvarkit bamstats05 exon ({args.sample_id})', cmd)    

    #bwa 결과에 read group 추가한 bam 파일에 mark duplicates 추가
    sample_dict["gatk_RG_markdup_bam"] = f'{sample_dict["sample_preprocess_dir"]}/{args.sample_id}.bwa.RG.markdup.bam'
    sample_dict["gatk_RG_markdup_metrics"] = f'{sample_dict["sample_preprocess_dir"]}/{args.sample_id}.bwa.RG.markdup.txt'
    cmd = f'docker run -u {os.getuid()}:{os.getgid()} --rm -v {sample_dict["sample_dir"]}:{sample_dict["sample_dir"]} {args.gatk} gatk --java-options "-Xmx8g" MarkDuplicates \
        --INPUT {sample_dict["gatk_RG_bam"]} \
        --OUTPUT {sample_dict["gatk_RG_markdup_bam"]} \
        --CREATE_INDEX true \
        --METRICS_FILE {sample_dict["gatk_RG_markdup_metrics"]}'
    log_subprocess_output(f'gatk MarkDuplicates ({args.sample_id})', cmd)
    # mark duplicates에 multiple metrics 계산 결과를 추가
    sample_dict["gatk_RG_markdup_multiple_metrics"] = f'{sample_dict["sample_preprocess_dir"]}/{args.sample_id}.bwa.RG.markdup.bam.multiplemetrics'
    cmd = f'docker run -u {os.getuid()}:{os.getgid()} --rm -d -v /data/ref/hg19:/data/ref/hg19 -v {sample_dict["sample_dir"]}:{sample_dict["sample_dir"]} {args.gatk} gatk --java-options "-Xmx8g" CollectMultipleMetrics  \
        -I {sample_dict["gatk_RG_markdup_bam"]} \
        -O {sample_dict["gatk_RG_markdup_multiple_metrics"]} \
        --PROGRAM CollectQualityYieldMetrics'
    if args.interval != None:
        cmd = cmd + f' --INTERVALS {args.interval}'
    log_subprocess_output(f'gatk CollectMultipleMetrics ({args.sample_id})', cmd)
    # mark duplicates에 hs metrics 계산 결과를 추가
    sample_dict["gatk_RG_markdup_hs_metrics"] = f'{sample_dict["sample_preprocess_dir"]}/{args.sample_id}.bwa.RG.markdup.bam.hsmetrics'
    cmd = f'docker run -u {os.getuid()}:{os.getgid()} --rm -d -v /data/ref/hg19:/data/ref/hg19 -v {sample_dict["sample_dir"]}:{sample_dict["sample_dir"]} {args.gatk} gatk --java-options "-Xmx8g" CollectHsMetrics  \
        -I {sample_dict["gatk_RG_markdup_bam"]} \
        -O {sample_dict["gatk_RG_markdup_hs_metrics"]}'
    if args.interval != None:
        cmd = cmd + f' -BI {args.interval} -TI {args.interval}'
         
    log_subprocess_output(f'gatk CollectHsMetrics ({args.sample_id})', cmd)
    # bamstats05로 유전자 별 depth, coverage 등을 계산
    sample_dict["gatk_RG_markdup_tx_bamstats05"] = f'{sample_dict["sample_preprocess_dir"]}/{args.sample_id}.bwa.RG.markdup.bam.tx_bamstats05.txt'
    cmd = f'docker run -u {os.getuid()}:{os.getgid()} --rm -d -v /data/ref/hg19:/data/ref/hg19 -v {sample_dict["sample_dir"]}:{sample_dict["sample_dir"]} {args.jvarkit} java -jar /opt/jvarkit/dist/jvarkit.jar bamstats05 \
        -o {sample_dict["gatk_RG_markdup_tx_bamstats05"]} \
        -B {args.gene_bed} \
        -m 10 \
        {sample_dict["gatk_RG_markdup_bam"]} --merge'
    log_subprocess_output(f'jvarkit bamstats05 tx ({args.sample_id})', cmd) 
    # bamstats05로 엑손 별 depth, coverage 등을 계산
    sample_dict["gatk_RG_markdup_exon_bamstats05"] = f'{sample_dict["sample_preprocess_dir"]}/{args.sample_id}.bwa.RG.markdup.bam.exon_bamstats05.txt'
    cmd = f'docker run -u {os.getuid()}:{os.getgid()} --rm -v /data/ref/hg19:/data/ref/hg19 -v {sample_dict["sample_dir"]}:{sample_dict["sample_dir"]} {args.jvarkit} java -jar /opt/jvarkit/dist/jvarkit.jar bamstats05 \
        -o {sample_dict["gatk_RG_markdup_exon_bamstats05"]} \
        -B {args.region_bed} \
        -m 10 \
        {sample_dict["gatk_RG_markdup_bam"]}'
    log_subprocess_output(f'jvarkit bamstats05 exon ({args.sample_id})', cmd)    

    #remove decoy read. The decoy mapped reads can be seen in BAM file viewer as mapped on the human chromosome. To prevent this problem, I will erase this.
    sample_dict["gatk_RG_markdup_deldecoy_bam"] = f'{sample_dict["sample_preprocess_dir"]}/{args.sample_id}.bwa.RG.markdup.deldecoy.bam'
    cmd = f'docker run -u {os.getuid()}:{os.getgid()} --rm -v {sample_dict["sample_dir"]}:{sample_dict["sample_dir"]} {args.samtools} samtools view \
        -@ {args.threads} \
        -bq3 \
        {sample_dict["gatk_RG_markdup_bam"]} \
        1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT > \
        {sample_dict["gatk_RG_markdup_deldecoy_bam"]}'
    log_subprocess_output(f'samtools view ({args.sample_id})', cmd)

    sample_dict["gatk_RG_markdup_deldecoy_bai"] = f'{sample_dict["sample_preprocess_dir"]}/{args.sample_id}.bwa.RG.markdup.deldecoy.bai'
    cmd = f'docker run -u {os.getuid()}:{os.getgid()} --rm -v {sample_dict["sample_dir"]}:{sample_dict["sample_dir"]} {args.samtools} samtools index \
        -@ {args.threads} \
        {sample_dict["gatk_RG_markdup_deldecoy_bam"]} \
        {sample_dict["gatk_RG_markdup_deldecoy_bai"]}'
    log_subprocess_output(f'samtools index ({args.sample_id})', cmd)
   
def run_gatk_variantcall():
    #gatk haplotypecaller로 vcf 파일 생성
    sample_dict["gatk_RG_markdup_deldecoy_bam"] = f'{sample_dict["sample_preprocess_dir"]}/{args.sample_id}.bwa.RG.markdup.deldecoy.bam'
    sample_dict["HaplotypeCaller_vcf"] = f'{sample_dict["sample_preprocess_dir"]}/{args.sample_id}.bwa.RG.markdup.deldecoy.hc.vcf.gz'
    cmd = f'docker run -u {os.getuid()}:{os.getgid()} --rm -v {sample_dict["sample_dir"]}:{sample_dict["sample_dir"]} -v /data/ref/hg19:/data/ref/hg19 {args.gatk} gatk --java-options "-Xmx8g" HaplotypeCaller \
        -R {hg19} \
        -I {sample_dict["gatk_RG_markdup_deldecoy_bam"]} \
        -O {sample_dict["HaplotypeCaller_vcf"]} \
        --native-pair-hmm-threads {args.threads} '
    log_subprocess_output(f'gatk HaplotypeCaller ({args.sample_id})', cmd) 
    #vcf에서 snp과 indel을 분리
    sample_dict["SelectVariants_snp_vcf"] = f'{sample_dict["sample_preprocess_dir"]}/{args.sample_id}.bwa.RG.markdup.deldecoy.hc.snp.vcf.gz'
    cmd = f'docker run -u {os.getuid()}:{os.getgid()} --rm -v {sample_dict["sample_dir"]}:{sample_dict["sample_dir"]} {args.gatk} gatk SelectVariants \
        -V {sample_dict["HaplotypeCaller_vcf"]} \
        -select-type SNP \
        -O {sample_dict["SelectVariants_snp_vcf"]}'
    log_subprocess_output(f'gatk snp SelectVariants ({args.sample_id})', cmd)
    #snp과 indel을 분리하여 각각에 상황에 맞는 hard filter를 적용
    sample_dict["VariantFiltration_snp_vcf"] = f'{sample_dict["sample_preprocess_dir"]}/{args.sample_id}.bwa.RG.markdup.deldecoy.hc.snp.filtered.vcf.gz'
    cmd = f'docker run -u {os.getuid()}:{os.getgid()} --rm -v {sample_dict["sample_dir"]}:{sample_dict["sample_dir"]} {args.gatk} gatk VariantFiltration \
        -V {sample_dict["SelectVariants_snp_vcf"]} \
        -filter "QD < 2.0" --filter-name "QD2" \
        -filter "QUAL < 30.0" --filter-name "QUAL30" \
        -filter "SOR > 12.0" --filter-name "SOR12" \
        -filter "FS > 60.0" --filter-name "FS60" \
        -filter "MQ < 40.0" --filter-name "MQ40" \
        -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
        -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
        -filter "DP<10" --filter-name "DP" \
        -O {sample_dict["VariantFiltration_snp_vcf"]}'
    log_subprocess_output(f'gatk snp VariantFiltration ({args.sample_id})', cmd)

    sample_dict["SelectVariants_indel_vcf"] = f'{sample_dict["sample_preprocess_dir"]}/{args.sample_id}.bwa.RG.markdup.deldecoy.hc.indel.vcf.gz'
    cmd = f'docker run -u {os.getuid()}:{os.getgid()} --rm -v {sample_dict["sample_dir"]}:{sample_dict["sample_dir"]} {args.gatk} gatk SelectVariants \
        -V {sample_dict["HaplotypeCaller_vcf"]} \
        -select-type INDEL \
        -O {sample_dict["SelectVariants_indel_vcf"]}'
    log_subprocess_output(f'gatk indel SelectVariants ({args.sample_id})', cmd)

    sample_dict["VariantFiltration_indel_vcf"] = f'{sample_dict["sample_preprocess_dir"]}/{args.sample_id}.bwa.RG.markdup.deldecoy.hc.indel.filtered.vcf.gz'
    cmd = f'docker run -u {os.getuid()}:{os.getgid()} --rm -v {sample_dict["sample_dir"]}:{sample_dict["sample_dir"]} {args.gatk} gatk VariantFiltration \
        -V {sample_dict["SelectVariants_indel_vcf"]} \
        -filter "QD<7.0" --filter-name "QD7" \
        -filter "QUAL<30.0" --filter-name "QUAL30" \
        -filter "FS>200.0" --filter-name "FS200" \
        -filter "DP<10" --filter-name "DP"  \
        -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
        -O {sample_dict["VariantFiltration_indel_vcf"]}'
    log_subprocess_output(f'gatk indel VariantFiltration ({args.sample_id})', cmd)
    #각각 필터를 적용한 snp과 indel 파일을 다시 합침
    sample_dict["merged_vcf"] = f'{sample_dict["sample_preprocess_dir"]}/{args.sample_id}.bwa.RG.markdup.deldecoy.hc.merged.vcf.gz'
    cmd = f'docker run -u {os.getuid()}:{os.getgid()} --rm -v {sample_dict["sample_dir"]}:{sample_dict["sample_dir"]} {args.gatk} gatk MergeVcfs \
        -I {sample_dict["VariantFiltration_snp_vcf"]}  \
        -I {sample_dict["VariantFiltration_indel_vcf"]} \
        -O {sample_dict["merged_vcf"]}'
    log_subprocess_output(f'gatk MergeVcfs ({args.sample_id})', cmd)
    #merged된 vcf 파일에 left align and trim variants 적용
    sample_dict["LeftAlignAndTrimVariants_vcf"] = f'{sample_dict["sample_preprocess_dir"]}/{args.sample_id}.bwa.RG.markdup.deldecoy.hc.merged.LATV.vcf.gz'
    cmd = f'docker run -u {os.getuid()}:{os.getgid()} --rm -v /data/ref/hg19:/data/ref/hg19 -v {sample_dict["sample_dir"]}:{sample_dict["sample_dir"]} {args.gatk} gatk LeftAlignAndTrimVariants  \
        -R {hg19} \
        -V {sample_dict["merged_vcf"]} \
        -O {sample_dict["LeftAlignAndTrimVariants_vcf"]} \
        --split-multi-allelics'
    log_subprocess_output(f'gatk LeftAlignAndTrimVariants ({args.sample_id})', cmd)
    #left align and trim variants 적용된 vcf 파일에 contig mapping 적용
    sample_dict["annotate_vcf"] = f'{sample_dict["sample_preprocess_dir"]}/{args.sample_id}.bwa.RG.markdup.deldecoy.hc.merged.LATV.chr.vcf.gz'
    cmd = f'docker run -u {os.getuid()}:{os.getgid()} --rm -v /data/ref/hg19:/data/ref/hg19 -v {sample_dict["sample_dir"]}:{sample_dict["sample_dir"]} {args.bcftools} bcftools annotate  \
        --rename-chrs {contig_mapping} \
        {sample_dict["LeftAlignAndTrimVariants_vcf"]} \
        -Oz \
        -o {sample_dict["annotate_vcf"]}'
    log_subprocess_output(f'bcftools annotate ({args.sample_id})', cmd)
    #contig mapping 적용된 vcf 파일에 tabix 적용
    cmd = f'docker run -u {os.getuid()}:{os.getgid()} --rm -v {sample_dict["sample_dir"]}:{sample_dict["sample_dir"]} {args.tabix} tabix \
        -p \
        vcf {sample_dict["annotate_vcf"]}'
    log_subprocess_output(f'tabix ({args.sample_id})', cmd)

def run_vep():
    #vep로 vcf 파일에 annotation 
    sample_dict["annotate_vcf"] = f'{sample_dict["sample_preprocess_dir"]}/{args.sample_id}.bwa.RG.markdup.deldecoy.hc.merged.LATV.chr.vcf.gz'
    sample_dict["vep_vcf"] = f'{sample_dict["sample_preprocess_dir"]}/{args.sample_id}.bwa.RG.markdup.deldecoy.hc.merged.LATV.chr.vep.vcf'
    cmd = f'docker run -u {os.getuid()}:{os.getgid()} --rm -v /data/ref/hg19:/data/ref/hg19 -v {sample_dict["sample_dir"]}:{sample_dict["sample_dir"]} {args.vep} vep \
        --fork {args.threads}  \
        --fasta /data/ref/hg19/hg19.fa \
        --force_overwrite \
        --everything \
        --offline \
        --refseq \
        --vcf \
        --minimal \
        --cache \
        --dir_cache /data/ref/hg19/hereditary_disease/vep \
        -i {sample_dict["annotate_vcf"]} \
        -o {sample_dict["vep_vcf"]} \
        --dir_plugins /data/ref/hg19/hereditary_disease/Plugins \
        --custom /data/ref/hg19/hereditary_disease/vep/clinvar/clinvar_20230826_xml_update.vcf.gz,ClinVar,vcf,exact,0,ALLELEID,CLNREVSTAT,CLNSIG,variation_id,variation_type,variation_name,number_of_submissions,number_of_submitters,datelastupdated,gene,allele_id,gene_list,chromosome,strand,cytogenetic_location,grch38_hgvs_g,grch37_hgvs_g,hgvs_c,hgvs_p,dbsnp_id,max_af,review_status,interpretation,trait_list,clinical_assertion_localkey,clingen,omim_id,dbvar_id,uniprot_id \
        --plugin CADD,/data/ref/hg19/hereditary_disease/vep/CADD/whole_genome_SNVs.tsv.gz,/data/ref/hg19/hereditary_disease/vep/CADD/InDels.tsv.gz \
        --plugin dbNSFP,/data/ref/hg19/hereditary_disease/vep/dbNSFP4/dbNSFP4.3c_grch37.gz,MutationTaster_score,MutationTaster_pred,MutationAssessor_score,MutationAssessor_pred,FATHMM_score,FATHMM_pred,PROVEAN_score,PROVEAN_pred,MetaSVM_score,MetaSVM_pred,GERP++_RS,phyloP100way_vertebrate,phyloP30way_mammalian,phastCons100way_vertebrate,phastCons30way_mammalian,SiPhy_29way_pi,SiPhy_29way_logOdds \
        --plugin SpliceAI,snv=/data/ref/hg19/hereditary_disease/vep/spliceai/spliceai_scores.raw.snv.hg19.vcf.gz,indel=/data/ref/hg19/hereditary_disease/vep/spliceai/spliceai_scores.raw.indel.hg19.vcf.gz \
        --plugin pLI,/data/ref/hg19/hereditary_disease/vep/pLI/pLI_values.txt'
    log_subprocess_output(f'vep ({args.sample_id})', cmd)        

    omim_df = pd.read_csv(f'/data/script/lims/Hereditary_Disease/omim_gene.csv')
    omim_df.fillna('', inplace=True)
    omim_dict = {}
    for index, row in omim_df.iterrows():
        if not row['Gene'] in omim_dict :
            omim_dict[row['Gene']] = []    
        omim_dict[row['Gene']].append('('+row['Phenotype']+'_'+row['Inheritance']+')')

    for gene_id in omim_dict:
        omim_dict[gene_id] = '/'.join(omim_dict[gene_id])
    
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
                csq_header = csq_header.split('"')[1].split(':')[-1].strip().split('|')+['CSQ_OTH']
                csq_df = pd.DataFrame(columns=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'GENOTPYPE', 'DEPTH', 'VAF', 'INFO'] +csq_header+ ['OMIM', 'FORMAT', 'SAMPLE'])
                vep_parsing_output_file.write('\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'GENOTPYPE', 'DEPTH', 'VAF', 'INFO']+csq_header+['OMIM', 'FORMAT', 'SAMPLE'])+'\n')
                max_af_index = csq_header.index('MAX_AF')
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
                        new_csq_list = []
                        for csq in csq_list :
                            csq_sp = csq.split('|')
                            gene_id = csq_sp[3]
                            transcript_id = csq_sp[6]
                            if csq_sp[max_af_index] == '' :
                                csq_sp[max_af_index] = '0'
                            # exon에서 /가 들어가면 excel에서 오류가 발생하므로 /를 _로 대체
                            csq_sp[8] = csq_sp[8].replace('/','_')
                            csq_sp[9] = csq_sp[9].replace('/','_')
                            new_csq = '|'.join(csq_sp)
                            if gene_id in summary_dict and transcript_id in summary_dict[gene_id] :                                    
                                new_csq_list.insert(0, new_csq)
                            else :
                                new_csq_list.append(new_csq)
                    else :
                        new_info.append(i)
                new_csq = new_csq_list[0].split('|')
                if len(new_csq_list) > 2 :
                    csq_oth = [','.join(new_csq_list[1:])]
                elif len(new_csq_list) == 2 :
                    csq_oth = [new_csq_list[1]]
                else :
                    csq_oth = ['']
                gene_id = new_csq[3]
                if gene_id in omim_dict :
                    omim_info = [omim_dict[gene_id]]
                else :
                    omim_info = ['']
                # IMPACT가 MODIFIER라면 제외
                if 'MODIFIER' in new_csq :
                    continue
                new_line = '\t'.join([chrom, pos, id, ref, alt, qual, filter, genotype, depth, vaf, ';'.join(new_info)] + new_csq + csq_oth + omim_info + [format, sample])
                vep_parsing_output_file.write(new_line+'\n')
    vep_parsing_output_file.close()

def run_qc():
    qc_dict = {}
    #multiqc로 qc 결과를 한 파일로 합침
    os.system(f'multiqc --force -o {sample_dict["sample_preprocess_dir"]} {sample_dict["sample_preprocess_dir"]}')
    multiqc_data_json = f'{sample_dict["sample_preprocess_dir"]}/multiqc_data/multiqc_data.json'
    sample_dict['sample_qc'] = f'{sample_dict["sample_preprocess_dir"]}/{args.sample_id}.qc.txt'
    with open(multiqc_data_json) as f:
        data = json.load(f)
    #qc 결과를 sample 별로 분리하여 저장
    column_list = ['TOTAL_READS', 'TOTAL_BASES', 'Q20_BASES', 'Q30_BASES']
    for sample_id in data['report_general_stats_data'][0] :
        for key in column_list:        
            qc_dict[key] = data['report_general_stats_data'][0][sample_id][key]
        qc_dict['Q20_RATE'] = qc_dict['Q20_BASES'] / qc_dict['TOTAL_BASES']
        qc_dict['Q30_RATE'] = qc_dict['Q30_BASES'] / qc_dict['TOTAL_BASES']
    #qc 결과를 파일로 출력
    column_list = ['TOTAL_READS', 'PCT_PF_READS_ALIGNED', 'PF_ALIGNED_BASES', 'ON_BAIT_BASES', 'NEAR_BAIT_BASES', 'OFF_BAIT_BASES', 'PCT_SELECTED_BASES', 'ON_BAIT_VS_SELECTED', 'PCT_USABLE_BASES_ON_BAIT', 'MEAN_BAIT_COVERAGE', 'FOLD_80_BASE_PENALTY', 'PCT_TARGET_BASES_10X', 'PCT_TARGET_BASES_20X', 'PCT_TARGET_BASES_30X', 'PCT_TARGET_BASES_50X', 'PCT_TARGET_BASES_100X', 'summed_mean', 'PERCENT_DUPLICATION']
    for sample_id in data['report_general_stats_data'][1] :
        for key in column_list:        
            qc_dict[key] = data['report_general_stats_data'][1][sample_id][key]
        qc_df = pd.DataFrame(qc_dict, index=[sample_id])
        #qc 결과를 파일로 출력
        qc_df.to_csv(f'{sample_dict["sample_qc"]}', sep='\t')

def run_make_excel():
    # 파일 경로를 저장할 딕셔너리 생성
    sample_dict['sample_qc'] = f'{sample_dict["sample_preprocess_dir"]}/{args.sample_id}.qc.txt'
    sample_dict["vep_vcf"] = f'{sample_dict["sample_preprocess_dir"]}/{args.sample_id}.bwa.RG.markdup.deldecoy.hc.merged.LATV.chr.vep.vcf'
    sample_dict["gatk_RG_markdup_exon_bamstats05"] = f'{sample_dict["sample_preprocess_dir"]}/{args.sample_id}.bwa.RG.markdup.bam.exon_bamstats05.txt'
    sample_dict["gatk_RG_markdup_tx_bamstats05"] = f'{sample_dict["sample_preprocess_dir"]}/{args.sample_id}.bwa.RG.markdup.bam.tx_bamstats05.txt'
    # 합칠 파일 목록
    file_list = [sample_dict['sample_qc'], sample_dict['gatk_RG_markdup_tx_bamstats05'], sample_dict['gatk_RG_markdup_exon_bamstats05'], sample_dict["vep_vcf"]+'.txt']
    # 엑셀 파일 생성 (또는 기존 엑셀 파일 열기)
    sample_dict["final_excel"] = f'{sample_dict["sample_preprocess_dir"]}/{args.sample_id}.xlsx'
    excel_file = pd.ExcelWriter(f'{sample_dict["final_excel"]}', engine='xlsxwriter')

    # 파일을 각각의 시트로 추가
    for file_name in file_list:
        # 파일을 DataFrame으로 읽어오기
        df = pd.read_csv(f'{file_name}', sep='\t')
        
        # 파일 이름에서 확장자 제거하여 시트 이름 설정
        sheet_name = file_name.split('.')[-2]

        # DataFrame을 시트로 추가
        df.to_excel(excel_file, sheet_name=sheet_name, index=False)
    lims_df = pd.DataFrame(columns=['Class', 'Gene', 'HGVSc', 'Zygosity', 'Disease', 'Inheritance', 'Interpretation'])
    lims_df.to_excel(excel_file, sheet_name='LIMS', index=False)

    # 엑셀 파일 저장
    excel_file.close()
    os.system(f'cp {sample_dict["sample_preprocess_dir"]}/{args.sample_id}.xlsx {sample_dict["sample_result_dir"]}/{args.sample_id}.xlsx')

def main():
    init_env_check()
    run_bwa()
    run_samtools()
    run_gatk_preprocess()
    run_gatk_variantcall()
    run_vep()
    run_qc()
    run_make_excel()
    finish_summary()

if __name__ == "__main__":
    main()

