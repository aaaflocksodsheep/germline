import shutil
import os
import argparse
import hashlib
import pandas as pd
import matplotlib.pyplot as plt
import json
import re

parser = argparse.ArgumentParser(description='QC Report')
parser.add_argument('run_id', help='(example: 221102_M70544_0051_000000000-KSLM_BAM)')
args = parser.parse_args()

docker_dict = {}
docker_dict['fastp'] = 'quay.io/biocontainers/fastp:0.23.2--h79da9fb_0'
docker_dict['multiqc'] = 'ewels/multiqc:v1.14'

raw_fastq_path = f'/data/raw_data/{args.run_id}/'
analysis_dir = f'/data/analysis/project/230412_qc_report_jieun/{args.run_id}/'
fastq_dir = f"{analysis_dir}0.fastq/"
mutiqc_dir = f"{analysis_dir}1.fastp_mutiqc/"
result_dir = f"{analysis_dir}2.result/"
output_md5sum_file_path = f'{result_dir}md5sum.csv'
json_file_path = f"{mutiqc_dir}multiqc_data/multiqc_data.json"
output_raw_file_path = f"{result_dir}raw_data_statistics.csv"

dirs_to_create = [fastq_dir, mutiqc_dir,result_dir]
for dir_path in dirs_to_create:
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

def copy_raw_fastq_files(raw_fastq_path, fastq_dir):
    for file in os.listdir(raw_fastq_path):
        if file.endswith(".fastq.gz"):
            #continue
            shutil.copy(os.path.join(raw_fastq_path, file), fastq_dir)
    fastq_files = [os.path.join(fastq_dir, f) for f in os.listdir(fastq_dir) if f.endswith('.fastq.gz')]
    return fastq_files

def run_fastp(fastq_files):
    for fastq_file in fastq_files:
        output_file = os.path.join(mutiqc_dir, os.path.basename(fastq_file).replace('.fastq.gz', '_fastp.fastq.gz'))
        cmd = f'docker run -u {os.getuid()}:{os.getgid()} --rm \
            -v /data/:/data/ \
            {docker_dict["fastp"]} fastp -i {fastq_file} -o {output_file} \
            --json {output_file.replace(".fastq.gz", ".json")} \
            --html {output_file.replace(".fastq.gz", ".html")}'
        print(cmd)
        os.system(cmd)

def run_multiqc():
    cmd = f'docker run -u {os.getuid()}:{os.getgid()} --rm \
        -v /data/:/data/  \
        {docker_dict["multiqc"]} multiqc {mutiqc_dir} -o {mutiqc_dir}'
    print(cmd)
    os.system(cmd)

def generate_md5sum_file(fastq_dir, output_md5sum_file_path):
    with open(output_md5sum_file_path, 'w') as output_file:
        output_file.write("Filename\tFilesize (bytes)\tMD5sum\n")
    
    df = pd.DataFrame(columns=["Filename", "Filesize (bytes)", "MD5sum"])
    
    for root, dirs, files in os.walk(fastq_dir):
        for file in files:
            if file.endswith(".fastq.gz"):
                file_path = os.path.join(root, file)
                with open(file_path, 'rb') as input_file:
                    md5sum = hashlib.md5(input_file.read()).hexdigest()
                    file_size = os.path.getsize(file_path)
                    file_name=file_path.split("/")[-1]
                    with open(output_md5sum_file_path, 'a') as output_file:
                        output_file.write(f"{file_name}\t{file_size}\t{md5sum}\n")
                    df = df.append({"Filename": file_name, "Filesize (bytes)": file_size, "MD5sum": md5sum}, ignore_index=True)
                    #print(df.head(3))
    return df

def process_multiqc_fastp(json_file_path, output_raw_file_path):
    with open(json_file_path, 'r') as f:
        data = json.load(f)    
    samples = {}
    for sample in data['report_saved_raw_data']['multiqc_fastp']:
        sample_name = os.path.basename(sample).split('_')[0]
        if sample_name not in samples:
            samples[sample_name] = {}
        if "R1" in sample:
            samples[sample_name]['before_filtering_r1'] = data['report_saved_raw_data']['multiqc_fastp'][sample]['summary']['before_filtering']
        elif "R2" in sample:
            samples[sample_name]['before_filtering_r2'] = data['report_saved_raw_data']['multiqc_fastp'][sample]['summary']['before_filtering']
    
    total_reads = {}
    total_bases = {}
    q20_bases = {}
    q30_bases = {}
    gc_content = {}
    at_content = {}
    at_percentages = {}
    for k, v in samples.items():
        for i in v:
            if i.startswith('before_filtering'):
                if 'total_reads' in v[i]:
                    if k not in total_reads:
                        total_reads[k] = 0
                    total_reads[k] += v[i]['total_reads']
                if 'total_bases' in v[i]:
                    if k not in total_bases:
                        total_bases[k] = 0
                    total_bases[k] += v[i]['total_bases']
                if 'q20_rate' in v[i]:
                    if k not in q20_bases:
                        q20_bases[k] = 0
                    q20_bases[k] += v[i]['q20_bases']
                if 'q30_rate' in v[i]:
                    if k not in q30_bases:
                        q30_bases[k] = 0
                    q30_bases[k] += v[i]['q30_bases']
                if 'gc_content' in v[i]:
                    if k not in gc_content:
                        gc_content[k] = 0
                    gc_content[k] += v[i]['gc_content'] * v[i]['total_bases']
                    at_content = 100 - v[i]['gc_content']
                    if k not in at_percentages:
                        at_percentages[k] = 0
                    at_percentages[k] += at_content * v[i]['total_bases']

    results = {'Sample': [], 'Total reads': [], 'Total bases': [], 'Q20 (%)': [], 'Q30 (%)': [], 'GC content': [], 'AT content': []}
    for k, v in samples.items():
        results['Sample'].append(k)
        results['Total reads'].append(total_reads[k])
        results['Total bases'].append(f"{total_bases[k]} ({total_bases[k] / 1000000000:.2f}GB)")
        if k in q20_bases and total_bases[k] > 0:
            q20_ratio = q20_bases[k] / total_bases[k] * 100
            results['Q20 (%)'].append(f"{q20_bases[k] / 1000000:.2f}M ({q20_ratio:.2f}%)")
        if k in q30_bases and total_bases[k] > 0:
            q30_ratio = q30_bases[k] / total_bases[k] * 100
            results['Q30 (%)'].append(f"{q30_bases[k] / 1000000:.2f}M ({q30_ratio:.2f}%)")
        if k in gc_content and total_bases[k] > 0:
            gc_ratio = gc_content[k] / total_bases[k] * 100
            results['GC content'].append(f"{gc_ratio:.2f}%")
            at_ratio = 100 - gc_ratio
            results['AT content'].append(f"{at_ratio:.2f}%")    

    df = pd.DataFrame(results)
    df_sorted = df.sort_values(by=['Sample'], ascending=  False)
    df_sorted = df_sorted[df_sorted['Sample'] != 'Undetermined']
    return df_sorted
    #with open(output_raw_file_path, 'w') as f:
    #    f.write(df_sorted.to_string(index=False))
    
def plot_total_bases(df_sorted):
    import matplotlib.pyplot as plt
    import re
    from matplotlib.patches import Rectangle

    df=df_sorted
    bases = df['Total bases'].apply(lambda x: float(re.findall(r'\((.*?)\)', x)[0].replace('GB','')))
    fig, ax = plt.subplots(figsize=(10,15))
    ax.barh(df['Sample'], bases, label='Total Bases(Gb)', color='#FFB6C1')

    x_values = list(range(0, int(max(bases))+2, 1))
    ax.set_xticks(x_values)
    ax.set_ylabel('')
    ax.set_xlabel('')
    ax.set_xlim(0, max(bases)*1.1)
    ax.set_title('Total bases')
    ax.set_xticklabels(x_values, color='gray')
    
    for val in x_values:
        ax.axvline(x=val, color='gray', linestyle='dotted')

    ax.bar_label(ax.containers[0], fmt='%.2f', color='#FFB6C1', fontsize=10)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.margins(y=0)
    ax.legend(loc='upper right', bbox_to_anchor=(1.2, 1), fontsize=10)
    plt.subplots_adjust(right=1)
    plt.savefig(f'{result_dir}total_bases.png', dpi=300)   
    plt.close()

def create_gc_content_graph(df_sorted):
    gc_percentages = df_sorted['GC content'].apply(lambda x: float(x.strip('%')))
    at_percentages = df_sorted['AT content'].apply(lambda x: float(x.strip('%')))
    fig, ax = plt.subplots(figsize=(10,20))
    
    ax.barh(df_sorted['Sample'], gc_percentages, label='GC Content', color='#87CEFA')
    ax.barh(df_sorted['Sample'], at_percentages, left=gc_percentages, label='AT Content', color='#FFC0CB')
    ax.set_ylabel('')
    ax.set_xlabel('')
    ax.set_xlim(0, 100)

    x_values = list(range(0, 101, 10))
    ax.set_xticks(x_values)
    ax.set_xticklabels(x_values, color='gray')

    for i, v in enumerate(gc_percentages):
        ax.text(v/2, i, f"{v:.2f}%", color='white', ha='center', va='center',fontsize=10)
        ax.text(100-v/2, i, f"{100-v:.2f}%", color='white', ha='center', va='center',fontsize=10)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.margins(y=0)
    ax.legend(loc='upper right', bbox_to_anchor=(1.2, 1), fontsize=10)
    plt.subplots_adjust(right=1)
    plt.savefig(f'{result_dir}gc_content.png')
    plt.close()

def create_q20_q30_graph(df_sorted):
    samples = list(df_sorted['Sample'])

    q20_perc = [float(val.split()[1].replace('(', '').replace('%)', '')) for val in df_sorted['Q20 (%)']]
    q30_perc = [float(val.split()[1].replace('(', '').replace('%)', '')) for val in df_sorted['Q30 (%)']]

    fig, ax = plt.subplots(figsize=(15,20))
    index = range(len(samples))
    bar_width = 0.35
    opacity = 0.8

    q20_bars = ax.barh(index, q20_perc, bar_width,
                    alpha=opacity, color='#ADD8E6', label='Q20')
    q30_bars = ax.barh([i + bar_width for i in index], q30_perc, bar_width,
                    alpha=opacity, color='#7FB3D5', label='Q30')

    ax.set_ylabel('')
    ax.set_xlabel('')
    ax.set_title('Q20 / Q30 (%)')
    ax.set_yticks([i + bar_width / 2 for i in index])
    ax.set_yticklabels(samples)
    
    x_values = list(range(0, 101, 10))
    ax.set_xticks(x_values)
    ax.set_xticklabels(x_values, color='gray')
    
    for val in x_values:
        ax.axvline(x=val, color='gray', linestyle='dotted')
    for i, v in enumerate(q20_perc):
        ax.text(v + 1, i - 0.1, str(v)+'%', color='#ADD8E6', fontsize=10)
    for i, v in enumerate(q30_perc):
        ax.text(v + 1, i + bar_width - 0.1, str(v)+'%', color='#7FB3D5', fontsize=10)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.margins(y=0)
    ax.legend(loc='upper right', bbox_to_anchor=(1.2, 1), fontsize=10)
    plt.subplots_adjust(right=0.7)
    ax.legend()
    fig.tight_layout()
    plt.savefig(f'{result_dir}q20_q30.png')
    plt.close()
    
def main():
    fastq_files = copy_raw_fastq_files(raw_fastq_path, fastq_dir)
    run_fastp(fastq_files)
    run_multiqc()
    md5sum_df = generate_md5sum_file(fastq_dir, output_md5sum_file_path)
    df_sorted =process_multiqc_fastp(json_file_path, output_raw_file_path)
    plot_total_bases(df_sorted)
    create_gc_content_graph(df_sorted)
    create_q20_q30_graph(df_sorted)
main()


