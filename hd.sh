#!/bin/bash

run_id=$1
sample_id=$2
datetime=`date +%Y%m%d_%H%M%S`
python_path=$(which python)

sbatch <<EOT
#! /bin/bash
#SBATCH -J ${sample_id}
#SBATCH --nodelist=com001
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=6GB
#SBATCH --ntasks=1
#SBATCH --output=/data/home/lims/slurm_log/${datetime}_${sample_id}_HD.out
#SBATCH --error=/data/home/lims/slurm_log/${datetime}_${sample_id}_HD.err

source /data/tools/miniconda3/bin/activate
conda activate tsoe

echo "sample id : ${sample_id}"
which python

python /data/script/lims/Hereditary_Disease/variant_call_script.py -o /data/analysis/service/hd/ ${run_id} ${sample_id}
EOT

