#!/bin/bash

main_dir=..
raw_dir=${main_dir}/raw_data
clean_dir=${main_dir}/clean_data

samples=(CTR1 CTR2 CTR3 GKD1 GKD2 GKD3 PKD1 PKD2 PKD3)
for sample in ${samples[@]}
do
  echo "------ ${sample} ------"
  fastp -i ${raw_dir}/${sample}_R1.fq.gz -o ${clean_dir}/${sample}_R1.fq.gz -I ${raw_dir}/${sample}_R2.fq.gz -O ${clean_dir}/${sample}_R2.fq.gz --compression 6 --report_title ${sample} --json ${clean_dir}/${sample}_fastp.json --html ${clean_dir}/${sample}.html --detect_adapter_for_pe
done
