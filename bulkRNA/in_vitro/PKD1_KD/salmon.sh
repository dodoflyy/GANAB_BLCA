#!/bin/bash

main_dir=/home/pengguoyu/Workshop/MultiOmics/GANAB2020/20220119UMUC3_PKD1_GANAB
ref_dir=/home/pengguoyu/Database/GRCh38/GENCODE
clean_dir=${main_dir}/clean_data
salmon_dir=${main_dir}/salmon
index_path=${ref_dir}/transcripts_index_salmon
anno_path=${ref_dir}/gencode.v35.annotation.gtf

Samples=(CTR1 CTR2 CTR3 GKD1 GKD2 GKD3 PKD1 PKD2 PKD3)
for Sample in ${Samples[@]}
do
  echo "++++++ ${Sample} ++++++"
  salmon quant -l IU -i ${index_path} -1 ${clean_dir}/${Sample}_R1.fq.gz -2 ${clean_dir}/${Sample}_R2.fq.gz -p 4 -g ${anno_path} -o ${salmon_dir}/${Sample}
done
