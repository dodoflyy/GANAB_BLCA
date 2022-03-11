#! /bin/bash
#$ -N Salmon
#$ -q all_el7.q
#$ -P el7_project
#$ -l vf=84G
#$ -wd /datapool/pengguoyu/ATACseq/20200916ATAC/20201110RNA/Log
#$ -e Salmon.e
#$ -o Salmon.o
# conda activate RNA-seq

CleanDir=/datapool/pengguoyu/ATACseq/20200916ATAC/20201110RNA/CleanData
SalmonDir=/datapool/pengguoyu/ATACseq/20200916ATAC/20201110RNA/Salmon
SampleList=("WT4" "WT5" "WT6" "KO4" "KO5" "KO6")
for Sample in ${SampleList[@]}
do
  salmon quant --libType ISR --index /share/database/openData/GRCh38_hg38/hsa_transcripts_index -p 8 -g /share/database/openData/GRCh38_hg38/gencode.v33.annotation.gtf --validateMappings -1 ${CleanDir}/${Sample}_R1.fq.gz -2 ${CleanDir}/${Sample}_R2.fq.gz -o ${SalmonDir}/${Sample}
done
