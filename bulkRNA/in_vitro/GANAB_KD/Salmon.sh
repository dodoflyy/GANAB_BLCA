#!/bin/bash

CleanDir=~/Workshop/MultiOmics/GANAB2020/20210813UMUC3_KD/RNAseq/CleanData
SalmonDir=~/Workshop/MultiOmics/GANAB2020/20210813UMUC3_KD/RNAseq/Salmon
IndexPath=~/Database/GRCh38/GENCODE/transcripts_index_salmon
GRCh38GTF=~/Database/GRCh38/GENCODE/gencode.v35.annotation.gtf

Samples=("CTR_1" "CTR_2" "CTR_3" "KD_1" "KD_2" "KD_3")
for Sample in ${Samples[@]}
do
  echo "++++++ Sample: ${Sample} ++++++"
  salmon quant -l IU -i ${IndexPath} -1 ${CleanDir}/${Sample}_R1.fq.gz -2 ${CleanDir}/${Sample}_R2.fq.gz -p 4 -g ${GRCh38GTF} -o ${SalmonDir}/${Sample}
done
