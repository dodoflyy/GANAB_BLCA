#!/bin/bash

CleanDir=../CleanData
SalmonDir=../Salmon
IndexPath=../GRCh38/transcripts_index_salmon
GRCh38GTF=../GRCh38/gencode.v35.annotation.gtf

Samples=("CTR_1" "CTR_2" "CTR_3" "KD_1" "KD_2" "KD_3")
for Sample in ${Samples[@]}
do
  echo "++++++ Sample: ${Sample} ++++++"
  salmon quant -l IU -i ${IndexPath} -1 ${CleanDir}/${Sample}_R1.fq.gz -2 ${CleanDir}/${Sample}_R2.fq.gz -p 4 -g ${GRCh38GTF} -o ${SalmonDir}/${Sample}
done
