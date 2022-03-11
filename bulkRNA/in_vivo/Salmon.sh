#! /bin/bash

CleanDir=../CleanData
SalmonDir=../Salmon
RefDir=../GRCh38
SampleList=("WT4" "WT5" "WT6" "KO4" "KO5" "KO6")

for Sample in ${SampleList[@]}
do
  salmon quant --libType ISR --index ${RefDir}/hsa_transcripts_index -p 8 -g ${RefDir}/gencode.v33.annotation.gtf --validateMappings -1 ${CleanDir}/${Sample}_R1.fq.gz -2 ${CleanDir}/${Sample}_R2.fq.gz -o ${SalmonDir}/${Sample}
done
