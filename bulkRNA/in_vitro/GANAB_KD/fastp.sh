#!/bin/bash

RawDir=~/Workshop/MultiOmics/GANAB2020/20210813UMUC3_KD/RNAseq/Rawdata
CleanDir=~/Workshop/MultiOmics/GANAB2020/20210813UMUC3_KD/RNAseq/CleanData

FileNames=("UMUC3-ctr1_FKDL210219094-1a" "UMUC3-ctr2_FKDL210219095-1a" "UMUC3-ctr3_FKDL210219096-1a" "UMUC3-ko1_FKDL210219097-1a" "UMUC3-ko2_FKDL210219098-1a" "UMUC3-ko3_FKDL210219099-1a")
SampleNames=("CTR_1" "CTR_2" "CTR_3" "KD_1" "KD_2" "KD_3")

for i in {0..5};
do
  fastp -i ${RawDir}/${FileNames[i]}_1.fq.gz -o ${CleanDir}/${SampleNames[i]}_R1.fq.gz -I ${RawDir}/${FileNames[i]}_2.fq.gz -O ${CleanDir}/${SampleNames[i]}_R2.fq.gz --compression 6 --thread 4 --report_title ${SampleNames[i]} --json ${CleanDir}/${SampleNames[i]}_fastp.json --html ${CleanDir}/${SampleNames[i]}_fastp.html --detect_adapter_for_pe
done
