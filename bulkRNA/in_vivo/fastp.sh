#! /bin/bash
#$ -N fastp
#$ -q all_el7.q
#$ -P el7_project
#$ -l vf=64G
#$ -wd /datapool/pengguoyu/ATACseq/20200916ATAC/Log
#$ -e fastp_2.e
#$ -o fastp_2.o

# conda activate DNA
RawDir=/share/service04/pengguoyu/ATAC/20200916ATAC/20201110RNA
CleanDir=/datapool/pengguoyu/ATACseq/20200916ATAC/20201110RNA/CleanData
fastp -i ${RawDir}/KO_L1_341X41.R1.fastq.gz -o ${CleanDir}/KO4_R1.fq.gz -I ${RawDir}/KO_L1_341X41.R2.fastq.gz -O ${CleanDir}/KO4_R2.fq.gz --detect_adapter_for_pe --cut_front --cut_tail --json ${CleanDir}/KO4.json --html ${CleanDir}/KO4.html --report_title KO4 
fastp -i ${RawDir}/KO-SNT_L1_342X42.R1.fastq.gz -o ${CleanDir}/KO5_R1.fq.gz -I ${RawDir}/KO-SNT_L1_342X42.R2.fastq.gz -O ${CleanDir}/KO5_R2.fq.gz --detect_adapter_for_pe --cut_front --cut_tail --json ${CleanDir}/KO5.json --html ${CleanDir}/KO5.html --report_title KO5 
fastp -i ${RawDir}/KO-SNT_L1_343X43.R1.fastq.gz -o ${CleanDir}/KO6_R1.fq.gz -I ${RawDir}/KO-SNT_L1_343X43.R2.fastq.gz -O ${CleanDir}/KO6_R2.fq.gz --detect_adapter_for_pe --cut_front --cut_tail --json ${CleanDir}/KO6.json --html ${CleanDir}/KO6.html --report_title KO6 
fastp -i ${RawDir}/WT1_L1_336X36.R1.fastq.gz -o ${CleanDir}/WT4_R1.fq.gz -I ${RawDir}/WT1_L1_336X36.R2.fastq.gz -O ${CleanDir}/WT4_R2.fq.gz --detect_adapter_for_pe --cut_front --cut_tail --json ${CleanDir}/WT4.json --html ${CleanDir}/WT4.html --report_title WT4
fastp -i ${RawDir}/WT2_L1_337X37.R1.fastq.gz -o ${CleanDir}/WT5_R1.fq.gz -I ${RawDir}/WT2_L1_337X37.R2.fastq.gz -O ${CleanDir}/WT5_R2.fq.gz --detect_adapter_for_pe --cut_front --cut_tail --json ${CleanDir}/WT5.json --html ${CleanDir}/WT5.html --report_title WT5
fastp -i ${RawDir}/WT3_L1_339X39.R1.fastq.gz -o ${CleanDir}/WT6_R1.fq.gz -I ${RawDir}/WT3_L1_339X39.R2.fastq.gz -O ${CleanDir}/WT6_R2.fq.gz --detect_adapter_for_pe --cut_front --cut_tail --json ${CleanDir}/WT6.json --html ${CleanDir}/WT6.html --report_title WT6
