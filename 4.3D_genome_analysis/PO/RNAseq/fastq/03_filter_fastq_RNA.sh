#!/bin/bash
#Filtering RNA-seq fastq files in parallel.
#Dependnecy: fastp v0.23.2

FASTQ_DIR=raw
OUT_DIR=filtered

files=($(ls $FASTQ_DIR/*_R1.fastq.gz))

id=()
for i in ${files[@]};do
	id+=($(basename $i _R1.fastq.gz))
done

mkdir $OUT_DIR

for i in ${id[@]};do
	fastp -i ${FASTQ_DIR}/${i}_R1.fastq.gz -o ${OUT_DIR}/${i}_filtered_1P.fq.gz -I ${FASTQ_DIR}/${i}_R2.fastq.gz -O ${OUT_DIR}/${i}_filtered_2P.fq.gz \
	--detect_adapter_for_pe -w 16 -h ${OUT_DIR}/${i}_fastp.html -j ${OUT_DIR}/${i}_fastp.json
done
