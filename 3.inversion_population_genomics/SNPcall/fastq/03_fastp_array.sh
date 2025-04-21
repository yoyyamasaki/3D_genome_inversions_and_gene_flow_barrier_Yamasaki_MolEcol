#!/bin/bash
#Filtering WGS fastq files in parallel.
#Dependnecy: fastp v0.23.2

DIR=./raw
THREADS=1

fastq=(XXX $(ls ${DIR}/*_1.fastq.gz))
id=($(basename ${fastq[${1}]} _1.fastq.gz))

mkdir -p dedup
mkdir -p filtered

#Remove PCR duplicates
fastp -i ${DIR}/${id}_1.fastq.gz -o ./dedup/${id}_dedup_1P.fq.gz \
-I ${DIR}/${id}_2.fastq.gz -O ./dedup/${id}_dedup_2P.fq.gz \
-w ${THREADS} -D \
-h ./dedup/${id}_fastp.html -j ./dedup/${id}_fastp.json

#Filtering by sliding window
fastp -i ./dedup/${id}_dedup_1P.fq.gz -o ./filtered/${id}_trimmed_1P.fq.gz \
-I ./dedup/${id}_dedup_2P.fq.gz -O ./filtered/${id}_trimmed_2P.fq.gz \
--detect_adapter_for_pe --cut_right --cut_window_size 4 --cut_mean_quality 20 -l 35 -w ${THREADS} \
-h ./filtered/${id}_fastp.html -j ./filtered/${id}_fastp.json
