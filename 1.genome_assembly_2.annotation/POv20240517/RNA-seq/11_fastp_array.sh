#! /bin/bash
#Run fastp to clean RNA-seq reads
#Dependency: fastp v0.23.4

DIR=./raw

THREADS=1

fastq=(XXX $(ls ${DIR}/*_R1.fastq.gz))
id=($(basename ${fastq[${1}]} _R1.fastq.gz))

mkdir -p filtered

fastp -i ${DIR}/${id}_R1.fastq.gz -o ./filtered/${id}_trimmed_1P.fq.gz \
-I ${DIR}/${id}_R2.fastq.gz  -O ./filtered/${id}_trimmed_2P.fq.gz \
--detect_adapter_for_pe --cut_right --cut_window_size 4 --cut_mean_quality 20 -l 35 -w ${THREADS} \
-h ./filtered/${id}_fastp.html -j ./filtered/${id}_fastp.json
