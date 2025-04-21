#!/bin/bash
#Clean Omnic reads
#Dependency: fastp v0.23.4

mkdir -p filtered

fastq=(XXX $(ls ./raw/*_L001_R1_001.fastq.gz))

id=$(basename ${fastq[1]} _L001_R1_001.fastq.gz)

fastp -i ./raw/${id}_L001_R1_001.fastq.gz -o ./filtered/JS_male_${id}.trimmed_R1.fastq.gz \
-I ./raw/${id}_L001_R2_001.fastq.gz -O ./filtered/JS_male_${id}.trimmed_R2.fastq.gz \
--detect_adapter_for_pe --qualified_quality_phred 30 --unqualified_percent_limit 80 -w 20 \
-h ./filtered/${id}_30_80_fastp.html -j ./filtered/${id}_30_80_fastp.json


fastq=(XXX $(ls ./raw/*_L001_R1_001.fastq.gz))

id=$(basename ${fastq[2]} _L001_R1_001.fastq.gz)

fastp -i ./raw/${id}_L001_R1_001.fastq.gz -o ./filtered/JS_male_${id}.trimmed_R1.fastq.gz \
-I ./raw/${id}_L001_R2_001.fastq.gz -O ./filtered/JS_male_${id}.trimmed_R2.fastq.gz \
--detect_adapter_for_pe --qualified_quality_phred 30 --unqualified_percent_limit 80 -w 20 \
-h ./filtered/${id}_30_80_fastp.html -j ./filtered/${id}_30_80_fastp.json
