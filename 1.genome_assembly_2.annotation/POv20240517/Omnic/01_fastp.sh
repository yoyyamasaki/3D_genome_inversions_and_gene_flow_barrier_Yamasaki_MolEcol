#!/bin/bash
#Clean Omnic reads
#Dependency: fastp v0.23.4

mkdir -p filtered

FASTQ=(XXX $(ls ./raw/*_L001_R1_001.fastq.gz))

ID=$(basename ${FASTQ[1]} _L001_R1_001.fastq.gz)

fastp -i ./raw/${ID}_L001_R1_001.fastq.gz -o ./filtered/PO_male_${ID}.trimmed_R1.fastq.gz \
-I ./raw/${ID}_L001_R2_001.fastq.gz -O ./filtered/PO_male_${ID}.trimmed_R2.fastq.gz \
--detect_adapter_for_pe --qualified_quality_phred 30 --unqualified_percent_limit 80 -w 20 \
-h ./filtered/${ID}_30_80_fastp.html -j ./filtered/${ID}_30_80_fastp.json


ID=$(basename ${FASTQ[2]} _L001_R1_001.fastq.gz)

fastp -i ./raw/${ID}_L001_R1_001.fastq.gz -o ./filtered/PO_male_${ID}.trimmed_R1.fastq.gz \
-I ./raw/${ID}_L001_R2_001.fastq.gz -O ./filtered/PO_male_${ID}.trimmed_R2.fastq.gz \
--detect_adapter_for_pe --qualified_quality_phred 30 --unqualified_percent_limit 80 -w 20 \
-h ./filtered/${ID}_30_80_fastp.html -j ./filtered/${ID}_30_80_fastp.json
