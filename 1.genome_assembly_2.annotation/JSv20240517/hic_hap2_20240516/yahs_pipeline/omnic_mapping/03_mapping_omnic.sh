#!/bin/bash
#Mapping omnic read to the contigs for the following yahs scaffolding.
#Dependency: bwa-mem2 v2.2.1, pairtools v1.0.3, samtools v1.19.2

#Contig file path
REF=/home/yo/HDD/yo_analysis/HDD4_20230223/stickleback/pacbio/JS_male_HiFi/hifiasm_DC_OmniCtrimmed_enzymemix_20240515/hifiasm_result20240515/JS_male_DC_OmniCtrimmed.hic.hap2.p_ctg.fa
#Omnic read path
HIC1_F=/home/yo/HDD/yo_analysis/HDD4_20230223/stickleback/pacbio/JS_male_HiFi/OmniC/enzyme_mix/filtered/JS_male_OmniC-G-nipponicus-2_S2.trimmed_R1.fastq.gz
HIC1_R=/home/yo/HDD/yo_analysis/HDD4_20230223/stickleback/pacbio/JS_male_HiFi/OmniC/enzyme_mix/filtered/JS_male_OmniC-G-nipponicus-2_S2.trimmed_R2.fastq.gz
HIC2_F=/home/yo/HDD/yo_analysis/HDD4_20230223/stickleback/pacbio/JS_male_HiFi/OmniC/enzyme_mix/filtered/JS_male_OmniC-G-nipponicus-2_S9.trimmed_R1.fastq.gz
HIC2_R=/home/yo/HDD/yo_analysis/HDD4_20230223/stickleback/pacbio/JS_male_HiFi/OmniC/enzyme_mix/filtered/JS_male_OmniC-G-nipponicus-2_S9.trimmed_R2.fastq.gz 

TMP_DIR=/home/yo/HDD/yo_analysis/HDD4_20230223/stickleback/pacbio/JS_male_HiFi/hifiasm_DC_OmniCtrimmed_enzymemix_20240515/hic_hap2_20240516/yahs_pipeline/omnic_mapping/tmp
OUT_HEADER=JS_male_hic_hap2_20240515_omnic
THREAD=40

mkdir -p ${TMP_DIR}

#Mapping omnic read to assembled contigs and extract valid omnic reads
bwa-mem2 index ${REF}

bwa-mem2 mem -5SP -T0 -t ${THREAD} ${REF} ${HIC1_F} ${HIC1_R}| \
pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in ${THREAD} --nproc-out ${THREAD} --chroms-path ${REF} | \
pairtools sort --tmpdir=${TMP_DIR} --nproc ${THREAD}|\
pairtools dedup --nproc-in ${THREAD} --nproc-out ${THREAD} --mark-dups --output-stats stats_1.txt|
pairtools split --nproc-in ${THREAD} --nproc-out ${THREAD} --output-pairs ${OUT_HEADER}_1.pairs --output-sam -|\
samtools view -bS -@${THREAD} | \
samtools sort -@${THREAD} -o ./${OUT_HEADER}_1.bam;samtools index ./${OUT_HEADER}_1.bam

bwa-mem2 mem -5SP -T0 -t ${THREAD} ${REF} ${HIC2_F} ${HIC2_R}| \
pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in ${THREAD} --nproc-out ${THREAD} --chroms-path ${REF} | \
pairtools sort --tmpdir=${TMP_DIR} --nproc ${THREAD}|\
pairtools dedup --nproc-in ${THREAD} --nproc-out ${THREAD} --mark-dups --output-stats stats_2.txt|
pairtools split --nproc-in ${THREAD} --nproc-out ${THREAD} --output-pairs ${OUT_HEADER}_2.pairs --output-sam -|\
samtools view -bS -@${THREAD} | \
samtools sort -@${THREAD} -o ./${OUT_HEADER}_2.bam;samtools index ./${OUT_HEADER}_2.bam

#Merge two omnic mapping results
samtools merge ./${OUT_HEADER}.bam ./${OUT_HEADER}_1.bam ./${OUT_HEADER}_2.bam

rm -r ${TMP_DIR}
