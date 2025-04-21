#!/bin/bash
#Conduct haplotype phased assembly by hifiasm.
#Dependency: hifiasm_0.19.9-r616

#Pacbio Hifi read path
HIFI=/home/yo/HDD/yo_analysis/HDD4_20230223/stickleback/pacbio/PO_male_HiFi/HiFi_DeepConsensus/m64521e_240509_120656.hifi_reads.fastq.gz
#Omnic read path
HIC1_F=/home/yo/HDD/yo_analysis/HDD4_20230223/stickleback/pacbio/PO_male_HiFi/OmniC/mnase/filtered/JS_male_OmniC-G-acureatus_S10.trimmed_R1.fastq.gz
HIC1_R=/home/yo/HDD/yo_analysis/HDD4_20230223/stickleback/pacbio/PO_male_HiFi/OmniC/mnase/filtered/JS_male_OmniC-G-acureatus_S10.trimmed_R2.fastq.gz
HIC2_F=/home/yo/HDD/yo_analysis/HDD4_20230223/stickleback/pacbio/PO_male_HiFi/OmniC/mnase/filtered/JS_male_OmniC-G-acureatus_S3.trimmed_R1.fastq.gz
HIC2_R=/home/yo/HDD/yo_analysis/HDD4_20230223/stickleback/pacbio/PO_male_HiFi/OmniC/mnase/filtered/JS_male_OmniC-G-acureatus_S3.trimmed_R2.fastq.gz

OUT=PO_male_DC_OmniCtrimmed
OUTDIR=hifiasm_result20240515
THREADS=60

mkdir -p ${OUTDIR}

#Run hifiasm
/home/yo/software/hifiasm_0.19.9-r616/hifiasm -o ./${OUTDIR}/${OUT} -t ${THREADS} --h1 ${HIC1_F},${HIC2_F} --h2 ${HIC1_R},${HIC2_R} ${HIFI}

#Make fasta from gfa
awk '/^S/{print ">"$2;print $3}' ./${OUTDIR}/${OUT}.hic.hap1.p_ctg.gfa > ./${OUTDIR}/${OUT}.hic.hap1.p_ctg.fa
awk '/^S/{print ">"$2;print $3}' ./${OUTDIR}/${OUT}.hic.hap2.p_ctg.gfa > ./${OUTDIR}/${OUT}.hic.hap2.p_ctg.fa
awk '/^S/{print ">"$2;print $3}' ./${OUTDIR}/${OUT}.hic.p_utg.gfa > ./${OUTDIR}/${OUT}.hic.p_utg.fa

