#!/bin/bash
#Mapping iconHi-C reads to the reference genome. Then extract valid Hi-C read pairs.
#Dependency: bwa-mem2 v2.2.1, pairtools v1.0.3, samtools v1.19.2, juicertools v1.22.01, hicexplorer v3.7.3

THREAD=48
REF=~/db/JS_male_merged_mask_chrXII_chrXXI_chrXIY_20240517/JS_male_merged.v20240517.maskXII_XXI_Y.fa
FASTQ_R1=/home/yo/HDD/yo_analysis/HDD4_20230223/stickleback/HiC/chromatin_analysis/JS/merged/P606_02_1_merged.R1.fq.gz
FASTQ_R2=/home/yo/HDD/yo_analysis/HDD4_20230223/stickleback/HiC/chromatin_analysis/JS/merged/P606_02_1_merged.R2.fq.gz
OUT_HEADER=iconhic_JS_refJS_merged_maskXII_XXI_Y_pairtools
DIR=./

CHR=($(grep ">chr" ${REF}|perl -pe 's/>//g'))

mkdir -p ${DIR}

tmp_dir=${DIR}/tmp

mkdir -p ${tmp_dir}

#Mapping omnic read to assembled contigs and extract valid omnic reads
bwa-mem2 mem -5SP -T0 -t ${THREAD} ${REF} ${FASTQ_R1} ${FASTQ_R2}| \
pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in ${THREAD} --nproc-out ${THREAD} --chroms-path ${REF} | \
pairtools sort --tmpdir=${tmp_dir} --nproc ${THREAD}|\
pairtools dedup --nproc-in ${THREAD} --nproc-out ${THREAD} --mark-dups --output-stats ${DIR}/stats.txt|\
pairtools split --nproc-in ${THREAD} --nproc-out ${THREAD} --output-pairs ${DIR}/${OUT_HEADER}.pairs --output-sam -|\
samtools view -bS -@${THREAD} | \
samtools sort -@${THREAD} -o ${DIR}/${OUT_HEADER}.bam;samtools index ${DIR}/${OUT_HEADER}.bam

bgzip ${DIR}/${OUT_HEADER}.pairs

pairix ${DIR}/${OUT_HEADER}.pairs.gz

#Check mapping statistics
~/software/Omni-C-main/get_qc.py -p ${DIR}/stats.txt > ${DIR}/Hi-c_stats.txt

#Create .hic format contact matrix 
java -Xmx48000m -Djava.awt.headless=true -jar ~/software/juicer_tools_1.22.01.jar pre ${DIR}/${OUT_HEADER}.pairs.gz ${DIR}/${OUT_HEADER}.hic ${REF}.genome

#Convert .hic format to .cool format
hicConvertFormat -m ${DIR}/${OUT_HEADER}.hic --inputFormat hic --outputFormat cool -o ${DIR}/${OUT_HEADER}.hic.mcool 

