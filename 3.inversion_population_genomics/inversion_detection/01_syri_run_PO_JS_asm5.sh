#!/bin/bash
#Detect structural variation from pairwise alignment of the two genomes.
#Dependency: minimap2 v2.28-r1209, samtools v1.21, syri v1.7, protsr v1.1.1
REF=/home/yo/HDD/yo_analysis/HDD4_20230223/stickleback/pacbio/PO_male_HiFi/syri/PO_male_merged.syri.v20240517.fa
ASSEMBLY=/home/yo/HDD/yo_analysis/HDD4_20230223/stickleback/pacbio/PO_male_HiFi/syri/JS_male_merged.syri.v20240517.fa
DIR=syri20240917_PO_JS_asm5
OUT_PREFIX=syri20240917_PO_JS_asm5
NUMTH=5
GENOMES=genomes_PO_JS.txt

mkdir -p $DIR
cd $DIR

#Pairwise alignment of the two genomes
minimap2 -a -x asm5 --eqx -t $NUMTH $REF $ASSEMBLY > ${OUT_PREFIX}.sam
samtools sort -m 8G -@ $NUMTH ${OUT_PREFIX}.sam > ${OUT_PREFIX}.sort.bam
samtools index ${OUT_PREFIX}.sort.bam
rm ${OUT_PREFIX}.sam

#Detect SVs and plot
syri -c ${OUT_PREFIX}.sort.bam -r $REF -q $ASSEMBLY -F B --prefix ${OUT_PREFIX}  --nc $NUMTH　
plotsr -s 1000 -o ${OUT_PREFIX}syri.pdf --sr ${OUT_PREFIX}syri.out --genomes ../${GENOMES}

