#!/bin/bash
#Write respective chromosomes in different files. Also count gap number.
#Dependency: samtools v 1.18, assembly-stats v1.0.1

REF=JS_male_DC_OmniCtrimmed.enzymemix.hic.hap1.yahs.min1000_JBAT.FINAL.rename.fa
HEADER=$(basename ${REF} .fa)
DIR=chr

CHR=($(grep ">" ${REF}|grep chr|perl -pe 's/>//g'))

mkdir -p ${DIR}

rm ${HEADER}.chr_gaps.txt
touch ${HEADER}.chr_gaps.txt
for i in ${CHR[@]};do
 #separate chr fasta
 samtools faidx ${REF} ${i} > ${DIR}/${HEADER}.${i}.fa
 
 #count gap number
 gap=$(assembly-stats ${DIR}/${HEADER}.${i}.fa |grep Gaps|perl -pe 's/Gaps = //g')
 echo -e ${i}"\t"${gap} >> ${HEADER}.chr_gaps.txt
  
done
