#!/bin/bash
#Convert scaffold name to chromosome number, and convert direction of chromosome to that of reference.
#chromosome name and direction was identified by dgenies.
#Make table:1st column=scaffold name, 2nd column=chromosome name, 3rd column=direction of chromosome (+ or -)
#Dependency: seqkit v2.8.0

REF=JS_male_DC_OmniCtrimmed.enzymemix.hic.hap1.yahs.min1000_JBAT.FINAL.fa
CONVERT_TABLE=scaffold_directions.txt
OUT_HEADER=JS_male_DC_OmniCtrimmed.enzymemix.hic.hap1.yahs.min1000_JBAT.FINAL.rename
OUTDIR=renamed

mkdir -p ${OUTDIR}

rm ${OUTDIR}/${OUT_HEADER}.fa
touch ${OUTDIR}/${OUT_HEADER}.fa

while read line; do
 data=(${line})
 scaffold=${data[0]}
 chr=${data[1]}
 direction=${data[2]}
 
 if test "${direction}" = "+"; then
 
  seqkit grep -n -p ${scaffold} ${REF}|perl -pe "s/${scaffold}/${chr}/g" >> ${OUTDIR}/${OUT_HEADER}.fa
 
 elif test "${direction}" = "-"; then
 
  seqkit grep -n -p ${scaffold} ${REF}|seqkit seq -p -r -v -t dna|perl -pe "s/${scaffold}/${chr}/g" >> ${OUTDIR}/${OUT_HEADER}.fa
 
 else
 
  seqkit grep -n -p ${scaffold} ${REF} >> ${OUTDIR}/${OUT_HEADER}.fa
   
 fi
done < ${CONVERT_TABLE}

seqkit sort -n ${OUTDIR}/${OUT_HEADER}.fa > ${OUTDIR}/${OUT_HEADER}.sorted.fa 
