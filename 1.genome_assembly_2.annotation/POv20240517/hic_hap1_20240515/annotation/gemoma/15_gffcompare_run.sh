#!/bin/bash
#Compare gemoma and braker2 results, and add blaker2 unique result to gemoma result.
#Dependency: GeMoMa v1.9, gffcompare v0.12.6, gffread v0.12.7

TARGET=../PO_male_hap1.v20240517.fa
GEMOMA_RESULT=final_annotation.gff
BRAKER_RESULT=../braker3_prot/braker.gtf

#Extract genes which originally predicetd by Braker2
ln -s ${BRAKER_RESULT} ./
gffcompare  -r ${GEMOMA_RESULT} -o ./gffcmp $(basename ${BRAKER_RESULT})

cat gffcmp.braker.gtf.tmap | awk '$3=="u" || $3=="i" || $3=="y" { print $4 }' > braker.keep.list

cat ${BRAKER_RESULT} | grep -wf braker.keep.list > braker.keep.gtf

#Add blaker2 result to Gemoma result
java -jar ~/software/GeMoMa-1.9/GeMoMa-1.9.jar CLI AnnotationEvidence a=braker.keep.gtf g=${TARGET}

java -jar ~/software/GeMoMa-1.9/GeMoMa-1.9.jar CLI GAF g=${GEMOMA_RESULT} g=annotation_with_attributes.gff tf=true

#Add gene names
java -jar ~/software/GeMoMa-1.9/GeMoMa-1.9.jar CLI AnnotationFinalizer g=${TARGET} \
 a=filtered_predictions_1.gff \
 u=YES i=denoised_introns.gff \
 tf=true \
 rename=COMPOSED \
 p=POhap1 \
 infix=G \
 s="" \
 csp="" \
 crp="" \
 d=6

#Make transcript, CDS, and protein sequence based on the annotation
gffread -g ${TARGET} -w PO_male_hap1.v20240517.transcript.fasta -x PO_male_hap1.v20240517.CDS.fasta -y PO_male_hap1.v20240517.protein.fasta final_annotation_1.gff

