#!/bin/bash
#Compare gemoma and braker2 results, and add blaker2 unique result to gemoma result.
#Dependency: GeMoMa v1.9, gffcompare v0.12.6, gffread v0.12.7
#We did not use blaker2 result for the annotation of JS_male_hap2.v20240517 because of the error of braker2.

TARGET=../JS_male_hap2.v20240517.softmasked.fa
GEMOMA_RESULT=final_annotation.gff

#Add gene names
java -jar ~/software/GeMoMa-1.9/GeMoMa-1.9.jar CLI AnnotationFinalizer g=${TARGET} \
 a=final_annotation.gff \
 u=YES i=denoised_introns.gff \
 tf=true \
 rename=COMPOSED \
 p=JShap2 \
 infix=G \
 s="" \
 csp="" \
 crp="" \
 d=6
 
#Make transcript, CDS, and protein sequence based on the annotation
gffread -g ${TARGET} -w JS_male_hap2.v20240517.transcript.fasta -x JS_male_hap2.v20240517.CDS.fasta -y JS_male_hap2.v20240517.protein.fasta final_annotation_1.gff



