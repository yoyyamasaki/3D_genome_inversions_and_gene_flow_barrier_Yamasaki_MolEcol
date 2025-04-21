#!/bin/bash
#Visualize Hi-C heat map with TADs and inversions.
#Dependency: HiGlass v1.13.4

#File path
REF=~/db/PO_merged_maskXII_XXI_Y_20240517/PO_male_merged.v20240517.maskXII_XXI_Y.fa

GENOME=PO_male_merged.v20240517.maskXII_XXI_Y.fa.genome

ASSEMBLY=POmask.XII.XXI.Y.v20240517

GTF_INGEST=~/db/PO_merged_maskXII_XXI_Y_20240517/PO_male_merged.v20240517.mod.maskXII_XXI_Y.hgbed
GTF_INGEST_FILE=~/HDD/yo_analysis/HDD4_20230223/stickleback/HiC/chromatin_analysis/PO_hifi/iconic_PO_ref_PO_merge_maskXII_XXI_Y_paitrools/PO_male_merged.v20240517.mod.maskXII_XXI_Y.hgbed

PROJECT=PO_refPO

MCOOL=~/HDD/yo_analysis/HDD4_20230223/stickleback/HiC/chromatin_analysis/PO_hifi/iconic_PO_ref_PO_merge_maskXII_XXI_Y_paitrools/iconic_PO_refPO_merged_maskXII_XXI_Y_pairtools.hic.mcool

TAD1=~/HDD/yo_analysis/HDD4_20230223/stickleback/HiC/chromatin_analysis/PO_hifi/iconic_PO_ref_PO_merge_maskXII_XXI_Y_paitrools/iconic_PO_refPO_merged_maskXII_XXI_Y_pairtools_KR/spectraltad/iconic_PO_refPO_merged_maskXII_XXI_Y_pairtools_spectraltad_25kb_Level1.bed
TAD2=~/HDD/yo_analysis/HDD4_20230223/stickleback/HiC/chromatin_analysis/PO_hifi/iconic_PO_ref_PO_merge_maskXII_XXI_Y_paitrools/iconic_PO_refPO_merged_maskXII_XXI_Y_pairtools_KR/spectraltad/iconic_PO_refPO_merged_maskXII_XXI_Y_pairtools_spectraltad_25kb_Level2.bed

INVERSION=/home/yo/HDD/yo_analysis/HDD4_20230223/stickleback/pacbio/PO_male_HiFi/syri/syri20240917_PO_JS_asm5/syri20240917_PO_JS_asm5syri.INV.bed
INVERSION2=/home/yo/HDD/yo_analysis/HDD4_20230223/stickleback/pacbio/PO_male_HiFi/syri/syri20241127_PO_v5/syri20241127_PO_v5syri.INV.bed
INVERSION3=/home/yo/HDD/yo_analysis/HDD4_20230223/stickleback/pacbio/PO_male_HiFi/syri/syri20240917_PO_JS_asm5/syri20240917_PO_JS_asm5syri.INV.bed

cp ${REF}.genome ./${GENOME}

#higlass-manage
higlass-manage start  --temp-dir ${PWD} --data-dir ${PWD}/higlass

#Ingest mcool
higlass-manage ingest \
                      ${MCOOL} \
                      --project-name ${PROJECT} \
                      --assembly ${ASSEMBLY} \
                      --chromsizes-filename ${REF}.genome

#Ingest Chromosome size
higlass-manage ingest \
    ${PWD}/${GENOME} \
    --filetype chromsizes-tsv \
    --datatype chromsizes \
    --project-name ${PROJECT}  \
    --assembly ${ASSEMBLY} 
    
#Ingest gtf
clodius aggregate bedfile --chromsizes-filename ${PWD}/${GENOME} ${GTF_INGEST}

cp ${GTF_INGEST}.beddb ./
higlass-manage ingest \
                      ${GTF_INGEST_FILE}.beddb \
                      --project-name ${PROJECT} \
                      --filetype beddb \
                      --datatype gene-annotation  \
                      --assembly ${ASSEMBLY} \
                      --chromsizes-filename ${PWD}/${GENOME}
                     
#Ingest TAD
#TAD1
clodius aggregate bedpe \
    --chromsizes-filename ${PWD}/${GENOME} \
    --chr1-col 1 --from1-col 2 --to1-col 3 \
    --chr2-col 1 --from2-col 2 --to2-col 3 \
    --assembly ${ASSEMBLY} \
    --output-file ${TAD1}.bedpedb \
    ${TAD1}
    
higlass-manage ingest \
    --filetype bed2ddb \
    --datatype 2d-rectangle-domains \
    --project-name ${PROJECT} \
    --assembly ${ASSEMBLY} \
    --chromsizes-filename ${PWD}/${GENOME} \
    ${TAD1}.bedpedb

#TAD2
clodius aggregate bedpe \
    --chromsizes-filename ${PWD}/${GENOME} \
    --chr1-col 1 --from1-col 2 --to1-col 3 \
    --chr2-col 1 --from2-col 2 --to2-col 3 \
    --assembly ${ASSEMBLY} \
    --output-file ${TAD2}.bedpedb \
    ${TAD2}
    
higlass-manage ingest \
    --filetype bed2ddb \
    --datatype 2d-rectangle-domains \
    --project-name ${PROJECT} \
    --assembly ${ASSEMBLY} \
    --chromsizes-filename ${PWD}/${GENOME} \
    ${TAD2}.bedpedb


#Inversion
clodius aggregate bedfile \
                          --chromsizes-filename ${PWD}/${GENOME} \
                          --delimiter $'\t' \
                          --assembly ${ASSEMBLY} \
                          ${INVERSION}
                          
higlass-manage ingest \
                      ${INVERSION}.beddb \
                      --project-name ${PROJECT} \
                      --datatype bedlike  \
                      --assembly ${ASSEMBLY} \
                      --chromsizes-filename ${PWD}/${GENOME}
                      
clodius aggregate bedfile \
                          --chromsizes-filename ${PWD}/${GENOME} \
                          --delimiter $'\t' \
                          --assembly ${ASSEMBLY} \
                          ${INVERSION2}
                          
higlass-manage ingest \
                      ${INVERSION2}.beddb \
                      --project-name ${PROJECT} \
                      --datatype bedlike  \
                      --assembly ${ASSEMBLY} \
                      --chromsizes-filename ${PWD}/${GENOME}
                      
clodius aggregate bedfile \
                          --chromsizes-filename ${PWD}/${GENOME} \
                          --delimiter $'\t' \
                          --assembly ${ASSEMBLY} \
                          ${INVERSION3}
                          
higlass-manage ingest \
                      ${INVERSION3}.beddb \
                      --project-name ${PROJECT} \
                      --datatype bedlike  \
                      --assembly ${ASSEMBLY} \
                      --chromsizes-filename ${PWD}/${GENOME}
                      
#URL: http://localhost:8989/app

higlass-manage stop
higlass-manage shell

higlass-manage create superuser

#To remove data, login from the following URL
#URL: http://localhost:8989/admin/

