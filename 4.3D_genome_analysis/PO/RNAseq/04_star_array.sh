#!/bin/bash
#Conduct RNA-seq read mapping by star in parallel.
#Dependnecy: STAR v.2.7.11b, samtools v1.6, featurecount v2.0.6

DIR=../fastq/filtered
REF=/home/yo/db/PO_merged_maskXII_XXI_Y_20240517
GTF=${REF}/PO_male_merged.v20240517.mod.maskXII_XXI_Y.gtf

mkdir -p mapping

fastq=(XXX $(ls ${DIR}/*_filtered_1P.fq.gz))
id=($(basename ${fastq[${1}]} _filtered_1P.fq.gz))

STAR_2.7.11 --runMode alignReads --genomeDir ${REF} \
--readFilesIn ${DIR}/${id}_filtered_1P.fq.gz ${DIR}/${id}_filtered_2P.fq.gz \
--sjdbGTFfile ${GTF} --readFilesCommand zcat --runThreadN 1 --outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix ./mapping/${id} --limitBAMsortRAM 40000000000 --outSAMunmapped Within --outSAMattributes Standard --outFilterMultimapNmax 1

samtools index ./mapping/${id}Aligned.sortedByCoord.out.bam
