#!/bin/bash
#Conduct RNA-seq read mapping by star in parallel, then count mapped reads by featurecount.
#Dependnecy: STAR v.2.7.11b, samtools v1.6, featurecount v2.0.6

GTF=/home/yo/db/PO_merged_20240517/PO_male_merged.v20240517.mod.gtf

seq 1 42|parallel --jobs 10 bash ./star_array.sh

featureCounts -p -T 48 -t exon -g gene_id -a ${GTF} -o PO_count.PO_male_merged.v20240517.mod.maskXII_XXI_Y.txt ./mapping/*.bam
