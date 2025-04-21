#!/bin/bash
#Conduct RNA-seq read mapping by star in parallel, then count mapped reads by featurecount.
#Dependnecy: STAR v.2.7.11b, samtools v1.6, featurecount v2.0.6

GTF=/home/yo/db/JS_male_merged_20240517/JS_male_merged.v20240517.mod.gtf

seq 1 42|parallel --jobs 10 bash ./star_array.sh

featureCounts -p -T 48 -t exon,CDS -g gene_id -a ${GTF} -o JS_count.merged.maskXII_XXI_Y.txt ./mapping/*.bam
