#!/bin/bash
#Select reliable annotations

ENTAP_OUT=./entap_outfiles/final_results/entap_results.tsv 
INGFF=../gemoma/final_annotation_1.gff
OUT_TRANSCRIPT_LIST=/home/yo/HDD/yo_analysis/HDD4_20230223/stickleback/pacbio/JS_male_HiFi/hifiasm_DC_OmniCtrimmed_enzymemix_20240515/hic_hap2_20240516/annotation/entap_gemoma_only/entap_outfiles/final_results/filtered_transcript_list.txt
OUT_GENE_LIST=/home/yo/HDD/yo_analysis/HDD4_20230223/stickleback/pacbio/JS_male_HiFi/hifiasm_DC_OmniCtrimmed_enzymemix_20240515/hic_hap2_20240516/annotation/entap_gemoma_only/entap_outfiles/final_results/filtered_gene_list.txt
OUT_GFF=JS_male_hap2.v20240517.filtered.20220612.gff

#Select genes which have annotations of "Fishes", "Aves", "Mammals", "Animals", or "Eukaryotes".
cat ${ENTAP_OUT} | grep -e "Fishes" -e "Aves" -e "Mammals" -e "Animals" -e "Eukaryotes" | cut -f 1 > ${OUT_TRANSCRIPT_LIST}

#Extract gene name
cat ${INGFF} | grep -f ${OUT_TRANSCRIPT_LIST} |awk '$3=="mRNA" {print $9}'|awk -F ";" '{print $1}'|perl -pe 's/Name=//g'|perl -pe 's/\.[1234567890]//g'|uniq > ${OUT_GENE_LIST}

#Make annotation entry list to be extracted
cat ${OUT_TRANSCRIPT_LIST} ${OUT_GENE_LIST} > grep_list.txt

#Extract annotations
cat ${INGFF} | grep -f grep_list.txt -e "#" > ${OUT_GFF}
