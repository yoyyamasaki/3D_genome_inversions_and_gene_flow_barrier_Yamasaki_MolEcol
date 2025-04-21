#!/bin/bash
#Tandem repeat analysis by TRASH
#Dependency: TRASH v1.2

FASTA_DIR=/home/yo/HDD/yo_analysis/HDD4_20230223/stickleback/pacbio/PO_male_HiFi/hifiasm_DC_OmniCtrimmed20240514/hic_hap2_20240515/yahs_pipeline/min1000/renamed/chr
OUTDIR=/home/yo/HDD/yo_analysis/HDD4_20230223/stickleback/pacbio/PO_male_HiFi/hifiasm_DC_OmniCtrimmed20240514/hic_hap2_20240515/yahs_pipeline/min1000/renamed/trash_results

fasta_list=(xxx $(ls ${FASTA_DIR}/*.fa))
fasta_file=${fasta_list[$1]}

mkdir -p ${OUTDIR}
~/software/TRASH/TRASH_run.sh ${fasta_file} --o ${OUTDIR}/
