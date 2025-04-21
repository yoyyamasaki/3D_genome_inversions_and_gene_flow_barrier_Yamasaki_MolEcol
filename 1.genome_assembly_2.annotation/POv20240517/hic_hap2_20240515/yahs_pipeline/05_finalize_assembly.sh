#!/bin/bash
#Apply the result of manual correction the juicebox to the fasta file.
#Dependency: juicer v1.2

#Contig file path
REF=/home/yo/HDD/yo_analysis/HDD4_20230223/stickleback/pacbio/PO_male_HiFi/hifiasm_DC_OmniCtrimmed20240514/hifiasm_result20240515/PO_male_DC_OmniCtrimmed.hic.hap2.p_ctg.fa
MIN_CONTIG=1000
OUT_HEADER=PO_male_DC_OmniCtrimmed.hic.hap2.yahs.min${MIN_CONTIG}

#Generate final assembly modified by juicebox
~/software/yahs/juicer post -o ./min${MIN_CONTIG}/${OUT_HEADER}_JBAT ./min${MIN_CONTIG}/${OUT_HEADER}_JBAT.review.assembly ./min${MIN_CONTIG}/${OUT_HEADER}_JBAT.liftover.agp ${REF}
