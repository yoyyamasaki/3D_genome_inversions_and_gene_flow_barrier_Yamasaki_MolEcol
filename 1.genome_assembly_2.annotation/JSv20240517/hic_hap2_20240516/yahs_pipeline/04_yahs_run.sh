#!/bin/bash
#Run yahs scaffolding.
#After finissing the scaffolding, check and correct scaffolding manually by juicebox.
#Dependency: yahs v1.2, juicer v1.2, juicer v1.9.9, PretextMap v0.1.9, PretextSnapshot v0.0.4

#Contig file path
REF=/home/yo/HDD/yo_analysis/HDD4_20230223/stickleback/pacbio/JS_male_HiFi/hifiasm_DC_OmniCtrimmed_enzymemix_20240515/hifiasm_result20240515/JS_male_DC_OmniCtrimmed.hic.hap2.p_ctg.fa
#Mapped omnic reads
BAM=/home/yo/HDD/yo_analysis/HDD4_20230223/stickleback/pacbio/JS_male_HiFi/hifiasm_DC_OmniCtrimmed_enzymemix_20240515/hic_hap2_20240516/yahs_pipeline/omnic_mapping/JS_male_hic_hap2_20240515_omnic.bam
MIN_CONTIG=1000
OUT_HEADER=JS_male_DC_OmniCtrimmed.enzymemix.hic.hap2.yahs.min${MIN_CONTIG}

mkdir -p min${MIN_CONTIG}

samtools faidx ${REF}
~/software/yahs/yahs -o ./min${MIN_CONTIG}/${OUT_HEADER} -l ${MIN_CONTIG}  ${REF} ${BAM}


#### this is to generate input file for juicer_tools - non-assembly mode or for PretextMap
## here we use 8 CPUs and 32Gb memory for sorting - adjust it according to your device
(~/software/yahs/juicer pre ./min${MIN_CONTIG}/${OUT_HEADER}.bin ./min${MIN_CONTIG}/${OUT_HEADER}_scaffolds_final.agp ${REF}.fai 2>./min${MIN_CONTIG}/tmp_juicer_pre.log | LC_ALL=C sort -k2,2d -k6,6d -T ./min${MIN_CONTIG} --parallel=8 -S32G | awk 'NF' > ./min${MIN_CONTIG}/alignments_sorted.txt.part) && (mv ./min${MIN_CONTIG}/alignments_sorted.txt.part ./min${MIN_CONTIG}/alignments_sorted.txt)

## prepare chromosome size file from samtools index file
 samtools faidx ./min${MIN_CONTIG}/${OUT_HEADER}_scaffolds_final.fa
 cut -f1-2 ./min${MIN_CONTIG}/${OUT_HEADER}_scaffolds_final.fa.fai >./min${MIN_CONTIG}/${OUT_HEADER}_scaffolds_final.chrom.sizes

## do juicer hic map
(java -Xmx32G -jar ~/software/juicer_tools.1.9.9_jcuda.0.8.jar pre ./min${MIN_CONTIG}/alignments_sorted.txt ./min${MIN_CONTIG}/${OUT_HEADER}.hic.part ./min${MIN_CONTIG}/${OUT_HEADER}_scaffolds_final.chrom.sizes) && (mv ./min${MIN_CONTIG}/${OUT_HEADER}.hic.part ./min${MIN_CONTIG}/${OUT_HEADER}.hic)

## do Pretext hic map
(awk 'BEGIN{print "## pairs format v1.0"} {print "#chromsize:\t"$1"\t"$2} END {print "#columns:\treadID\tchr1\tpos1\tchr2\tpos2\tstrand1\tstrand2"}' ./min${MIN_CONTIG}/${OUT_HEADER}_scaffolds_final.chrom.sizes; awk '{print ".\t"$2"\t"$3"\t"$6"\t"$7"\t.\t."}' ./min${MIN_CONTIG}/alignments_sorted.txt) | PretextMap -o ./min${MIN_CONTIG}/${OUT_HEADER}.pretext

# and a pretext snapshot
PretextSnapshot -m ./min${MIN_CONTIG}/${OUT_HEADER}.pretext --sequences "=full" -o ./min${MIN_CONTIG}

#### this is to generate input file for juicer_tools - assembly (JBAT) mode (-a)
~/software/yahs/juicer pre -a -o ./min${MIN_CONTIG}/${OUT_HEADER}_JBAT ./min${MIN_CONTIG}/${OUT_HEADER}.bin ./min${MIN_CONTIG}/${OUT_HEADER}_scaffolds_final.agp ${REF}.fai 2>./min${MIN_CONTIG}/tmp_juicer_pre_JBAT.log
cat ./min${MIN_CONTIG}/tmp_juicer_pre_JBAT.log | grep "PRE_C_SIZE" | cut -d' ' -f2- >./min${MIN_CONTIG}/${OUT_HEADER}_JBAT.chrom.sizes
(java -Xmx32G -jar ~/software/juicer_tools.1.9.9_jcuda.0.8.jar pre ./min${MIN_CONTIG}/${OUT_HEADER}_JBAT.txt ./min${MIN_CONTIG}/${OUT_HEADER}_JBAT.hic.part ./min${MIN_CONTIG}/${OUT_HEADER}_JBAT.chrom.sizes) && (mv ./min${MIN_CONTIG}/${OUT_HEADER}_JBAT.hic.part ./min${MIN_CONTIG}/${OUT_HEADER}_JBAT.hic)

#After the yahs scaffolding, check Hi-C heatmap and manually correct scaffolding by juicebox
