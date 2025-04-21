#! /bin/bash
#Obtain RNA-seq mapping data for genome annotation.
#RNA-seq reads were cleaned by mapping for assembled transcriptome sequence. Only correctly mapped reads were used for the annotation.
#Dependency: hisat2 v2.2.1

#Map RNA-seq data to trinity assembled RNA-seq contigs. Then, extract uniquely mapped reads. 
hisat2-build -p 8 trinity_out_dir.Trinity.fasta INDEX

hisat2 -x INDEX -1 ../filtered/brain_PP1_m_trimmed_1P.fq.gz,../filtered/brain_PP2_f_trimmed_1P.fq.gz,../filtered/brain_PP3_m_trimmed_1P.fq.gz,../filtered/brain_PP4_f_trimmed_1P.fq.gz,../filtered/brain_PP5_m_trimmed_1P.fq.gz,../filtered/brain_PP6_f_trimmed_1P.fq.gz,../filtered/brain_PP7_m_trimmed_1P.fq.gz,../filtered/gill_PP1_m_trimmed_1P.fq.gz,../filtered/gill_PP2_f_trimmed_1P.fq.gz,../filtered/gill_PP3_m_trimmed_1P.fq.gz,../filtered/gill_PP4_f_trimmed_1P.fq.gz,../filtered/gill_PP5_m_trimmed_1P.fq.gz,../filtered/gill_PP6_f_trimmed_1P.fq.gz,../filtered/gill_PP7_m_trimmed_1P.fq.gz,../filtered/gonad_PP1_m_trimmed_1P.fq.gz,../filtered/gonad_PP2_f_trimmed_1P.fq.gz,../filtered/gonad_PP3_m_trimmed_1P.fq.gz,../filtered/gonad_PP4_f_trimmed_1P.fq.gz,../filtered/gonad_PP5_m_trimmed_1P.fq.gz,../filtered/gonad_PP6_f_trimmed_1P.fq.gz,../filtered/gonad_PP7_m_trimmed_1P.fq.gz,../filtered/liver_PP1_m_trimmed_1P.fq.gz,../filtered/liver_PP2_f_trimmed_1P.fq.gz,../filtered/liver_PP3_m_trimmed_1P.fq.gz,../filtered/liver_PP4_f_trimmed_1P.fq.gz,../filtered/liver_PP5_m_trimmed_1P.fq.gz,../filtered/liver_PP6_f_trimmed_1P.fq.gz,../filtered/liver_PP7_m_trimmed_1P.fq.gz,../filtered/pec_musAD_PP1_m_trimmed_1P.fq.gz,../filtered/pec_musAD_PP2_f_trimmed_1P.fq.gz,../filtered/pec_musAD_PP3_m_trimmed_1P.fq.gz,../filtered/pec_musAD_PP4_f_trimmed_1P.fq.gz,../filtered/pec_musAD_PP5_m_trimmed_1P.fq.gz,../filtered/pec_musAD_PP6_f_trimmed_1P.fq.gz,../filtered/pec_musAD_PP7_m_trimmed_1P.fq.gz,../filtered/pit_PP1_m_trimmed_1P.fq.gz,../filtered/pit_PP2_f_trimmed_1P.fq.gz,../filtered/pit_PP3_m_trimmed_1P.fq.gz,../filtered/pit_PP4_f_trimmed_1P.fq.gz,../filtered/pit_PP5_trimmed_1P.fq.gz,../filtered/pit_PP6_f_trimmed_1P.fq.gz,../filtered/pit_PP7_m_trimmed_1P.fq.gz \
-2 ../filtered/brain_PP1_m_trimmed_2P.fq.gz,../filtered/brain_PP2_f_trimmed_2P.fq.gz,../filtered/brain_PP3_m_trimmed_2P.fq.gz,../filtered/brain_PP4_f_trimmed_2P.fq.gz,../filtered/brain_PP5_m_trimmed_2P.fq.gz,../filtered/brain_PP6_f_trimmed_2P.fq.gz,../filtered/brain_PP7_m_trimmed_2P.fq.gz,../filtered/gill_PP1_m_trimmed_2P.fq.gz,../filtered/gill_PP2_f_trimmed_2P.fq.gz,../filtered/gill_PP3_m_trimmed_2P.fq.gz,../filtered/gill_PP4_f_trimmed_2P.fq.gz,../filtered/gill_PP5_m_trimmed_2P.fq.gz,../filtered/gill_PP6_f_trimmed_2P.fq.gz,../filtered/gill_PP7_m_trimmed_2P.fq.gz,../filtered/gonad_PP1_m_trimmed_2P.fq.gz,../filtered/gonad_PP2_f_trimmed_2P.fq.gz,../filtered/gonad_PP3_m_trimmed_2P.fq.gz,../filtered/gonad_PP4_f_trimmed_2P.fq.gz,../filtered/gonad_PP5_m_trimmed_2P.fq.gz,../filtered/gonad_PP6_f_trimmed_2P.fq.gz,../filtered/gonad_PP7_m_trimmed_2P.fq.gz,../filtered/liver_PP1_m_trimmed_2P.fq.gz,../filtered/liver_PP2_f_trimmed_2P.fq.gz,../filtered/liver_PP3_m_trimmed_2P.fq.gz,../filtered/liver_PP4_f_trimmed_2P.fq.gz,../filtered/liver_PP5_m_trimmed_2P.fq.gz,../filtered/liver_PP6_f_trimmed_2P.fq.gz,../filtered/liver_PP7_m_trimmed_2P.fq.gz,../filtered/pec_musAD_PP1_m_trimmed_2P.fq.gz,../filtered/pec_musAD_PP2_f_trimmed_2P.fq.gz,../filtered/pec_musAD_PP3_m_trimmed_2P.fq.gz,../filtered/pec_musAD_PP4_f_trimmed_2P.fq.gz,../filtered/pec_musAD_PP5_m_trimmed_2P.fq.gz,../filtered/pec_musAD_PP6_f_trimmed_2P.fq.gz,../filtered/pec_musAD_PP7_m_trimmed_2P.fq.gz,../filtered/pit_PP1_m_trimmed_2P.fq.gz,../filtered/pit_PP2_f_trimmed_2P.fq.gz,../filtered/pit_PP3_m_trimmed_2P.fq.gz,../filtered/pit_PP4_f_trimmed_2P.fq.gz,../filtered/pit_PP5_trimmed_2P.fq.gz,../filtered/pit_PP6_f_trimmed_2P.fq.gz,../filtered/pit_PP7_m_trimmed_2P.fq.gz -p 24 | \
samtools view -f 0x2 -bh | \
samtools sort -@24 -o Trinity.bam && samtools index -@24 Trinity.bam

#Convert bam to fastq
samtools bam2fq -n -1 CLEAN_1.fq.gz -2 CLEAN_2.fq.gz Trinity.bam

#Map clean reads to hardmasked hap1 genome
ln -s ~/HDD/yo_analysis/HDD4_20230223/stickleback/pacbio/JS_male_HiFi/hifiasm_DC_OmniCtrimmed_enzymemix_20240515/hic_hap1_20240516/annotation/repeatmask/JS_male_hap1.v20240517.hardmasked.fa ./
hisat2-build -p 8 JS_male_hap1.v20240517.hardmasked.fa JS_male_hap1.v20240517.hardmasked
hisat2 -x JS_male_hap1.v20240517.hardmasked -1 CLEAN_1.fq.gz -2 CLEAN_2.fq.gz -p 60 --max-intronlen 100000 | samtools sort -@60 -o JS_male_hap1.v20240517.hardmasked.bam && samtools index -@60 JS_male_hap1.v20240517.hardmasked.bam

#Map clean reads to hardmasked hap2 genome
ln -s ~/HDD/yo_analysis/HDD4_20230223/stickleback/pacbio/JS_male_HiFi/hifiasm_DC_OmniCtrimmed_enzymemix_20240515/hic_hap2_20240516/annotation/repeatmask/JS_male_hap2.v20240517.hardmasked.fa ./
hisat2-build -p 8 JS_male_hap2.v20240517.hardmasked.fa JS_male_hap2.v20240517.hardmasked
hisat2 -x JS_male_hap2.v20240517.hardmasked -1 CLEAN_1.fq.gz -2 CLEAN_2.fq.gz -p 40 --max-intronlen 100000 | samtools sort -@40 -o JS_male_hap2.v20240517.hardmasked.bam && samtools index -@40 JS_male_hap2.v20240517.hardmasked.bam

