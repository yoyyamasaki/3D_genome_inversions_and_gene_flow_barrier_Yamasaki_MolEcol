#! /bin/bash
#Obtain RNA-seq mapping data for genome annotation.
#RNA-seq reads were cleaned by mapping for assembled transcriptome sequence. Only correctly mapped reads were used for the annotation.
#Dependency: hisat2 v2.2.1

#Map RNA-seq data to trinity assembled RNA-seq contigs. Then, extract uniquely mapped reads. 
hisat2-build -p 8 trinity_out_dir.Trinity.fasta INDEX

hisat2 -x INDEX -1 ../filtered/brain_JJ1_m_trimmed_1P.fq.gz,../filtered/brain_JJ2_f_trimmed_1P.fq.gz,../filtered/brain_JJ3_m_trimmed_1P.fq.gz,../filtered/brain_JJ4_f_trimmed_1P.fq.gz,../filtered/brain_JJ5_m_trimmed_1P.fq.gz,../filtered/brain_JJ6_f_trimmed_1P.fq.gz,../filtered/brain_JJ7_m_trimmed_1P.fq.gz,../filtered/gill_JJ1_m_trimmed_1P.fq.gz,../filtered/gill_JJ2_f_trimmed_1P.fq.gz,../filtered/gill_JJ3_m_trimmed_1P.fq.gz,../filtered/gill_JJ4_f_trimmed_1P.fq.gz,../filtered/gill_JJ5_m_trimmed_1P.fq.gz,../filtered/gill_JJ6_f_trimmed_1P.fq.gz,../filtered/gill_JJ7_m_trimmed_1P.fq.gz,../filtered/gonad_JJ1_m_trimmed_1P.fq.gz,../filtered/gonad_JJ2_f_trimmed_1P.fq.gz,../filtered/gonad_JJ3_m_trimmed_1P.fq.gz,../filtered/gonad_JJ4_f_trimmed_1P.fq.gz,../filtered/gonad_JJ5_m_trimmed_1P.fq.gz,../filtered/gonad_JJ6_f_trimmed_1P.fq.gz,../filtered/gonad_JJ7_m_trimmed_1P.fq.gz,../filtered/liver_JJ1_m_trimmed_1P.fq.gz,../filtered/liver_JJ2_f_trimmed_1P.fq.gz,../filtered/liver_JJ3_m_trimmed_1P.fq.gz,../filtered/liver_JJ4_f_trimmed_1P.fq.gz,../filtered/liver_JJ5_m_trimmed_1P.fq.gz,../filtered/liver_JJ6_f_trimmed_1P.fq.gz,../filtered/liver_JJ7_m_trimmed_1P.fq.gz,../filtered/pec_musAD_JJ2_f_trimmed_1P.fq.gz,../filtered/pec_musAD_JJ3_m_trimmed_1P.fq.gz,../filtered/pec_musAD_JJ4_f_trimmed_1P.fq.gz,../filtered/pec_musAD_JJ5_m_trimmed_1P.fq.gz,../filtered/pec_musAD_JJ6_f_trimmed_1P.fq.gz,../filtered/pec_musAD_JJ7_m_trimmed_1P.fq.gz,../filtered/pec_musAD__trimmed_1P.fq.gz,../filtered/pit_JJ1_m_trimmed_1P.fq.gz,../filtered/pit_JJ2_f_trimmed_1P.fq.gz,../filtered/pit_JJ3_m_trimmed_1P.fq.gz,../filtered/pit_JJ4_f_trimmed_1P.fq.gz,../filtered/pit_JJ5_m_trimmed_1P.fq.gz,../filtered/pit_JJ6_f_trimmed_1P.fq.gz,../filtered/pit_JJ7_m_trimmed_1P.fq.gz \
-2 ../filtered/brain_JJ1_m_trimmed_2P.fq.gz,../filtered/brain_JJ2_f_trimmed_2P.fq.gz,../filtered/brain_JJ3_m_trimmed_2P.fq.gz,../filtered/brain_JJ4_f_trimmed_2P.fq.gz,../filtered/brain_JJ5_m_trimmed_2P.fq.gz,../filtered/brain_JJ6_f_trimmed_2P.fq.gz,../filtered/brain_JJ7_m_trimmed_2P.fq.gz,../filtered/gill_JJ1_m_trimmed_2P.fq.gz,../filtered/gill_JJ2_f_trimmed_2P.fq.gz,../filtered/gill_JJ3_m_trimmed_2P.fq.gz,../filtered/gill_JJ4_f_trimmed_2P.fq.gz,../filtered/gill_JJ5_m_trimmed_2P.fq.gz,../filtered/gill_JJ6_f_trimmed_2P.fq.gz,../filtered/gill_JJ7_m_trimmed_2P.fq.gz,../filtered/gonad_JJ1_m_trimmed_2P.fq.gz,../filtered/gonad_JJ2_f_trimmed_2P.fq.gz,../filtered/gonad_JJ3_m_trimmed_2P.fq.gz,../filtered/gonad_JJ4_f_trimmed_2P.fq.gz,../filtered/gonad_JJ5_m_trimmed_2P.fq.gz,../filtered/gonad_JJ6_f_trimmed_2P.fq.gz,../filtered/gonad_JJ7_m_trimmed_2P.fq.gz,../filtered/liver_JJ1_m_trimmed_2P.fq.gz,../filtered/liver_JJ2_f_trimmed_2P.fq.gz,../filtered/liver_JJ3_m_trimmed_2P.fq.gz,../filtered/liver_JJ4_f_trimmed_2P.fq.gz,../filtered/liver_JJ5_m_trimmed_2P.fq.gz,../filtered/liver_JJ6_f_trimmed_2P.fq.gz,../filtered/liver_JJ7_m_trimmed_2P.fq.gz,../filtered/pec_musAD_JJ2_f_trimmed_2P.fq.gz,../filtered/pec_musAD_JJ3_m_trimmed_2P.fq.gz,../filtered/pec_musAD_JJ4_f_trimmed_2P.fq.gz,../filtered/pec_musAD_JJ5_m_trimmed_2P.fq.gz,../filtered/pec_musAD_JJ6_f_trimmed_2P.fq.gz,../filtered/pec_musAD_JJ7_m_trimmed_2P.fq.gz,../filtered/pec_musAD__trimmed_2P.fq.gz,../filtered/pit_JJ1_m_trimmed_2P.fq.gz,../filtered/pit_JJ2_f_trimmed_2P.fq.gz,../filtered/pit_JJ3_m_trimmed_2P.fq.gz,../filtered/pit_JJ4_f_trimmed_2P.fq.gz,../filtered/pit_JJ5_m_trimmed_2P.fq.gz,../filtered/pit_JJ6_f_trimmed_2P.fq.gz,../filtered/pit_JJ7_m_trimmed_2P.fq.gz -p 24 | \
samtools view -f 0x2 -bh | \
samtools sort -@24 -o Trinity.bam && samtools index -@24 Trinity.bam

#Convert bam to fastq
samtools bam2fq -n -1 CLEAN_1.fq.gz -2 CLEAN_2.fq.gz Trinity.bam

#Map clean reads to hardmasked hap1 genome
ln -s ln -s ~/HDD/yo_analysis/HDD4_20230223/stickleback/pacbio/PO_male_HiFi/hifiasm_DC_OmniCtrimmed20240514/hic_hap1_20240515/annotation/repeatmask/PO_male_hap1.v20240517.hardmasked.fa ./
hisat2-build -p 8 PO_male_hap1.v20240517.hardmasked.fa PO_male_hap1.v20240517.hardmasked
hisat2 -x PO_male_hap1.v20240517.hardmasked -1 CLEAN_1.fq.gz -2 CLEAN_2.fq.gz -p 60 --max-intronlen 100000 | samtools sort -@60 -o PO_male_hap1.v20240517.hardmasked.bam && samtools index -@60 PO_male_hap1.v20240517.hardmasked.bam

#Map clean reads to hardmasked hap2 genome
ln -s ~/HDD/yo_analysis/HDD4_20230223/stickleback/pacbio/PO_male_HiFi/hifiasm_DC_OmniCtrimmed20240514/hic_hap2_20240515/annotation/repeatmask/PO_male_hap2.v20240517.hardmasked.fa ./
hisat2-build -p 8 PO_male_hap2.v20240517.hardmasked.fa PO_male_hap2.v20240517.hardmasked
hisat2 -x PO_male_hap2.v20240517.hardmasked -1 CLEAN_1.fq.gz -2 CLEAN_2.fq.gz -p 40 --max-intronlen 100000 | samtools sort -@40 -o PO_male_hap2.v20240517.hardmasked.bam && samtools index -@40 PO_male_hap2.v20240517.hardmasked.bam

