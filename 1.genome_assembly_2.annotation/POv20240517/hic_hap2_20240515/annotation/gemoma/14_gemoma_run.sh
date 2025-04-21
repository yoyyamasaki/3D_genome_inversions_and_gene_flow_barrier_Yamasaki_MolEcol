#!/bin/bash
#Run gemoma to annotate assembled genome using hiquality annotaions of other species
#Dependency: GeMoMa v1.9, MMseqs2 v.13.45111

RNA=~/HDD/yo_analysis/HDD4_20230223/stickleback/pacbio/PO_male_HiFi/RNAseq_annotation/trinity_out_dir/PO_male_hap2.v20240517.hardmasked.bam
WORKING_DIR=/raid/yo_analysis/HDD4_20230223/stickleback/pacbio/PO_male_HiFi/hifiasm_DC_OmniCtrimmed20240514/hic_hap2_20240515/annotation/gemoma
TARGET=../PO_male_hap2.v20240517.fa
THREAD=10

#Extract introns
java -jar ~/software/GeMoMa-1.9/GeMoMa-1.9.jar CLI ERE s=FR_UNSTRANDED m=${RNA} v=LENIENT c=true mmq=40 mc=1 e=0 

#Check introns
java -jar ~/software/GeMoMa-1.9/GeMoMa-1.9.jar CLI CheckIntrons t=${TARGET} i=introns.gff > intron_dist.txt

#Denoise 
java -jar ~/software/GeMoMa-1.9/GeMoMa-1.9.jar CLI DenoiseIntrons i=introns.gff c=UNSTRANDED coverage_unstranded=coverage.bedgraph m=500000

#Search positions of reference transcripts in target genome
#make mmseqs2 database for target genome
mkdir TARGET
mmseqs createdb ${TARGET} ./TARGET/GenomeDB -v 2

#####################################################################################################################################################################
##BluefinTuna
REF_SPECIES=BluefinTuna
REF_FASTA=/home/yo/db/BluefinTuna_fThuMac1.1/GCF_910596095.1/GCF_910596095.1_fThuMac1.1_genomic.fna
REF_GFF=/home/yo/db/BluefinTuna_fThuMac1.1/GCF_910596095.1/genomic.gff
OUTDIR=${WORKING_DIR}/${REF_SPECIES}

mkdir -p ${OUTDIR}
java -jar ~/software/GeMoMa-1.9/GeMoMa-1.9.jar CLI Extractor a=${REF_GFF} g=${REF_FASTA} outdir=${OUTDIR}

##Run mmseqs2
mkdir ${OUTDIR}/QUERY
mmseqs createdb ${OUTDIR}/cds-parts.fasta ${OUTDIR}/QUERY/cdsDB -v 2

mmseqs search ${OUTDIR}/QUERY/cdsDB ./TARGET/GenomeDB ${OUTDIR}/${REF_SPECIES}_mmseqs.out ${OUTDIR}/mmseqs_tmp -e 100.0 --threads ${THREAD} -s 8.5 -a --comp-bias-corr 0 --max-seqs 500 --mask 0 --orf-start-mode 1 -v 2

mmseqs convertalis ${OUTDIR}/QUERY/cdsDB ./TARGET/GenomeDB ${OUTDIR}/${REF_SPECIES}_mmseqs.out ${OUTDIR}/${REF_SPECIES}_search.txt --threads ${THREAD} --format-output "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,empty,raw,nident,empty,empty,empty,qframe,tframe,qaln,taln,qlen,tlen" -v 2

##Run GeMoMa
java -Xmx50G -jar ~/software/GeMoMa-1.9/GeMoMa-1.9.jar CLI GeMoMa s=${OUTDIR}/${REF_SPECIES}_search.txt \
  c=${OUTDIR}/cds-parts.fasta a=${OUTDIR}/assignment.tabular \
  t=${TARGET} sort=true Score=ReAlign i=denoised_introns.gff \
  coverage=UNSTRANDED coverage_unstranded=coverage.bedgraph \
  outdir=${OUTDIR}

#####################################################################################################################################################################
##Clownfish
REF_SPECIES=Clownfish
REF_FASTA=~/db/Clownfish_ASM2253959v1/GCF_022539595.1/GCF_022539595.1_ASM2253959v1_genomic.fna
REF_GFF=~/db/Clownfish_ASM2253959v1/GCF_022539595.1/genomic.gff
OUTDIR=${WORKING_DIR}/${REF_SPECIES}

mkdir -p ${OUTDIR}
java -jar ~/software/GeMoMa-1.9/GeMoMa-1.9.jar CLI Extractor a=${REF_GFF} g=${REF_FASTA} outdir=${OUTDIR}

##Run mmseqs2
mkdir ${OUTDIR}/QUERY
mmseqs createdb ${OUTDIR}/cds-parts.fasta ${OUTDIR}/QUERY/cdsDB -v 2

mmseqs search ${OUTDIR}/QUERY/cdsDB ./TARGET/GenomeDB ${OUTDIR}/mmseqs.out ${OUTDIR}/mmseqs_tmp -e 100.0 --threads ${THREAD} -s 8.5 -a --comp-bias-corr 0 --max-seqs 500 --mask 0 --orf-start-mode 1 -v 2

mmseqs convertalis ${OUTDIR}/QUERY/cdsDB ./TARGET/GenomeDB ${OUTDIR}/mmseqs.out ${OUTDIR}/${REF_SPECIES}_search.txt --threads ${THREAD} --format-output "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,empty,raw,nident,empty,empty,empty,qframe,tframe,qaln,taln,qlen,tlen" -v 2

##Run GeMoMa
java -Xmx50G -jar ~/software/GeMoMa-1.9/GeMoMa-1.9.jar CLI GeMoMa s=${OUTDIR}/${REF_SPECIES}_search.txt \
  c=${OUTDIR}/cds-parts.fasta a=${OUTDIR}/assignment.tabular \
  t=${TARGET} sort=true Score=ReAlign i=denoised_introns.gff \
  coverage=UNSTRANDED coverage_unstranded=coverage.bedgraph \
  outdir=${OUTDIR}

#####################################################################################################################################################################
##Medaka
REF_SPECIES=Medaka
REF_FASTA=~/db/medaka_ASM223467v1/GCF_002234675.1/GCF_002234675.1_ASM223467v1_genomic.fna
REF_GFF=~/db/medaka_ASM223467v1/GCF_002234675.1/genomic.gff
OUTDIR=${WORKING_DIR}/${REF_SPECIES}

mkdir -p ${OUTDIR}
java -jar ~/software/GeMoMa-1.9/GeMoMa-1.9.jar CLI Extractor a=${REF_GFF} g=${REF_FASTA} outdir=${OUTDIR}

##Run mmseqs2
mkdir ${OUTDIR}/QUERY
mmseqs createdb ${OUTDIR}/cds-parts.fasta ${OUTDIR}/QUERY/cdsDB -v 2

mmseqs search ${OUTDIR}/QUERY/cdsDB ./TARGET/GenomeDB ${OUTDIR}/mmseqs.out ${OUTDIR}/mmseqs_tmp -e 100.0 --threads ${THREAD} -s 8.5 -a --comp-bias-corr 0 --max-seqs 500 --mask 0 --orf-start-mode 1 -v 2

mmseqs convertalis ${OUTDIR}/QUERY/cdsDB ./TARGET/GenomeDB ${OUTDIR}/mmseqs.out ${OUTDIR}/${REF_SPECIES}_search.txt --threads ${THREAD} --format-output "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,empty,raw,nident,empty,empty,empty,qframe,tframe,qaln,taln,qlen,tlen" -v 2

##Run GeMoMa
java -Xmx50G -jar ~/software/GeMoMa-1.9/GeMoMa-1.9.jar CLI GeMoMa s=${OUTDIR}/${REF_SPECIES}_search.txt \
  c=${OUTDIR}/cds-parts.fasta a=${OUTDIR}/assignment.tabular \
  t=${TARGET} sort=true Score=ReAlign i=denoised_introns.gff \
  coverage=UNSTRANDED coverage_unstranded=coverage.bedgraph \
  outdir=${OUTDIR}
  
#####################################################################################################################################################################
##Stickleback
REF_SPECIES=Stickleback
REF_FASTA=/home/yo/db/stickleback_v5/ncbi_dataset/ncbi_dataset/data/GCF_016920845.1/GCF_016920845.1_GAculeatus_UGA_version5_genomic.fna
REF_GFF=/home/yo/db/stickleback_v5/ncbi_dataset/ncbi_dataset/data/GCF_016920845.1/genomic.gff
OUTDIR=${WORKING_DIR}/${REF_SPECIES}

mkdir -p ${OUTDIR}
java -jar ~/software/GeMoMa-1.9/GeMoMa-1.9.jar CLI Extractor a=${REF_GFF} g=${REF_FASTA} outdir=${OUTDIR}

##Run mmseqs2
mkdir ${OUTDIR}/QUERY
mmseqs createdb ${OUTDIR}/cds-parts.fasta ${OUTDIR}/QUERY/cdsDB -v 2

mmseqs search ${OUTDIR}/QUERY/cdsDB ./TARGET/GenomeDB ${OUTDIR}/mmseqs.out ${OUTDIR}/mmseqs_tmp -e 100.0 --threads ${THREAD} -s 8.5 -a --comp-bias-corr 0 --max-seqs 500 --mask 0 --orf-start-mode 1 -v 2

mmseqs convertalis ${OUTDIR}/QUERY/cdsDB ./TARGET/GenomeDB ${OUTDIR}/mmseqs.out ${OUTDIR}/${REF_SPECIES}_search.txt --threads ${THREAD} --format-output "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,empty,raw,nident,empty,empty,empty,qframe,tframe,qaln,taln,qlen,tlen" -v 2

##Run GeMoMa
java -Xmx50G -jar ~/software/GeMoMa-1.9/GeMoMa-1.9.jar CLI GeMoMa s=${OUTDIR}/${REF_SPECIES}_search.txt \
  c=${OUTDIR}/cds-parts.fasta a=${OUTDIR}/assignment.tabular \
  t=${TARGET} sort=true Score=ReAlign i=denoised_introns.gff \
  coverage=UNSTRANDED coverage_unstranded=coverage.bedgraph \
  outdir=${OUTDIR}

#####################################################################################################################################################################
##Tilapia
REF_SPECIES=Tilapia
REF_FASTA=~/db/tilapia_O_niloticus_UMD_NMBU/GCF_001858045.2/GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fna
REF_GFF=~/db/tilapia_O_niloticus_UMD_NMBU/GCF_001858045.2/genomic.gff
OUTDIR=${WORKING_DIR}/${REF_SPECIES}

mkdir -p ${OUTDIR}
java -jar ~/software/GeMoMa-1.9/GeMoMa-1.9.jar CLI Extractor a=${REF_GFF} g=${REF_FASTA} outdir=${OUTDIR}

##Run mmseqs2
mkdir ${OUTDIR}/QUERY
mmseqs createdb ${OUTDIR}/cds-parts.fasta ${OUTDIR}/QUERY/cdsDB -v 2

mmseqs search ${OUTDIR}/QUERY/cdsDB ./TARGET/GenomeDB ${OUTDIR}/mmseqs.out ${OUTDIR}/mmseqs_tmp -e 100.0 --threads ${THREAD} -s 8.5 -a --comp-bias-corr 0 --max-seqs 500 --mask 0 --orf-start-mode 1 -v 2

mmseqs convertalis ${OUTDIR}/QUERY/cdsDB ./TARGET/GenomeDB ${OUTDIR}/mmseqs.out ${OUTDIR}/${REF_SPECIES}_search.txt --threads ${THREAD} --format-output "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,empty,raw,nident,empty,empty,empty,qframe,tframe,qaln,taln,qlen,tlen" -v 2

##Run GeMoMa
java -Xmx50G -jar ~/software/GeMoMa-1.9/GeMoMa-1.9.jar CLI GeMoMa s=${OUTDIR}/${REF_SPECIES}_search.txt \
  c=${OUTDIR}/cds-parts.fasta a=${OUTDIR}/assignment.tabular \
  t=${TARGET} sort=true Score=ReAlign i=denoised_introns.gff \
  coverage=UNSTRANDED coverage_unstranded=coverage.bedgraph \
  outdir=${OUTDIR}

#####################################################################################################################################################################
#Aggregate GeMoMa results
java -jar ~/software/GeMoMa-1.9/GeMoMa-1.9.jar CLI GAF \
p=BluefinTuna_fThuMac11 g=BluefinTuna/predicted_annotation.gff \
p=Clownfish_ASM2253959v1 g=Clownfish/predicted_annotation.gff \
p=medaka_ASM223467v1 g=Medaka/predicted_annotation.gff \
p=stickleback_v5 g=Stickleback/predicted_annotation.gff \
p=tilapia_O_niloticus_UMD_NMBU g=Tilapia/predicted_annotation.gff

#Add UTR information
java -jar ~/software/GeMoMa-1.9/GeMoMa-1.9.jar CLI AnnotationFinalizer g=${TARGET} a=filtered_predictions.gff u=YES i=denoised_introns.gff c=UNSTRANDED coverage_unstranded=coverage.bedgraph rename=NO






