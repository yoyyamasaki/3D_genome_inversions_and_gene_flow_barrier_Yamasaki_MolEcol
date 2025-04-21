# README for inversion detection and population genomic analysis

## General information
- Author Information
  - Script writing\
    Name: Yo Yamasaki\
    Institution: National Institute of Genetics\
    Address: Mishima, Shizuoka, Japan\
    Email: [yamasaki@nig.ac.jp](mailto:yamasaki@nig.ac.jp)
- Information about funding sources that supported the collection of the data
    - JSPS Kakenhi (22H04983, 20J01503, 21H02542, and 22KK0105)
    - JST CREST (JPMJCR20S2)
    
## Overview
- This directory includes the scripts for inversion detection and population genomic analysis.
- We first conducted inversion detection using Syri. Then we conducted Fst and dxy calculation and D-statistics analysis.
- Each script has a number at the head of the script name. This is the order in which the scripts are executed.
- If you want to use these scripts for your datasets, please set an appropriate path in your environment.
- The software used and its versions are written in each script.

## Sequencing data
### Genome
- Gasterosteus nipponicus genome (JSv20240517)
- G. aculeatus genome(POv20240517)

### Short read
- G. nipponicus (Accession: DRP001192)
- G. aculeatus Japanese population (Accession: DRP001192)
- G. aculeatus Quebec population (Accession: PRJNA935671)
- G. wheatlandi (Accession: DRP001133)

## Brief explanation for scripts
- 01_syri_run_PO_JS_asm5.sh: Detect structural variation from pairwise alignment of the two genomes.
- 02_make_bed_inversion_colinear.sh: Extract inversion and colinear regions in bed format.
- 03_exec_fastp_array.sh: Filtering WGS fastq files in parallel.
- 03_fastp_array.sh: Filtering WGS fastq files in parallel.
- 04_exec_bcftools_call_array.sh: SNP call by bcftools. SNP call was conducted for each chromosome in parallel, then they were merged.
- 04_bcftools_call_array.sh: SNP call by bcftools. SNP call was conducted for each chromosome in parallel, then they were merged.
- 05_vcf_filter_bcftools.sh: Filtering of the vcf file for pixy and dsuite.
- 06_pixy_run.sh: Run pixy to calculate Fst and dxy.
- 07_plot_fst_dxy.r: Plot Fst and dxy, and conduct statistical test.
- 08_dsuite_run.sh: Conduct D statistics analysis by dsuite.

## Directory structure
```
.
├── README.md
├── SNPcall
│   ├── 04_bcftools_call_array.sh
│   ├── 04_exec_bcftools_call_array.sh
│   ├── bcftools_call
│   │   ├── 05_vcf_filter_bcftools.sh
│   │   ├── dsuite_inversion
│   │   │   ├── 08_dsuite_run.sh
│   │   │   ├── popmap_all.txt
│   │   │   └── tree_dsuite_tree_Quebec_PO_JS_WLD.tre
│   │   └── pixy_inversion
│   │       ├── 06_pixy_run.sh
│   │       └── 07_plot_fst_dxy.r
│   └── fastq
│       ├── 03_exec_fastp_array.sh
│       └── 03_fastp_array.sh
└── inversion_detection
    ├── 01_syri_run_PO_JS_asm5.sh
    └── 02_make_bed_inversion_colinear.sh
```

## Sharing/Access information

1. Licenses/restrictions placed on the data: CC0 1.0 Universal (CC0 1.0) Public Domain
2. Links to publications that cite or use the data:\
   Yo Y. Yamasaki, Atsushi Toyoda, Mitsutaka Kadota, Shigehiro Kuraku, Jun Kitano
   3D genome constrains breakpoints of inversions that can act as barriers to gene flow in the stickleback
