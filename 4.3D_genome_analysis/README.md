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
- This directory includes the scripts for 3D genome analysis using iconHi-C reads
- We first mapped iconHi-C reads. Then, we conducted AB compartment analysis and TAD analysis. Finally, we visualized them using HiGlass.
- Each script has a number at the head of the script name. This is the order in which the scripts are executed.
- If you want to use these scripts for your datasets, please set an appropriate path in your environment.
- The software used and its versions are written in each script.

## Sequencing data
### Genome
- Gasterosteus nipponicus genome (JSv20240517)
- G. aculeatus genome(POv20240517)

### iconHi-C reads
- G. nipponicus (Accession: PRJDB19958)
- G. aculeatus Japanese population (Accession:PRJDB19958)

## Brief explanation for scripts
- 01_mapping_pairtools.sh: Mapping iconHi-C reads to the reference genome. Then extract valid Hi-C read pairs.
- 02_fanc_run.sh: Calculate A/B compartment and make inputs for spectralTAD by fanc.
- 03_filter_fastq_RNA.sh: Filtering RNA-seq fastq files in parallel.
- 04_exec_star_parallel_featurecount.sh: Conduct RNA-seq read mapping by star in parallel, then count mapped reads by featurecount.
- 04_star_array.sh: Conduct RNA-seq read mapping by star in parallel.
- 05_calc_tpm.r: Conduct TPM normalization
- 06_asseign_genes_to_compartment.r: Assign compartment status to each gene.
- 07_plot_expression.r: Plot expression level for each compartment and conduct statistical test.
- 08_calc_compartment_stats.r: Calculate basic statistics of AB compartment.
- 09_spectraltads_run.r: Run SpectralTAD.
- 10_extract_unique_boundaries_from_bed.r: Caluculate unique TAD boundary number and TAD size.
- 11_higlass_forpublish.sh: Draw Hi-C heat map with TADs and inversions.

## Directory structure
```
.
├── JS
│   ├── 01_mapping_pairtools.sh
│   ├── 02_fanc_run.sh
│   ├── 11_higlass_run.sh
│   ├── RNAseq
│   │   ├── 04_exec_star_parallel_featurecount.sh
│   │   ├── 04_star_array.sh
│   │   ├── RNAseq_analysis
│   │   │   └── 05_calc_tpm.r
│   │   └── fastq
│   │       └── 03_filter_fastq_RNA.sh
│   └── iconhic_JS_refJS_merged_maskXII_XXI_Y_pairtools_KR
│       ├── compartments
│       │   └── 100kb
│       │       └── plot_expression
│       │           ├── 06_asseign_genes_to_compartment.r
│       │           ├── 07_plot_expression.r
│       │           └── 08_calc_compartment_stats.r
│       └── spectraltads
│           ├── 09_spectraltads_run.r
│           └── 10_extract_unique_boundaries_from_bed.r
├── PO
│   ├── 01_mapping_pairtools.sh
│   ├── 02_fanc_run.sh
│   ├── 11_higlass.sh
│   ├── RNAseq
│   │   ├── 04_exec_star_parallel_featurecount.sh
│   │   ├── 04_star_array.sh
│   │   ├── RNAseq_analysis
│   │   │   └── 05_calc_tpm.r
│   │   └── fastq
│   │       └── 03_filter_fastq_RNA.sh
│   └── iconic_PO_refPO_merged_maskXII_XXI_Y_pairtools_KR
│       ├── compartments
│       │   └── 100kb
│       │       └── plot_expression
│       │           ├── 06_asseign_genes_to_compartment.r
│       │           ├── 07_plot_expression.r
│       │           └── 08_calc_compartment_stats.r
│       └── spectraltad
│           ├── 09_spectraltads_run.r
│           └── 10_extract_unique_boundaries_from_bed.r
└── README.md
```

## Sharing/Access information

1. Licenses/restrictions placed on the data: CC0 1.0 Universal (CC0 1.0) Public Domain
2. Links to publications that cite or use the data:\
   Yo Y. Yamasaki, Atsushi Toyoda, Mitsutaka Kadota, Shigehiro Kuraku, Jun Kitano
   3D genome constrains breakpoints of inversions that can act as barriers to gene flow in the stickleback
