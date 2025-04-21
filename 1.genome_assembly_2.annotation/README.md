# README for genome assembly and annotation

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
- This directory includes the scripts for genome assembly and its annotation.
- We first conducted haplotype-phased assembly for Gasterosteus nipponicus (JS) and G. aculeatus (PO), then obtained two haplotypes for each species.
- Omni-C scaffolding and annotation were conducted for each haplotype (in total, four haplotypes).
- Each script has a number at the head of the script name. This is the order in which the scripts are executed.
- If you want to use these scripts for your datasets, please set an appropriate path in your environment.
- The software used and its versions are written in each script.

## Sequencing data
### Gasterosteus nipponicus
- Pacbio HiFi reads (Accession: DRR629085)
- Omni-C reads (Accession: DRR629086, DRR629087)
- RNA-seq (Accession: DRA007146, use G. nipponicus data only)

### Gasterosteus aculeatus
- Pacbio HiFi reads (Accession: DRR631041)
- Omni-C reads (Accession: DRR631042, DRR631043)
- RNA-seq (Accession: DRA007146, use G. aculeatus data only)

## Brief explanation for scripts
- 01_fastp.sh: Clean Omnic reads.
- 02_hifiasm_run.sh: Conduct haplotype phased assembly by hifiasm.
- 03_mapping_omnic.sh: Mapping omnic read to the contigs for the following yahs scaffolding.
- 04_yahs_run.sh: Run yahs scaffolding.
- 05_finalize_assembly.sh: Apply the result of manual correction the juicebox to the fasta file.
- 06_convert_name_directions.sh: Convert scaffold name to chromosome number, and convert direction of chromosome to that of reference.
- 07_separate_chr.sh: Write respective chromosomes in different files. Also count gap number.
- 08_exec_TRASH_array.sh: Tandem repeat analysis by TRASH.
- 08_TRASH_array.sh: Tandem repeat analysis by TRASH.
- 09_repeatmodeler_repeatmasker_docker.sh: Run repeatmodeler and repeatmasker.
- 10_braker_run.sh: Run braker2 annotation pipeline.
- 11_exec_fastp_array.sh: Run fastp to clean RNA-seq reads.
- 11_fastp_array.sh: Run fastp to clean RNA-seq reads.
- 12_trinity_singularity.sh: Run trinity to assemble transcriptome.
- 13_hisat2_mapping_JS.sh: Obtain RNA-seq mapping data for genome annotation.
- 14_gemoma_run.sh: Run gemoma to annotate assembled genome using hiquality annotaions of other species.
- 15_gffcompare_run.sh: Compare gemoma and braker2 results, and add blaker2 unique result to gemoma result.
- 16_entap_run.sh: Conduct functional annotation by entap.
- 17_select_annotated_genes.sh: Select reliable annotations.

## Directory structure
```
.
├── JSv20240517
│   ├── 02_hifiasm_run.sh
│   ├── Omnic
│   │   └── 01_fastp.sh
│   ├── RNA-seq
│   │   ├── 11_exec_fastp_array.sh
│   │   ├── 11_fastp_array.sh
│   │   ├── 12_trinity_singularity.sh
│   │   └── 13_hisat2_mapping_JS.sh
│   ├── hic_hap1_20240516
│   │   ├── annotation
│   │   │   ├── braker_prot
│   │   │   │   └── 10_braker_run.sh
│   │   │   ├── entap
│   │   │   │   ├── 16_entap_run.sh
│   │   │   │   └── 17_select_annotated_genes.sh
│   │   │   ├── gemoma
│   │   │   │   ├── 14_gemoma_run.sh
│   │   │   │   └── 15_gffcompare_run.sh
│   │   │   └── repeatmask
│   │   │       └── 09_repeatmodeler_repeatmasker_docker.sh
│   │   └── yahs_pipeline
│   │       ├── 04_yahs_run.sh
│   │       ├── 05_finalize_assembly.sh
│   │       ├── min1000
│   │       │   ├── 06_convert_name_directions.sh
│   │       │   └── renamed
│   │       │       ├── 07_separate_chr.sh
│   │       │       ├── 08_TRASH_array.sh
│   │       │       └── 08_exec_TRASH_array.sh
│   │       └── omnic_mapping
│   │           └── 03_mapping_omnic.sh
│   └── hic_hap2_20240516
│       ├── annotation
│       │   ├── entap_gemoma_only
│       │   │   ├── 16_entap_run.sh
│       │   │   └── 17_select_annotated_genes.sh
│       │   ├── gemoma
│       │   │   ├── 14_gemoma_run.sh
│       │   │   └── 15_gffcompare_run.sh
│       │   └── repeatmask
│       │       └── 09_repeatmodeler_repeatmasker_docker.sh
│       └── yahs_pipeline
│           ├── 04_yahs_run.sh
│           ├── 05_finalize_assembly.sh
│           ├── min1000
│           │   ├── 06_convert_name_directions.sh
│           │   └── renamed
│           │       ├── 07_separate_chr.sh
│           │       ├── 08_TRASH_array.sh
│           │       └── 08_exec_TRASH_array.sh
│           └── omnic_mapping
│               └── 03_mapping_omnic.sh
├── POv20240517
│   ├── 02_hifiasm_run.sh
│   ├── Omnic
│   │   └── 01_fastp.sh
│   ├── RNA-seq
│   │   ├── 11_exec_fastp_array.sh
│   │   ├── 11_fastp_array.sh
│   │   ├── 12_trinity_singularity.sh
│   │   └── 13_hisat2_mapping_PO.sh
│   ├── hic_hap1_20240515
│   │   ├── annotation
│   │   │   ├── braker_prot
│   │   │   │   └── 10_braker_run.sh
│   │   │   ├── entap
│   │   │   │   ├── 16_entap_run.sh
│   │   │   │   └── 17_select_annotated_genes.sh
│   │   │   ├── gemoma
│   │   │   │   ├── 14_gemoma_run.sh
│   │   │   │   └── 15_gffcompare_run.sh
│   │   │   └── repeatmask
│   │   │       └── 09_repeatmodeler_repeatmasker_docker.sh
│   │   └── yahs_pipeline
│   │       ├── 04_yahs_run.sh
│   │       ├── 05_finalize_assembly.sh
│   │       ├── min1000
│   │       │   ├── 06_convert_name_directions.sh
│   │       │   └── renamed
│   │       │       ├── 07_separate_chr.sh
│   │       │       ├── 08_TRASH_array.sh
│   │       │       └── 08_exec_TRASH_array.sh
│   │       └── omnic_mapping
│   │           └── 03_mapping_omnic.sh
│   └── hic_hap2_20240515
│       ├── annotation
│       │   ├── braker_prot
│       │   │   └── 10_braker_run.sh
│       │   ├── entap
│       │   │   ├── 16_entap_run.sh
│       │   │   └── 17_select_annotated_genes.sh
│       │   ├── gemoma
│       │   │   ├── 14_gemoma_run.sh
│       │   │   └── 15_gffcompare_run.sh
│       │   └── repeatmask
│       │       └── 09_repeatmodeler_repeatmasker_docker.sh
│       └── yahs_pipeline
│           ├── 04_yahs_run.sh
│           ├── 05_finalize_assembly.sh
│           ├── min1000
│           │   ├── 06_convert_name_directions.sh
│           │   └── renamed
│           │       ├── 07_separate_chr.sh
│           │       ├── 08_TRASH_array.sh
│           │       └── 08_exec_TRASH_array.sh
│           └── omnic_mapping
│               └── 03_mapping_omnic.sh
└── README.md
```
## Sharing/Access information

1. Licenses/restrictions placed on the data: CC0 1.0 Universal (CC0 1.0) Public Domain
2. Links to publications that cite or use the data:\
   Yo Y. Yamasaki, Atsushi Toyoda, Mitsutaka Kadota, Shigehiro Kuraku, Jun Kitano
   3D genome constrains breakpoints of inversions that can act as barriers to gene flow in the stickleback
