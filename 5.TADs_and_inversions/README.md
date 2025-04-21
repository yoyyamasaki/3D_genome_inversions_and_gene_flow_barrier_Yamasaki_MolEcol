# README for TADs and inversion analysis

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
- This directory includes the scripts for testing the association between TADs and inversions
- We first tested the association between TADs and inversions in G. nipponicus and G. aculeatus.
- We conducted the same analysis for the inversions between G. nipponicus and the stickleback v5 genome and inversions between G. aculeatus and the stickleback v5 genome.
- We also tested the association between AB compartment boundaries and inversion breakpoints.
- Each script has a number at the head of the script name. This is the order in which the scripts are executed.
- If you want to use these scripts for your datasets, please set an appropriate path in your environment.
- The software used and its versions are written in each script.

## Data
- Inversion data estimated in section 3
- TAD estimated in section 4
- AB compartment estimated in section 4

## Brief explanation for scripts
- 01_select_inversion.r: Make inversion bed files. They were based on G. aculeatus and G. nipponicus genome coordinate. Also make 25kb < inversion bed files.
- 02_calc_distance_inversion_TAD.r: Make table and plot which describe the distance between inversion breakpoints and the nearest TAD boundaries.
- 03_count_INVbreakpoint_overlapped_bins.r: Test association betwen inversion breakpoints and TAD boundaries.
- 04_syri_run_asm5.sh: Detect structural variations from pairwise alignment of the two genomes.
- 05_make_bed_inversion.sh: Extract inversion regions in bed format.
- 06_calc_distance_inversion_TAD.r: Make table and plot which describe the distance between inversion breakpoints and the nearest TAD boundaries.
- 07_calc_distance_inversion_AB.r: Make table and plot which describe the distance between inversion breakpoints and the nearest AB compartment boundaries.
- 08_count_AB_associted_boundary_rates.r: Test association betwen inversion breakpoints and AB boundaries.

## Directory structure
```
.
├── JS_v5
│   ├── 04_syri_run_asm5.sh
│   └── syri20241127_JS_v5_asm5
│       ├── 05_make_bed_inversion.sh
│       └── TAD_SV
│           └── 06_calc_distance_inversion_TAD.r
├── PO_JS
│   ├── 01_select_inversion.r
│   ├── JS
│   │   ├── SV_ABcompartment
│   │   │   ├── 07_calc_distance_inversion_AB.r
│   │   │   └── 08_count_AB_associted_boundary_rates.r
│   │   └── SV_TAD
│   │       ├── 02_calc_distance_inversion_TAD.r
│   │       └── 03_count_INVbreakpoint_overlapped_bins.r
│   └── PO
│       ├── SV_ABcompartment
│       │   ├── 07_calc_distance_inversion_AB.r
│       │   └── 08_count_AB_associted_boundary_rates.r
│       └── SV_TAD
│           ├── 02_calc_distance_inversion_TAD.r
│           └── 03_count_INVbreakpoint_overlapped_bins.r
├── PO_v5
│   ├── 04_syri_run_asm5.sh
│   └── syri20241127_PO_v5_asm5
│       ├── 05_make_bed_inversion.sh
│       └── SV_TAD
│           └── 06_calc_distance_inversion_TAD.r
└── README.md
```

## Sharing/Access information

1. Licenses/restrictions placed on the data: CC0 1.0 Universal (CC0 1.0) Public Domain
2. Links to publications that cite or use the data:\
   Yo Y. Yamasaki, Atsushi Toyoda, Mitsutaka Kadota, Shigehiro Kuraku, Jun Kitano
   3D genome constrains breakpoints of inversions that can act as barriers to gene flow in the stickleback
