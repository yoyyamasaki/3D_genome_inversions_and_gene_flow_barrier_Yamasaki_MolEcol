# Scripts used in "3D genome constrains breakpoints of inversions that can act as barriers to gene flow in the stickleback"

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
- This is the scripts used in "Yamasaki et al. 3D genome constrains breakpoints of inversions that can act as barriers to gene flow in the stickleback."
- Each directory has a number at the head of the name. They correspond to the sections of Materials & Methods in the paper.
- Scripts for 1.genome assembly and 2.annotation were in one directory together.
- If you want to use these scripts for your datasets, please set an appropriate path in your environment.
- The software used and its versions are written in each script.

## Directories
- 1.genome_assembly_2.annotation  
This directory includes the scripts for genome assembly and annotation. We first conducted haplotype-phased assembly for each species and obtained two haplotypes per species. Then, we conducted annotation for each haplotype.

- 3.inversion_population_genomics  
This directory includes the scripts for inversion detection and population genomic analysis using short-read sequencing. We first detected inversions by whole genome alignment. Then, we calculated some statistics within inversions and colinear regions.

- 4.3D_genome_analysis  
This directory includes the script for 3D genome analysis. We first made Hi-C contact matrices for each species. Then, we detected A/B compartment structure and topologically associating domain (TAD) structure.

- 5.TADs_and_inversions  
This directory includes the script for testing the association between 3D genome structure and inversions. First, we tested the association between TAD structures and inversions between Gasterosteus nipponicus (JS) and G. aculeatus (PO) comparison (PO_JS). Second, we tested the relationship between the stcikeback version 5 reference genome and G. nipponicus (JS_v5) or G. aculeatus (PO_v5). Third, we tested the association between A/B compartment structure and inversions (PO_JS).

## Sharing/Access information
1. Licenses/restrictions placed on the data: CC0 1.0 Universal (CC0 1.0) Public Domain
2. Links to publications that cite or use the data:\
   Yo Y. Yamasaki, Atsushi Toyoda, Mitsutaka Kadota, Shigehiro Kuraku, Jun Kitano
   3D genome constrains breakpoints of inversions that can act as barriers to gene flow in the stickleback
