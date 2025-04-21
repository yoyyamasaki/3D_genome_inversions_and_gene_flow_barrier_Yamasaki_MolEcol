#! /bin/bash
#Run fastp to clean RNA-seq reads
#Dependency: fastp v0.23.4

parallel -j 10 bash fastp_array.sh ::: {1..42}
