#!/bin/bash
#Filtering WGS fastq files in parallel.
#Dependnecy: fastp v0.23.2

parallel --jobs 11 bash ./fastp_array.sh ::: {1..11}

