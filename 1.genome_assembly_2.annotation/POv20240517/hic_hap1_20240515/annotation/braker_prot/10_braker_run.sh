#! /bin/bash
#Run braker2 annotation pipeline
#Dependency: BRAKER2

BRAKER_PATH=/home/yo/singularity_container/braker3.sif
GENOME=genome.fasta
PROTEIN=/home/yo/db/OrthoDB11_Vertebrata/Vertebrata.fa

#Protein data only
singularity exec --bind ${PWD}:/mnt --bind /home:/home ${BRAKER_PATH} braker.pl --species PO_male_hap1 --genome /mnt/${GENOME} --threads 20 --gff3 --prot_seq ${PROTEIN} --workingdir /mnt


