#! /bin/bash
#Run braker2 annotation pipeline
#Dependency: BRAKER2

BRAKER_PATH=/home/yo/singularity_container/braker3.sif
GENOME=PO_male_hap2.v20240517.softmasked.fa
PROTEIN=/home/yo/db/OrthoDB11_Vertebrata/Vertebrata.fa

#Protein data only
singularity exec --bind ${PWD}:/mnt --bind /home:/home --bind /raid:/raid ${BRAKER_PATH} braker.pl --species PO_male_hap2 --genome /mnt/${GENOME} --threads 20 --gff3 --prot_seq ${PROTEIN} --workingdir /mnt


