#!/bin/bash
#Tandem repeat analysis by TRASH
#This script is for executing 07_TRASH_array.sh in parallel.
#Dependency: TRASH v1.2

parallel -j 10 bash 07_TRASH_array.sh ::: {1..20}
