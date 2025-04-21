#!/bin/bash
#Calculate A/B compartment and make inputs for spectralTAD by fanc.
#Dependency: fanc v0.9.1

REF=/home/yo/db/PO_merged_maskXII_XXI_Y_20240517/PO_male_merged.v20240517.maskXII_XXI_Y.fa
HEADER=iconic_PO_refPO_merged_maskXII_XXI_Y_pairtools
NORM_METHOD=KR
THREAD=20
BINS=(2.5mb 1mb 500kb 250kb 100kb 50kb 25kb 10kb)          

CHR=($(grep ">chr" ${REF}|perl -pe 's/>//g'))

#Normalization: KR
#Generating contact matrix
mkdir -p ${HEADER}_${NORM_METHOD}/hic/binned

for i in ${BINS[@]};do
	fanc hic --deepcopy --norm-method ${NORM_METHOD} ${HEADER}.hic@${i} ${HEADER}_${NORM_METHOD}/hic/binned/${HEADER}_${i}.hic
done 

#Make input files for SpectralTAD
for i in 10kb 25kb 50kb 100kb 250kb 500kb 1mb 2mb 5mb;do

	BIN=${i}
        fanc dump --uncorrected --only-intra ${HEADER}_${NORM_METHOD}/hic/binned/${HEADER}_${BIN}.hic ${HEADER}_${NORM_METHOD}/hic/binned/${HEADER}_uncorrected_${BIN}.txt
done

for i in 10kb 25kb 50kb 100kb 250kb 500kb 1mb 2mb 5mb;do

	BIN=${i}
        fanc dump --only-intra ${HEADER}_${NORM_METHOD}/hic/binned/${HEADER}_${BIN}.hic ${HEADER}_${NORM_METHOD}/hic/binned/${HEADER}_${NORM_METHOD}_${BIN}.txt
done
         
#AB compartment analysis
for i in ${BINS[@]};do
	BIN=${i} 
	mkdir -p ${HEADER}_${NORM_METHOD}/compartments/${BIN}
	
	fanc compartments --domains ${HEADER}_${NORM_METHOD}/compartments/${BIN}/${HEADER}_${NORM_METHOD}_${BIN}_domains.bed \
        	          --eigenvector ${HEADER}_${NORM_METHOD}/compartments/${BIN}/${HEADER}_${NORM_METHOD}_${BIN}_eigen.bed \
        	          --enrichment-profile ${HEADER}_${NORM_METHOD}/compartments/${BIN}/${HEADER}_${NORM_METHOD}_${BIN}_profile.png \
        	          --genome ${REF} \
        	          -f \
        	          ${HEADER}_${NORM_METHOD}/hic/binned/${HEADER}_${BIN}.hic ${HEADER}_${NORM_METHOD}/compartments/${BIN}/${HEADER}_${NORM_METHOD}_${BIN}.ab
        	          
	for j in ${CHR[@]};do
			REGION=${j}
			fancplot -o ${HEADER}_${NORM_METHOD}/compartments/${BIN}/${HEADER}_${NORM_METHOD}_${BIN}.ab.${REGION}.png ${REGION} \
				 -p square ${HEADER}_${NORM_METHOD}/compartments/${BIN}/${HEADER}_${NORM_METHOD}_${BIN}.ab \
     				 -vmin -0.75 -vmax 0.75 -c RdBu_r \
     				 -p line ${HEADER}_${NORM_METHOD}/compartments/${BIN}/${HEADER}_${NORM_METHOD}_${BIN}_eigen.bed
     	done
     	
done	 

