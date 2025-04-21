#!/bin/bash
#Extract inversion and colinear regions in bed format.

VCF=syri20240917_PO_JS_asm5syri.vcf
OUT_INV=syri20240917_PO_JS_asm5syri.INV.bed
OUT_BEDPE=syri20240917_PO_JS_asm5syri.inversion.bedpe
OUT_SYNAL=syri20240917_PO_JS_asm5syri.SYNAL.bed

#Extract inversion region
cat ${VCF}|awk -v 'OFS=\t' '$3 ~ /^INV[1234567890]/ {print $1, $2-1, $8}' ${VCF}|perl -pe 's/END=//g'|perl -pe 's/;.*$//g' > ${OUT_INV}

#Extract inversion regions in bedpe format
cat ${VCF}|awk -v 'OFS=\t' '$3 ~ /^INV[1234567890]/ {print $1, $2-1, $8}' ${VCF}|perl -pe 's/;Parent.*$//g'|perl -pe 's/;/\t/g'|\
perl -pe 's/(END=|ChrB=|StartB=|EndB=)/\t/g'|awk -v 'OFS=\t' '{print $1, $2, $3, $4, $5-1, $6}' > ${OUT_BEDPE}

#Extract colinear region
cat ${VCF}|awk -v 'OFS=\t' '$3 ~ /^SYNAL[1234567890]/ {print $1, $2-1, $8}' ${VCF}|perl -pe 's/END=//g'|perl -pe 's/;.*$//g' > ${OUT_SYNAL}
