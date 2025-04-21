#!/bin/bash
#Extract inversion regions in bed format.

VCF=syri20241127_PO_v5_asm5syri.vcf
OUT=syri20241127_PO_v5_asm5syri.INV.bed
OUT_SVs=syri20241127_PO_v5_asm5syri.SVs.bed
OUT_BEDPE=syri20241127_PO_v5_asm5syri.inversion.bedpe

cat ${VCF}|awk -v 'OFS=\t' '$3 ~ /^INV[1234567890]/ {print $1, $2-1, $8}' ${VCF}|perl -pe 's/END=//g'|perl -pe 's/;.*$//g' > ${OUT}

cat ${VCF}|awk -v 'OFS=\t' '$3 ~ /^(CPG|CPL|TDM|INV|TRANS|INVTR|DUP|INVDP|INS|DEL)[1234567890]/ {print $1, $2-1, $8, $3}' |perl -pe 's/END=//g'|perl -pe 's/;.*\t/\t/g' > ${OUT_SVs}

cat ${VCF}|awk -v 'OFS=\t' '$3 ~ /^INV[1234567890]/ {print $1, $2-1, $8}' ${VCF}|perl -pe 's/;Parent.*$//g'|perl -pe 's/;/\t/g'|\
perl -pe 's/(END=|ChrB=|StartB=|EndB=)/\t/g'|awk -v 'OFS=\t' '{print $1, $2, $3, $4, $5-1, $6}' > ${OUT_BEDPE}

