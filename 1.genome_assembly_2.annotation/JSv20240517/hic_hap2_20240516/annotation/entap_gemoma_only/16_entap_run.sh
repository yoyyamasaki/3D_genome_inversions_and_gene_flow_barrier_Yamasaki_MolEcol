#!/bin/bash
#Conduct functional annotation by entap. We also conduct functional annotation by EggNOG mapper v2.1.12 in a web page. 
#Dependency: entap v1.0.1

ln -s ../gemoma/JS_male_hap2.v20240517.protein.fasta ./

#Configuration
~/software/EnTAP-1.0.1/EnTAP --config -d /home/yo/db/entap_DB/uniprot_sprot.fasta -d /home/yo/db/entap_DB/uniprot_trembl.fasta -d /home/yo/db/entap_DB/vertebrate_other.protein.faa -t 20 --ini ~/software/EnTAP-1.0.1/entap_config.ini --out-dir /home/yo/db/entap_DB

#Exectute entap. Check NCBI database first, then check uniprot.
~/software/EnTAP-1.0.1/EnTAP --runP -i ./JS_male_hap2.v20240517.protein.fasta -d /home/yo/db/entap_DB/bin/vertebrate_other.dmnd -d /home/yo/db/entap_DB/bin/uniprot_sprot.dmnd -d /home/yo/db/entap_DB/bin/uniprot_trembl.dmnd -t 20 --ini /home/yo/software/EnTAP-1.0.1/entap_config.ini

