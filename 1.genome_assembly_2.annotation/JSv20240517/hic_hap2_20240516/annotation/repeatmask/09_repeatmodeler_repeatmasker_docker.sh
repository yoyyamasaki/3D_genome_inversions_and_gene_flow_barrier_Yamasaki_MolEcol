#! /bin/bash
#Run repeatmodeler and repeatmasker
#Dependency: repeatmodeler v2.0.5, repeatmasker v4.1.6

#Activate docker image
docker run -v ${PWD}:/opt/src -v /home/yo/db/repeatmasker_lib/Libraries:/opt/RepeatMasker/Libraries -v /home:/home -v /raid:/raid -it --rm dfam/tetools:latest

BuildDatabase -name JS_male_hap2.v20240517.fa JS_male_hap2.v20240517.fa

RepeatModeler -database JS_male_hap2.v20240517.fa -LTRStruct -threads 20

exit

docker run -v ${PWD}:/opt/src -v /home/yo/db/repeatmasker_lib/Libraries:/opt/RepeatMasker/Libraries -v /home:/home -v /raid:/raid -it --rm dfam/tetools:latest 

export LIBDIR=/opt/RepeatMasker/Libraries

RepeatMasker -lib JS_male_hap2.v20240517.fa-families.fa -par 48  -html -gff -xsmall /opt/src/JS_male_hap2.v20240517.fa

docker run -v ${PWD}:/opt/src -v /home/yo/db/repeatmasker_lib/Libraries:/opt/RepeatMasker/Libraries -v /home:/home -v /raid:/raid -it --rm dfam/tetools:latest  ProcessRepeats -maskSource /opt/src/JS_male_hap2.v20240517.fa  JS_male_hap2.v20240517.fa.cat.gz

