#!/bin/bash

rm -rf external_med_fix.dat
for JOB in {1..100}
do
	for PT in {20,30,40,50,60,70,80,90,110,135,165,200}
	do
	FILE=/dcache/alice/jlomker/Jewel/med/analysisResults/JetToyHIResultTimeClusJewel_${PT}_${JOB}.root     
	if [ -f $FILE ]; then
   		echo "File $FILE exists."
	else
   		echo "${PT}, ${JOB}" >> external_med_fix.dat
	fi

	done	
done
