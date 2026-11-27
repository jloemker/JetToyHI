#!/bin/bash

bins=("20" "30" "40" "50" "60" "70" "80" "90" "110" "135" "165" "200")
j="200"
job=1
for JOB in {1..10};
do
    for i in "${bins[@]}";
    do
        #echo "JetToyHIResultTimeClusJewel_${j}_${JOB} JetToyHIResultTimeClusJewel_${i}_${JOB}"
        sed -i "s/JetToyHIResultTimeClusJewel_${j}_${job}/JetToyHIResultTimeClusJewel_${i}_${JOB}/" runJetTreeSig.C jetTreeSig.h
        j=${i}
        job=${JOB}
        root -b -q runJetTreeSig.C
	mv JetToyHIResultTimeClusJewel_${i}_${JOB}_Histos.root jewelVacMinPt60/JetToyHIResultTimeClusJewel_${i}_${JOB}_Histos.root 
    done
done
