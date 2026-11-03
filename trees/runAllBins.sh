#!/bin/bash

bins=("20" "30" "40" "50" "60" "70" "80" "90" "110" "135" "165" "200")
j="200"

for i in "${bins[@]}";
do
    #sed -i "s/JetToyHIResultTimeClusBkg_13p6_$j/JetToyHIResultTimeClusBkg_13p6_$i/" runJetTreeSig.C jetTreeSig.h
    sed -i "s/JetToyHIResultTimeClusBkgJewel_$j/JetToyHIResultTimeClusJewel_$i/" runJetTreeSig.C jetTreeSig.h
    j=$i
    root -b -q runJetTreeSig.C
done

