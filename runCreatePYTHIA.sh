#!/bin/bash

ptHat=(
	"200"
	"20"
	"30"
	"40"
	"50"
	"60"	
	"70"
	"80"
	"90"
	"110"
	"135"
	"165"
)

for pthat in "${ptHat[@]}"; do
       echo "file number: $pthat"
       ./runCreatePythiaEvents -pthat $pthat -nev 1000000 > stats$pthat.txt
       mv PythiaEventsTune14PtHat$pthat.pu14 /dcache/alice/jlomker/PYTHIA/LHC_13p6TeV/PythiaEventsTune14PtHat$pthat.pu14
       mv stats$pthat.txt /dcache/alice/jlomker/PYTHIA/LHC_13p6TeV/
done

echo "all done!"

