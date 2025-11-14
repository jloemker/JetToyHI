#!/bin/bash

ptHat=(
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
       ./runCreatePythiaEvents -pthat $pthat -nev 100000
       mv PythiaEventsTune14PtHat$pthat.pu14 /dcache/alice/jlomker/PYTHIA/LHC_5p36TeV/PythiaEventsTune14PtHat$pthat.pu14
done

echo "all done!"

