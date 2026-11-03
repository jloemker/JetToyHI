#!/bin/bash

#$1=$fileNumber

ptHat=(
	"200"
	"50"
	"60"
	"70"
	"80"
	"90"
	"110"
	"135"
	"165"
	"40"
	"30"
	"20"
)

source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_105 x86_64-el9-gcc12-opt

#for fileNumber in "${ptHat[@]}"; do
cd $TMPDIR
mkdir $1
cd $1
cp /data/alice/jlomker/jetToyHI/JetToyHI/runTimeClusBkg .
mkdir include/
cd include/
cp /data/alice/jlomker/jetToyHI/JetToyHI/include/tr_eff.root .
cd ../
echo "file number: $1"

./runTimeClusBkg -hard /dcache/alice/jlomker/PYTHIA/LHC_13p6TeV/PythiaEventsTune14PtHat$1.pu14 -nev 500000
mv JetToyHIResultTimeClusBkg.root /dcache/alice/jlomker/PYTHIA/LHC_13p6TeV/analysisResults/JetToyHIResultTimeClusBkg_13p6_$1.root

echo "done with: $1"

#echo "all done!"
