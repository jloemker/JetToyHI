#!/bin/bash

source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_105 x86_64-el9-gcc12-opt

cd $TMPDIR
mkdir $1_$2
cd $1_$2
cp /data/alice/jlomker/jetToyHI/JetToyHI/runTimeClusJewel .

echo "pTbin: $1"
echo "file Number: $2"
./runTimeClusJewel -hard /dcache/alice/jlomker/Jewel/vac/2.2.0_$1_$2.hepmc -hardtype HepMC2 -nev 100000
mv JetToyHIResultTimeClusJewel.root /dcache/alice/jlomker/Jewel/vac/analysisResults/JetToyHIResultTimeClusJewel_$1_$2.root

echo "done with pt: $1 and file Number: $2"

#echo "all done!"
