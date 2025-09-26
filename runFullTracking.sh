#!/bin/bash

files=(
	"10"
	"11"
	"12"
	"13"
	"14"
	"15"
	"16"
	"17"
	"18"
	"19"
	"20"
        "21"
        "22"
        "23"
        "24"
        "25"
        "26"
        "27"
        "28"
        "29"
        "30"
        "31"
        "32"
        "33"
        "34"
        "35"
        "36"
        "37"
        "38"
        "39"
        "40"
)

for fileNumber in "${files[@]}"; do
       echo "file number: $fileNumber"
        ./runTimeClus_tracking -hard /dcache/alice/bashof/Pythia/5TeV_pu14/PythiaEventsTune14PtHat132_$fileNumber.pu14 -hardtype PU14 -nev 10000
	#./runTimeClus_MC -hard input/vectorTree_LHC25a2b.root -hardtype ROOT -hardtreename o2TracksFormatted -hardvarname particle_truth -reco input/vectorTree_LHC25a2b.root -hardtype ROOT -hardtreename o2TracksFormatted -hardvarname particle_reco -nev 10000
       mv JetToyHIResultTimeClus_tracking.root /dcache/alice/jlomker/resultsTracking/JetToyHIResultTimeClus_tracking_$fileNumber.root
done

echo "all done!"

