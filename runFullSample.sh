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

data=(
	"0"
	"1"
	"2"
	"3"
	"4"	
	"5"
	"6"
	"7"
	"8"
	"9"
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
)

for fileNumber in "${data[@]}"; do
       echo "file number: $fileNumber"
       ./runTimeClus_data -hard /dcache/alice/jlomker/LHC22o_pass7_minBias/526641/vectorTree_LHC22o_pass7_$fileNumber.root -hardtype ROOT -hardtreename o2Tracks -hardvarname particle_data -nev 50000
       mv JetToyHIResultTimeClus_data.root /dcache/alice/jlomker/LHC22o_pass7_minBias/526641/results/JetToyHIResultTimeClus_data_$fileNumber.root
 #./runTimeClus_MC -hard input/vectorTree_LHC25a2b.root -hardtype ROOT -hardtreename o2TracksFormatted -hardvarname particle_truth -reco input/vectorTree_LHC25a2b.root -hardtype ROOT -hardtreename o2TracksFormatted -hardvarname particle_reco -nev 10000
       # first we run over 'all' pythia 
      # ./runTimeClus_modeldependence -hard /dcache/alice/bashof/Pythia/5TeV_pu14/PythiaEventsTune14PtHat132_$fileNumber.pu14 -nev 10000
      # for hepMC:
      #./runTimeClus_modeldependence -hard /dcache/alice/jlomker/Herwig/LHC_5p36TeV/LHC-Matchbox-08.hepmc -hardtype HepMC2 -nev 100
       #mv JetToyHIResultTimeClus_model.root /dcache/alice/jlomker/resultsPYTHIA/JetToyHIResultTimeClus_modelPYTHIA_$fileNumber.root
       
       #./runTimeClus_modeldependence -hard /dcache/alice/bashof/Herwig/5TeV_pu14/Herwig_5TeV_pt_9_job_$fileNumber.pu14 -nev 10000
       #mv JetToyHIResultTimeClus_model.root /dcache/alice/jlomker/resultsHERWIG/JetToyHIResultTimeClus_modelHERWIG_$fileNumber.root
done

echo "all done!"

