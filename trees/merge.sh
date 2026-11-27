#!/bin/bash

rm -rf mergeList
mkdir mergeList

for PT in {20,30,40,50,60,70,80,90,110,135,165,200}
do
	for JOB in {1..100}
	do
		#echo /dcache/alice/jlomker/Jewel/vac/analysisResults/JetToyHIResultTimeClusJewel_${PT}_${JOB}.root >> hadd_${PT}_this.list		
		echo /data/alice/jlomker/jetToyHI/JetToyHI/trees/jewelVacMinPt60/JetToyHIResultTimeClusJewel_${PT}_${JOB}_Histos.root >> hadd_${PT}_this.list		
	done
	mv hadd_${PT}_this.list mergeList/hadd_${PT}_this.list
done

echo "created merge lists"

for PT in {20,30,40,50,60,70,80,90,110,135,165,200}
#for PT in {135,165,200}
do
	echo "merging for pt bin ${PT}"
	#hadd -k -f /dcache/alice/jlomker/Jewel/vac/analysisResults/JetToyHIResultTimeClusJewel_${PT}.root @mergeList/hadd_${PT}_this.list
	hadd -k -f /data/alice/jlomker/jetToyHI/JetToyHI/trees/jewelVacMinPt60/JetToyHIResultTimeClusJewel_${PT}_Histos.root @mergeList/hadd_${PT}_this.list
done

echo "done merging things"

for PT in {20,30,40,50,60,70,80,90,110,135,165,200}
do
        for JOB in {1..100}
        do
                mv /data/alice/jlomker/jetToyHI/JetToyHI/trees/jewelVacMinPt60/JetToyHIResultTimeClusJewel_${PT}_${JOB}_Histos.root /dcache/alice/jlomker/Jewel/vac/analysisResults/histograms/JetToyHIResultTimeClusJewel_${PT}_${JOB}.root                                  
        done
done

echo "done moving intermediate histogram files to dcache"
