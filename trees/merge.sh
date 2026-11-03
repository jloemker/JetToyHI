#!/bin/bash

rm -rf mergeList
mkdir mergeList

for PT in {20,30,40,50,60,70,80,90,110,135,165,200}
do
	for JOB in {1..100}
	do
		echo /dcache/alice/jlomker/Jewel/vac/analysisResults/JetToyHIResultTimeClusJewel_${PT}_${JOB}.root >> hadd_${PT}_this.list		
	done
	mv hadd_${PT}_this.list mergeList/hadd_${PT}_this.list
done

echo "created merge lists"

for PT in {20,30,40,50,60,70,80,90,110,135,165,200}
do
	echo "testing"
	#hadd -k -f /dcache/alice/jlomker/Jewel/vac/analysisResults/JetToyHIResultTimeClusJewel_${PT}.root mergeLists/@hadd_${PT}_this.list
done

echo "done merging things"
