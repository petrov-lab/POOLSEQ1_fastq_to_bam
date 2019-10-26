#!/bin/bash

myDir=${1}

## PCT TRIMMED
paste <(ls $myDir/trimmed/*.log | sed 's/.*\///' | sed 's/R1_001.fastq-trimmed.log//') \
<(cat $myDir/trimmed/*.log | grep ") trimmed" | cut -f2 -d'(' | cut -f1 -d')') \
 > $myDir/allSamps.pctTrimmed

## PCT MERGED
paste <(ls $myDir/jobScripts/merge.*.output | sed 's/.*\///' | cut -f1-3 -d'_' | sed 's/merge.//' | sed 's/.output//') \
<(cat $myDir/jobScripts/merge.*.output | grep 'Assembled reads \.' | cut -f2 -d'(' | cut -f1 -d')') \
 > $myDir/allSamps.pctMerged

## PCT DUP
paste <(ls $myDir/mapped/*.metrics | sed 's/.*\///' | sed 's/.dedup.bam.metrics//') \
<(cat $myDir/mapped/*.metrics | grep 'Unknown' | cut -f8) \
 > $myDir/allSamps.pctDup

##CHROM COV (#after running mosdepth from biodap)
paste <(ls $myDir/mapped/*.dedup.realign.bam.mosdepth.regions.bed | sed 's/.*\///' | sed 's/.dedup.realign.bam.mosdepth.regions.bed//' ) \
<(cat $myDir/mapped/*.dedup.realign.bam.mosdepth.regions.bed | grep ^[23X] | cut -f4 | awk '{
	covString=covString","$1;if((NR%5)==0){print covString; covString=""}
}') | sed 's/\t//' > $myDir/allSamps.chromCov

#paste <(ls $myDir/mapped/*.dedup.bqsr.bam.mosdepth.regions.bed | sed 's/.*\///' | sed 's/.dedup.bqsr.bam.mosdepth.regions.bed//' ) \
<(cat $myDir/mapped/*.dedup.bqsr.bam.mosdepth.regions.bed | grep ^[23X] | cut -f4 | awk '{
	covString=covString","$1;if((NR%5)==0){print covString; covString=""}
}') | sed 's/\t//' > $myDir/allSamps.chromCov

## PCT MAPPED
paste <(ls $myDir/mapped/*.flagstat | sed 's/.*\///' | sed 's/.dedup.realign.bam.flagstat//' ) \
<(cat $myDir/mapped/*.flagstat | grep "mapped (" | cut -f2 -d'(' | cut -f1 -d' ') \
> $myDir/allSamps.pctMapped
