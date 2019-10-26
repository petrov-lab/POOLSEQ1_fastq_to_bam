#!/bin/bash

bam=$1
outDir=${2:-$(dirname $bam)}

dedupbam=$outDir/$(basename $bam | sed 's/.bam/.dedup.bam/')
picardDir=$PI_HOME/SOFTWARE/picard-tools
#picardDir=software/picard-tools-1.140/ #for biodap20

java -Xmx16g -jar $picardDir/picard.jar MarkDuplicates ASSUME_SORTED=true \
		INPUT=$bam \
		OUTPUT=$dedupbam \
		METRICS_FILE=${dedupbam}.metrics \
		VALIDATION_STRINGENCY=LENIENT \
		REMOVE_DUPLICATES=TRUE
samtools index $dedupbam
