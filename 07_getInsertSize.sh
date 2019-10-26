#!/bin/bash

bam=$1
outDir=${2:-$(dirname $bam)}

metricsFile=$outDir/$(basename $bam | sed 's/.bam/.insertMetrics.txt/')
picardDir=$PI_HOME/SOFTWARE/picard-tools
#picardDir=software/picard-tools-1.140/ #for biodap20

java -jar $picardDir/picard.jar CollectInsertSizeMetrics I=$bam O=$metricsFile H=hist.pdf M=.5
