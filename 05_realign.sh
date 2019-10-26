#!/bin/bash

bam=$1
intervals=$2
ref=${3:-dmel}
refFasta=$($(dirname $0)/../get_ref_fasta.sh $ref)
if [ -e $refFasta ]; then echo "using $refFasta as reference"; else echo "reference fasta not found; exiting"; exit 1; fi
gatkDir=/home/groups/dpetrov/SOFTWARE/bin/
#gatkDir=software  #for biodap20
picardDir=/home/groups/dpetrov/SOFTWARE/picard-tools

# make ref dictionary if not already there
refDict=$(echo $refFasta | sed 's/.fasta$/.dict/')
if [ ! -e $refDict ]; then 
	java -jar $picardDir/picard.jar CreateSequenceDictionary \
	REFERENCE=$refFasta \
	OUTPUT=$refDict
fi

#re-align bam
newbam=$(echo $bam | sed 's/.bam/.realign.bam/')
java -jar $gatkDir/GenomeAnalysisTK.jar -T IndelRealigner -R $refFasta -I $bam -targetIntervals $intervals -o $newbam --filter_mismatching_base_and_quals
samtools index $newbam


