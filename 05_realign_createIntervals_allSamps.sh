#!/bin/bash

bamDir=$1
threads=$2
ref=${3:-dmel}
intervalFile=${4:-$bamDir/allSamps.intervals}

refFasta=$($(dirname $0)/get_ref_fasta.sh $ref)
if [ -e $refFasta ]; then echo "using $refFasta as reference"; else echo "reference fasta not found; exiting"; exit 1; fi
gatkDir=/home/groups/dpetrov/SOFTWARE/bin/
picardDir=/home/groups/dpetrov/SOFTWARE/picard-tools

# make ref dictionary if not already there
refDict=$(echo $refFasta | sed 's/.fasta$/.dict/')
if [ ! -e $refDict ]; then 
	java -jar $picardDir/picard.jar CreateSequenceDictionary \
	REFERENCE=$refFasta \
	OUTPUT=$refDict
fi

# make list of input bams
inputString=""
for bam in $bamDir/*.dedup.bam; do
	inputString=$inputString" -I "$bam
done

#make target indel file
java -jar $gatkDir/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $refFasta -nt $threads \
 -o $intervalFile \
$inputString 

