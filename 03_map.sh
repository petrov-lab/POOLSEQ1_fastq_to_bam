#!/bin/bash

#### map trimmed reads to the d mel ref genome

## SET VARS
R1=$1
R2=$2
samp=$3
bamDir=$4
threads=${5}
ref=${6:-dmel}
refFasta=$($HOME/scripts/PIPELINE1_fastq_to_bam/get_ref_fasta.sh $ref)
if [ -e $refFasta ]; then echo "using $refFasta as reference"; else echo "reference fasta not found for $ref; exiting"; exit 1; fi

picardDir=$PI_HOME/SOFTWARE/picard-tools
ml biology samtools java bwa; 

echo "running $samp"


## MAP PAIRED READS WITH PAIRED END BWA	
echo "---mapping paired end reads..."
	bwa mem -t $threads -R "@RG\tID:${samp}.pe\tSM:${samp}" \
	$refFasta \
	$R1 \
	$R2 | \
	samtools view -Suh - | \
	samtools sort -@10 - > ${samp}.pe.bam


echo "---indexing bam file..."
	samtools index ${bamDir}/${samp}.bam

echo "---finished."





	
