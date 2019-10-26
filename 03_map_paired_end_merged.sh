#!/bin/bash

#### map trimmed reads to the d mel ref genome

## SET VARS

samp=${1}
fastqDir=${2}
bamDir=${3}
threads=${4}
ref=${5:-dmel}
refFasta=$($HOME/scripts/PIPELINE1_fastq_to_bam/get_ref_fasta.sh $ref)
if [ -e $refFasta ]; then echo "using $refFasta as reference"; else echo "reference fasta not found for $ref; exiting"; exit 1; fi

picardDir=$PI_HOME/SOFTWARE/picard-tools

echo "running $samp"

## MAP ASSEMBLED AND UNPAIRED READS WITH SINGLE END BWA
	echo "---mapping single end reads..."
	cat ${fastqDir}/${samp}.assembled.fastq.gz | \
	bwa mem -t $threads -R "@RG\tID:${samp}.se\tSM:${samp}" \
	$refFasta \
	- | \
	samtools view -Suh - | \
	samtools sort -@10 - > ${samp}.se.bam

## MAP PAIRED READS WITH PAIRED END BWA	
	echo "---mapping paired end reads..."
	bwa mem -t $threads -R "@RG\tID:${samp}.pe\tSM:${samp}" \
	$refFasta \
	${fastqDir}/${samp}.unassembled.forward.fastq.gz \
	${fastqDir}/${samp}.unassembled.reverse.fastq.gz | \
	samtools view -Suh - | \
	samtools sort -@10 - > ${samp}.pe.bam

## MERGE BAM FILES AND INDEX 
	echo "---merging bam files..."
	java -jar $picardDir/picard.jar MergeSamFiles \
	I=${samp}.se.bam \
	I=${samp}.pe.bam \
	O=${bamDir}/${samp}.bam
	
	echo "---indexing bam file..."
	samtools index ${bamDir}/${samp}.bam
	
	echo "---cleaning up..."
	rm ${samp}.se.bam ${samp}.pe.bam
	
	echo "---finished."







	
