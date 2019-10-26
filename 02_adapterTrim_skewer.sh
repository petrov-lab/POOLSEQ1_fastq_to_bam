#!/bin/bash

#### perform adapter trimming with skewer
R1=$1
R2=$2
outDir=$3
threads=${4:-4}
qualMin=${5:-20}
adapterType=${6:-nextera}   
adapterFile=$GROUP_HOME/poolseq_pipelines/PIPELINE1_fastq_to_bam/adapters/adapter_${adapterType}R
skewerDir=$PI_HOME/SOFTWARE/bin

$skewerDir/skewer -t $threads -q $qualMin -m pe --compress -x ${adapterFile}1.fa -y ${adapterFile}2.fa $R1 $R2

mv $(echo $R1 | sed 's/.gz/-trimmed*/') $outDir
