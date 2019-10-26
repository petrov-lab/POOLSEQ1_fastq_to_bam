#!/bin/bash

rootDir=/scratch/users/greensi/orchard2017/
jobscriptDir=/scratch/users/greensi/orchard2017/jobScripts
if [ ! -e $rootDir ]; then mkdir -p $rootDir; fi
if [ ! -e $jobscriptDir ]; then mkdir -p $jobscriptDir; fi
picardDir=$PI_HOME/SOFTWARE/picard-tools
gdrive="petrov-td:Dmel_cageExperiments_SGreenblum/2017-orchard-dmel"
sampKey=$SCRATCH/18057-02-SampleKey.txt
create_batch_script=/home/groups/dpetrov/poolseq_pipelines/submit_sbatch_jobs.sh

doGETFASTQS=0
doRENAME=0
doFASTQC=0
doTRIM=1
doMERGE_READS=1
doMAP=1
doDEDUP=1
doREALIGN=1
doCOV=1
doINSERTSIZE=0
doBACKUP=1

ml biology; ml devel; ml samtools; ml java; ml bwa; ml bedtools; ml system; ml parallel; ml rclone
chroms=(2L 2R 3L 3R X)

## TRANSFER FASTQ FILES FROM GDRIVE TO SHERLOCK SCRATCH
if [ $doGETFASTQS -gt 0 ]; then
	rclone copy $gdrive/01_fastqs $rootDir/raw --include=*.fastq.gz
fi

#  RENAME
if [ $doRENAME -gt 0 ]; then
	for f in $rootDir/raw/*.gz; do
		ID=$(basename $f | cut -f2,3 -d'-' | cut -f1 -d'_');
		samp=$(cat $sampKey | awk -v id="$ID" 'match($1, id"$")' | cut -f2);
		if [ $samp ]; then		
			new=$samp"_"$(basename $f | cut -f4,5 -d'_'); 
			mv $f $rootDir/raw/$new
		else echo "no matching sample for $(basename $f) $ID"		
		fi		
	done
fi


# GET SAMPS
samples=$(ls $rootDir/raw | grep gz$ | cut -f1 -d'R' | sort | uniq )
for samp in $samples; do
	echo $samp ; 
	depend=none


## FASTQC
if [ $doFASTQC -gt 0 ]; then
	fastqcDir=$PI_HOME/SOFTWARE/bin/
	cmd="$fastqcDir/fastqc --threads 2 $rootDir/raw/${samp}R1_001.fastq.gz   $rootDir/raw/${samp}R2_001.fastq.gz"
	$create_batch_script "$cmd" $samp fastqc  1:00:00 2000 1 1 $jobscriptDir	
fi


## TRIM
if [ $doTRIM -gt 0 ]; then
	echo "trimming"
	if [ ! -e $rootDir/trimmed ]; then mkdir $rootDir/trimmed; fi
	R1=$rootDir/raw/${samp}R1_001.fastq.gz  
	R2=$rootDir/raw/${samp}R2_001.fastq.gz
	qualMin=20; adapterType=truseq; threads=5; time="2:00:00"; mem=4000
	if [ ! -e $rootDir/trimmed/${samp}R1_001.fastq-trimmed-pair1.fastq.gz ]; then
		cmd="$HOME/scripts/PIPELINE1_fastq_to_bam/02_adapterTrim_skewer.sh $R1 $R2 $rootDir/trimmed $threads $qualMin $adapterType" ; 
		#echo $cmd
		depend=$($create_batch_script "$cmd" $samp trim $time $mem 1 $threads $depend $jobscriptDir) 
	fi
fi

## MERGE READS
if [ $doMERGE_READS -gt 0 ]; then
	if [ ! -e $rootDir/merged ]; then mkdir $rootDir/merged; fi
	threads=5; time="5:00:00"; mem=6000
	if [ ! -e $rootDir/merged/${samp}.assembled.fastq.gz ]; then
		echo "merging paired end reads"
		R1=$rootDir/trimmed/${samp}R1_001.fastq-trimmed-pair1.fastq.gz
		R2=$rootDir/trimmed/${samp}R1_001.fastq-trimmed-pair2.fastq.gz
		cmd="$HOME/scripts/PIPELINE1_fastq_to_bam/01_pairedEndMerge.sh $R1 $R2 $rootDir/merged/$samp $threads"
		#echo $cmd
		depend=$($create_batch_script "$cmd" $samp merge  $time $mem 1 $threads $depend $jobscriptDir)
	fi
fi

# MAP
if [ $doMAP -gt 0 ]; then
	if [ ! -e $rootDir/mapped ]; then mkdir $rootDir/mapped; fi
	if [ ! -e $rootDir/mapped/${samp}.bam ]; then
		echo "mapping"
		time="04:00:00"; mem=12000; threads=10    
		cmd="$HOME/scripts/PIPELINE1_fastq_to_bam/03_map_paired_end_merged.sh  $samp $rootDir/merged $rootDir/mapped $threads"; 
		#echo $cmd
		depend=$($create_batch_script "$cmd" $samp map $time $mem 1 $threads $depend $jobscriptDir)
	fi
fi

## DEDUP
if [ $doDEDUP -gt 0 ]; then
	bam=$rootDir//mapped/$samp.bam
	if [ ! -e $rootDir//mapped/$samp.dedup.bam ]; then
		echo "deduplicating"
		mem=16000; time="02:00:00"; threads=4
		cmd="$HOME/scripts/PIPELINE1_fastq_to_bam/04_dedup.sh $bam" 	
		depend=$($create_batch_script "$cmd" $samp dedup  $time $mem 1 $threads $depend $jobscriptDir)
	fi
fi

## RE-ALIGN
if [ $doREALIGN -gt 0 ]; then
	echo "re-aligning"
	# make intervals file
	intervals=$SCRATCH/orchard2017_microbiome_fastqs/allSamps.intervals
	bam=$rootDir//mapped/$samp.dedup.bam
	if [ ! -e $rootDir//mapped/$samp.dedup.realign.bam ]; then
		cmd="$HOME/scripts/PIPELINE1_fastq_to_bam/05_realign.sh $bam $intervals"; 
		#echo $cmd
		time="1:00:00"; mem=16000; threads=10
		depend=$($create_batch_script "$cmd" $samp realign  $time $mem 1 $threads $depend $jobscriptDir)	
	fi	
fi

## TRANSFER TO TEAM DRIVE
if [ $doBACKUP -gt 0 ]; then
	echo "copying bam to team drive"
	cmd="rclone copy $rootDir/mapped $gdrive/bams --include=${samp}.dedup.realign.bam"
	time="00:30:00"; mem=2000; threads=1
	depend=$($create_batch_script "$cmd" $samp backup  $time $mem 1 $threads $depend $jobscriptDir)	
fi

# GET COVERAGE
if [ $doCOV -gt 0 ]; then
	bam=$rootDir/mapped/$samp.dedup.realign.bam
	bed=$PI_SCRATCH/REFERENCE/dgrp/snps/all_lines/freeze2.allSNPS.bed.2
	#note that this is the list of segregating DGRP sites - not necessarily the sites that vary in this dataset, but a good sampling to get a sense of coverage before calling actual snps
	echo "calculating coverage"
	cmd="$HOME/scripts/PIPELINE1_fastq_to_bam/06_getCov.sh $bam $bed; "	
	mem=12000; time="1:00:00"; threads=3
	$create_batch_script "$cmd" $samp cov  $time $mem 1 $threads $depend $jobscriptDir
fi

# GET INSERT SIZE DISTRIBUTION
if [ $doINSERTSIZE -gt 0 ]; then
	echo "getting insert size distribution"
	ml java; ml R
	for samp in $samples; do
		bam=$rootDir//mapped_dmel_and_bacteria/$samp.dedup.bam
		if [ -e $rootDir//mapped/$samp.dedup.bam ] && [ ! -e $rootDir//mapped/$bam.insertMetrics.txt ]; then
			echo $bam
			cmd="$HOME/scripts/PIPELINE1_fastq_to_bam/07_getInsertSize.sh $bam"		
			mem=8000; time="2:00:00"		
			$create_batch_script "$cmd" $samp insertsize  $time $mem 1 1 $jobscriptDir
			
		fi	
	done
fi
done
