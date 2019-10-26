#!/bin/bash

#### perform paired end read merging

R1=${1}
R2=${2}
outfile=${3}
threads=${4:-5}
standardBaseQual=${5:-33}
	
	echo "##running pear..."
	pear-0.9.6-bin-64 -b $standardBaseQual -j $threads -f $R1 -r $R2 -o $outfile 
	echo "##zipping..."
	gzip $outfile.*
	echo "##finished "$(basename $outfile)

	
