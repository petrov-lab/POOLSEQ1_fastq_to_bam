#!/bin/bash

bam=$1
bed=$2
samp=$(basename $bam | sed 's/.bam//')

#coverageBed -abam $bam -b $bed -mean | grep -v "\-1" | cut -f1-3 | \
#coverageBed -a $bed -b - -mean | \
#awk -v samp="$samp" '{
#	print samp"\t"$0
#}' > ${bam}.cov

#cat $bam.cov | awk -v samp="$samp" '
#{ cov[$2]=cov[$2]+$5; sites[$2]=sites[$2]+1} 
#END { print samp"\t"cov["2L"]/sites["2L"]"\t"cov["2R"]/sites["2R"]"\t"cov["3L"]/sites["3L"]"\t"cov["3R"]/sites["3R"]"\t"cov["X"]/sites["X"]}
#' > $bam.cov.chromMeans

samtools depth -b $bed $bam > ${bam}.cov
cat $bam.cov | awk -v samp="$samp" '
{ cov[$1]=cov[$1]+$3; sites[$1]++} 
END { print samp"\t"cov["2L"]/sites["2L"]"\t"cov["2R"]/sites["2R"]"\t"cov["3L"]/sites["3L"]"\t"cov["3R"]/sites["3R"]"\t"cov["X"]/sites["X"]}
' > $bam.cov.chromMeans