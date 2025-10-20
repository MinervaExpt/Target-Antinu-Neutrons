#!/bin/bash

indir=$1
outName=$2
fileNum=0

rm $outName\_*.txt

for file in `ls -d $indir/*.root`;
do
    echo $file | awk -F"pnfs" '{print "root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr"$2}' > $outName\_$fileNum.txt
    fileNum=$((fileNum+1))
done
