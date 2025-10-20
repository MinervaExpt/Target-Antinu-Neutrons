#!/bin/bash

file=$1
outName=$2
fileNum=0

while IFS= read -r line;
do
    echo $line > $outName\_$fileNum.txt
    fileNum=$((fileNum+1))
done < $file
