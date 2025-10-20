#!/bin/bash

file=$1
outName=$2
max=$3

counter=1
fileNum=0
total=0

while IFS= read -r line;
do
    if [ $counter -eq 1 ];
    then
	echo $line > $outName\_$fileNum.txt
	counter=$((counter+1))
    elif [ $counter -eq $max ];
    then
	counter=1
	echo $line >> $outName\_$fileNum.txt
	fileNum=$((fileNum+1))
    else
	echo $line >> $outName\_$fileNum.txt
	counter=$((counter+1))
    fi
    total=$((total+1))
done < $file
