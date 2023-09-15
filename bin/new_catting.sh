#!/bin/bash
#run in processed directory
#first part finds each filename and puts it into $f
find $1 -type d | while read f ;
do

f2=$(echo "$f" | awk -F/ '{print $3}');
cat $f/*.fastq.gz > $f/$f2.fq.gz;

#this will remove the original fastq file
rm $f/*.fastq.gz;

printf "\n"; done
