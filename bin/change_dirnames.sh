#!/bin/bash
# Find dirnames and put it into $f
find $1 -type d | while read f ;
do

f2=$(echo "$f" | awk -F- '{print $1}');
mv $f $f2;

printf "\n"; 
done