#!/bin/bash

# Initialize variables
ffull=$1;
fbase=$(basename -- "$FILE");
f2=$(echo "$f" | awk -F/ '{print $3}');

# For forward read
## include reads that are 2nd in a pair (128);
## exclude reads that are mapped to the reverse strand (16)
samtools view -b -f 128 -F 16 $ffull > a.fwd1.bam

## exclude reads that are mapped to the reverse strand (16) and
## first in a pair (64): 64 + 16 = 80
samtools view -b -f 80 $ffull > $fbase.fwd2.bam

## combine the temporary files
samtools merge -f fwd.bam a.fwd1.bam a.fwd2.bam

## index the filtered BAM file
samtools index fwd.bam

## run bamCoverage
bamCoverage -b fwd.bam -o a.fwd.bw

## remove the temporary files
rm a.fwd*.bam

# For reverse read
## include reads that map to the reverse strand (128)
## and are second in a pair (16): 128 + 16 = 144
samtools view -b -f 144 results/QUO1006522/bwa_aln/wt_39C_1.bam > a.rev1.bam

## include reads that are first in a pair (64), but
## exclude those ones that map to the reverse strand (16)
samtools view -b -f 64 -F 16 results/QUO1006522/bwa_aln/wt_39C_1.bam > a.rev2.bam

## merge the temporary files
samtools merge -f rev.bam a.rev1.bam a.rev2.bam

## index the merged, filtered BAM file
samtools index rev.bam

## run bamCoverage
bamCoverage -b rev.bam -o a.rev.bw

## remove temporary files
rm a.rev*.bam