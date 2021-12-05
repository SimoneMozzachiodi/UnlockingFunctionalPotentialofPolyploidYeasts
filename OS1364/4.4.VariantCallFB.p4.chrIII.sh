#!/bin/bash

MyRef="S288c"
sample=#here add all the samples

genomes="/home/smozzachiodi/Batch2RTG-OS1364/index"
InpDir="/home/smozzachiodi/Batch2RTG-OS1364/mappResuStnd/Deduplicated"
OutDir="/home/smozzachiodi/Batch2RTG-OS1364/mappResuStnd"
wd="/home/smozzachiodi/Batch2RTG-OS1364"

for ref in $MyRef
do
   for indsample in $sample
    do
  /opt/tvs/bin/freebayes -f $genomes/$ref.genome.chrref.fa -@ $wd/"OS1364S288c.Filtvarp4Pos.vcf.gz" -p4 -m 30 -q 20 -i -X -u  $InpDir/$indsample/$indsample.srt.rmd.bam > $OutDir/$indsample/$indsample$ref/$indsample$ref.var.p4.2.vcf 

   done

done

# m:mapping quality, q:quality variant, -i:no indels, -X:no mnps, -u:no complex
