#!/bin/bash
wd="/home/smozzachiodi/Batch2RTG-OS1364"
myRef="S288c"
mySample="OS1364"
###Ploidy p3
ploidy="p3"#then add p4

for indRef in $myRef
 do 
 for indSample in $mySample
 do
  for indploidy in $ploidy
   do
    bcftools filter -i 'TYPE="snp" && QUAL>20 && FORMAT/DP>10' $wd/$indSample$indRef".var.$indploidy.vcf" > $wd/$indSample$indRef".Filtvar$indploidy.vcf"
# | grep -w "TYPE=snp" > $wd/$indSample$indRef".Filtvar$indploidy.vcf"
    done

   done
  done
