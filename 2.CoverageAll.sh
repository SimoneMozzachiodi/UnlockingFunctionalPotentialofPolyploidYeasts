#!/bin/bash
#Coverage for chromosome

WorkDir="/home/smozzachiodi/Batch2RTG-OS1364/mappResuStnd/Deduplicated"
MySample=#here add all samples
allRef="S288c"
myChr="chrI chrII chrIII chrIV chrV chrVI chrVII chrVIII chrIX chrX chrXI chrXII chrXIII chrXIV chrXV chrXVI"
OutDir="/home/smozzachiodi/Batch2RTG-OS1364/mappResuStnd/Deduplicated"
chrDir="ChrCov"

for myRef in $allRef
 do
 for sample in $MySample

 do 
 OutDirSam=$OutDir/$sample
  #mkdir OutDirSam
   OutDirCov=$OutDir/$sample/$chrDir
    mkdir $OutDirCov

   for i in $myChr 
    do

    cat $OutDir/$sample/$sample.srt.rmd.bam.dat | grep -w $i > $OutDirCov/$sample$myRef$i.cov.txt

    done


 done
done
