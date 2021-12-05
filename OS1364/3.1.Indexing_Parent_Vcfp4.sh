#!/bin/bash

wd="/home/smozzachiodi/Batch2RTG-OS1364"

#Parent vcf added below for example, remember to set the correct path
#Generate the vcf file to use for call markers status in RTG derived from the parent

#Filter the original vcf of the parent for quality and depth, plus select only Snp
bcftools filter -i 'TYPE="snp" && QUAL>20 && FORMAT/DP>10' $wd/OS1364S288c.var.p4.vcf > $wd/OS1364S288c.Filtvarp4.vcf 

#select the first 7 fields from the filtered vcf
cut -f -7 $wd/OS1364S288c.Filtvarp4.vcf > $wd/OS1364S288c.Filtvarp4Pos.vcf

#make the gz file con bgzip needed for tabiz index
bgzip -c $wd/OS1364S288c.Filtvarp4Pos.vcf > $wd/OS1364S288c.Filtvarp4Pos.vcf.gz

#make the tabix index
tabix -p vcf $wd/OS1364S288c.Filtvarp4Pos.vcf.gz


#Used the index made in professorx folder BeerRTG
