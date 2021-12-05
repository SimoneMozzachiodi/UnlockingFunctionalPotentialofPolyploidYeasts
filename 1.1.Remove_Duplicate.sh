#!/bin/bash

OutDir="Deduplicated"
BamPath="/home/smozzachiodi/Batch2RTG-OS1364/mappResuStnd"
OutPath="/home/smozzachiodi/Batch2RTG-OS1364/mappResuStnd/$OutDir"

mkdir $OutPath

MySample="D536R91 D415R86 D167R80 D167R74 D167R75 D167R76 D167R77 D167R78 D167R79 D178R80 D287R097 D287R098 D287R099 D287R100 D287R101 D287R102 D287R104 D287R112 D287R113 
D287R114 D287R115 D287R116
D287R117 D287R118 D287R119 D287R120 D287R121 D287R122 D287R123 D287R124 D287R125 D287R126 D287R127 D287R128 D287R129
D287R130 D287R131 D287R132 D287R133 D287R134 D287R135 D287R136 D287R137 D287R138 D287R139 D287R140 D287R141 D287R142 D287R143
D287R144 D287R145 D287R146 D287R147 D287R148 D287R149 D287R150 D287R151 D287R152 D287R153 D287R154 D287R155 D287R156 D287R157 D287R158 D287R159 D287R160 D287R161 D287R162 D287R163 D287R164 D287R165 D287R166 D287R169 D287R170 D167R80 D415R85 
D415R86 D415R88 D415R89 D415R90 D415R91 D415R92 D415R93 D415R94 D415R96"


Nthreads=30

for IndSamp in $MySample

do   
   mkdir $OutPath/$IndSamp
    samtools sort -@ 20 -n -o $BamPath/$IndSamp/$IndSamp"S288c"/$IndSamp"S288c.srt.name.bam" $BamPath/$IndSamp/$IndSamp"S288c"/$IndSamp"S288c.srt.bam"
    samtools fixmate -@ 20 -r -p -m -O bam $BamPath/$IndSamp/$IndSamp"S288c"/$IndSamp"S288c.srt.name.bam" $BamPath/$IndSamp/$IndSamp"S288c"/$IndSamp"S288c.srt.fix.bam"
    samtools sort -@ 20 -o $BamPath/$IndSamp/$IndSamp"S288c"/$IndSamp"S288c.srt.fix.srt.bam" $BamPath/$IndSamp/$IndSamp"S288c"/$IndSamp"S288c.srt.fix.bam" 
    samtools markdup -@ 20 -r $BamPath/$IndSamp/$IndSamp"S288c"/$IndSamp"S288c.srt.fix.srt.bam" $OutPath/$IndSamp/$IndSamp.srt.rmd.bam 
    rm $BamPath/$IndSamp/$IndSamp"S288c"/$IndSamp"S288c.srt.name.bam"
    rm $BamPath/$IndSamp/$IndSamp"S288c"/$IndSamp"S288c.srt.fix.bam"
    rm $BamPath/$IndSamp/$IndSamp"S288c"/$IndSamp"S288c.srt.fix.srt.bam"
   samtools index $OutPath/$IndSamp/$IndSamp.srt.rmd.bam
   samtools depth -a $OutPath/$IndSamp/$IndSamp.srt.rmd.bam > $OutPath/$IndSamp/$IndSamp.srt.rmd.bam.dat

done
