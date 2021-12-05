WorkDir="/home/smozzachiodi/Batch2RTG-OS1364/RTG-Beer3"
MySample="D287R131 D287R132 D287R133 D287R134 D287R135 D287R136 D287R137 D287R138 D287R139 D287R140 D287R141 D287R142 D287R143
D287R144 D287R145 D287R146 D287R147 D287R148 D287R149 D287R150 D287R151 D287R152 D287R153 D287R154 D287R155 D287R156 D287R157
D287R158 D287R159 D287R160 D287R161 D287R162 D287R163 D287R164 D287R165 D287R166 D287R169 D287R170 D415R92 
D415R93 D415R94 D415R96"

beginrefseq="/home/smozzachiodi/Batch2RTG-OS1364/index"
OutDir="/home/smozzachiodi/Batch2RTG-OS1364/mappResuStnd"
MyRef="S288c"
#Ref S288c
mkdir $OutDir

for ref in $MyRef
do
#ad chrref if you use reference S288C
refseq=$beginrefseq/$ref.genome.chrref.fa

for sample in $MySample

do
OutDirAll=$OutDir/$sample
mkdir $OutDirAll
OutDirSam=$OutDir/$sample/$sample$ref
mkdir $OutDirSam

 
 bwa mem -M -t 20 $refseq $WorkDir/$sample/$sample".R1.fastq.gz" $WorkDir/$sample/$sample".R2.fastq.gz" > $OutDirSam/$sample$ref.sam
 samtools sort --threads 6 -o $OutDirSam/$sample$ref.srt.bam $OutDirSam/$sample$ref.sam
 samtools index $OutDirSam/$sample$ref.srt.bam
 samtools depth -a $OutDirSam/$sample$ref.srt.bam > $OutDirSam/$sample$ref.bam.dat
  rm $OutDirSam/$sample$ref.sam
done

done
