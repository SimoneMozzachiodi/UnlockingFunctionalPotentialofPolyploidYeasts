WorkDir="/home/smozzachiodi/Batch2RTG-OS1364/RTG-Beer3"
MySample="D536R91 D287R097 D287R098 D287R099 D287R100 D287R101 D287R102 D287R104 D287R112 D287R113 D287R114 D287R115 D287R116
D287R117 D287R118 D287R119 D287R120 D287R121 D287R122 D287R123 D287R124 D287R125 D287R126 D287R127 D287R128 
D287R129 D287R130 D415R85 D415R88 D415R89 D415R90 D415R91"
beginrefseq="/home/smozzachiodi/Batch2RTG-OS1364/index"
OutDir="/home/smozzachiodi/Batch2RTG-OS1364/mappResuStnd"
MyRef="S288c"

mkdir $OutDir

for ref in $MyRef
do

refseq=$beginrefseq/$ref.genome.chrref.fa

for sample in $MySample

do
OutDirAll=$OutDir/$sample
mkdir $OutDirAll
OutDirSam=$OutDir/$sample/$sample$ref
mkdir $OutDirSam

 
 bwa mem -M -t 11 $refseq $WorkDir/$sample/$sample".R1.fastq.gz" $WorkDir/$sample/$sample".R2.fastq.gz" > $OutDirSam/$sample$ref.sam
 samtools sort --threads 11 -o $OutDirSam/$sample$ref.srt.bam $OutDirSam/$sample$ref.sam
 samtools index $OutDirSam/$sample$ref.srt.bam
 samtools depth -a $OutDirSam/$sample$ref.srt.bam > $OutDirSam/$sample$ref.bam.dat
  rm $OutDirSam/$sample$ref.sam
done

done
