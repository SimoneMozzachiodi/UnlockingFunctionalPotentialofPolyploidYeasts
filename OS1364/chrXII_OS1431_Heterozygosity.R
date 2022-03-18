#### check heterozygosity OS1431 chromosome XII ####

pathVcf <- "/home/simone/Scrivania/Scrivania/NewBeerStrains/OS1431/OS1431.HTmarkers.p4.vcf"

library(vcfR)

vcfOS1431 <- read.table(pathVcf)

chrXII_OS1431 <- vcfOS1431[which(vcfOS1431$V1 == "chrXII"),]

AFos1431_chrXII <- sapply(strsplit(as.character(chrXII_OS1431$V8), split = c(";")), "[[", 4)
posOS1431_chrXII <- vcfOS1431[which(vcfOS1431$V1 == "chrXII"), 2]

plotAF <- as.numeric(as.character(sapply(strsplit(as.character(AFos1431_chrXII), split = c("=")), "[[", 2)))


plotAFpos_chrXII <- data.frame(cbind(posOS1431_chrXII, plotAF))

colnames(plotAFpos_chrXII) <- c("Position", "AF")

library(ggplot2)

ggplot(plotAFpos_chrXII, aes(x=Position, y=AF))+geom_point()+theme_bw()+xlim(340000, 900000)+ylim(1,0)+
  geom_vline(xintercept = 482000)+geom_vline(xintercept = 628000)+geom_vline(xintercept = 693000)

