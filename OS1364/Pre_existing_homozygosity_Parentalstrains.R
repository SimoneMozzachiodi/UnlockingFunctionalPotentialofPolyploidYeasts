#Heterozygosity levels genome wide----

pathHet1364_3n <- "/home/simone/Scrivania/Scrivania/NewBeerStrains/FinalFilesPaper/Figure1/OS1364S288c.Filtvarp3.vcf"
pathHet1364_3nplot <- "/home/simone/Scrivania/Scrivania/NewBeerStrains/FinalFilesPaper/Figure1/OS1364S288c.Filtvarp3_FiltHomo.vcf"
pathHet1364_4n <- "/home/simone/Scrivania/Scrivania/NewBeerStrains/FinalFilesPaper/Figure1/OS1364S288c.Filtvarp4.vcf"
pathHet1364_4nplot <- "/home/simone/Scrivania/Scrivania/NewBeerStrains/FinalFilesPaper/Figure1/OS1364S288c.Filtvarp4_FiltHomo.vcf"
pathHet1431_4n <- "/home/simone/Scrivania/Scrivania/NewBeerStrains/FinalFilesPaper/Figure1/OS1431S288c.Filtvarp4.vcf"
pathHet1431_4nplot <- "/home/simone/Scrivania/Scrivania/NewBeerStrains/FinalFilesPaper/Figure1/OS1431S288c.Filtvarp4_FiltHomo.vcf"


Het1364_3n <- read.table(pathHet1364_3n)
Het1364_3nplot <- read.table((pathHet1364_3nplot))
Het1364_4n <- read.table(pathHet1364_4n)
Het1364_4nplot <- read.table((pathHet1364_4nplot))
Het1431_4n <- read.table(pathHet1431_4n)
Het1431_4nplot <- read.table((pathHet1431_4nplot))

#remove chromosme III from 3n 1364-----

Het1364_3n <- Het1364_3n[which(Het1364_3n$V1 != "chrIII"),]
Het1364_3nplot <- Het1364_3nplot[which(Het1364_3nplot$V1 != "chrIII"),]

Het1364_4n <- Het1364_4n[which(Het1364_4n$V1 == "chrIII"),]
Het1364_4nplot <- Het1364_4nplot[which(Het1364_4nplot$V1 == "chrIII"),]

AFin_1364_3n <- table(sapply(strsplit(as.character(Het1364_3n$V8), split = c(";") ), "[[" , 4))

AFin_1364_4n <- table(sapply(strsplit(as.character(Het1364_4n$V8), split = c(";") ), "[[" , 4))

AFin_1431_4n <- table(sapply(strsplit(as.character(Het1431_4n$V8), split = c(";") ), "[[" , 4))



#### Count number of LOH regions ###########

Chrom <- c("chrI","chrII","chrIII","chrIV","chrV","chrVI","chrVII","chrVIII",
           "chrIX","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI")
Start <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
chrScer <- rev(c(954457,1091343,777615,930506,1075542,666862,751611,440036,
                 581049,1091538,271539,583092,1566853,341580,813597,219929))

flag <- 1
allChrom <- c()

mergeHet1364_3nplot <- rbind(Het1364_3nplot, Het1364_4nplot)

for (indChrom in Chrom) {
  
  subVar <- mergeHet1364_3nplot[which(mergeHet1364_3nplot$V1 == indChrom),]
  
  indexCovSC <- seq(from=50000, to=chrScer[flag], by=50000)
  
  #loop for 1000bin 
  for(indLowSc in indexCovSC){
    nMarkersregion <- nrow(subVar[which(subVar$V2 < indLowSc & subVar$V2 > (indLowSc-50000)),])
    allChrom <- rbind(allChrom,(cbind(indChrom, indLowSc, nMarkersregion)))
  }
  
  flag <- flag+1
}

dfAllChr <- data.frame(allChrom)
LOH <- dfAllChr[which(as.numeric(as.character(dfAllChr$nMarkersregion)) < 10),]
cent <- c(148940,232532,111477,451488,152858,101306, 490686,94814,
          351037,427100,440359,151953, 250657,611292,327894,545074)

dfAllCent <- data.frame(cent, Chrom, chrScer)
library(ggplot2)

cust_theme <- theme_bw() + theme(legend.position="none", 
                                 axis.title = element_blank(),
                                 strip.background = element_blank(), panel.margin = unit(0, "lines"), 
                                 panel.border = element_rect(size = 0.25, color = "black"), 
                                 panel.grid = element_blank())

dfAllChr$indChrom <- factor(dfAllChr$indChrom, levels = Chrom)

pdf("/home/simone/Scrivania/Scrivania/NewBeerStrains/FinalFilesPaper/Sfig1/Het_OS1364.pdf", height = 1, width = 4.7)
ggplot(dfAllChr, aes(x=as.numeric(as.character(indLowSc)), y=as.numeric(as.character(nMarkersregion))))+ylab("")+
  geom_area()+theme_bw()+facet_grid(~indChrom,  scales = "free", space='free')+cust_theme
dev.off()  

pdf("/home/simone/Scrivania/Scrivania/NewBeerStrains/FinalFilesPaper/Sfig1/HetCent_OS1364.pdf", height = 1, width = 4.7)
ggplot(dfAllCent, aes(x=chrScer, y=400))+ geom_vline(data = dfAllCent, aes(xintercept = cent))+facet_grid(~Chrom,scales = "fixed", space='free')+theme_bw()+ylab("")+cust_theme
dev.off()


####OS1431#####

#### Count number of LOH regions ###########

Chrom <- c("chrI","chrII","chrIII","chrIV","chrV","chrVI","chrVII","chrVIII",
           "chrIX","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI")
Start <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
chrScer <-rev(c(954457,1091343,777615,930506,1075542,666862,751611,440036,
                581049,1091538,271539,583092,1566853,341580,813597,219929))

flag <- 1
allChrom <- c()

for (indChrom in Chrom) {
  
  subVar <- Het1431_4nplot[which(Het1431_4nplot$V1 == indChrom),]
  
  indexCovSC <- seq(from=50000, to=chrScer[flag], by=50000)
  
  #loop for 1000bin 
  for(indLowSc in indexCovSC){
    nMarkersregion <- nrow(subVar[which(subVar$V2 < indLowSc & subVar$V2 > (indLowSc-50000)),])
    allChrom <- rbind(allChrom,(cbind(indChrom, indLowSc, nMarkersregion)))
  }
  
  flag <- flag+1
}

dfAllChr <- data.frame(allChrom)
LOH <- dfAllChr[which(as.numeric(as.character(dfAllChr$nMarkersregion)) < 10),]
cent <- c(148940,232532,111477,451488,152858,101306, 490686,94814,
          351037,427100,440359,151953, 250657,611292,327894,545074)

dfAllCent <- data.frame(cent, Chrom, chrScer)
library(ggplot2)

cust_theme <- theme_bw() + theme(legend.position="none", 
                                 axis.title = element_blank(),
                                 strip.background = element_blank(), panel.margin = unit(0, "lines"), 
                                 panel.border = element_rect(size = 0.25, color = "black"), 
                                 panel.grid = element_blank())

dfAllChr$indChrom <- factor(dfAllChr$indChrom, levels = Chrom)

pdf("/home/simone/Scrivania/Scrivania/NewBeerStrains/FinalFilesPaper/Sfig1/Het_OS1431.pdf", height = 1, width = 4.7)
ggplot(dfAllChr, aes(x=as.numeric(as.character(indLowSc)), y=as.numeric(as.character(nMarkersregion))))+ylab("")+
  geom_area()+theme_bw()+facet_grid(~indChrom,  scales = "free", space='free')+cust_theme
dev.off()  

pdf("/home/simone/Scrivania/Scrivania/NewBeerStrains/FinalFilesPaper/Sfig1/HetCent_OS1431.pdf", height = 1, width = 4.7)
ggplot(dfAllCent, aes(x=chrScer, y=400))+ geom_vline(data = dfAllCent, aes(xintercept = cent))+facet_grid(~Chrom,scales = "free", space='free')+theme_bw()+ylab("")+cust_theme
dev.off()

