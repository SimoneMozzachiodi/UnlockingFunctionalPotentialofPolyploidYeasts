
# check region CNV----

pathCov <- "/home/simone/Scrivania/Scrivania/NewBeerStrains/OS1431/CovLys2Ura"

#covT0 <- read.csv(pathCov, sep="\t")

#calculate median chromosome

#medianCov <- median(covT0$X2)

myChrSc <- c("S288cchrI","S288cchrII","S288cchrIII","S288cchrIV","S288cchrV","S288cchrVI","S288cchrVII","S288cchrVIII",
             "S288cchrIX","S288cchrX","S288cchrXI","S288cchrXII","S288cchrXIII","S288cchrXIV","S288cchrXV","S288cchrXVI")

myChr <- c("I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV","XVI")

RTG <- c("D287R103")

library(data.table)

for (indRTG in RTG) {
  n <-1
  medianMat <- c()
  covFileSCmean <- c()
  allSCmean <- c()
  medianCov <- c()
  #loop 10000 bp window
  for (i in myChrSc) {
    
    covFileSC <- fread(file=paste(pathCov, indRTG, "ChrCov", paste(i, "cov.txt", sep="."), sep="/"))
    indexCovSC <- seq(from=10000, to=nrow(covFileSC), by=10000)
    medianCov <- rbind(medianCov, covFileSC)
    #loop for bin 
    for(indLowSc in indexCovSC){
      covFileSCmean <- rbind(covFileSCmean,cbind((covFileSC[indLowSc ,]$V2), mean(covFileSC[(indLowSc-9999):indLowSc ,]$V3)))
      
    }
    covFileSCmean <- as.data.frame(covFileSCmean)
    #added after
    allSCmean <- rbind(allSCmean, cbind(covFileSCmean, rep(x=myChr[n], time=nrow(covFileSCmean)) ) )
    
    covFileSCmean <- c()
    covFileSEmean <- c()
    n <- n+1
    
  }
  
  medianCovall <- median(medianCov$V3)
  colnames(allSCmean) <- c("Pos","CovN","Chr")
  allSCmean$CovN <- log2(allSCmean$CovN/medianCovall)
  
  library(tidyverse)
  chr_median <- allSCmean %>% 
    group_by(Chr) %>% 
    summarize(median_val = median(CovN))
  #print(chr_median)
  
  library(ggplot2)
  
  
  a <-  ggplot(allSCmean, aes(x=Pos, y=CovN, fill=Chr) )+ geom_point(size=1, col="purple1")+facet_grid(~ Chr,scales = "free_x",space = "free_x")+
    theme(legend.position="none")+theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
    guides(fill=guide_legend(nrow=1,byrow=TRUE))+ylab("Coverage")+theme_bw()+ylim(-2,2)+
    xlab("Genomic Coordinates") + geom_hline(data= chr_median, aes(yintercept=median_val))+
    theme(axis.text=element_blank(),axis.title=element_text(size=14,face="bold"))+ 
    theme(strip.text.x = element_blank())+
    theme(legend.position="none")
  
  
  #dir.create(file.path(paste(pathCov, sep="/" )))
  ggsave(filename=paste(indRTG,"coverage.pdf", sep="."),device = "pdf" ,path=paste(pathCov,sep="/" ), plot=a, width=8, height=4) 
  chr_median <- c()}