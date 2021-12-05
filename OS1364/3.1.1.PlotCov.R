#!/usr/bin/env Rscript
#nedded directories
pathCov46 <- "/home/smozzachiodi/Batch2RTG-OS1364/mappResuStnd/Deduplicated"

myWcol <- c("D415R86")
#"D167R74","D167R75","D167R76","D167R77","D167R78","D167R79","D167R80","D287R143","D287R144","D287R145","D287R146","D287R147","D287R148","D287R149","D287R150",
#            "D287R151","D287R152","D287R153","D287R154","D287R155","D287R156","D287R157","D287R158",
#            "D287R159","D287R160","D287R161","D287R162","D287R163","D287R164","D287R165","D287R166",
#            "D287R169","D287R170","D287R097","D287R098","D287R099","D287R100","D287R101","D287R102",
#            "D287R104","D287R112","D287R113","D287R114","D287R115","D287R116","D287R117","D287R118",
#            "D287R119","D287R120","D287R121","D287R122","D287R123","D287R124","D287R125","D287R126",
#            "D287R127","D287R128","D287R129","D287R130","D287R131","D287R132","D287R133","D287R134",
#            "D287R135","D287R136","D287R137","D287R138","D287R139","D287R140","D287R141","D287R142",
#             "D415R85","D415R88","D415R89","D415R92","D415R91","D415R93","D415R94","D415R96")


#compmapp-----
myChrSc <- c("S288cchrI","S288cchrII","S288cchrIII","S288cchrIV","S288cchrV","S288cchrVI","S288cchrVII","S288cchrVIII","S288cchrIX","S288cchrX",
             "S288cchrXI","S288cchrXII","S288cchrXIII","S288cchrXIV","S288cchrXV","S288cchrXVI")

myChr <- c("I","II","III","IV","V","VI","VII","VIII","IX","X",
           "XI","XII","XIII","XIV","XV","XVI")

StatAllSampl <- c()

for (indexW in myWcol) {
  n <-1
  medianMat <- c()
  covFileSCmean <- c()
  allSCmean <- c()
  allSampCov <- c()
  #ho due for per Sc e Se e in ogni loop ho un altro for per fare la media ogni 1000pb
  for (i in myChrSc) {
    
    covFileSC <- read.table(file=paste(pathCov46, indexW, "ChrCov", paste(i, "cov.txt", sep="."), sep="/"))
    indexCovSC <- seq(from=10000, to=nrow(covFileSC), by=10000)
    
    #loop for 100bin 
    for(indLowSc in indexCovSC){
      covFileSCmean <- rbind(covFileSCmean,cbind((covFileSC[indLowSc ,]$V2), mean(covFileSC[(indLowSc-9999):indLowSc ,]$V3)))
      
    }
    covFileSCmean <- as.data.frame(covFileSCmean)
    #added after
    allSCmean <- rbind(allSCmean, cbind(covFileSCmean, rep(x=myChr[n], time=nrow(covFileSCmean)) ) )
    
    allSampCov <- rbind(allSampCov,matrix(covFileSC$V3, ncol=1))  
    covFileSCmean <- c()
    covFileSEmean <- c()
    n <- n+1
    
  }
  colnames(allSCmean) <- c("Pos","Cov","Chr")
  
  library(ggplot2)
  
  #Statistics whole sample----
  
  StatAllSampl <- rbind(StatAllSampl, cbind(mean(allSampCov[, 1]), sd(allSampCov[, 1])) )   
  
  p <- ggplot(allSCmean, aes(x=Pos, y=Cov, fill=Chr) )+ geom_point(size=1, col="steelblue")+facet_grid(~ Chr,scales = "free_x",space = "free_x")+theme(legend.position="none")+theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+guides(fill=guide_legend(nrow=1,byrow=TRUE))+ylab("Coverage")+xlab("Genomic Coordinates")+ylim(0,400) + theme(axis.text=element_text(size=12,face="bold"),axis.title=element_text(size=14,face="bold"))+ theme(strip.text.x = element_text(size=9, face = "bold"))
  
  
  a <- p + theme(panel.background = element_rect(fill = "white" ,colour = "black",size = 0.5, linetype = "solid"   ),  panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "white"), panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "white"))
  ggsave(filename=paste(indexW,"coverage.pdf", sep="."),device = "pdf" ,path=paste(pathCov46, "CoverageResults", sep="/" ), plot=a, width = 8, height = 4) 
  
  
}

colnames(StatAllSampl) <- c("Mean_Coverage", "SD_Coverage")

write.table(StatAllSampl, file="/home/smozzachiodi/Batch2RTG-OS1364/mappResuStnd/Deduplicated/Statistcs_coverage.csv", row.names = T, col.names = T)


