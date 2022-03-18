#Detection microsites deletion or amplification-----

pathCov <- "/home/simone/Scrivania/Scrivania/NewBeerStrains/OS1364/CoverageSeptAll/ChrCov"


myChrSc <- c("S288cchrI","S288cchrII","S288cchrIII","S288cchrIV","S288cchrV","S288cchrVI","S288cchrVII","S288cchrVIII",
             "S288cchrIX","S288cchrX","S288cchrXI","S288cchrXII","S288cchrXIII","S288cchrXIV","S288cchrXV","S288cchrXVI")

myChr <- c("I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV","XVI")

RTG <- c("OS1364")

allSCmeanAll <- c()

library(data.table)

for (indRTG in RTG) {
  n <-1
  medianMat <- c()
  covFileSCmean <- c()
  allSCmean <- c()
  medianCov <- c()
  #ho due for per Sc e Se e in ogni loop ho un altro for per fare la media ogni 10000pb
  for (i in myChrSc) {
    
    covFileSC <- fread(file=paste(pathCov, paste( indRTG, i, ".cov.txt", sep=""), sep="/"))
    indexCovSC <- seq(from=10000, to=nrow(covFileSC), by=10000)
    medianCov <- rbind(medianCov, covFileSC)
    #loop for 100bin 
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
allSCmeanAll <- rbind(allSCmeanAll, cbind(allSCmean, rep.int(indRTG, times = nrow(allSCmean))))
}

# select windows with potential amp/del
# log2 0.4 variation is more or less then 3/4 copies
allSCmean_loss <- allSCmeanAll[which(allSCmeanAll$CovN < -0.58),]
allSCmean_gain <- allSCmeanAll[which(allSCmeanAll$CovN > 0.4),]

allSCmean_loss <- cbind(allSCmean_loss$Pos-10000, allSCmean_loss)
allSCmean_gain <- cbind(allSCmean_gain$Pos-10000, allSCmean_gain)

colnames(allSCmean_loss) <- c("start","end","CovN","Chr","Sample")
colnames(allSCmean_gain) <- c("start","end","CovN","Chr","Sample")

#cycle through the windows with 100bp

subset_loss <- c()
allSCmeanAll2 <- c()
for (indRTG in RTG) {
  n <-1
  medianMat <- c()
  covFileSCmean <- c()
  allSCmean <- c()
  medianCov <- c()
  #ho due for per Sc e Se e in ogni loop ho un altro for per fare la media ogni 10000pb
  for (i in myChrSc) {
    
    covFileSC <- read.table(file=paste(pathCov, paste( indRTG, i, ".cov.txt", sep=""), sep="/"))
    indexCovSC <- seq(from=1000, to=nrow(covFileSC), by=1000)
    medianCov <- rbind(medianCov, covFileSC)
    #loop for 100bin 
    for(indLowSc in indexCovSC){
      covFileSCmean <- rbind(covFileSCmean,cbind((covFileSC[indLowSc ,]$V2), mean(covFileSC[(indLowSc-999):indLowSc ,]$V3)))
      
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
allSCmeanAll2 <- rbind(allSCmeanAll2, cbind(allSCmean, rep.int(indRTG, times = nrow(allSCmean))))
}
## Subselect interval of 1000 bp
colnames(allSCmeanAll2) <- c("Pos","CovN","Chr","Sample")
sub_loss_1000 <- c()

for (indRow in 1:nrow(allSCmean_loss)) {
  
  subregions_1000 <- allSCmeanAll2[which(allSCmeanAll2$Pos >= allSCmean_loss$start[indRow] & 
                                           allSCmeanAll2$Pos <= allSCmean_loss$end[indRow] &
                                           allSCmeanAll2$Chr == allSCmean_loss$Chr[indRow] &
                                           allSCmeanAll2$Sample== allSCmean_loss$Sample[indRow]),]
  
  subregions_1000 <- subregions_1000[which(subregions_1000$CovN < -0.4),]
  
  sub_loss_1000 <- rbind(sub_loss_1000, subregions_1000)
  
}

#gain
sub_gain_1000 <- c()

for (indRow in 1:nrow(allSCmean_gain)) {
  
  subregions_1000 <- allSCmeanAll2[which(allSCmeanAll2$Pos >= allSCmean_gain$start[indRow] & 
                                           allSCmeanAll2$Pos <= allSCmean_gain$end[indRow] &
                                           allSCmeanAll2$Chr == allSCmean_gain$Chr[indRow]&
                                           allSCmeanAll2$Sample== allSCmean_loss$Sample[indRow]),]
  
  subregions_1000 <- subregions_1000[which(
    subregions_1000$CovN > 0.4),]
  
  sub_gain_1000 <- rbind(sub_gain_1000, subregions_1000)
  
}


#Filter for subtelomeric regions--------

subtelL <- c(31606, 9957, 11686, 22362, 26495, 33994, 10876, 19197,
             30617, 31885, 3107, 20776, 9849, 11117, 34288, 26457)
subtelR <- c(197027, 802713, 314566, 1528575, 566584, 266753, 1069643, 
             540634, 425000, 733194, 642915, 1061277, 921042, 765393, 1070680, 924020)

myChr <- c("I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV","XVI")

subtel_Regions_gain <- c()
subtel_Regions_loss <- c()
subcore_Regions_loss <- c()
subcore_Regions_gain <- c()

for (indRowLoss in 1:length(subtelL)) {
  
  subgain_1000_subtel <- sub_gain_1000[which((sub_gain_1000$Pos >= subtelR[indRowLoss] | 
                                                sub_gain_1000$Pos <= subtelL[indRowLoss]) &
                                               sub_gain_1000$Chr == myChr[indRowLoss]),]
  subgain_1000_core <- sub_gain_1000[which((sub_gain_1000$Pos < subtelR[indRowLoss] & 
                                              sub_gain_1000$Pos > subtelL[indRowLoss]) &
                                             sub_gain_1000$Chr == myChr[indRowLoss]),]
  
  subloss_1000_subtel <- sub_loss_1000[which((sub_loss_1000$Pos >= subtelR[indRowLoss] | 
                                                sub_loss_1000$Pos <= subtelL[indRowLoss]) &
                                               sub_loss_1000$Chr == myChr[indRowLoss]),]
  subloss_1000_core <- sub_loss_1000[which((sub_loss_1000$Pos < subtelR[indRowLoss] & 
                                              sub_loss_1000$Pos > subtelL[indRowLoss]) &
                                             sub_loss_1000$Chr == myChr[indRowLoss]),]
  
  subtel_Regions_gain <- rbind(subtel_Regions_gain, subgain_1000_subtel)
  subtel_Regions_loss <- rbind(subtel_Regions_loss, subloss_1000_subtel)
  subcore_Regions_loss <- rbind(subcore_Regions_loss, subloss_1000_core)
  subcore_Regions_gain <- rbind(subcore_Regions_gain, subgain_1000_core)
  
  
}

subtel_Regions_gain <- cbind(subtel_Regions_gain, rep.int("SubtelGain", times = nrow(subtel_Regions_gain)))
subtel_Regions_loss <- cbind(subtel_Regions_loss, rep.int("SubtelLoss", times = nrow(subtel_Regions_loss))) 
subcore_Regions_gain <- cbind(subcore_Regions_gain, rep.int("SubcoreGain", times = nrow(subcore_Regions_gain)))
subcore_Regions_loss <- cbind(subcore_Regions_loss, rep.int("SubcoreLoss", times = nrow(subcore_Regions_loss)))

colnames(subtel_Regions_gain) <- c("Pos","CovN","Chr","Sample","Status")
colnames(subtel_Regions_loss) <- c("Pos","CovN","Chr","Sample","Status")
colnames(subcore_Regions_gain) <- c("Pos","CovN","Chr","Sample","Status")
colnames(subcore_Regions_loss) <- c("Pos","CovN","Chr","Sample","Status")

allCNV1364RTG <- rbind(subtel_Regions_gain, subtel_Regions_loss, subcore_Regions_gain, subcore_Regions_loss)

write.table(allCNV1364RTG, "/home/simone/Scrivania/Scrivania/NewBeerStrains/FinalFilesPaper/Figure1/CNVOS1364RTG.csv")

library(ggplot2)


truelengthNA <- rev(c(954457,1091343,777615,930506,1075542,666862,751611,440036,
                      581049,1091538,271539,583092,1566853,341580,813597,219929))

myChr <- rev(c("XVI","XV","XIV","XIII","XII","XI","X",
               "IX", "VIII", "VII", "VI", "V", "IV", "III",
               "II", "I"))




# plot white region and all the chromosome per samples
for(indRTG in RTG) {
  
  allCNV1364 <- allCNV1364RTG[which(allCNV1364RTG$Sample == indRTG),]
pdf(file=paste("/home/simone/Scrivania/Scrivania/NewBeerStrains/OS1364/CNV", paste(indRTG, ".pdf", sep=""), sep="/"), height = 6, width = 12)
plot(c(0,0),c(0,0),xlim=c(0,sum(truelengthNA)),ylim=c(0,50), yaxt = 'n', xaxt='n',ann=FALSE, bty = 'n')
n <-0


cLL <- c()
cUL <- c()
cUR <- c()
cLR <- c()
chrCum <- 0
for (chr in truelengthNA) {
  
  cLL <- 0
  cUL <- 50
  cUR <- 50
  cLR <- 0
  
  polygon(c(0+chrCum,0+chrCum,chr+chrCum,chr+chrCum), c(cLL,cUL,cUR,cLR), col="white", border = "black")
  chrCum <- chrCum+chr
}


n <- 0

cLL <- c()
cUL <- c()
cUR <- c()
cLR <- c()
n <-1
cumChr <- 0
for (indChr in myChr) {
  
  subCNV <- allCNV1364[which(allCNV1364$Chr == indChr),]
  subCNVloss <- subCNV[which(subCNV$Status == "SubtelLoss" | subCNV$Status == "SubcoreLoss"),]
  subCNVgain <- subCNV[which(subCNV$Status == "SubtelGain" | subCNV$Status == "SubcoreGain"),]
  
  
  cLL <- 0
  cUL <- 25
  cUR <- 25
  cLR <- 0
  #Qui controllo che ci siano LOH per quel chr dato un campione, se non c'è nulla non plotto nulla
  if (nrow(subCNVloss) > 1) {
    for ( indexCNV in (1:nrow(subCNVloss))){
      polygon(c(subCNVloss$Pos[indexCNV] - 1000+cumChr,
                subCNVloss$Pos[indexCNV] - 1000+cumChr,
                subCNVloss$Pos[indexCNV]+cumChr,
                subCNVloss$Pos[indexCNV]+cumChr), 
              c(cLL,cUL,cUR,cLR), col="steelblue", border=NA)
    } }
  
  if(nrow(subCNVgain) > 1) {
    for ( indexCNV in (1:nrow(subCNVgain))){
      polygon(c(subCNVgain$Pos[indexCNV] - 1000+cumChr,
                subCNVgain$Pos[indexCNV] - 1000+cumChr,
                subCNVgain$Pos[indexCNV]+cumChr,
                subCNVgain$Pos[indexCNV]+cumChr), 
              c(cLL+25,cUL+25,cUR+25,cLR+25), col="tomato2", border=NA)
    } 
  }
  cumChr <- cumChr+truelengthNA[n]
  n <- n+1
}
#centMA <- c(139082, 249895, 105556, 448555, 152546, 156803, 370013, 97364, 345199,
#            802499, 453762, 144257, 251046, 605144, 317642, 538828)
#centNA <- c(148940,232532,111477,451488,152858,101306, 490686,94814,
#             351037,427100,440359,151953, 250657,611292,327894,545074)

# from 16  to 1
centNA <- rev(c(545074,327894,611292,250657,151953,440359,427100,351037,
                94814,490686,101306,152858,451488,111477,232532,148940 ))

cLL <- c()
cUL <- c()
cUR <- c()
cLR <- c()
cumChr <- 0
n <- 1
for (indCent in centNA) {
  
  segments(centNA[n]+cumChr,
           0,
           centNA[n]+cumChr,
           50, lty = 2)
  
  cumChr <- cumChr+truelengthNA[n]
  n <- n+1
}


dev.off()

#pdf("/home/simone/Scrivania/Scrivania/NewBeerStrains/OS1364/CNV/partitionCNV_OS1364.pdf")
#ggplot(allCNV1364, aes(x=Status))+geom_bar(Stat="count")+theme_bw()
#dev.off()
}
