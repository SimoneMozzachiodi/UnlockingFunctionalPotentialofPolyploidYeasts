#Filter markers called as variation in AF but falling in AF shift regions/LOH----
#down to samples markers area designation----

dfLOHmarkpath <- "/home/simone/Scrivania/Scrivania/NewBeerStrains/OS1431/LOHanalysis/LastCallLOH/MarkerStatusOS1431_RTG_5marks_16052020.csv"

dfLOHmark <- read.csv(dfLOHmarkpath)
tableDfLOH <- read.table("/home/simone/Scrivania/Scrivania/NewBeerStrains/OS1431/LOHanalysis/LastCallLOH/Filt_wtLOH_CNV_OS1431.csv")

dfreciprocalF5 <- tableDfLOH
  
mySample <- c("D287R167","D287R168") ###e.g. for RTG 

myChr <- c("chrI","chrII","chrIII","chrIV","chrV","chrVI","chrVII","chrVIII",
           "chrIX","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI")

dfLOHmarkerUpd <- c()

for (indSample in mySample[1]) {
  
  subsetLOHsample <- dfreciprocalF5[which(dfreciprocalF5$Sample == indSample),] 
  subsetdfLOHmark <- dfLOHmark[which(dfLOHmark$indSample == indSample),]
  for (indChr in myChr) {
    #down to chromosome
    subsetLOHchr <- subsetLOHsample[which(subsetLOHsample$Chromosome == indChr),]
    subsetdfLOHmarkChr <- subsetdfLOHmark[which(subsetdfLOHmark$indChr == indChr),]
    
    for (indLOHlocal in 1:nrow(subsetLOHchr)) {
      #down to LOH single
      minimumLOH <- subsetLOHchr[indLOHlocal ,]
      localdfLOHMark <- subsetdfLOHmarkChr[which(subsetdfLOHmarkChr$pos >= minimumLOH$Start &
                                                   subsetdfLOHmarkChr$pos <= minimumLOH$End),]
      
      
      for (indMarkerlocal in 1:nrow(localdfLOHMark)) {
        #find group of consecutive LOH marker shift
        #if first marker-----
        if (indMarkerlocal == 1)
        {
          leftMark <- localdfLOHMark$pos[indMarkerlocal]
          rightMark <- (localdfLOHMark$pos[indMarkerlocal]+localdfLOHMark$pos[indMarkerlocal+1])/2
          dfLOHmarkerUpd <- rbind(dfLOHmarkerUpd, cbind(localdfLOHMark[indMarkerlocal ,], leftMark, rightMark))
        }
        
        #if marker in between-----
        if (indMarkerlocal > 1 & indMarkerlocal < nrow(localdfLOHMark))
        {
          leftMark <- (localdfLOHMark$pos[indMarkerlocal-1]+localdfLOHMark$pos[indMarkerlocal])/2
          rightMark <- (localdfLOHMark$pos[indMarkerlocal]+localdfLOHMark$pos[indMarkerlocal+1])/2
          dfLOHmarkerUpd <- rbind(dfLOHmarkerUpd, cbind(localdfLOHMark[indMarkerlocal ,], leftMark, rightMark))
          
        }
        
        #if marker at the end-----
        if (indMarkerlocal == nrow(localdfLOHMark))
        {
          leftMark <- (localdfLOHMark$pos[indMarkerlocal-1]+localdfLOHMark$pos[indMarkerlocal])/2
          rightMark <- localdfLOHMark$pos[indMarkerlocal]
          dfLOHmarkerUpd <- rbind(dfLOHmarkerUpd, cbind(localdfLOHMark[indMarkerlocal ,], leftMark, rightMark))
          
        }
      } 
      
    }  
  }
  
}




#Plotting LOH regions-------

truelengthNA <- c(954457,1091343,777615,930506,1075542,666862,751611,440036,
                  581049,1091538,271539,583092,1566853,341580,813597,219929)

myChr <- c("chrXVI","chrXV","chrXIV","chrXIII","chrXII","chrXI","chrX",
           "chrIX", "chrVIII", "chrVII", "chrVI", "chrV", "chrIV", "chrIII",
           "chrII", "chrI")

library(berryFunctions)

### Load marker status and then filter it on the new LOH regions-----
for (indSample in mySample) {
  
  n <- 0
  
  cLL <- c()
  cUL <- c()
  cUR <- c()
  cLR <- c()
  
  # plot white region and all the chromosome per samples
  pdf(paste("/home/simone/Scrivania/Scrivania/NewBeerStrains/OS1431/LOHanalysis/LastCallLOH/OS1431RTGWT/",paste(indSample, "invCNV.pdf",sep=""), sep ="/"), height = 12, width = 12)
  plot(c(0,0),c(0,0),xlim=c(0,max(truelengthNA)),ylim=c(0,1700), yaxt = 'n', xaxt='n',ann=FALSE, bty = 'n')
  n <-0
  for (chr in truelengthNA) {
    n <- n + 100
    cLL <- n - 29
    cUL <- n - 2
    cUR <- n - 2
    cLR <- n - 29
    
    polygon(c(0,0,chr,chr), c(cLL,cUL,cUR,cLR), col="grey", border = NA)
    
  }
  
  
  n <- 0
  
  cLL <- c()
  cUL <- c()
  cUR <- c()
  cLR <- c()
  
  subLOH <- dfLOHmarkerUpd[which(dfLOHmarkerUpd$indSample == indSample),]
  for (indChr in myChr) {
    
    SubdfLOHmarkerUpd <- subLOH[which(subLOH$indChr == indChr),]
    subLOHref <- SubdfLOHmarkerUpd[which(as.numeric(as.character(SubdfLOHmarkerUpd$AFshift)) > 0),]
    subLOHalt <- SubdfLOHmarkerUpd[which(as.numeric(as.character(SubdfLOHmarkerUpd$AFshift)) < 0),]
    
    n <- n + 100
    cLL <- n - 29
    cUL <- n - 2
    cUR <- n - 2
    cLR <- n - 29
    #Check if there are LOH/AF shift for that sample
    
    if (nrow(subLOHalt) >0) {
      for ( indexLOH in (1:nrow(subLOHalt))){
        polygon(c(subLOHalt$leftMark[indexLOH], subLOHalt$leftMark[indexLOH], 
                  subLOHalt$rightMark[indexLOH], subLOHalt$rightMark[indexLOH] ), 
                c(cLL,cUL,cUR,cLR), col="tomato2", border=NA)
      } }
    if (nrow(subLOHref) >0) {
      for ( indexLOH in (1:nrow(subLOHref))){
        polygon(c(subLOHref$leftMark[indexLOH], subLOHref$leftMark[indexLOH], 
                  subLOHref$rightMark[indexLOH], subLOHref$rightMark[indexLOH] ), 
                c(cLL,cUL,cUR,cLR), col="steelblue", border=NA)
      } }
    
  }
  
  
  # from 16  to 1
  centNA <- c(545074,327894,611292,250657,151953,440359,427100,351037,
              94814,490686,101306,152858,451488,111477,232532,148940 )
  n <- 0
  
  cLL <- c()
  cUL <- c()
  cUR <- c()
  cLR <- c()
  
  for (indCent in centNA) {
    
    n <- n + 100
    cLL <- n - 29
    cUL <- n - 2
    cUR <- n - 2
    cLR <- n - 29 
    circle(indCent,n-14, r=c(10000,7),border="white",col="white",lty=1,lwd=1)
  }
  dev.off()
  
}
