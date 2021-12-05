#Build set of variants haplotype-----
pathVcfparents <- "/home/smozzachiodi/Batch2RTG-OS1364/"



refSetVar1364 <- read.table(file= paste(pathVcfparents, paste("OS1364S288c.Filtvarp4.vcf"), sep="")) 
  
myChr <- c("chrIII")  

#blockVar1364 <- c()
library(stringr)
  #Select AF---
  info1364 <- data.frame(refSetVar1364$V1, refSetVar1364$V2,  unlist(lapply(strsplit(as.character(refSetVar1364$V8), ";"), '[[', 4)))
  colnames(info1364) <- c("chr", "pos", "AF")

#Load RTG vcf filtered and isec----------------

pathVcf <- "/home/smozzachiodi/Batch2RTG-OS1364/mappResuStnd"

mySample <- c("D287R097","D287R098","D287R099",
              "D287R100","D287R101","D415R90","D287R112","D415R92","D415R88","D287R114","D415R89","D287R115","D287R116","D287R117",
              "D287R118","D287R119","D287R120","D287R121","D287R122","D287R123","D287R124","D287R125","D287R126","D415R93","D287R127",
              "D415R85","D287R129","D287R130","D287R131","D287R146","D287R147","D287R148","D287R149",
              
              "D287R104","D415R91","D287R132","D287R133","D287R134","D287R135","D287R136","D287R137","D287R138",
              "D287R139","D287R157","D287R158","D415R86","D287R170","D287R163","D287R164","D287R165","D287R166",
              
              "D287R102","D287R140","D287R141","D287R142","D415R96","D287R143","D287R144","D415R94","D287R145",
              "D287R150","D287R151","D287R152","D287R153","D287R154","D287R155","D287R156",
              
              "D167R74","D167R75","D167R76","D167R77","D167R78","D167R79","D167R80","D536R91"
)

LOHsummary <- c()
LOHListStatus <- c()
dfLOHmark <- c()

for(indSample in mySample) {

n <- 1 #flag for subtel filtering
myChr <- c("chrI","chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII",
           "chrIX", "chrX", "chrXI", "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI")

chrScer <- c(230218,813184,316620,1531933,576874,270161,1090940,562643
             ,439888,745751,666816,1078177,924431,784333,1091291,948066)

vcfSample <- read.table(file= paste(pathVcf,indSample,paste(indSample,"S288c", sep=""), paste(indSample, ".p4.isec.var.vcf", sep=""), sep="/")) 
infoVcfSample <- data.frame(vcfSample$V1, vcfSample$V2,  unlist(lapply(strsplit(as.character(vcfSample$V8), ";"), '[[', 4)))
colnames(infoVcfSample) <- c("chr", "pos", "AF")

require(dplyr) 
library(data.table)
library(rlist)

for(indChr in myChr)  {
  
  info1364chr <- info1364[which (info1364$chr == indChr),]
  infoVcfSampleChr <- infoVcfSample[which (infoVcfSample$chr == indChr),]

  #Filterig for subtelomere regions---------
  
  subtelL <- c(31606, 9957, 11686, 22362, 26495, 33994, 10876, 19197, 30617, 31885,
               3107, 20776, 9849, 11117, 34288, 26457)
  subtelR <- c(197027, 802713, 314566, 1528575, 566584, 266753, 1069643, 540634, 425000,
               733194, 642915, 1061277, 921042, 765393, 1070680, 924020)
  
  
  info1364chr <- info1364chr[which(info1364chr$pos > subtelL[n] & info1364chr$pos < subtelR[n]),]
  infoVcfSampleChr <- infoVcfSampleChr[which(infoVcfSampleChr$pos > subtelL[n] & infoVcfSampleChr$pos < subtelR[n]),]
  
  
#remove variants with AF1 shared by sample and ref, they are not changing
AFref1 <- info1364chr[which(info1364chr$AF == "AF=1"),]
#AFsample1 <- infoVcfSampleChr[which(infoVcfSampleChr$AF == "AF=1"),] not nedded
commonAF1 <- data.frame(AFref1$pos)#intersect 
colnames(commonAF1) <- c("common")

info1364chr <- info1364chr[!info1364chr$pos %in% commonAF1$common,]
infoVcfSampleChr <- infoVcfSampleChr[!infoVcfSampleChr$pos %in% commonAF1$common ,]
#select the variants present in the sample but not ref for AF or pos
prova <- anti_join(infoVcfSampleChr , info1364chr)
posizioni <- prova$pos

#select variants which are present in both but with different AF because I select only position that are 
#shared between samples and control
differentAF <- info1364chr[info1364chr$pos %in% posizioni ,]

#Select variants present in ref but not sample, could be shifted to AF 0
shiftedto01 <- anti_join(info1364chr, infoVcfSampleChr)
posizioniShifted <- shiftedto01[which(shiftedto01$AF != "AF=1" ),]$pos
#select variants that are not present in the sample
differentAFshifted <- info1364chr[info1364chr$pos %in% posizioniShifted ,]


#Added 15052020 shift in AF------
#afShiftAll <- c()

#for (indShift in 1:nrow(differentAF)) {
  
#  afRef <- info1364chr[which(info1364chr$pos == differentAF$pos[indShift]),]$AF
#  afRef <- as.numeric(strsplit(as.character(afRef), split="=")[[1]][2])
#  afAlt <- infoVcfSampleChr[which(infoVcfSampleChr$pos == differentAF$pos[indShift]),]$AF
#  afAlt <- as.numeric(strsplit(as.character(afAlt), split="=")[[1]][2])
  
#  afShift <- afRef -afAlt
#  afShiftAll <- as.data.frame(rbind(afShiftAll, cbind(afShift, differentAF$pos[indShift], indChr, indSample) ) )
#}

#indexDiscard <- row.names(afShiftAll[which(afShiftAll$afShift == 0.666667 | afShiftAll$afShift == -0.666667),])

#if(length(indexDiscard) > 0)
  
#{differentAF <- differentAF[-c(as.numeric(indexDiscard)),] 
#    }

#afShiftAll01 <- c()
#for (indShift in 1:nrow(differentAFshifted)) {
  
#  afRef <- info1364chr[which(info1364chr$pos == differentAFshifted$pos[indShift]),]$AF
#  afRef <- as.numeric(strsplit(as.character(afRef), split="=")[[1]][2])
#  afAlt <- infoVcfSampleChr[which(infoVcfSampleChr$pos == differentAFshifted$pos[indShift]),]$AF
  #added if to check AF ==0 in the sample
#  if(length(afAlt)==0){afAlt <- 0}
  
#  else(afAlt <- as.numeric(strsplit(as.character(afAlt), split="=")[[1]][2]) )

#    afShift <- afRef -afAlt
#  afShiftAll01 <- as.data.frame(rbind(afShiftAll01, cbind(afShift, differentAFshifted$pos[indShift], indChr, indSample) ) )
#}
#remove from differentAFshifted the markers with AFshift = 0.6667 for 3n

#indexDiscard <- row.names(afShiftAll01[which(afShiftAll01$afShift == 0.666667 | afShiftAll01$afShift == -0.666667),])

#if(length(indexDiscard) > 0)
  
#{
#differentAFshifted <- differentAFshifted[-c(as.numeric(indexDiscard)),]
#}

#-----------------------

mergeAFvar <- rbind(differentAF, differentAFshifted)
mergeAFvar <- mergeAFvar[order(mergeAFvar$pos),]
mergeAFvar <- mergeAFvar[which(mergeAFvar$chr != "NA"),]
indexAFvar <- unique(match(mergeAFvar$pos, info1364chr$pos))

Breaks <- c(0, which(diff((indexAFvar)) != 1), length((indexAFvar)) )
#calculate index in each subgroups----
groupMark <- sapply(seq(length(Breaks) - 1),
                    function(i) indexAFvar[(Breaks[i] + 1):Breaks[i+1]]) 
indexLOHseg <- which(lengths(groupMark) >= 10)#here modify from 5 to 10 all 5 where likely false positive so went to 10 to be more stringent
groupMarkLOH <- as.numeric(unlist(groupMark[indexLOHseg]))


#trying to group for AF shifting towards ref or alternative-----
for (indIndexLOHseg in indexLOHseg) {
  #take index in ref file of each marker---
  markersPos <- groupMark[indIndexLOHseg]
  
  #retrivering AF Ref and RTG sample-----
  for (indMarkFlag in 1:length(unlist(markersPos))) {
    
    indMarkerLOH <- unlist(markersPos)[indMarkFlag]
    positionRef <- info1364chr[indMarkerLOH ,][2]
    retriveredAFsample <- infoVcfSampleChr[which(infoVcfSampleChr$pos == as.numeric(positionRef)),][3]
    if (nrow(retriveredAFsample) == 0) 
    {
      AFsample <- 0
    }
    if (nrow(retriveredAFsample) == 1) 
    {
    AFsample <- as.numeric(strsplit(as.character(retriveredAFsample[1, 1]), split = "=")[[1]][2])
    }
    
    AFref <- as.numeric(strsplit(as.character(info1364chr[indMarkerLOH ,][1, 3]), split = "=")[[1]][2])
    AFshift <- AFref - AFsample 
    dfLOHmark <- rbind(dfLOHmark, cbind(indSample, indChr, positionRef, AFref, AFsample, AFshift) )
      
  }
  #LOHListStatus <- list.append(LOHListStatus , as.list(as.data.frame((dfLOHmark)))) 
  
}



#isolate index of fst and last marker to find back the position
fstmarker <- unlist(lapply(groupMark[indexLOHseg], `[[`, 1))
lastmarker <- unlist(lapply(groupMark[indexLOHseg], tail, n = 1L))
n <- n+1 #flag for next subtel

if (length(fstmarker) >= 1) {
for (indMarker in 1:length(fstmarker)){
  #I take also the info from the vcf of the RTG sample------
  lengthLOH <- (as.numeric(info1364chr[lastmarker[indMarker] ,][2])-as.numeric(info1364chr[fstmarker[indMarker] ,][2]))
  regionLOH <- unlist(c(info1364chr[fstmarker[indMarker] ,][2] , info1364chr[lastmarker[indMarker] ,][2]))
  nmarkerLOH <- as.numeric(length(subset(info1364chr$pos, subset = info1364chr$pos >= regionLOH[1] & info1364chr$pos <= regionLOH[2]) ))
  LOHsummary <- rbind(LOHsummary, cbind(indSample,indChr,info1364chr[fstmarker[indMarker] ,][2-3], info1364chr[lastmarker[indMarker] ,][2-3], lengthLOH, nmarkerLOH))
}
}
   else{next}
}

Breaks <- c()
fstmarker <- c()
lastmarker <- c()
posizioni <- c()
prova <- c()
shiftedto01 <- c()
posizioniShifted <- c()
groupMark <- c()
} #not sure if the last parentesis is correct

colnames(LOHsummary) <- c("Sample", "Chromosome", "Start", "StartAF", "End", "EndAF", "LengthLOH", "nMarkers")
sumdata <- c()
LOHsummary <- LOHsummary[order(LOHsummary$Start),]
LOHsummary <- LOHsummary[order(LOHsummary$Chromosome),]
#LOHsummary <- LOHsummary[which(LOHsummary$Chromosome != "chrIII"),]
LOHsummary <- LOHsummary[which(LOHsummary$Length != "NA"),]
LOHsummary10k <- LOHsummary[which(LOHsummary$Length > 10000),]
LOHsummary1k <- LOHsummary[which(LOHsummary$Length > 1000),]

write.table(x=LOHsummary, file="/home/smozzachiodi/Batch2RTG-OS1364/mappResuStnd/LOH_summaryp4_5marks_11022021.csv")


#Filtering for shared region------------

#---------------------
#down to samples markers area designation----


dfLOHmarkerUpd <- c()

for (indSample in mySample) {
  
  subsetLOHsample <- LOHsummary[which(LOHsummary$Sample == indSample),] 
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

#Until filtering shared regions not done-----


#Rescaling length --------------------

chrScerCumSum <- c(0, 230218, 1043402, 1360022, 2891955, 3468829, 3738990, 4829930, 5392573, 
                   5832461, 6578212, 7245028, 8323205,9247636, 10031969, 11123260)

sumdata <- c()

  for (n in 1:16) {
    
    subChr <- LOHsummary[which(LOHsummary$Chromosome == myChr[n]),]
    sumleng <- (chrScerCumSum[n])
    totposStart <- (subChr$Start + sumleng)
    totposEnd <- (subChr$End + sumleng)
    if(length(subChr)>0) {
      sumdata <- rbind(sumdata, cbind(totposStart, totposEnd))
    }
  }


LOHsummaryPlot <- cbind(LOHsummary, sumdata)


#Rescaling length marker wise--------------------

chrScerCumSum <- c(0, 230218, 1043402, 1360022, 2891955, 3468829, 3738990, 4829930, 5392573, 
                   5832461, 6578212, 7245028, 8323205,9247636, 10031969, 11123260)

dfLOHmarkerUpd <- data.frame(dfLOHmarkerUpd)
dfLOHmarkerUpd <- dfLOHmarkerUpd[order(dfLOHmarkerUpd$pos),]
dfLOHmarkerUpd$indChr <- factor(dfLOHmarkerUpd$indChr, levels=c("chrI","chrII","chrIII","chrIV","chrV","chrVI","chrVII","chrVIII","chrIX","chrX",
                                                                "chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI"))
dfLOHmarkerUpd <- dfLOHmarkerUpd[order(dfLOHmarkerUpd$indChr),]
#dfLOHmarkerUpd <- dfLOHmarkerUpd[which(dfLOHmarkerUpd$Chromosome != "chrIII"),]
dfLOHmarkerUpd <- na.omit(dfLOHmarkerUpd)
#added filter AFshift too big---------
dfLOHmarkerUpd <- dfLOHmarkerUpd[which(dfLOHmarkerUpd$AFshift != 0.75 ),]
dfLOHmarkerUpd <- dfLOHmarkerUpd[which(dfLOHmarkerUpd$AFshift != -0.75 ),]





sumdataMarkerLOH <- c()
for (n in 1:16) {
  
  subChr <- c()
  sumleng <- c()
  totposStart <- c()
  totposEnd <- c()
  
  subChr <- dfLOHmarkerUpd[which(dfLOHmarkerUpd$indChr == myChr[n]),]
  sumleng <- (chrScerCumSum[n])
  totposStart <- (as.numeric(subChr$leftMark) + as.numeric(sumleng))
  totposEnd <- (as.numeric(subChr$rightMark) + as.numeric(sumleng))
  
  if(length(subChr)>0) {
    sumdataMarkerLOH <- rbind(sumdataMarkerLOH, cbind(totposStart, totposEnd))
  }
}


dfLOHmarkerUpdPlot <- cbind(dfLOHmarkerUpd, sumdataMarkerLOH)
write.table(dfLOHmarkerUpdPlot, file="/home/smozzachiodi/Batch2RTG-OS1364/mappResuStnd/MarkerStatusOS1364_RTG_5marksp4_11022021.csv",sep = ",")

myChr <- c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII",
           "chrIX", "chrX", "chrXI", "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI")

pdf("/home/smozzachiodi/Batch2RTG-OS1364/mappResuStnd/OS1364.LOH.5markers.p4.pdf", height = 6, width = 12)
plot(c(0,0),c(0,0),xlim=c(0,12071326),ylim=c(0,1500), yaxt = 'n', xaxt='n',ann=FALSE, bty = 'n')
n <-0
for (sample in mySample) {
  n <- n + 40
  cLL <- n - 19
  cUL <- n - 2
  cUR <- n - 2
  cLR <- n - 19
  
  polygon(c(0,0,12071326,12071326), c(cLL,cUL,cUR,cLR), col="lightgrey")
}  
  
n <- 0

cLL <- c()
cUL <- c()
cUR <- c()
cLR <- c()

for (indSample in mySample) {
  dfLOHmarkerUpdPlotSub <- c()
  n <- n + 40
  cLL <- n - 19
  cUL <- n - 2
  cUR <- n - 2
  cLR <- n - 19
        dfLOHmarkerUpdPlotSub <- dfLOHmarkerUpdPlot[which(dfLOHmarkerUpdPlot$indSample == indSample),]
  #Qui controllo che ci siano LOH per quel chr dato un campione, se non c'è nulla non plotto nulla

    for ( indexLOH in (1:nrow(dfLOHmarkerUpdPlotSub))){
      
      if (as.numeric(dfLOHmarkerUpdPlotSub$AFshift[indexLOH])>0)
        
      {polygon(c(dfLOHmarkerUpdPlotSub$totposStart[indexLOH], dfLOHmarkerUpdPlotSub$totposStart[indexLOH], dfLOHmarkerUpdPlotSub$totposEnd[indexLOH], dfLOHmarkerUpdPlotSub$totposEnd[indexLOH] ), c(cLL,cUL,cUR,cLR), col="steelblue", border=NA)}
    
      if (as.numeric(dfLOHmarkerUpdPlotSub$AFshift[indexLOH])<0)
      
      {polygon(c(dfLOHmarkerUpdPlotSub$totposStart[indexLOH], dfLOHmarkerUpdPlotSub$totposStart[indexLOH], dfLOHmarkerUpdPlotSub$totposEnd[indexLOH], dfLOHmarkerUpdPlotSub$totposEnd[indexLOH] ), c(cLL,cUL,cUR,cLR), col="firebrick1", border=NA)}

      } 
       
}  
i <- 1

n <- 0

for (sample in mySample) {
  
  n <- n + 40
  y0 <- n - 39
  y1 <- n
  
  
  i <- 1
  for (allChr in myChr) {
    
    segments(chrScerCumSum[i],y0, x1 =chrScerCumSum[i], y1, lwd=2)
    i <- i +1
    
  }}
dev.off()

