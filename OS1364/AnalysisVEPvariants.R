#Analysis VEP mutation triploid OS1364------

pathVEPoutput <- "/home/simone/Scrivania/Scrivania/NewBeerStrains/GenomWideVariantsLOF/OS1364_VEP.txt"
pathVEPoutputChrIII <- "/home/simone/Scrivania/Scrivania/NewBeerStrains/GenomWideVariantsLOF/OS1364_p4_VEP.txt"

pathVcfParent <- "/home/simone/Scrivania/Scrivania/NewBeerStrains/GenomWideVariantsLOF/OS1364SGD_2010.var.p3.vcf"
pathVcfParentp4 <- "/home/simone/Scrivania/Scrivania/NewBeerStrains/GenomWideVariantsLOF/OS1364SGD_2010.var.p4.vcf"

vepOutput <- read.table(pathVEPoutput)
vepOutputChrIII <- read.table(pathVEPoutputChrIII)
vcfParent <- read.table(pathVcfParent)
vcfParentp4 <- read.table(pathVcfParentp4)

library(tidyverse)
#filter for chrIII p3 variants called:
filtVepOutput <- vepOutput[which( sapply(strsplit(as.character(vepOutput$V2), split = c(":") ), "[[", 1) != "III" ),]

#remove frameshift variants derived from large indel
filtVepOutput <- filtVepOutput[which(#nchar(as.character(filtVepOutput$V3)) == 1 & 
                                       filtVepOutput$V4 != "upstream_gene_variant" &
                                       filtVepOutput$V4 != "downstream_gene_variant" & 
                                       filtVepOutput$V4 != "frameshift_variant" #& 
                                       #filtVepOutput$V4 != "frameshift_variant,stop_lost" &
                                       #filtVepOutput$V4 != "frameshift_variant,stop_gain"),]
),]
                                         
filtVepOutputChrIII <- vepOutputChrIII[which(#nchar(as.character(vepOutputChrIII$V3)) == 1 & 
                                               vepOutputChrIII$V4 != "upstream_gene_variant" &
                                               vepOutputChrIII$V4 != "downstream_gene_variant" & 
                                               vepOutputChrIII$V4 != "frameshift_variant" #& 
                                              # vepOutputChrIII$V4 != "frameshift_variant,stop_lost" &
                                              # vepOutputChrIII$V4 != "frameshift_variant,stop_gain"),]
),]
#Do not merge since saved with different columns elements
highImpactVariants <- filtVepOutput[which(filtVepOutput$V5 == "HIGH"),]
moderateImpactVariants <- filtVepOutput[which(filtVepOutput$V5 == "MODERATE"),]
lowImpactVariants <- filtVepOutput[which(filtVepOutput$V5 == "LOW"),]

highImpactVariantsChrIII <- filtVepOutputChrIII[which(filtVepOutputChrIII$V5 == "HIGH"),]
moderateImpactVariantsChrIII <- filtVepOutputChrIII[which(filtVepOutputChrIII$V5 == "MODERATE"),]
lowImpactVariantsChrIII <- filtVepOutputChrIII[which(filtVepOutputChrIII$V5 == "LOW"),]

#Analysis of deleterious variants----
#Remove variants in subtelomeres----

subtelL <- c(31606, 9957, 11686, 22362, 26495, 33994, 10876, 19197, 30617, 31885,
             3107, 20776, 9849, 11117, 34288, 26457)
subtelR <- c(197027, 802713, 314566, 1528575, 566584, 266753, 1069643, 540634, 425000,
             733194, 642915, 1061277, 921042, 765393, 1070680, 924020)

Chrom <- c("chrI","chrII","chrIII","chrIV","chrV","chrVI","chrVII","chrVIII",
           "chrIX","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI")
#This has been skipped to make the data comparable with those of the 1011
highImpactVariantsNoSubtel <- c()
nChrom <-1

for (indChr in Chrom) {
  
  subHighImpact <- highImpactVariants[which(paste("chr", (sapply(strsplit(as.character(highImpactVariants$V2), split = ":"),"[[",1)), sep="") == indChr),] 

  
  highImpactVariantsNoSubtel <- rbind(highImpactVariantsNoSubtel,
                                      subHighImpact[which(as.numeric(sapply(strsplit(as.character(subHighImpact$V2), split = "-"),"[[",2))  > subtelL[nChrom] &
                                                         as.numeric(sapply(strsplit(as.character(subHighImpact$V2), split = "-"),"[[",2)) < subtelR[nChrom]),]  )
  nChrom <- nChrom+1
}  
#chromosome III
Chrom <- c("chrI","chrII","chrIII","chrIV","chrV","chrVI","chrVII","chrVIII",
           "chrIX","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI")

highImpactVariantsNoSubtelIII <- c()
nChrom <-3

for (indChr in Chrom[3]) {
  
  subHighImpactIII <- highImpactVariantsChrIII[which(paste("chr", (sapply(strsplit(as.character(highImpactVariantsChrIII$V2), split = ":"),"[[",1)), sep="") == indChr),] 
  
  
  highImpactVariantsNoSubtelIII <- rbind(highImpactVariantsNoSubtelIII,
                                      subHighImpactIII[which(as.numeric(sapply(strsplit(as.character(subHighImpactIII$V2), split = "-"),"[[",2))  > subtelL[nChrom] &
                                                            as.numeric(sapply(strsplit(as.character(subHighImpactIII$V2), split = "-"),"[[",2)) < subtelR[nChrom]),]  )
  
}  


#Search for GO terms enrichment----
genesHighImpactNoSubtel <- unique(highImpactVariantsNoSubtel$V7)
classesHighImpactNoSubtel <- table(highImpactVariantsNoSubtel$V4)
genesHighImpact_IIINoSubtel <- unique(highImpactVariantsNoSubtelIII$V7)
classesHighImpact_IIINoSubtel <- table(highImpactVariantsNoSubtelIII$V4)
allgenesHighImpactNoSubtel <- as.character(genesHighImpactNoSubtel,genesHighImpact_IIINoSubtel)

#Here did the same including frameshift and subtelomeres
genesHighImpact <- unique(highImpactVariants$V7)
classesHighImpact <- table(highImpactVariants$V4)
genesHighImpact_III <- unique(highImpactVariantsChrIII[which(paste("chr", (sapply(strsplit(as.character(highImpactVariantsChrIII$V2), split = ":"),"[[",1)), sep="") == "chrIII"), 7])
classesHighImpact_III <- table(highImpactVariantsChrIII[which(paste("chr", (sapply(strsplit(as.character(highImpactVariantsChrIII$V2), split = ":"),"[[",1)), sep="") == "chrIII"), 4])
allgenesHighImpact <- c(as.character(genesHighImpact),as.character(genesHighImpact_III))


write.csv(allgenesHighImpact, file ="/home/simone/Scrivania/Scrivania/NewBeerStrains/GenomWideVariantsLOF/High_impact_genes_OS1364_subtel.txt", row.names = F)

write.csv(highImpactVariants, file="/home/simone/Scrivania/Scrivania/NewBeerStrains/GenomWideVariantsLOF/OS1364Final_files/OS1364_HighImp.csv")
write.csv(highImpactVariantsChrIII, file="/home/simone/Scrivania/Scrivania/NewBeerStrains/GenomWideVariantsLOF/OS1364Final_files/OS1364_HighImpIII.csv")

#Compare with essential genes list

essentialPath <- "/home/simone/Scrivania/Scrivania/NewBeerStrains/GenomWideVariantsLOF/EssentialGenesList.csv"

essentialGenes <- read.csv(essentialPath, sep = "\t", header = T)

essentialGenesMut1364 <- highImpactVariants[highImpactVariants$V7 %in% essentialGenes$Systematic.name ,]

#Find heterozygosity level for essential/non essential highly impactful mutations
#the datasets are: highImpactVariants (including subtel and frameshifts) and highImpactVariantsNoSubtel (where frameshift and subtelomeres are not included)
posHighImpact1364 <- cbind(sapply(strsplit(
                           sapply(strsplit(as.character(highImpactVariants$V2), split = "-"),"[[",1),
                           split =":"), "[[",2),
                           paste("chr", (sapply(strsplit(as.character(highImpactVariants$V2), split = ":"),"[[",1)), sep=""),
                           as.character(highImpactVariants$V7))

statusVariant <- c()
for (indHigh in 1:nrow(posHighImpact1364)) {
  
  variantSubChr <- vcfParent[which(vcfParent$V1 == posHighImpact1364[indHigh, 2]),] 
  
  variantSub1364 <- variantSubChr[which(variantSubChr$V2 == posHighImpact1364[indHigh, 1]),]
  
  AFvariant <- sapply(strsplit(as.character(variantSub1364$V10), split = ":"), "[[", 1)
  
  statusVariant <- rbind(statusVariant, cbind(AFvariant, posHighImpact1364[indHigh, 3]))
}
#Same for chromosome III
highImpactVariantsChrIIIFilt <- highImpactVariantsChrIII[which(paste("chr", (sapply(strsplit(as.character(highImpactVariantsChrIII$V2), split = ":"),"[[",1)), sep="") == "chrIII"), ]
posHighImpact1364III <- cbind(sapply(strsplit(
                           sapply(strsplit(as.character(highImpactVariantsChrIIIFilt$V2), split = "-"),"[[",1),
                           split=":"), "[[",2),
                           paste("chr", (sapply(strsplit(as.character(highImpactVariantsChrIIIFilt$V2), split = ":"),"[[",1)), sep=""),
                           as.character(highImpactVariantsChrIIIFilt$V7))

statusVariantIII <- c()

for (indHigh in 1:nrow(posHighImpact1364III)) {
  
  variantSubChrIII <- vcfParentp4[which(vcfParentp4$V1 == posHighImpact1364III[indHigh, 2]),] 
  
  variantSub1364III <- variantSubChrIII[which(variantSubChrIII$V2 == posHighImpact1364III[indHigh, 1]),]
  
  AFvariantIII <- sapply(strsplit(as.character(variantSub1364III$V10), split = ":"), "[[", 1)
  
  statusVariantIII <- rbind(statusVariantIII, cbind(AFvariantIII, posHighImpact1364III[indHigh, 3]))
}

#Check heterozygosity level of high impact function

statusVariantAll <- data.frame(rbind(statusVariant, statusVariantIII))

essential <- c()
status <- c()

for (indVariant in 1:nrow(statusVariantAll)) {
  
  checkEssential <- statusVariantAll[statusVariantAll$V2[indVariant] %in% essentialGenesMut1364$V7][indVariant ,]
  
  if (ncol(checkEssential) > 1) {essential <- rbind(essential, "Essential")}  
  if (ncol(checkEssential) < 1) {essential <- rbind(essential, "NonEssential")} 
  colnames(checkEssential) <- c()
  
  if (statusVariantAll$AFvariant[indVariant] == "1/1/1" || statusVariantAll$AFvariant[indVariant] == "1/1/1/1")
  { status <- rbind(status, "Homozygous")}
  
  if (statusVariantAll$AFvariant[indVariant] != "1/1/1" && statusVariantAll$AFvariant[indVariant] != "1/1/1/1")
  { status <- rbind(status, "Heterozygous")}  
  
}

statusVariantAllInfo1364 <- cbind(statusVariantAll, essential, status)

#Check synonymous vs non synonymous-----
#filt for subtelomeres
subtelL <- c(31606, 9957, 11686, 22362, 26495, 33994, 10876, 19197, 30617, 31885,
             3107, 20776, 9849, 11117, 34288, 26457)
subtelR <- c(197027, 802713, 314566, 1528575, 566584, 266753, 1069643, 540634, 425000,
             733194, 642915, 1061277, 921042, 765393, 1070680, 924020)

Chrom <- c("chrI","chrII","chrIII","chrIV","chrV","chrVI","chrVII","chrVIII",
           "chrIX","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI")

filtVepOutputNoSubtel <- c()
nChrom <-1

for (indChr in Chrom) {
  
  subHighImpact <- filtVepOutput[which(paste("chr", (sapply(strsplit(as.character(filtVepOutput$V2), split = ":"),"[[",1)), sep="") == indChr),] 
  
  
  filtVepOutputNoSubtel <- rbind(filtVepOutputNoSubtel,
                                      subHighImpact[which(as.numeric(sapply(strsplit(as.character(subHighImpact$V2), split = "-"),"[[",2))  > subtelL[nChrom] &
                                                            as.numeric(sapply(strsplit(as.character(subHighImpact$V2), split = "-"),"[[",2)) < subtelR[nChrom]),]  )
  nChrom <- nChrom+1
}  
#chromosome III
Chrom <- c("chrI","chrII","chrIII","chrIV","chrV","chrVI","chrVII","chrVIII",
           "chrIX","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI")

filtVepOutputNoSubtelIII <- c()
nChrom <-3

for (indChr in Chrom[3]) {
  
  subHighImpactIII <- filtVepOutputChrIII[which(paste("chr", (sapply(strsplit(as.character(filtVepOutputChrIII$V2), split = ":"),"[[",1)), sep="") == indChr),] 
  
  
  filtVepOutputNoSubtelIII <- rbind(filtVepOutputNoSubtelIII,
                                         subHighImpactIII[which(as.numeric(sapply(strsplit(as.character(subHighImpactIII$V2), split = "-"),"[[",2))  > subtelL[nChrom] &
                                                                  as.numeric(sapply(strsplit(as.character(subHighImpactIII$V2), split = "-"),"[[",2)) < subtelR[nChrom]),]  )
  
}  

select_missense_variants <- filtVepOutputNoSubtel[which(filtVepOutputNoSubtel$V4 =="missense_variant"),]
select_missense_variantsIII <- filtVepOutputNoSubtelIII[which(filtVepOutputNoSubtelIII$V4 =="missense_variant"),]

#Find heterozygosity of missense variants
select_missense_variants <- select_missense_variants[-9210 ,]
posHighImpact1364 <- cbind(sapply(strsplit(as.character(select_missense_variants$V2), split = "-"),"[[",2),
                           paste("chr", (sapply(strsplit(as.character(select_missense_variants$V2), split = ":"),"[[",1)), sep=""),
                           as.character(select_missense_variants$V7))

statusVariantMissense <- c()
for (indHigh in 1:nrow(posHighImpact1364)) {
  
  variantSubChr <- vcfParent[which(vcfParent$V1 == posHighImpact1364[indHigh, 2]),] 
  
  variantSub1364 <- variantSubChr[which(variantSubChr$V2 == posHighImpact1364[indHigh, 1]),]
  
  AFvariant <- sapply(strsplit(as.character(variantSub1364$V10), split = ":"), "[[", 1)
  
  statusVariantMissense <- rbind(statusVariantMissense, cbind(AFvariant, posHighImpact1364[indHigh, 3]))
}
#Same for chromosome III
posHighImpact1364III <- cbind(sapply(strsplit(as.character(select_missense_variantsIII$V2), split = "-"),"[[",2),
                              paste("chr", (sapply(strsplit(as.character(select_missense_variantsIII$V2), split = ":"),"[[",1)), sep=""),
                              as.character(select_missense_variantsIII$V7))

statusVariantIII <- c()

for (indHigh in 1:nrow(posHighImpact1364III)) {
  
  variantSubChrIII <- vcfParentp4[which(vcfParentp4$V1 == posHighImpact1364III[indHigh, 2]),] 
  
  variantSub1364III <- variantSubChrIII[which(variantSubChrIII$V2 == posHighImpact1364III[indHigh, 1]),]
  
  AFvariantIII <- sapply(strsplit(as.character(variantSub1364III$V10), split = ":"), "[[", 1)
  
  statusVariantIII <- rbind(statusVariantIII, cbind(AFvariantIII, posHighImpact1364III[indHigh, 3]))
}

library(ggplot2)

dfstatusVariants <- data.frame(rbind(statusVariantIII, statusVariantMissense))

write.csv(dfstatusVariants$V2, file = "/home/simone/Scrivania/Scrivania/NewBeerStrains/OS1364/Variants_VEP_Prediction/Gene_missense_OS1364.csv")

pdf("/home/simone/Scrivania/Scrivania/NewBeerStrains/FinalFilesPaper/Sfig1/OS1364_missense")
ggplot(dfstatusVariants, aes(x=AFvariantIII))+geom_bar(stat = "count")+theme_bw()
dev.off


#Analysis VEP mutation tetraploid OS1431------

pathVEPoutput1431 <- "/home/simone/Scrivania/Scrivania/NewBeerStrains/GenomWideVariantsLOF/OS1431vep/OS1431_VEP.txt"

pathVcfParent1431 <- "/home/simone/Scrivania/Scrivania/NewBeerStrains/GenomWideVariantsLOF/OS1431vep/OS1431SGD_2010.var.p4.vcf"

vepOutput1431 <- read.table(pathVEPoutput1431)

vcfParent1431 <- read.table(pathVcfParent1431)

library(tidyverse)

#remove frameshift variants derived from large indel
filtVepOutput1431 <- vepOutput1431[which(#nchar(as.character(vepOutput1431$V3)) == 1 & 
                                       vepOutput1431$V4 != "frameshift_variant" & 
                                      vepOutput1431$V4 != "upstream_gene_variant" &
                                       vepOutput1431$V4 != "downstream_gene_variant"# & 
                                       #vepOutput1431$V4 != "frameshift_variant" & 
                                      # vepOutput1431$V4 != "frameshift_variant,stop_lost" &
                                       #vepOutput1431$V4 != "frameshift_variant,stop_gain"),]
),]
#Do not merge since saved with different columns elements
highImpactVariants1431 <- filtVepOutput1431[which(filtVepOutput1431$V5 == "HIGH"),]
moderateImpactVariants1431 <- filtVepOutput1431[which(filtVepOutput1431$V5 == "MODERATE"),]
lowImpactVariants1431 <- filtVepOutput1431[which(filtVepOutput1431$V5 == "LOW"),]

#Analysis of deleterious variants----
#Remove variants in subtelomeres----

subtelL <- c(31606, 9957, 11686, 22362, 26495, 33994, 10876, 19197, 30617, 31885,
             3107, 20776, 9849, 11117, 34288, 26457)
subtelR <- c(197027, 802713, 314566, 1528575, 566584, 266753, 1069643, 540634, 425000,
             733194, 642915, 1061277, 921042, 765393, 1070680, 924020)

Chrom <- c("chrI","chrII","chrIII","chrIV","chrV","chrVI","chrVII","chrVIII",
           "chrIX","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI")

highImpactVariantsNoSubtel1431 <- c()
nChrom <-1

for (indChr in Chrom) {
  
  subHighImpact1431 <- highImpactVariants1431[which(paste("chr", (sapply(strsplit(as.character(highImpactVariants1431$V2), split = ":"),"[[",1)), sep="") == indChr),] 
  
  
  highImpactVariantsNoSubtel1431 <- rbind(highImpactVariantsNoSubtel1431,
                                      subHighImpact1431[which(as.numeric(sapply(strsplit(as.character(subHighImpact1431$V2), split = "-"),"[[",2))  > subtelL[nChrom] &
                                                            as.numeric(sapply(strsplit(as.character(subHighImpact1431$V2), split = "-"),"[[",2)) < subtelR[nChrom]),]  )
  nChrom <- nChrom+1
}  

#Search for GO terms enrichment----
genesHighImpact1431 <- unique(highImpactVariantsNoSubtel1431$V7)
classesHighImpact1431 <- table(highImpactVariantsNoSubtel1431$V4)
allgenesHighImpact1431 <- as.character(genesHighImpact1431)

#Search for GO terms enrichment subtel and frameshift----
genesHighImpact1431 <- unique(highImpactVariants1431$V7)
classesHighImpact1431 <- table(highImpactVariants1431$V4)
allgenesHighImpact1431 <- as.character(genesHighImpact1431)

write.csv(allgenesHighImpact1431, file ="/home/simone/Scrivania/Scrivania/NewBeerStrains/GenomWideVariantsLOF/High_impact_genes_OS1431_subtelFrame.txt", row.names = F)
write.csv(highImpactVariants1431, file="/home/simone/Scrivania/Scrivania/NewBeerStrains/GenomWideVariantsLOF/OS1431Final_files/OS1431_HighImpact.csv")
#Compare with essential genes list

essentialPath <- "/home/simone/Scrivania/Scrivania/NewBeerStrains/GenomWideVariantsLOF/EssentialGenesList.csv"

essentialGenes <- read.csv(essentialPath, sep = "\t", header = T)

essentialGenesMut1431 <- highImpactVariants1431[highImpactVariants1431$V7 %in% essentialGenes$Systematic.name ,]

#Find heterozygosity level for essential/non essential highly impactful mutations

posHighImpact1431 <- cbind(sapply(strsplit(
                           sapply(strsplit(as.character(highImpactVariants1431$V2), split = "-"),"[[",1),
                           split = ":"), "[[", 2),
                           paste("chr", (sapply(strsplit(as.character(highImpactVariants1431$V2), split = ":"),"[[",1)), sep=""),
                           as.character(highImpactVariants1431$V7))

statusVariant1431 <- c()
for (indHigh in 1:nrow(posHighImpact1431)) {
  
  variantSubChr1431 <- vcfParent1431[which(vcfParent1431$V1 == posHighImpact1431[indHigh, 2]),] 
  
  variantSub1431 <- variantSubChr1431[which(variantSubChr1431$V2 == posHighImpact1431[indHigh, 1]),]
  
  AFvariant1431 <- sapply(strsplit(as.character(variantSub1431$V10), split = ":"), "[[", 1)
  
  statusVariant1431 <- rbind(statusVariant1431, cbind(AFvariant1431, posHighImpact1431[indHigh, 3]))
}

#Check heterozygosity level of high impact function

statusVariantAll1431 <- data.frame(statusVariant1431)

essential1431 <- c()
status1431 <- c()

for (indVariant in 1:nrow(statusVariantAll1431)) {
  
  checkEssential1431 <- statusVariantAll1431[statusVariantAll1431$V2[indVariant] %in% essentialGenesMut1431$V7][indVariant ,]
  
  if (ncol(checkEssential1431) > 1) {essential1431 <- rbind(essential1431, "Essential")}  
  if (ncol(checkEssential1431) < 1) {essential1431 <- rbind(essential1431, "NonEssential")} 
  colnames(checkEssential1431) <- c()
  
  if (statusVariantAll1431$AFvariant[indVariant] == "1/1/1/1")
  { status1431 <- rbind(status1431, "Homozygous")}
  
  if (statusVariantAll1431$AFvariant[indVariant] != "1/1/1/1")
  { status1431 <- rbind(status1431, "Heterozygous")}  
  
}

statusVariantAllInfo1431 <- cbind(statusVariantAll1431, essential1431, status1431)

statusVariantAllInfo1364_plot <- cbind(statusVariantAllInfo1364, rep.int("OS1364", times = nrow(statusVariantAllInfo1364)))
statusVariantAllInfo1431 <- cbind(statusVariantAllInfo1431, rep.int("OS1431", times = nrow(statusVariantAllInfo1431)))

colnames(statusVariantAllInfo1364_plot) <- c("AF","Gene","Essential","status","Strain")
colnames(statusVariantAllInfo1431) <- c("AF","Gene","Essential","status","Strain")

statusVariantAllStrainsInfo <- rbind(statusVariantAllInfo1364_plot,
                                     statusVariantAllInfo1431)

PlotInfo <- c()
for (indElement in 1:nrow(statusVariantAllStrainsInfo)) {
  
  VariantInfo <- paste(statusVariantAllStrainsInfo$Essential[indElement], statusVariantAllStrainsInfo$status[indElement], sep="_")
  PlotInfo <- rbind(PlotInfo, VariantInfo)
}

#Plot figure 1C---------
statusVariantPlot <- cbind(statusVariantAllStrainsInfo, PlotInfo)
Info1364 <- table(statusVariantPlot[which(statusVariantPlot$Strain == "OS1364"),]$PlotInfo)
Info1431 <- table(statusVariantPlot[which(statusVariantPlot$Strain == "OS1431"),]$PlotInfo)

plotVariant <- data.frame(rbind(cbind(Info1364, rep.int("OS1364", times = nrow(Info1364))),
                     cbind(Info1431, rep.int("OS1431", times = nrow(Info1431)))),
                     rep.int(rownames(Info1364), times = 2))

colnames(plotVariant) <- c("Number", "Strain", "Type")  
library(ggplot2)

pdf("/home/simone/Scrivania/Scrivania/NewBeerStrains/FinalFilesPaper/Figure1/PlotVariants_v2.pdf", height = 1.84, width = 2.636)
ggplot(plotVariant, aes(fill=Type, y=as.numeric(as.character(Number)), x=Strain)) + 
  geom_bar(position="stack", stat="identity")+theme_bw()+ylab("Number of variants")+
  ylim(0,350)+scale_fill_manual(values = c("#f05d4d","#80372c",
                                           "#bab8b7", "#7a7776"))+
  theme(legend.text = element_text(size=8),  axis.title.x = element_text(size=8), axis.text.y  = element_text(size=5), axis.text.x= element_text(size=8, angle=45, vjust = 0.5),axis.title.y=element_text(size=8))
dev.off()


#Check synonymous vs non synonymous 1431-----
#filt for subtelomeres
subtelL <- c(31606, 9957, 11686, 22362, 26495, 33994, 10876, 19197, 30617, 31885,
             3107, 20776, 9849, 11117, 34288, 26457)
subtelR <- c(197027, 802713, 314566, 1528575, 566584, 266753, 1069643, 540634, 425000,
             733194, 642915, 1061277, 921042, 765393, 1070680, 924020)

Chrom <- c("chrI","chrII","chrIII","chrIV","chrV","chrVI","chrVII","chrVIII",
           "chrIX","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI")

filtVepOutputNoSubtel1431 <- c()
nChrom <-1

for (indChr in Chrom) {
  
  subHighImpact1431 <- filtVepOutput1431[which(paste("chr", (sapply(strsplit(as.character(filtVepOutput1431$V2), split = ":"),"[[",1)), sep="") == indChr),] 
  
  
  filtVepOutputNoSubtel1431 <- rbind(filtVepOutputNoSubtel1431,
                                 subHighImpact1431[which(as.numeric(sapply(strsplit(as.character(subHighImpact1431$V2), split = "-"),"[[",2))  > subtelL[nChrom] &
                                                       as.numeric(sapply(strsplit(as.character(subHighImpact1431$V2), split = "-"),"[[",2)) < subtelR[nChrom]),]  )
  nChrom <- nChrom+1
}  

select_missense_variants1431 <- filtVepOutputNoSubtel1431[which(filtVepOutputNoSubtel1431$V4 =="missense_variant"),]

#Find heterozygosity of missense variants
select_missense_variants1431 <- select_missense_variants1431[-12952 , ]
posHighImpact1431 <- cbind(sapply(strsplit(as.character(select_missense_variants1431$V2), split = "-"),"[[",2),
                           paste("chr", (sapply(strsplit(as.character(select_missense_variants1431$V2), split = ":"),"[[",1)), sep=""),
                           as.character(select_missense_variants1431$V7))

statusVariantMissense1431 <- c()
for (indHigh in 1:nrow(posHighImpact1431)) {
  
  variantSubChr1431 <- vcfParent1431[which(vcfParent1431$V1 == posHighImpact1431[indHigh, 2]),] 
  
  variantSub1431 <- variantSubChr1431[which(variantSubChr1431$V2 == posHighImpact1431[indHigh, 1]),]
  
  AFvariant1431 <- sapply(strsplit(as.character(variantSub1431$V10), split = ":"), "[[", 1)
  
  statusVariantMissense1431 <- rbind(statusVariantMissense1431, cbind(AFvariant1431, posHighImpact1431[indHigh, 3]))
}


library(ggplot2)

dfstatusVariants <- data.frame(rbind(statusVariantMissense1431))

write.csv(dfstatusVariants$V2, file = "/home/simone/Scrivania/Scrivania/NewBeerStrains/OS1431/Gene_missense_OS1364.csv")

pdf("/home/simone/Scrivania/Scrivania/NewBeerStrains/FinalFilesPaper/Sfig1/OS1431_missense")
ggplot(dfstatusVariants, aes(x=AFvariant1431))+geom_bar(stat = "count")+theme_bw()+ylim(0,8000)
dev.off()

