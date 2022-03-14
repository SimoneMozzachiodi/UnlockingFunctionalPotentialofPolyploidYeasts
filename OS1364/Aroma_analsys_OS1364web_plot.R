###### Aroma analysis plot #####

aromaPath <- "/home/simone/Scrivania/Scrivania/NewBeerStrains/FinlandPhenotyping/Aroma_Analysis/AromaResults.csv"

aromaData <- read.csv(aromaPath, sep = ",")

library(fmsb)

colors_border=c( rgb(0.2,0.5,0.5,0.9), rgb(0.8,0.2,0.5,0.9), rgb(0.7,0.5,0.1,0.9) )
colors_in=c( rgb(0.2,0.5,0.5,0.2), rgb(0.8,0.2,0.5,0.2), rgb(0.7,0.5,0.1,0.2) )

commdata1364 <- rbind(aromaData[7 , 2:13], aromaData[1, 2:13], aromaData[3, 2:13])
commWeb1364 <- data.frame(rbind(rep(0,12), rep(100,12), commdata1364))

rownames(commWeb1364) <- c("1","2","Commercial","OS1364","RTG1")
pdf("/home/simone/Scrivania/Scrivania/NewBeerStrains/FinalFilesPaper/Sfig8/Aroma1364.pdf",width = 6, height = 6)
radarchart(commWeb1364, pcol = colors_border)+
  legend( legend = rownames(commWeb1364[-c(1,2),]), x=1, y=1, bty = "n", pch=20 , col=colors_in , text.col = "grey", cex=1.2, pt.cex=3)
dev.off()
