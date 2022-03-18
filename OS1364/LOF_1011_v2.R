#Counting LOF in 1011------


pathLOF <- "/home/simone/Scrivania/Scrivania/NewBeerStrains/LOF_1011/all_LOF.csv"
pathmeta <- "/home/simone/Scrivania/Scrivania/NewBeerStrains/LOF_1011/supplementary_data_S1.csv"

LOF_all <- read.csv(pathLOF, sep="")
metadata <- read.csv(pathmeta)

LOFwgenes <- LOF_all[1:6031, 8:1018]
LOFfind <- c()

for (indSample in 1:ncol(LOFwgenes)) {
  
  countLOFsamp <- table(LOFwgenes[, indSample])[2]
  LOFfind <- rbind(LOFfind, cbind(indSample, countLOFsamp))
  
  
}

LOFfinddf <- data.frame(LOFfind)
LOFfinddf2 <- cbind(LOFfinddf, metadata$Clade_class, metadata$Ploidy)
colnames(LOFfinddf2) <- c("indSample","countLOFsamp","Origin","Ploidy")
library(ggplot2)

pdf("/home/simone/Scrivania/Scrivania/NewBeerStrains/FinalFilesPaper/Sfig1/LOFall.pdf", height = 4, width = 6)
ggplot(LOFfinddf2, aes(x=countLOFsamp,  col=Origin))+geom_density()+
  geom_vline(xintercept = 308)+geom_vline(xintercept = 188,col=c("red"))+theme_bw()+xlim(0,600)#+facet_wrap(~Ploidy)
dev.off()

ggplot(LOFfinddf2, aes(x=countLOFsamp,  col=Origin))+geom_boxplot()+
   geom_vline(xintercept = 308)+geom_vline(xintercept = 188,col=c("red"))+theme_bw()+facet_wrap(~Ploidy, ncol = 5)+
  xlim(0,500)


dom <-LOFfinddf2[which(LOFfinddf2$Origin == "dom"),]

polypl <- dom[which(dom$Ploidy %in% c(3,4)),]
dip <- dom[which(dom$Ploidy %in% c(1,2)),]
median(dom$countLOFsamp)
mean(dom$countLOFsamp)
sd(dom$countLOFsamp)

median(polypl$countLOFsamp)
mean(polypl$countLOFsamp)
sd(polypl$countLOFsamp)
pdf("/home/simone/Scrivania/Scrivania/NewBeerStrains/FinalFilesPaper/Sfig1/LOFallpolyp.pdf", height = 4, width = 6)
ggplot(polypl, aes(x=countLOFsamp))+geom_density()+
  geom_vline(xintercept = 308)+geom_vline(xintercept = 188,col=c("red"))+theme_bw()+xlim(0,500)
dev.off()
pdf("/home/simone/Scrivania/Scrivania/NewBeerStrains/FinalFilesPaper/Sfig1/LOFallDip.pdf", height = 4, width = 6)
ggplot(dip, aes(x=countLOFsamp))+geom_density()+
  geom_vline(xintercept = 308)+geom_vline(xintercept = 188,col=c("red"))+theme_bw()+xlim(0,500)
dev.off()

median(LOFfinddf$countLOFsamp)
mean(LOFfinddf$countLOFsamp)
sd(LOFfinddf$countLOFsamp)




