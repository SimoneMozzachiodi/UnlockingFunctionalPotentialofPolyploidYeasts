#Fermentation kinetic 2L--------


path1364 <- "/home/simone/Scrivania/Scrivania/NewBeerStrains/FinlandPhenotyping/Sugar_consumption_Fermentation/2Lfermentation_OS1364.csv"

dataFerm2l1364 <- read.csv(path1364, sep="\t")

#-----SummarySE------
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

FermPlot1364 <- summarySE(dataFerm2l1364, measurevar="ABV.",  groupvars=c("Strain","Time") )

library(ggplot2)
pdf("/home/simone/Scrivania/Scrivania/NewBeerStrains/FinalFilesPaper/Sfig8/Fermentation2L_OS1364.pdf")
ggplot(FermPlot1364, aes(x=Time, y=ABV., fill=Strain, col=Strain))+ 
  geom_errorbar(aes(ymin=ABV.-se, ymax=ABV.+se), width=1, col="black") +
  geom_point(size=0.5,  position=position_dodge(width=0.5), aes(fill=Strain))+
  geom_line(position=position_dodge(width=1))+xlab("Time (hours)")+ 
  ylab("ABV(%)")+
  theme_bw()+
  theme(legend.text = element_text(size=8),  axis.title.x = element_text(size=8), axis.text.y  = element_text(size=5), axis.text.x= element_text(size=8, angle=45, vjust = 0.5),axis.title.y=element_text(size=8))+
  ylim(0,10)+xlim(-10,300)+
  theme(strip.text.x = element_text(size = 8))+scale_fill_manual(values=c('#c6eaae','#7b916c',"forestgreen","black"))
dev.off()

pdf("/home/simone/Scrivania/Scrivania/NewBeerStrains/FinalFilesPaper/Sfig8/Fermentation2L_OS1364_zoom.pdf")
ggplot(FermPlot1364, aes(x=Time, y=ABV., fill=Strain, col=Strain))+ 
  geom_errorbar(aes(ymin=ABV.-se, ymax=ABV.+se), width=1, col="black") +
  geom_point(size=0.5,  position=position_dodge(width=0.5), aes(fill=Strain))+
  geom_line(position=position_dodge(width=1))+xlab("Time (hours)")+ 
  ylab("ABV(%)")+
  theme_bw()+
  theme(legend.text = element_text(size=8),  axis.title.x = element_text(size=8), axis.text.y  = element_text(size=5), axis.text.x= element_text(size=8, angle=45, vjust = 0.5),axis.title.y=element_text(size=8))+
  ylim(6,10)+xlim(80,300)+
  theme(strip.text.x = element_text(size = 8))+scale_fill_manual(values=c('#c6eaae','#7b916c',"forestgreen","black"))
dev.off()


#Fermentation kinetic 2L 1431--------


path1431 <- "/home/simone/Scrivania/Scrivania/NewBeerStrains/FinlandPhenotyping/Sugar_consumption_Fermentation/2Lfermentation_OS1431.csv"

dataFerm2l1431 <- read.csv(path1431, sep=",")

#-----SummarySE------
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

FermPlot1431 <- summarySE(dataFerm2l1431, measurevar="ABV.",  groupvars=c("Strain","Time") )

library(ggplot2)
pdf("/home/simone/Scrivania/Scrivania/NewBeerStrains/FinalFilesPaper/Sfig8/Fermentation2L_OS1431.pdf")
ggplot(FermPlot1431, aes(x=Time, y=ABV., fill=Strain, col=Strain))+ 
  geom_errorbar(aes(ymin=ABV.-se, ymax=ABV.+se), width=1, col="black") +
  geom_point(size=0.5,  position=position_dodge(width=0.5), aes(fill=Strain))+
  geom_line(position=position_dodge(width=1))+xlab("Time (hours)")+ 
  ylab("ABV(%)")+
  theme_bw()+
  theme(legend.text = element_text(size=8),  axis.title.x = element_text(size=8), axis.text.y  = element_text(size=5), axis.text.x= element_text(size=8, angle=45, vjust = 0.5),axis.title.y=element_text(size=8))+
  ylim(0,10)+xlim(-10,300)+
  theme(strip.text.x = element_text(size = 8))+scale_fill_manual(values=c('#c6eaae','#7b916c',"forestgreen","black"))
dev.off()

pdf("/home/simone/Scrivania/Scrivania/NewBeerStrains/FinalFilesPaper/Sfig8/Fermentation2L_OS1431_zoom.pdf")
ggplot(FermPlot1431, aes(x=Time, y=ABV., fill=Strain, col=Strain))+ 
  geom_errorbar(aes(ymin=ABV.-se, ymax=ABV.+se), width=1, col="black") +
  geom_point(size=0.5,  position=position_dodge(width=0.5), aes(fill=Strain))+
  geom_line(position=position_dodge(width=1))+xlab("Time (hours)")+ 
  ylab("ABV(%)")+
  theme_bw()+
  theme(legend.text = element_text(size=8),  axis.title.x = element_text(size=8), axis.text.y  = element_text(size=5), axis.text.x= element_text(size=8, angle=45, vjust = 0.5),axis.title.y=element_text(size=8))+
  ylim(8,9.5)+xlim(150,300)+
  theme(strip.text.x = element_text(size = 8))+scale_fill_manual(values=c('#c6eaae','#7b916c',"forestgreen","black"))
dev.off()













