setwd("~/Documents/MPI/ClimateAndLanguage/Grollemund/analysis/Discrete_Fossilised//")

source("../bayestraitr_helper.R")
library(diagram)
library(plotrix)
library(phylobase)
#BiocManager::install("ggtree")
library("ggtree")
library(ggimage)
library(phangorn)


d_all = bt_read.log("../../results/BayesTraitsOutput/Discrete_Fossilised/FossilCorB.Log.txt")

plot(d_all$`q34`)
#d = d[100:nrow(d),]
hist(d_all$`q34`)

getLogMarginalLikelihoodFromStonesFile = function(f){
  a = readLines(f)
  return(as.numeric(strsplit(tail(a,1),"\t")[[1]][2]))
}

# indep
logMarginalLikelihood_a = getLogMarginalLikelihoodFromStonesFile("../../results/BayesTraitsOutput/Discrete_Fossilised/FossilCorA.Stones.txt")
# dependent
logMarginalLikelihood_b = getLogMarginalLikelihoodFromStonesFile("../../results/BayesTraitsOutput/Discrete_Fossilised/FossilCorB.Stones.txt")

LogBF = 2*(logMarginalLikelihood_b - logMarginalLikelihood_a)
LogBF



####

d=apply(d_all[,4:11],2,mean)
d = signif(d)
widths = (d/5)


pdf("../../results/BayesTraitsOutput/Discrete_Fossilised/DependentModel.pdf")
plot(c(0.5,2.5),c(0.5,2.5),type='n',xaxt='n',yaxt='n',xlab='',ylab='',bty='n')

# humidity / tones
LN = c(1,1)
HN = c(2,1)
LY = c(1,2)
HY = c(2,2)

# 00 -> 01
curvedarrow(LN,LY, curve = 0.2, lwd = widths["q12"],arr.width=0.5,lcol="green")
curvedarrow(LY,LN, curve = 0.2, lwd = widths["q21"],arr.width=0.5,lcol="red")
# 00 -> 10
curvedarrow(LN,HN, curve = 0.2, lwd = widths["q13"],arr.width=0.5)
curvedarrow(HN, LN, curve = 0.2, lwd = widths["q31"],arr.width=0.5)
# 01 -> 1,1
curvedarrow(LY, HY, curve = 0.2, lwd = widths["q24"],arr.width=0.5)
curvedarrow(HY, LY, curve = 0.2, lwd = widths["q42"],arr.width=0.5)
# 10 -> 1,1
curvedarrow(HN, HY, curve = 0.2, lwd = widths["q34"],arr.width=0.5,,lcol="green")
curvedarrow(HY, HN, curve = 0.2, lwd = widths["q43"],arr.width=0.5,lcol="red")

draw.circle(1,1,0.2,col='#ffdaa3')
draw.circle(2,1,0.2,col='#acffa3')
draw.circle(1,2,0.2,col='#ffdaa3')
draw.circle(2,2,0.2,col='#acffa3')
text(c(1,2,1,2),c(1,1,2,2),c("Non-complex\nTones","Non-complex\nTones","Complex\nTones","Complex\nTones"))
text(c(1,2),c(2.5,2.5),c("Dry","Humid"))

# text(1.5,0.7,rates["qNS"])
# text(1.5,1.3,rates["qSN"])
# text(2.5,0.7,rates["qSC"])
# text(2.5,1.3,rates["qCS"])
# text(2,0.2,rates["qNC"])
# text(2,1.8,rates["qCN"])
dev.off()

library(ggplot2)

dx = data.frame(
  Estimate = c(d_all$q12,d_all$q34),
  Transition = c(rep(c("GainTones(dry)","GainTones(humid)"),each=nrow(d_all)))
)

pdf("../../results/BayesTraitsOutput/Discrete_Fossilised/GainTones.pdf",
    width=4,height=3)
ggplot(mapping=aes(x=Estimate)) +
  geom_histogram(data=dx[dx$Transition=="GainTones(dry)",],fill="#99660080") +
  geom_histogram(data=dx[dx$Transition=="GainTones(humid)",],fill="#16a60080")
dev.off()


dx = data.frame(
  Estimate = c(d_all$q21,d_all$q43),
  Transition = c(rep(c("LoseTones(dry)","LoseTones(humid)"),each=nrow(d_all)))
)

pdf("../../results/BayesTraitsOutput/Discrete_Fossilised/LoseTones.pdf",
    width=4,height=3)
ggplot(mapping=aes(x=Estimate)) +
  geom_histogram(data=dx[dx$Transition=="LoseTones(dry)",],fill="#99660080") +
  geom_histogram(data=dx[dx$Transition=="LoseTones(humid)",],fill="#16a60080")
dev.off()

