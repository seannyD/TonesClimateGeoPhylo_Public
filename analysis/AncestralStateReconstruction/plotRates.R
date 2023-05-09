# Load node ancestra state probabilities and plot tree with tip data and edge humidity.
setwd("~/Documents/MPI/ClimateAndLanguage/Grollemund/analysis/AncestralStateReconstruction/")

source("../bayestraitr_helper.R")
library(diagram)
library(plotrix)
library(phylobase)
#BiocManager::install("ggtree")
library("ggtree")
library(ggimage)
library(phangorn)

d = bt_read.log("../../results/BayesTraitsOutput/AncestralStateReconstruction/Anc_Tones_3cat.Log.txt")

rates=apply(d[,4:9],2,mean)
rates = signif(rates,2)
widths = (rates/min(rates))

plot(c(0,4),c(0,3),type='n',xaxt='n',yaxt='n',xlab='',ylab='',bty='n')

N = c(1,1)
S = c(2,1)
C = c(3,1)

curvedarrow(N,S, curve = 0.2, lwd = widths["qNS"])
curvedarrow(S,C, curve = 0.2, lwd = widths["qSC"])
curvedarrow(N, C, curve = 0.5, lwd = widths["qNC"])
curvedarrow(S, N, curve = 0.2, lwd = widths["qSN"])
curvedarrow(C, S, curve = 0.2, lwd = widths["qCS"])
curvedarrow(C, N, curve = 0.5, lwd = widths["qCN"])

draw.circle(1,1,0.2,col='white')
draw.circle(2,1,0.2,col='white')
draw.circle(3,1,0.2,col='white')
text(1:3,c(1,1,1),c("N","S","C"))

text(1.5,0.7,rates["qNS"])
text(1.5,1.3,rates["qSN"])
text(2.5,0.7,rates["qSC"])
text(2.5,1.3,rates["qCS"])
text(2,0.2,rates["qNC"])
text(2,1.8,rates["qCN"])

load("../../data/processed/GTree_phylo4g_withClimateSimData.RDat")
tree.chosen = tree4d
tips.to.remove  = tree4d@label[is.na(tree4d@data[1:length(tree4d@label),]$Tones)]
tree.chosen = prune(tree.chosen,tips.exclude = tips.to.remove, trim.internal=T)
tree.chosen@data$TonesCat = cut(tree.chosen@data$Tones,c(-1,1,2,12),labels = c("N","S","C"))

anc = apply(d[,grepl("^A",names(d))],2,mean)
anc2 = matrix(anc,nrow=3)
colnames(anc2) = unique(sapply(names(anc),function(X){strsplit(X," ")[[1]][1]}))
rownames(anc2) = c("C","N","S")
mx = match(paste0("A",rownames(tree.chosen@data)),colnames(anc2))
tree.chosen@data$Tones.Anc.C = anc2["C",mx]
tree.chosen@data$Tones.Anc.N = anc2["N",mx]
tree.chosen@data$Tones.Anc.S = anc2["S",mx]
tree.chosen@data$node = as.numeric(rownames(tree.chosen@data))


# Cut humidity into 3 bins
humidityVar  = "specH.mean.sim.bilinear"
humBreaks = c(0,0.0125,0.0165,max(tree.chosen@data[,humidityVar]))
#humBreaks = seq(min(tree.chosen@data[,humidityVar],na.rm=T),max(tree.chosen@data[,humidityVar],na.rm=T),length.out=4)
#humBreaks = seq(min(tree.chosen@data$specH.mean.sim.bilinear,na.rm=T),
#                max(tree.chosen@data$specH.mean.sim.bilinear),length.out=4)
#plotmeans(tree.chosen@data$specH.mean.sim.bilinear~tree.chosen@data$TonesCat)
#abline(h=humBreaks)
hist(tree.chosen@data[,humidityVar])
abline(v=humBreaks,col=2)

tree.chosen@data$specH.mean.sim.bilinear.cat = 
  cut(tree.chosen@data[,humidityVar],
    breaks=humBreaks,
    include.lowest = T,labels = c("L","M","H"))

p <- ggtree(tree.chosen,aes(colour=specH.mean.sim.bilinear.cat)) 
p <- p + geom_tippoint(aes(colour=TonesCat))+ scale_color_manual(values=c(H="green", M="orange",L="red",C="green",S="orange",N="red"))
p


# Maximum likelihood ancestral tones
tree.chosen@data$Tones.Anc.Max = 
  factor(c("None","Simple","Complex")[apply(tree.chosen@data[,c("Tones.Anc.N","Tones.Anc.S","Tones.Anc.C")],
        1,function(X){which(X==max(X))[1]})],
        levels = c("None","Simple","Complex"))


p <- ggtree(tree.chosen,aes(colour=specH.mean.sim.bilinear.cat)) 
p <- p + geom_tippoint(aes(colour=TonesCat))+ scale_color_manual(values=c(H="green", M="orange",L="red",C="green",S="orange",N="red"))
pies <- nodepie(tree.chosen@data,c("Tones.Anc.N","Tones.Anc.S","Tones.Anc.C"),
                color=c(Tones.Anc.N="red",Tones.Anc.S="orange",Tones.Anc.C="green"), 
                alpha=.6)
pdf("../../results/BayesTraitsOutput/AncestralStateReconstruction/Phylogeny_Tones_Humidity_3cat_BayesTraits.pdf",height=15)
inset(p, pies,width=1,height=1)
dev.off()

# Just with max likelihood
dx = tree.chosen@data
for(i in 1:nrow(dx)){
  if(!is.na(dx[i,"Tones.Anc.N"])){
    X = dx[i,names(dx) %in% c("Tones.Anc.N","Tones.Anc.S","Tones.Anc.C")]
    if(max(X)>0.45){
      x = c(0,0,0)
      x[max.col(X,ties.method = "first")] = 1
      dx[i,names(dx) %in% c("Tones.Anc.N","Tones.Anc.S","Tones.Anc.C")] = x
    } else{
      dx[i,names(dx) %in% c("Tones.Anc.N","Tones.Anc.S","Tones.Anc.C")] = c(1/3,1/3,1/3)
    }
  }
}

p <- ggtree(tree.chosen,aes(colour=specH.mean.sim.bilinear.cat)) 
p <- p + geom_tippoint(aes(colour=TonesCat))+ scale_color_manual(values=c(H="green", M="orange",L="red",C="green",S="orange",N="red"))
pies <- nodepie(dx,c("Tones.Anc.N","Tones.Anc.S","Tones.Anc.C"),
                color=c(Tones.Anc.N="red",Tones.Anc.S="orange",Tones.Anc.C="green"), 
                alpha=.6)
pdf("../../results/BayesTraitsOutput/AncestralStateReconstruction/Phylogeny_Tones_Humidity_3cat_BayesTraits_MaxLikelihood.pdf",height=15)
inset(p, pies,width=1,height=1)
dev.off()


# Raw data
library(gplots)
plotmeans(specH.mean.sim.bilinear~Tones.Anc.Max,
          data = tree.chosen@data)

g = ggplot(tree.chosen@data,
       aes(y=specH.mean.sim.bilinear,
           x=Tones.Anc.N)) +
  geom_point(alpha=0.3,col="green") +
  geom_point(aes(y=specH.mean.sim.bilinear,
                 x=Tones.Anc.S),
             alpha=0.3,col="orange") +
  geom_point(aes(y=specH.mean.sim.bilinear,
                 x=Tones.Anc.C),
             alpha=0.3,col="red") +
  stat_smooth(colour="green") +
  stat_smooth(aes(y=specH.mean.sim.bilinear,
                  x=Tones.Anc.S),colour="orange") +
  stat_smooth(aes(y=specH.mean.sim.bilinear,
                  x=Tones.Anc.C),colour="red") +
  xlab("Ancestral state probability") +
  ylab("Humidity")

pdf("../../results/BayesTraitsOutput/AncestralStateReconstruction/ANC_vs_Hum.pdf", height=4.5,width=6)
g
dev.off()


# Rates of change
nodeData(tree.chosen)$change = NA
for(i in 1:nrow(nodeData(tree.chosen))){
  node = rownames(nodeData(tree.chosen))[i]
  anc = as.character(phylobase::ancestor(tree.chosen,as.numeric(node)))
  nodeProb = nodeData(tree.chosen)[i,c("Tones.Anc.N","Tones.Anc.S","Tones.Anc.C")]
  ancProb = nodeData(tree.chosen)[anc,c("Tones.Anc.N","Tones.Anc.S","Tones.Anc.C")]
  nodeData(tree.chosen)$change[i] = sum((nodeProb-ancProb)^2)
}

tree.chosen@data$change2 = 1-tree.chosen@data$change

bars <- nodepie(tree.chosen@data,c("change","change2"),alpha=.6)
pdf("../../results/BayesTraitsOutput/AncestralStateReconstruction/Phylogeny_Tones_Humidity_rate_BayesTraits.pdf",height=12)
inset(p, bars,width=0.005,height=20)
dev.off()

#----
# tx = as.phylo(tree.chosen)
# tx$tip.label = gsub("_","",tx$tip.label)
# dx = data.frame(taxa=gsub("_","",rownames(tipData(tree.chosen))),
#                 TonesCat= tipData(tree.chosen)[,"TonesCat"])
# 
# 
# pies <- nodepie(tree.chosen@data,c("Tones.Anc.N","Tones.Anc.S","Tones.Anc.C"),
#                 color=c("green","orange","red"), alpha=.6)
# 
# p <- ggtree(tx,aes(color=I(color)))
# p <- p %<+% dx +geom_tippoint(aes(color=TonesCat))  + scale_color_manual(values=c("green", "orange","red"))
# pdf("tmp.pdf",height=20)
# inset(p, pies,width=1,height=1)
# dev.off()