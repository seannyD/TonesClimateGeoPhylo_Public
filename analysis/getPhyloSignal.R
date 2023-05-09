library(phylobase)
library(ape)
library(gplots)
library(ggplot2)
#source("https://bioconductor.org/biocLite.R")
#biocLite("ggtree",type="source")
#library(ggtree)

library(phylosignal)
library(ggrepel)

setwd("~/OneDrive - Cardiff University/Research/MPI/ClimateAndLanguage/Grollemund/analysis")

load("../data/processed/GTree_phylo4g_withClimateSimData.RDat")
tree.chosen = tree4d
humidityVar = "specH.mean.sim.gam"

# Prune tree (remove where tones is NA)
#tips.to.remove  = tree.chosen@label[is.na(tree.chosen@data[1:length(tree.chosen@label),]$Tones)]
#tree.chosen = prune(tree.chosen,tips.exclude = tips.to.remove, trim.internal=T)

tree.chosen.min = tree.chosen
tree.chosen.min@data = tree.chosen.min@data[,c("Tones",humidityVar,"Vowels","Consonants","ASJPVowelRatio")]

pdf("../results/graphs/PhylogenyData.pdf", width=14, height=22)
dotplot.phylo4d(tree.chosen.min,grid.horizontal=F,center=F,scale=F)
dev.off()

dx = tree.chosen@data
dx$Tones[dx$Tones==1] = 0
dx = dx[!is.na(dx$Tones),]
pdf("../results/graphs/Tones_vs_Humidity.pdf",height=4,width=4.5)
ggplot(dx,aes(x=as.factor(Tones),y=specH.mean)) +
  geom_boxplot() +
  xlab("Tones") +
  ylab("Humidity")
dev.off()

set.seed(1289)
for(x in c("Tones",humidityVar,"Vowels","Consonants","ASJPVowelRatio")){
  tips.to.remove  = tree.chosen@label[is.na(tree.chosen@data[1:length(tree.chosen@label),x])]
  tx = prune(tree.chosen,tips.exclude = tips.to.remove, trim.internal=T)
  tx@data = tx@data[,c(x,humidityVar)]
  sig = phyloSignal(tx,reps=100000)
  print(x)
  print(sig)
}



# Same, but in phytools
p4 <- phylobase::extractTree(tree.chosen.min)
phy <- as(p4, "phylo")
phytools::phylosig(phy,tipData(tree.chosen.min)$specH.mean,test=T)
phytools::phylosig(phy,tipData(tree.chosen.min)$Tones,test=T)

# Phylogenetic signal for each node
sigINT.tones = phyloSignalINT(tree.chosen.min,method="K",trait = "Tones")
sigINT.tones@data$tipLabel = c(sigINT.tones@label,rep("",nrow(sigINT.tones@data) - length(sigINT.tones@label)))
p <- ggtree(sigINT.tones,layout = "circular") + 
  geom_nodepoint(aes(color=stat.K.Tones,size=stat.K.Tones)) +
  scale_color_gradientn(colors = RColorBrewer::brewer.pal(9, "YlOrRd")) +
  geom_tiplab(aes(label=Tones,color=Tones)) +
  geom_tiplab2(offset=0.003) 
p

sigINT.humidity = phyloSignalINT(tree.chosen.min,method="K",trait = humidityVar)
sigINT.humidity@data$K = sigINT.humidity@data[,paste0("stat.K.",humidityVar)]
p <- ggtree(sigINT.humidity,layout = "circular") + 
  geom_nodepoint(aes(color=K,size=K)) +
  scale_color_gradientn(colors = RColorBrewer::brewer.pal(9, "YlOrRd")) +
  geom_tippoint(aes(colour=Tones)) +
  geom_tiplab2(offset=0.003) 
p


p <- ggtree(sigINT.humidity,layout = "circular",mapping=aes(colour=specH.mean.sim.bilinear))
p


pc = phyloCorrelogram(tree.chosen.min)

# Local Indicators of Spatial Association
# One simple and well‐described LISA is the local Moran's I (Equation (3)), noted Ii (Anselin 1995), which can be used to detect hotspots of positive and negative autocorrelation. The same statistic can be applied to phylogenetic data to detect species with similar neighbors and species with different neighbors. In this context, we call these indicators Local Indicators of Phylogenetic Association (LIPA), for sake of consistency in terminology, although the statistic remains the same. 
local.i = lipaMoran(tree.chosen.min, trait = "Tones",prox.phylo = "nNodes", as.p4d = TRUE)
points.col = lipaMoran(tree.chosen.min, trait = "Tones", prox.phylo = "nNodes")$p.value
points.col = ifelse(points.col < 0.05, "red", "black")
dotplot.phylo4d(local.i, dot.col = points.col)


phyloCorrelogram(tree.chosen.min,c("Tones",humidityVar))


sig.int = phyloSignalINT(tree.chosen.min, trait = "Tones")
