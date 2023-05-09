try(setwd("~/OneDrive - Cardiff University/Research/MPI/ClimateAndLanguage/Grollemund/analysis/AncestralStateReconstruction/"))

source("../bayestraitr_helper.R")
library(diagram)
library(plotrix)
library(phylobase)
#BiocManager::install("ggtree")
library("ggtree")
library(ggimage)
library(phangorn)

load("../../data/processed/GTree_phylo4g_withClimateSimData.RDat")
tree.chosen = tree4d

# Prune tree
# This breaks if you give it too many things to trim at once, so do it the slow way:
tips.to.remove  = tree4d@label[is.na(tree4d@data[1:length(tree4d@label),]$Tones)]
tree.chosen = prune(tree.chosen,tips.exclude = tips.to.remove, trim.internal=T)

# Load BT output of ancestral state reconstruction
d = bt_read.log("../../results/BayesTraitsOutput/AncestralStateReconstruction/Anc_Tones_continuous_est.Log.txt")

# Add estimates to data
anc = apply(d[,grepl("^Est",names(d))],2,mean)
names(anc) = sapply(names(anc),function(X){strsplit(X," ")[[1]][2]})
mx = match(paste0("A",rownames(tree.chosen@data)),names(anc))
tree.chosen@data$Tones.Anc = anc[mx]
#tree.chosen@data$node = as.numeric(rownames(tree.chosen@data))

# Categorical humidity
humidityVar  = "specH.mean.sim.bilinear"
humBreaks = c(0,0.0125,0.0165,max(tree.chosen@data[,humidityVar]))
tree.chosen@data$specH.mean.sim.bilinear.cat = 
  cut(tree.chosen@data[,humidityVar],
      breaks=humBreaks,
      include.lowest = T,labels = c("L","M","H"))

# Break tones into categories
tree.chosen@data$Tones.Anc.cat = cut(tree.chosen@data$Tones.Anc,breaks=c(-1,0.5,1.5,2.5,3.5,100),labels = c(0,1,2,3,4))
tree.chosen@data$Tones.cat = cut(tree.chosen@data$Tones,breaks=c(-1,0.5,1.5,2.5,3.5,100),labels = c(0,1,2,3,4))

p <- ggtree(tree.chosen) +  
  geom_tippoint(aes(colour=Tones.cat)) +
  geom_nodepoint(aes(colour=Tones.Anc.cat)) +
  scale_color_manual(values=c("0"="blue", "1"="green","2"="orange","3"="red","4"="dark red"))
p


p2 <- ggtree(tree.chosen) +  
  geom_nodepoint(aes(colour=Tones.Anc)) +
  geom_tippoint(aes(colour=Tones)) +
  scale_colour_gradientn(colours = c(rev(terrain.colors(10)),"#00A600","#00A600","#00A600"))
p2


dV = bt_read.log("../../results/BayesTraitsOutput/AncestralStateReconstruction/Anc_VOWELS_continuous_est.Log.txt")
