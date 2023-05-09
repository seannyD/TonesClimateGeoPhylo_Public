library(phylobase)
library(ape)
library(gplots)
#source("https://bioconductor.org/biocLite.R")
#biocLite("ggtree",type="source")
#library(ggtree)


setwd("~/Documents/MPI/ClimateAndLanguage/Grollemund/analysis")

load("../data/processed/GTree_phylo4g_withClimateSimData.RDat")
tree.chosen = tree4d
humidityVar = "specH.mean.sim.bilinear"

# Prune tree
# This breaks if you give it too many things to trim at once, so do it the slow way:
tips.to.remove  = tree4d@label[is.na(tree4d@data[1:length(tree4d@label),]$Tones)]
tree.chosen = prune(tree.chosen,tips.exclude = tips.to.remove, trim.internal=T)
# for(t in tips.to.remove){
#   #print(t)
#   tree.chosen = phylobase::prune(tree.chosen,tips.exclude = t, trim.internal=F)
# }


pdf("../results/graphs/MissingGeoData.pdf", width = 10, height = 30)
plot(as(tree.chosen,'phylo'), show.tip.label=T)
nodelabels('x',which(is.na(tree.chosen@data[,humidityVar])))
#tiplabels()
dev.off()

ASE = ace(tipData(tree.chosen)$Tones, as(tree.chosen,'phylo'), type="continuous", method="REML")

nodeData(tree.chosen)$Tones = ASE$ace

###############
# Plots
##############

humidityVar = "specH.mean.sim.bilinear"

plotmeans(tree.chosen@data[,humidityVar]~round(tree.chosen@data$Tones),ylim=range(tree.chosen@data[,humidityVar],na.rm=T))
sapply(seq(min(tree.chosen@data[,humidityVar],na.rm=T),max(tree.chosen@data[,humidityVar],na.rm=T),length.out=4), function(X){abline(h=X)})

plotmeans(tree.chosen@data[,humidityVar]~cut(tree.chosen@data$Tones,c(0,2,3,8)),ylim=range(tree.chosen@data[,humidityVar],na.rm=T))

hist(tree.chosen@data[,humidityVar],breaks=20)
sapply(seq(min(tree.chosen@data[,humidityVar],na.rm=T),max(tree.chosen@data[,humidityVar],na.rm=T),length.out=4), function(X){abline(v=X,col=3)})

num.hum.colours = 10
hum.colour.scale = rev(terrain.colors(num.hum.colours))[2:num.hum.colours]
edge.specH.mean = tree.chosen@data[as.numeric(tree.chosen@edge[,2]),humidityVar]
edge.specH.mean.col = hum.colour.scale[as.numeric(cut(edge.specH.mean,length(hum.colour.scale)),include.lowest=T)]
edge.specH.mean.col[is.na(edge.specH.mean.col)] = "black"


ncat = min(tipData(tree.chosen)$Tones):max(tipData(tree.chosen)$Tones)
breaks = c(ncat,max(ncat)+1)


pdf.width = 10
pdf.height = 30

pdf("../results/graphs/Phylogeny_Tones_Humidity.pdf", pdf.width, pdf.height)
plot.phylo(as(tree.chosen,'phylo'), edge.color = edge.specH.mean.col, adj = 0.05)

tip.cols = heat.colors(max(ncat))[as.numeric(cut(tipData(tree.chosen)$Tones, breaks=breaks, include.lowest = T))]
tiplabels(pch=16, col= tip.cols, width = 0.3)

node.cols = heat.colors(max(ncat))[as.numeric(cut(ASE$ace, breaks=breaks, include.lowest = T))]
nodelabels(pch=16,col=node.cols, width = 0.3)

legend(0,130, c("low",'med','high'), col=heat.colors(3), pch=16, title="Tone")
legend(0.01,130, c("low",'med','high'), col=hum.colour.scale[c(1,ceiling(length(hum.colour.scale)/2),length(hum.colour.scale))], pch=16, title = "Humidity")
dev.off()


# 3 Categories
hum.colour.scale = c("brown",'orange','green')
edge.specH.mean = tree.chosen@data[as.character(tree.chosen@edge[,2]),humidityVar]
#edge.specH.mean.col = hum.colour.scale[cut(edge.specH.mean,length(hum.colour.scale))]
edge.specH.mean.col = hum.colour.scale[cut(edge.specH.mean,seq(min(tree.chosen@data[,humidityVar],na.rm=T),max(tree.chosen@data[,humidityVar],na.rm=T),length.out=4), include.lowest=T)]
edge.specH.mean.col[is.na(edge.specH.mean.col)] = "black"


breaks2 = c(-1,0,2,4,8)
pointCols = rev(rainbow(length(breaks2-1)))
pointCols =  c("brown",'orange','yellow','green')
pdf("../results/graphs/Phylogeny_Tones_Humidity_3cat.pdf", pdf.width, pdf.height)
plot.phylo(as(tree.chosen,'phylo'), edge.color = edge.specH.mean.col, edge.width = 2, adj = 0.05)
tip.cols = pointCols[as.numeric(cut(tipData(tree.chosen)$Tones, breaks=breaks2, include.lowest = T))]
tiplabels(pch=16, col= tip.cols, width = 0.3)
node.cols = pointCols[as.numeric(cut(ASE$ace, breaks=breaks2, include.lowest = T))]
nodelabels(pch=16,col=node.cols, width = 0.3)
legend(0,130, c("low",'med','high'), col=pointCols, pch=16, title="Tone")
legend(0.01,130, c("low",'med','high'), col=hum.colour.scale[c(1,ceiling(length(hum.colour.scale)/2),length(hum.colour.scale))], pch=16, title = "Humidity")
dev.off()


breaks2 = c(-1,0,2,8)
hum.colour.scale = c("brown",'orange','green')
pointCols = hum.colour.scale
pdf("../results/graphs/Phylogeny_Tones_Humidity_3cat_sameColour.pdf", pdf.width, pdf.height)
plot.phylo(as(tree.chosen,'phylo'), edge.color = edge.specH.mean.col, edge.width = 2, adj = 0.05)
tip.cols = pointCols[as.numeric(cut(tipData(tree.chosen)$Tones, breaks=breaks2, include.lowest = T),include.lowest=T)]
tiplabels(pch=16, col= tip.cols, width = 0.3)
node.cols = pointCols[as.numeric(cut(ASE$ace, breaks=breaks2, include.lowest = T))]
nodelabels(pch=16,col=node.cols, width = 0.3, cex=2)
legend(0,130, c("low",'med','high'), col=pointCols, pch=16, title="Tone")
legend(0.01,130, c("low",'med','high'), col=hum.colour.scale[c(1,ceiling(length(hum.colour.scale)/2),length(hum.colour.scale))], pch=16, title = "Humidity")
dev.off()

# Draw the same graph, but without colouring the edges
pdf("../results/graphs/Phylogeny_Tones_Humidity_3cat_sameColour_TonesOnly.pdf", pdf.width, pdf.height)
plot.phylo(as(tree.chosen,'phylo'), edge.color = 1, edge.width = 2, adj = 0.05)
tip.cols = pointCols[as.numeric(cut(tipData(tree.chosen)$Tones, breaks=breaks2, include.lowest = T))]
tiplabels(pch=16, col= tip.cols, width = 0.3)
node.cols = pointCols[as.numeric(cut(ASE$ace, breaks=breaks2, include.lowest = T))]
nodelabels(pch=16,col=node.cols, width = 0.3, cex=2)
legend(0,130, c("low",'med','high'), col=pointCols, pch=16, title="Tone")
legend(0.01,130, c("low",'med','high'), col=hum.colour.scale[c(1,ceiling(length(hum.colour.scale)/2),length(hum.colour.scale))], pch=16, title = "Humidity")
dev.off()

# Draw the same graph, but with no colours
pdf("../results/graphs/Phylogeny_Tones_Humidity_3cat_sameColour_BNW.pdf", pdf.width, pdf.height)
plot.phylo(as(tree.chosen,'phylo'), edge.color = 1, edge.width = 2, adj = 0.05)
tip.cols = pointCols[as.numeric(cut(tipData(tree.chosen)$Tones, breaks=breaks2, include.lowest = T))]
tiplabels(pch=16, col= tip.cols, width = 0.3, cex=2)
#node.cols = pointCols[as.numeric(cut(ASE$ace, breaks=breaks2, include.lowest = T))]
#nodelabels(pch=16,col=node.cols, width = 0.3)
legend(0,130, c("low",'med','high'), col=pointCols, pch=16, title="Tone")
legend(0.01,130, c("low",'med','high'), col=hum.colour.scale[c(1,ceiling(length(hum.colour.scale)/2),length(hum.colour.scale))], pch=16, title = "Humidity")
dev.off()

#ggtree(tree.chosen, aes(color=Tones)) +
#  scale_color_continuous(low='darkgreen', high='red') +
#  theme(legend.position="right")


## Draw geo tree


breaks2 = c(-1,0,2,8)
pointCols = c("brown",'orange','green')

tree.chosen@data$Tones.colour = pointCols[as.numeric(cut(tree.chosen@data$Tones, breaks=breaks2, include.lowest = T))]

edge.specH.mean = tree.chosen@data[as.character(tree.chosen@edge[,2]),humidityVar]
#edge.specH.mean.col = hum.colour.scale[cut(edge.specH.mean,length(hum.colour.scale))]
tree.chosen@data$edge.specH.mean.col = pointCols[cut(edge.specH.mean,seq(min(tree.chosen@data[,humidityVar],na.rm=T),max(tree.chosen@data[,humidityVar],na.rm=T),length.out=4), include.lowest=T)]



library(maps)
map(xlim = c(0,60),ylim=c(-40,20))
points(tree.chosen@data$Longitude,
       tree.chosen@data$Latitude,
       col=tree.chosen@data$Tones.colour,
       pch=16
       )
edges = phylobase::edges(tree.chosen)
segments(tree.chosen@data[edges[,1],]$Longitude,
       tree.chosen@data[edges[,1],]$Latitude,
       tree.chosen@data[edges[,2],]$Longitude,
       tree.chosen@data[edges[,2],]$Latitude,
       col=tree.chosen@data[edges[,1],]$edge.specH.mean.col)

