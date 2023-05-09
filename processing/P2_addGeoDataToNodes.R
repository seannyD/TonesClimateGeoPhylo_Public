# https://www.nature.com/articles/sdata20146


library(ape)
library(phytools)
#library(caper)
#library(nlme)
#library(geiger)
library(phylobase)
library(phangorn)


try(setwd("~/OneDrive - Cardiff University/Research/MPI/ClimateAndLanguage/Grollemund/analysis"))

tree = read.nexus("../data/trees/BP424_M1P_100_cv2_relaxed.trees")

tdata = read.delim("../data/trees/data.txt", stringsAsFactors = F)
tdata$Node.Name = gsub(" ",'',tdata$Node.Name)

#phylomorphospace
#http://blog.phytools.org/2012/02/mrca-for-set-of-taxa.html

parseNodeList = function(X){
  X  = substring(X,2,nchar(X)-1)
  X = gsub("\\'",'',X)
  strsplit(X,", ")
}

nodeList = parseNodeList(tdata$Node.List)

#nodeList = sapply(nodeList, function(X){
#  X[X %in% tree$tip.label]
#})

# TODO: The tree has more tips than the data file. 

x = tdata[tdata$Node.Group!="Int",]$Node.Name 
x[!x %in% tree$tip.label]
# These nodes are missing?
missing = tree$tip.label[!tree$tip.label %in% x]

# Unless we drop these tips, the tree is badly conformed
tree = drop.tip(tree, missing)

dim(tdata[tdata$Node.Group=="Int",])

internalNodes = unique(tree$edge[,1])
A = matrix(nrow=length(internalNodes),ncol=3)
rownames(A) = paste("Node",internalNodes,sep='-')

tip.data = matrix(nrow=length(tree$tip.label), ncol=3)
rownames(tip.data) = tree$tip.label

# all tips in DATA (some tree tips are not in data)
alltips = tdata[tdata$Node.Child.1==-2,]$Node.Name

mrca.tree = mrca(tree)

for(i in (length(nodeList):1)){
  thisNodeList = nodeList[[i]]
  if(length(thisNodeList)>1){
    # (some tips are not in data)
    thisNodeList = thisNodeList[thisNodeList %in% alltips]
    # xnode = findMRCA(tree,tips = thisNodeList)
    xnode = mrca.tree[thisNodeList,thisNodeList]
    xnode = min(xnode[upper.tri(xnode)])
    A[paste("Node",xnode,sep='-'),1] = tdata[i,]$Node.Long
    A[paste("Node",xnode,sep='-'),2] = tdata[i,]$Node.Lat
    A[paste("Node",xnode,sep='-'),3] = tdata[i,]$Node.Age
  } else{
    tip.data[nodeList[[i]][1],1] = tdata[i,]$Node.Long
    tip.data[nodeList[[i]][1],2] = tdata[i,]$Node.Lat
    tip.data[nodeList[[i]][1],3] = tdata[i,]$Node.Age
  }
}

A2 = A
rownames(A2) = as.character(sapply(rownames(A),function(X){strsplit(X,"-")[[1]][2]}))


###
# Patch together the missing geo data - find the closest match and set the node data to that

A3 = A2
missingGeoData = as.numeric(rownames(A2[is.na(A2[,1]),]))

for(i in missingGeoData){
  # get desendents of the missing node
  missing.node.children = tree$tip.label[Descendants(tree,i)[[1]]]
  # find a record in the data that matches this list best
  num.matches = sapply(nodeList, function(X){
    # (essentially number of edits)
    length(setdiff(missing.node.children,X)) + length(setdiff(X,missing.node.children))
  })
  # find the best matching node
  node = which(num.matches==min(num.matches))
  # if there are multiple, choose the largest node
  tx = tdata[node,]
  tx = tx[tx$Node.Count== max(tx$Node.Count),][1,]
  # Add the data
  A3[as.character(i),1] = tx$Node.Long
  A3[as.character(i),2] = tx$Node.Lat
  A3[as.character(i),3] = tx$Node.Age
}

############
# Make phylo4d object

tree4d = phylo4d(tree, tip.data = tip.data, node.data = A3)


#####################
# Add humidity data

# 1 degree data:
#http://iridl.ldeo.columbia.edu/SOURCES/.NOAA/.ESRL/.PSD/.reforecast2/.ensemble_mean/.spfh_2m/Y/%2820N%29%2840S%29RANGEEDGES/X/%280%29%2850E%29RANGEEDGES/ngridtable/5+ncoltable.html
# SOURCES .NOAA .ESRL .PSD .reforecast2 .ensemble_mean .spfh_2m
#Y (20N) (40S) RANGEEDGES
#X (0) (50E) RANGEEDGES

#http://iridl.ldeo.columbia.edu/SOURCES/.NOAA/.ESRL/.PSD/.reforecast2/.ensemble_mean/.spfh_2m/Y/%2820N%29%2840S%29RANGEEDGES/X/%280%29%2850E%29RANGEEDGES/L/%2812%29%2812%29RANGEEDGES/[X]data.tsv

#http://iridl.ldeo.columbia.edu/SOURCES/.NOAA/.NCEP-DOE/.Reanalysis-2/.Monthly/.flx/.flx_climo_7998/.spfh2m/

cols = 190
rows = 94
nrows = 3572
rotate <- function(x) t(apply(x, 2, rev))

datax = tree4d@data
names(datax)[1] = "Longitude"
names(datax)[2] = "Latitude"
names(datax)[3] = "Node.Age"
datax$Latitude.grid = as.numeric(cut(datax$Latitude,seq(90,-90,length.out=rows)))
datax$Longitude.grid = as.numeric(cut(datax$Longitude,seq(-180,180,length.out=cols)))  

h.mean = read.csv("../data/climate/MeanMonthlyMeans.csv", stringsAsFactors = F)
h.mean = as.matrix(h.mean[,2:ncol(h.mean)])

datax$specH.mean = apply(
  cbind(datax$Latitude.grid,datax$Longitude.grid),1,
  function(X){h.mean[X[2],X[1]]})

colx = heat.colors(10)[as.numeric(cut(datax$specH.mean,10))]
map(ylim=c(-35,10), xlim=c(5,50), interior = F)
points(datax$Longitude, datax$Latitude, col=colx, pch=16)

# put data back in tree

tree4d@data = datax

# Add language data
dx = read.csv("../data/processed/tones_combined.csv", stringsAsFactors = F)
# Remove langs with no tone and no vowel data
dx = dx[!(is.na(dx$Tones) & is.na(dx$Vowels) & is.na(dx$ASJPVowelRatio)),]
langData = dx[match(rownames(tipData(tree4d)), dx$Name_on_tree_tip2), c('Language','glottolog_id',"Tones","Phonemes","Consonants","Vowels","ASJPVowelRatio",'Lang.Source','InventoryID')]
tipData(tree4d) = cbind(tipData(tree4d), langData)

# Save tree

save(tree4d, file="../data/processed/GTree_phylo4g.RDat")
#load("../data/processed/GTree_phylo4g.RDat")

#Ax = A
#Ax[is.na(Ax[,1]),1] = mean(A[,1],na.rm=T)
#Ax[is.na(Ax[,2]),2] = mean(A[,2],na.rm=T)


# Prune tree to languages for which we have any language data
# (may still have missing data for tones or vowels etc.)
tips.to.remove  = tree4d@label[!tree4d@label %in% dx$Name_on_tree_tip2]

tree2 = tree4d

# This breaks if you give it too many things to trim at once, so do it the slow way:
for(t in tips.to.remove){
  print(t)
  tree2 = prune(tree2,tips.exclude = t, trim.internal=F)
}

# Add tone data

langData = dx[match(rownames(tipData(tree2)), dx$Name_on_tree_tip2), c('Language','glottolog_id',"Tones","Phonemes","Consonants","Vowels","ASJPVowelRatio",'Lang.Source','InventoryID')]

tipData(tree2) = cbind(tipData(tree2), langData)

save(tree2, file="../data/processed/GTree_phylo4g_combined.Rdat")

# Alt: Trim internal nodes

tree3 = prune(tree4d,tips.exclude = tips.to.remove, trim.internal=T)

langData3 = dx[match(rownames(tipData(tree3)), dx$Name_on_tree_tip2), c('Language','glottolog_id',"Tones","Phonemes","Consonants","Vowels","ASJPVowelRatio",'Lang.Source','InventoryID')]

tipData(tree3) = cbind(tipData(tree3), langData3)

save(tree3, file="../data/processed/GTree_phylo4g_combined_trimmed.Rdat")

######################
##### Plot the trees

plotGeoTree = function(tree4dX, tip.cols=rgb(0,0,0,0.3),add=F){
  
  n1.x = tree4dX@data[as.character(tree4dX@edge[,1]),1]
  n1.y = tree4dX@data[as.character(tree4dX@edge[,1]),2]
  n2.x = tree4dX@data[as.character(tree4dX@edge[,2]),1]
  n2.y = tree4dX@data[as.character(tree4dX@edge[,2]),2]
  
  tips = tree4dX@data[1:length(tree4dX@label),]
  par(mar=c(0,0,0,0))
  if(!add){
    map(ylim=c(-35,10), xlim=c(5,50), interior = F)
  }
  points(tips[,1], tips[,2], pch=16, 
         col= tip.cols)
  arrows(n1.x,n1.y,n2.x,n2.y, length=0.1)
}

plotGeoTree(tree4d)


plot(tree2)

tree2.phylo = as(tree2, "phylo")
plot(tree2.phylo)

dx$tone.colour = cm.colors(max(dx$Tones)+1, alpha = 0.8)[dx$Tones+1]

plotGeoTree(tree2, tip.cols = dx$tone.colour)


h.mean.x = seq(-180,180,length.out=cols)
h.mean.y = seq(-90,90,length.out=rows)
map(ylim=c(-35,10), xlim=c(5,50), interior = F)
map()
colx = gray.colors(10)[as.numeric(cut(h.mean,10))]
points(
  rep(h.mean.x,length(h.mean.y)),
  rep((h.mean.y),each=length(h.mean.x)), 
  col=colx, pch=15,cex=2)

map(ylim=c(-35,10), xlim=c(5,50), interior = F,add=T)
