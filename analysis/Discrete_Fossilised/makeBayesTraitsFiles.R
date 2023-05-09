library(phylobase)
library(ape)
library(gplots)
#source("https://bioconductor.org/biocLite.R")
#biocLite("ggtree",type="source")
#library(ggtree)


setwd("~/Documents/MPI/ClimateAndLanguage/Grollemund/analysis/Discrete_Fossilised/")

load("../../data/processed/GTree_phylo4g_withClimateSimData.RDat")
tree.chosen = tree4d

# Prune tree
# This breaks if you give it too many things to trim at once, so do it the slow way:
tips.to.remove  = tree4d@label[is.na(tree4d@data[1:length(tree4d@label),]$Tones)]
tree.chosen = prune(tree.chosen,tips.exclude = tips.to.remove, trim.internal=T)

# Tree file
write.nexus(as(tree.chosen,'phylo'),file = "../../data/processed/GTree_phylo4g_combined_trimmed.nex")

humidityVar = "specH.mean.sim.bilinear"
humidityVarCat = paste0(humidityVar,".cat")
tree.chosen@data[,humidityVarCat] = tree.chosen@data[,humidityVar]


# # cut humidity into 3 categories
# cuts = quantile(tree.chosen@data[,humidityVarCat],probs = seq(0,1,length.out=4))
# tree.chosen@data[,humidityVarCat] = cut(tree.chosen@data[,humidityVarCat],
#                                     cuts,labels=c("L","M","H"),include.lowest = T)
# # cut into 0 or 1 tones, 2 tones, 3 or more tones
# toneCuts = c(-1,1,2,8)
# tree.chosen@data$TonesCat = 
#   cut(tree.chosen@data$Tones,
#       toneCuts,labels = c("N","S","C"))

# cut humidity into 2 categories
cuts = c(0,0.0148,1)
tree.chosen@data[,humidityVarCat] = cut(tree.chosen@data[,humidityVarCat],
                                        cuts,labels=c("0","1"),include.lowest = T)
# cut into 0/1/2 vs 3 or more
toneCuts = c(-1,2,10)
tree.chosen@data$TonesCat = 
  cut(tree.chosen@data$Tones,
      toneCuts,labels = c("0","1"))
ProtoBantuValue = "0"

# Data file
write.table(tipData(tree.chosen[,c(humidityVarCat,"TonesCat")]),
            sep=" ",col.names=F, quote=F,
            file="../../data/processed/Discrete_Fossilised/GTree_Data_Tones_3cat.txt")

outCor1 = 
"2
2
LogFile ../../results/BayesTraitsOutput/Discrete_Fossilised/FossilCorA
Stones 100 1000
"

outCor2 = 
"3
2
LogFile ../../results/BayesTraitsOutput/Discrete_Fossilised/FossilCorB
Stones 100 1000
"


ProtoBantuNodeId = rownames(nodeData(tree.chosen))[which(nodeData(tree.chosen)$Node.Age==max(nodeData(tree.chosen)$Node.Age))]


outTag = ""
for(nodeId in rownames(nodeData(tree.chosen))){
  anc = names(descendants(tree.chosen,as.numeric(nodeId),type="tips"))
  tag = paste0("N",nodeId)
  nodeName = paste0("A",nodeId)
  fossilNodeName = paste0("F",nodeId)
  fossilCommand = ""
  # The '-' is added at the end because BayesTraits expects a value for each site
  #fossilCommand  =  paste("Fossil",fossilNodeName,tag,nodeData(tree.chosen)[nodeId,humidityVarCat],"-")
  # Fossilise proto-bantu to have 2 tones
  if(nodeId==ProtoBantuNodeId){
    fossilCommand  =  paste("Fossil",fossilNodeName,tag,nodeData(tree.chosen)[nodeId,humidityVarCat],ProtoBantuValue)
  }
  
  outTag = paste(outTag,
                 paste("AddTag",tag,paste(anc,collapse=" "),collapse=" "),
                 #paste("AddNode",nodeName,tag,collapse=" "),
                 fossilCommand,
                 sep="\n")
}

outCor1 = paste(outCor1,outTag,"Run","\n",sep="\n")
outCor2 = paste(outCor2,outTag,"Run","\n",sep="\n")

cat(outCor1,file="GTree_FossiliseCorrelCommand_A.txt")
cat(outCor2,file="GTree_FossiliseCorrelCommand_B.txt")


runFileA = "/Applications/BayesTraitsV3.0.1-OSX/BayesTraitsV3 ../../data/processed/GTree_phylo4g_combined_trimmed.nex ../../data/processed/Discrete_Fossilised/GTree_Data_Tones_3cat.txt < GTree_FossiliseCorrelCommand_A.txt"
runFileB = "/Applications/BayesTraitsV3.0.1-OSX/BayesTraitsV3 ../../data/processed/GTree_phylo4g_combined_trimmed.nex ../../data/processed/Discrete_Fossilised/GTree_Data_Tones_3cat.txt < GTree_FossiliseCorrelCommand_B.txt"

cat(paste(runFileA,runFileB,"\n",sep="\n"),file="runBayesTraitsFossilisedANC.sh")


###

d = read.table