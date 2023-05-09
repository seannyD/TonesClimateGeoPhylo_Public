# Make BayesTraits file for discrete evolution of tones

library(phylobase)
library(ape)
library(gplots)
#source("https://bioconductor.org/biocLite.R")
#biocLite("ggtree",type="source")
#library(ggtree)

setwd("~/Documents/MPI/ClimateAndLanguage/Grollemund/analysis/AncestralStateReconstruction/")

load("../../data/processed/GTree_phylo4g_withClimateSimData.RDat")
tree.chosen = tree4d

# Prune tree
# This breaks if you give it too many things to trim at once, so do it the slow way:
tips.to.remove  = tree4d@label[is.na(tree4d@data[1:length(tree4d@label),]$Tones)]
tree.chosen = prune(tree.chosen,tips.exclude = tips.to.remove, trim.internal=T)

# Tree file
write.nexus(as(tree.chosen,'phylo'),file = "../../data/processed/GTree_phylo4g_combined_trimmed.nex")

# cut into 0 or 1 tones, 2 tones, 3 or more tones
toneCuts = c(-1,1,2,8)
tree.chosen@data$TonesCat = 
  cut(tree.chosen@data$Tones,
      toneCuts,labels = c("N","S","C"))

# Data file
write.table(tipData(tree.chosen[,c("TonesCat")]),
            sep=" ",col.names=F, quote=F,
            file="../../data/processed/AncestralStateReconstruction/GTree_Data_Tones_3cat.txt")

# Command file (multistate MCMC)
outCor1 = 
"1
2
LogFile ../../results/BayesTraitsOutput/AncestralStateReconstruction/Anc_Tones_3cat
Burnin 50000
Iterations 1050000
Stones 100 1000
"


# Proto-Bantu had two tones:
# (we've already defined a tag for the top node above)
#https://books.google.co.uk/books?hl=en&lr=&id=M8cHBAAAQBAJ&oi=fnd&pg=PP1&dq=proto-bantu+tone&ots=BNCIXOsOl1&sig=Y4cLBrgPAnW_NHrWWLL-9wLeYD4#v=snippet&q=proto-bantu%20tone&f=false
#https://www.tandfonline.com/doi/pdf/10.1080/00437956.1948.11659343
#https://www.jstor.org/stable/pdf/1264405.pdf?casa_token=2UIv0j03Fk8AAAAA:L8YUTZwX6fqDWb-uLwBah07PjLrBxY1KHfisKrWAvOzlV87P-RM4ly4Roj9V1pMEBzxI5O5fge9Leh8Jwg-UgeYcvdouzLiCWVxipr1vel2olI31gQ
#https://www.persee.fr/docAsPDF/aflin_2033-8732_1967_num_3_1_873.pdf
ProtoBantuNodeId = rownames(nodeData(tree.chosen))[which(nodeData(tree.chosen)$Node.Age==max(nodeData(tree.chosen)$Node.Age))]


# Add tags to track prob at each node
outTag = ""
for(nodeId in rownames(nodeData(tree.chosen))){
  anc = names(descendants(tree.chosen,as.numeric(nodeId),type="tips"))
  tag = paste0("N",nodeId)
  nodeName = paste0("A",nodeId)
  outTag = paste(outTag,
    paste("AddTag",tag,paste(anc,collapse=" "),collapse=" "),
    paste("AddNode",nodeName,tag,collapse=" "),
    sep="\n")
  if(nodeId==ProtoBantuNodeId){
    fossilNodeName = paste0("F",nodeId)
    fossilCommand  =  paste("Fossil",fossilNodeName,tag,"S")
    outTag = paste(outTag,fossilCommand,sep="\n")
  }
}


outCor1 = paste(outCor1,outTag,"Run","\n",sep="\n")


cat(outCor1,file="GTree_Anc_Tones_3cat.txt")

# SH file
runFileA = "/Applications/BayesTraitsV3.0.1-OSX/BayesTraitsV3 ../../data/processed/GTree_phylo4g_combined_trimmed.nex ../../data/processed/AncestralStateReconstruction/GTree_Data_Tones_3cat.txt < GTree_Anc_Tones_3cat.txt"

cat(paste(runFileA,"\n",sep="\n"),file="runBayesTraits_Anc_Tones_3cat.sh")
