# Make BayesTraits file for discrete evolution of tones
# Using contious model A
# continuous variables cannot be fossilised for indepdent contrasts

library(phylobase)
library(ape)
library(gplots)
#source("https://bioconductor.org/biocLite.R")
#biocLite("ggtree",type="source")
#library(ggtree)

setwd("~/Documents/MPI/ClimateAndLanguage/Grollemund/analysis/ContinuousA/")

load("../../data/processed/GTree_phylo4g_withClimateSimData.RDat")
tree.chosen = tree4d
humidityVar = "specH.mean.sim.bilinear"

# Prune tree
# This breaks if you give it too many things to trim at once, so do it the slow way:
tips.to.remove  = tree4d@label[is.na(tree4d@data[1:length(tree4d@label),]$Tones)]
tree.chosen = prune(tree.chosen,tips.exclude = tips.to.remove, trim.internal=T)

# Tree file
write.nexus(as(tree.chosen,'phylo'),file = "../../data/processed/GTree_phylo4g_combined_trimmed.nex")

# Data file
write.table(tipData(tree.chosen[,c("Tones",humidityVar)]),
            sep=" ",col.names=F, quote=F,
            file="../../data/processed/ContinuousA/GTree_Data_Tones_SpecHSimBi.txt")

# Command file (multistate MCMC)
outCor1 = 
"4
2
LogFile ../../results/BayesTraitsOutput/ContinuousA/ContA_Tones_SpecHSimBi
Burnin 50000
Iterations 1050000
Stones 100 1000
"

outCor2 = 
"4
2
TestCorrel
Burnin 50000
Iterations 1050000
LogFile ../../results/BayesTraitsOutput/ContinuousA/ContA_Tones_SpecHSimBi_TestCorrel
Stones 100 1000"

# Add tags to track prob at each node
#  and fossilise humidity variable
#  and proto-bantu has two tones
#https://books.google.co.uk/books?hl=en&lr=&id=M8cHBAAAQBAJ&oi=fnd&pg=PP1&dq=proto-bantu+tone&ots=BNCIXOsOl1&sig=Y4cLBrgPAnW_NHrWWLL-9wLeYD4#v=snippet&q=proto-bantu%20tone&f=false
#https://www.tandfonline.com/doi/pdf/10.1080/00437956.1948.11659343
#https://www.jstor.org/stable/pdf/1264405.pdf?casa_token=2UIv0j03Fk8AAAAA:L8YUTZwX6fqDWb-uLwBah07PjLrBxY1KHfisKrWAvOzlV87P-RM4ly4Roj9V1pMEBzxI5O5fge9Leh8Jwg-UgeYcvdouzLiCWVxipr1vel2olI31gQ
#https://www.persee.fr/docAsPDF/aflin_2033-8732_1967_num_3_1_873.pdf

ProtoBantuNodeId = rownames(nodeData(tree.chosen))[which(nodeData(tree.chosen)$Node.Age==max(nodeData(tree.chosen)$Node.Age))]

outTag = ""
for(nodeId in rownames(nodeData(tree.chosen))){
  anc = names(descendants(tree.chosen,as.numeric(nodeId),type="tips"))
  tag = paste0("N",nodeId)
  nodeName = paste0("A",nodeId)
  fossilNodeName = paste0("F",nodeId)
  fossilCommand = ""
  # The '-' is added at the end because BayesTraits expects a value for each site
  fossilCommand  =  paste("Fossil",fossilNodeName,tag,"-",nodeData(tree.chosen)[nodeId,humidityVar])
  # Fossilise proto-bantu to have 2 tones
  if(nodeId==ProtoBantuNodeId){
    fossilCommand  =  paste("Fossil",fossilNodeName,tag,"2",nodeData(tree.chosen)[nodeId,humidityVar])
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

runFileA = "/Applications/BayesTraitsV3.0.1-OSX/BayesTraitsV3 ../../data/processed/GTree_phylo4g_combined_trimmed.nex ../../data/processed/ContinuousA/GTree_Data_Tones_SpecHSimBi.txt < GTree_FossiliseCorrelCommand_A.txt"
runFileB = "/Applications/BayesTraitsV3.0.1-OSX/BayesTraitsV3 ../../data/processed/GTree_phylo4g_combined_trimmed.nex ../../data/processed/ContinuousA/GTree_Data_Tones_SpecHSimBi.txt < GTree_FossiliseCorrelCommand_B.txt"

cat(paste(runFileA,
          runFileB,
          "\n",sep="\n"),file="runBayesTraitsFossilisedContinuousA.sh")

## Version without fossilising:

ProtoBantuFossilCommand = paste("Fossil F",ProtoBantuNodeId," N",ProtoBantuNodeId," 2 -",sep="")

removeNodeFossils = function(X){
  X = strsplit(X,"\n")[[1]]
  X = X[!grepl("Fossil",X)]
  X = paste(X,collapse="\n")
  X = gsub("\nRun",paste("\n",ProtoBantuFossilCommand,"\nRun",sep="\n"),X)
}

outCor1NoF = removeNodeFossils(outCor1)
outCor1NoF = gsub("LogFile ../../results/BayesTraitsOutput/ContinuousA/ContA_Tones_SpecHSimBi",
                  "LogFile ../../results/BayesTraitsOutput/ContinuousA/ContA_Tones_SpecHSimBi_NoNodeFossils",
                  outCor1NoF)
outCor2NoF = removeNodeFossils(outCor2)
outCor2NoF = gsub("LogFile ../../results/BayesTraitsOutput/ContinuousA/ContA_Tones_SpecHSimBi_TestCorrel",
                  "LogFile ../../results/BayesTraitsOutput/ContinuousA/ContA_Tones_SpecHSimBi_NoNodeFossils_TestCorrel",
                  outCor2NoF)

cat(outCor1NoF,file="GTree_FossiliseCorrelCommand_A_NoNodeFossils.txt")
cat(outCor2NoF,file="GTree_FossiliseCorrelCommand_B_NoNodeFossils.txt")

runFileANoF = "/Applications/BayesTraitsV3.0.1-OSX/BayesTraitsV3 ../../data/processed/GTree_phylo4g_combined_trimmed.nex ../../data/processed/ContinuousA/GTree_Data_Tones_SpecHSimBi.txt < GTree_FossiliseCorrelCommand_A_NoNodeFossils.txt"
runFileBNoF = "/Applications/BayesTraitsV3.0.1-OSX/BayesTraitsV3 ../../data/processed/GTree_phylo4g_combined_trimmed.nex ../../data/processed/ContinuousA/GTree_Data_Tones_SpecHSimBi.txt < GTree_FossiliseCorrelCommand_B_NoNodeFossils.txt"

cat(paste(runFileANoF,
          runFileBNoF,
          "\n",sep="\n"),file="runBayesTraitsFossilisedContinuousA_NoNodeFossils.sh")
