# Make BayesTraits file for continuous ancestral state reconstruction of tones,
# then also for vowels

library(phylobase)
library(ape)
library(gplots)
#source("https://bioconductor.org/biocLite.R")
#biocLite("ggtree",type="source")
#library(ggtree)

setwd("~/OneDrive - Cardiff University/Research/MPI/ClimateAndLanguage/Grollemund_public/analysis/AncestralStateReconstruction/")

load("../../data/processed/GTree_phylo4g_withClimateSimData_withToneSim.RDat")
tree.chosen = tree4d

# Prune tree
# This breaks if you give it too many things to trim at once, so do it the slow way:
tips.to.remove  = tree4d@label[is.na(tree4d@data[1:length(tree4d@label),]$Tones)]
tree.chosen = phylobase::prune(tree.chosen,tips.exclude = tips.to.remove,trim.internal=TRUE)

# Tree file
write.nexus(as(tree.chosen,'phylo'),file = "../../data/processed/GTree_phylo4g_combined_trimmed.nex")

# cut into 0 or 1 tones, 2 tones, 3 or more tones
#toneCuts = c(-1,1,2,8)
#tree.chosen@data$TonesCat = 
#  cut(tree.chosen@data$Tones,
#      toneCuts,labels = c("N","S","C"))

# Data file
write.table(tipData(tree.chosen[,c("Tones")]),
            sep=" ",col.names=F, quote=F,
            file="../../data/processed/AncestralStateReconstruction/GTree_Data_Tones_continuous.txt")

# Chain lengths (as string to avoid floating point conversion)
burnin = "1000000"
iterations = "2000000"
# First command file to create model (continuous MCMC)
outCor1 = 
"4
2
SaveModels ../../results/BayesTraitsOutput/AncestralStateReconstruction/Anc_Tones_continuous_models.bin
LogFile ../../results/BayesTraitsOutput/AncestralStateReconstruction/Anc_Tones_continuous_make
Burnin BURNINNUM
Iterations ITERATIONNUM
PriorAll uniform 0 12
Seed 9898
"

outCor1 = gsub("BURNINNUM",burnin,outCor1)
outCor1 = gsub("ITERATIONNUM",iterations,outCor1)

# Proto-Bantu had two tones:
# (we've already defined a tag for the top node above)
#https://books.google.co.uk/books?hl=en&lr=&id=M8cHBAAAQBAJ&oi=fnd&pg=PP1&dq=proto-bantu+tone&ots=BNCIXOsOl1&sig=Y4cLBrgPAnW_NHrWWLL-9wLeYD4#v=snippet&q=proto-bantu%20tone&f=false
#https://www.tandfonline.com/doi/pdf/10.1080/00437956.1948.11659343
#https://www.jstor.org/stable/pdf/1264405.pdf?casa_token=2UIv0j03Fk8AAAAA:L8YUTZwX6fqDWb-uLwBah07PjLrBxY1KHfisKrWAvOzlV87P-RM4ly4Roj9V1pMEBzxI5O5fge9Leh8Jwg-UgeYcvdouzLiCWVxipr1vel2olI31gQ
#https://www.persee.fr/docAsPDF/aflin_2033-8732_1967_num_3_1_873.pdf
ProtoBantuNodeId = rownames(nodeData(tree.chosen))[which(nodeData(tree.chosen)$Node.Age==max(nodeData(tree.chosen)$Node.Age))]
ProtoBantuAnc = names(descendants(tree.chosen,as.numeric(ProtoBantuNodeId),type="tips"))
ProtoBantuTag = "PBX"
nodeCommand = paste("AddTag",ProtoBantuTag,paste(ProtoBantuAnc,collapse=" "),collapse=" ")
fossilNodeName = "PBXF"
ProtoBantuNodeValue = "2"
fossilCommand  =  paste("Fossil",fossilNodeName,ProtoBantuTag,ProtoBantuNodeValue)
# Add fossil command to output
outCor1 = paste(outCor1,nodeCommand,sep="")
outCor1 = paste(outCor1,fossilCommand,sep="\n")
outCor1 = paste(outCor1,"\nRun",sep="")

# Write commands
cat(outCor1,file="GTree_Anc_Tones_continuous_makeModels.txt")

# Second model to estimate data
outCor2 =
"4
2
LoadModels ../../results/BayesTraitsOutput/AncestralStateReconstruction/Anc_Tones_continuous_models.bin
LogFile ../../results/BayesTraitsOutput/AncestralStateReconstruction/Anc_Tones_continuous_est
Burnin BURNINNUM
Iterations ITERATIONNUM
PriorAll uniform 0 12
Seed 1212
"

outCor2 = gsub("BURNINNUM",burnin,outCor2)
outCor2 = gsub("ITERATIONNUM",iterations,outCor2)


# Add tags to track prob at each node
outTag = ""
for(nodeId in rownames(nodeData(tree.chosen))){
  anc = names(descendants(tree.chosen,as.numeric(nodeId),type="tips"))
  tag = paste0("N",nodeId)
  nodeName = paste0("A",nodeId)
  outTag = paste(outTag,
    paste("AddTag",tag,paste(anc,collapse=" "),collapse=" "),
    sep="\n")
  if(nodeId==ProtoBantuNodeId){
    fossilCommand  =  paste("Fossil",fossilNodeName,tag,ProtoBantuNodeValue)
    outTag = paste(outTag,fossilCommand,sep="\n")
  } else{
    outTag = paste(outTag,
                   paste("AddMRCA",nodeName,tag,collapse=" "),
                   sep="\n")
  }
}


outCor2 = paste(outCor2,outTag,"Run","\n",sep="\n")


cat(outCor2,file="GTree_Anc_Tones_continuous_estimateAnc.txt")

# There's some problem with Bayes Traits V4 producing continuous estimates
#btVersion = "/Applications/BayesTraitsV4.0.0-OSX/BayesTraitsV4"
btVersion = "/Applications/BayesTraitsV3.0.5-OSX/BayesTraitsV3"

# SH file
runFileA = paste0(btVersion," ../../data/processed/GTree_phylo4g_combined_trimmed.nex ../../data/processed/AncestralStateReconstruction/GTree_Data_Tones_continuous.txt < GTree_Anc_Tones_continuous_makeModels.txt\n",
                  btVersion, " ../../data/processed/GTree_phylo4g_combined_trimmed.nex ../../data/processed/AncestralStateReconstruction/GTree_Data_Tones_continuous.txt < GTree_Anc_Tones_continuous_estimateAnc.txt\n")

cat(paste(runFileA,"\n",sep="\n"),file="p5_runBayesTraits_Anc_Tones_continuous.sh")

################
# Run ancestral states for simulated tones

fhSim = read.csv("../../data/reconstrutions/purelin_df_altsim.csv",stringsAsFactors = F)

tree.chosen@data$Tones.simFH = 8*(fhSim$Tones/max(fhSim$Tones))

tree.chosen@data[is.na(tree.chosen@data$Tones.sim),c("Tones.sim")] = 3
tree.chosen@data[is.na(tree.chosen@data$Tones.sim2),c("Tones.sim2")] = 3
tree.chosen@data[is.na(tree.chosen@data$Tones.sim10),c("Tones.sim10")] = 3

write.table(tipData(tree.chosen[,c("Tones.sim")]),
            sep=" ",col.names=F, quote=F,
            file="../../data/processed/AncestralStateReconstruction/GTree_Data_Tones_SIM1_continuous.txt")
write.table(tipData(tree.chosen[,c("Tones.sim2")]),
            sep=" ",col.names=F, quote=F,
            file="../../data/processed/AncestralStateReconstruction/GTree_Data_Tones_SIM2_continuous.txt")
write.table(tipData(tree.chosen[,c("Tones.sim10")]),
            sep=" ",col.names=F, quote=F,
            file="../../data/processed/AncestralStateReconstruction/GTree_Data_Tones_SIM10_continuous.txt")
write.table(round(tipData(tree.chosen[,c("Tones.simFH")]),3),
            sep=" ",col.names=F, quote=F,
            file="../../data/processed/AncestralStateReconstruction/GTree_Data_Tones_SIMFH_continuous.txt")

outCor1SIM = gsub("Anc_Tones_continuous","Anc_Tones_SIM1_continuous",outCor1)
cat(outCor1SIM,file="GTree_Anc_Tones_SIM1_continuous_makeModels.txt")
outCor2SIM = gsub("Anc_Tones_continuous","Anc_Tones_SIM1_continuous",outCor2)
cat(outCor2SIM,file="GTree_Anc_Tones_SIM1_continuous_estimateAnc.txt")

outCor1SIM = gsub("Anc_Tones_continuous","Anc_Tones_SIM2_continuous",outCor1)
cat(outCor1SIM,file="GTree_Anc_Tones_SIM2_continuous_makeModels.txt")
outCor2SIM = gsub("Anc_Tones_continuous","Anc_Tones_SIM2_continuous",outCor2)
cat(outCor2SIM,file="GTree_Anc_Tones_SIM2_continuous_estimateAnc.txt")

outCor1SIM = gsub("Anc_Tones_continuous","Anc_Tones_SIM10_continuous",outCor1)
cat(outCor1SIM,file="GTree_Anc_Tones_SIM10_continuous_makeModels.txt")
outCor2SIM = gsub("Anc_Tones_continuous","Anc_Tones_SIM10_continuous",outCor2)
cat(outCor2SIM,file="GTree_Anc_Tones_SIM10_continuous_estimateAnc.txt")

outCor1SIM = gsub("Anc_Tones_continuous","Anc_Tones_SIMFH_continuous",outCor1)
cat(outCor1SIM,file="GTree_Anc_Tones_SIMFH_continuous_makeModels.txt")
outCor2SIM = gsub("Anc_Tones_continuous","Anc_Tones_SIMFH_continuous",outCor2)
cat(outCor2SIM,file="GTree_Anc_Tones_SIMFH_continuous_estimateAnc.txt")

runFileASIM = paste0(btVersion," ../../data/processed/GTree_phylo4g_combined_trimmed.nex ../../data/processed/AncestralStateReconstruction/GTree_Data_Tones_SIM1_continuous.txt < GTree_Anc_Tones_SIM1_continuous_makeModels.txt\n",
                  btVersion, " ../../data/processed/GTree_phylo4g_combined_trimmed.nex ../../data/processed/AncestralStateReconstruction/GTree_Data_Tones_SIM1_continuous.txt < GTree_Anc_Tones_SIM1_continuous_estimateAnc.txt\n")

cat(paste(runFileASIM,"\n",sep="\n"),file="p5_runBayesTraits_Anc_Tones_SIM1_continuous.sh")

runFileASIM2 = gsub("SIM1","SIM2",runFileASIM)
cat(paste(runFileASIM2,"\n",sep="\n"),file="p5_runBayesTraits_Anc_Tones_SIM2_continuous.sh")

runFileASIM10 = gsub("SIM1","SIM10",runFileASIM)
cat(paste(runFileASIM10,"\n",sep="\n"),file="p5_runBayesTraits_Anc_Tones_SIM10_continuous.sh")

runFileASIMFH = gsub("SIM1","SIMFH",runFileASIM)
cat(paste(runFileASIMFH,"\n",sep="\n"),file="p5_runBayesTraits_Anc_Tones_SIMFH_continuous.sh")

##########3

# Switch back to Bayestraits v4
btVersion = "/Applications/BayesTraitsV4.0.0-OSX/BayesTraitsV4"

#########################
# Same process for Vowels (vowel ratio)

load("../../data/processed/GTree_phylo4g_withClimateSimData.RDat")
tree.chosen = tree4d
tree.chosen@data$VowelRatio = tree.chosen@data$Vowels/(tree.chosen@data$Vowels + tree.chosen@data$Consonants)

tips.to.remove  = tree.chosen@label[is.na(tree.chosen@data[1:length(tree.chosen@label),]$VowelRatio)]
tree.chosen = prune(tree.chosen,tips.exclude = tips.to.remove, trim.internal=T)

write.nexus(as(tree.chosen,'phylo'),file = "../../data/processed/GTree_phylo4g_combined_trimmed_VOWELS.nex")

write.table(tipData(tree.chosen[,c("VowelRatio")]),
            sep=" ",col.names=F, quote=F,
            file="../../data/processed/AncestralStateReconstruction/GTree_Data_VOWELS_continuous.txt")

write.table(tipData(tree.chosen[,c("Vowels")]),
            sep=" ",col.names=F, quote=F,
            file="../../data/processed/AncestralStateReconstruction/GTree_Data_VOWELSRAW_continuous.txt")

write.table(tipData(tree.chosen[,c("Consonants")]),
            sep=" ",col.names=F, quote=F,
            file="../../data/processed/AncestralStateReconstruction/GTree_Data_CONSRAW_continuous.txt")


# Chain lengths (as string to avoid floating point conversion)
burnin = "1000000"
iterations = "2000000"
# First command file to create model (continuous MCMC)
outCor1 = 
  "4
2
SaveModels ../../results/BayesTraitsOutput/AncestralStateReconstruction/Anc_VOWELS_continuous_models.bin
LogFile ../../results/BayesTraitsOutput/AncestralStateReconstruction/Anc_VOWELS_continuous_VOWELS_make
Burnin BURNINNUM
Iterations ITERATIONNUM
PriorAll gamma 5.8 19.7
Seed 1565
"

outCor1 = gsub("BURNINNUM",burnin,outCor1)
outCor1 = gsub("ITERATIONNUM",iterations,outCor1)

# (No fossilisation, because we don't have a reconstruction)
outCor1 = paste(outCor1,"\nRun",sep="")

# Write commands
cat(outCor1,file="GTree_Anc_VOWELS_continuous_makeModels.txt")

# Second model to estimate data
outCor2 =
  "4
2
LoadModels ../../results/BayesTraitsOutput/AncestralStateReconstruction/Anc_VOWELS_continuous_models.bin
LogFile ../../results/BayesTraitsOutput/AncestralStateReconstruction/Anc_VOWELS_continuous_est
Burnin BURNINNUM
Iterations ITERATIONNUM
PriorAll gamma 5.8 19.7
Seed 1873
"

outCor2 = gsub("BURNINNUM",burnin,outCor2)
outCor2 = gsub("ITERATIONNUM",iterations,outCor2)

# Add tags to track prob at each node
outTag = ""
for(nodeId in rownames(nodeData(tree.chosen))){
  anc = names(descendants(tree.chosen,as.numeric(nodeId),type="tips"))
  tag = paste0("N",nodeId)
  nodeName = paste0("A",nodeId)
  outTag = paste(outTag,
                 paste("AddTag",tag,paste(anc,collapse=" "),collapse=" "),
                 sep="\n")
  outTag = paste(outTag,
                   paste("AddMRCA",nodeName,tag,collapse=" "),
                   sep="\n")
  
}


outCor2 = paste(outCor2,outTag,"Run","\n",sep="\n")


cat(outCor2,file="GTree_Anc_VOWELS_continuous_estimateAnc.txt")

# SH file
runFileA = "
/Applications/BayesTraitsV4.0.0-OSX/BayesTraitsV4 ../../data/processed/GTree_phylo4g_combined_trimmed_VOWELS.nex ../../data/processed/AncestralStateReconstruction/GTree_Data_VOWELS_continuous.txt < GTree_Anc_VOWELS_continuous_makeModels.txt
/Applications/BayesTraitsV4.0.0-OSX/BayesTraitsV4 ../../data/processed/GTree_phylo4g_combined_trimmed_VOWELS.nex ../../data/processed/AncestralStateReconstruction/GTree_Data_VOWELS_continuous.txt < GTree_Anc_VOWELS_continuous_estimateAnc.txt"

cat(paste(runFileA,"\n",sep="\n"),file="p6_runBayesTraits_Anc_VOWELS_continuous.sh")


#####################
# Raw Vowels

# First command file to create model (continuous MCMC)
outCor1 = 
  "4
2
SaveModels ../../results/BayesTraitsOutput/AncestralStateReconstruction/Anc_VOWELSRAW_continuous_models.bin
LogFile ../../results/BayesTraitsOutput/AncestralStateReconstruction/Anc_VOWELSRAW_continuous_make
Burnin BURNINNUM
Iterations ITERATIONNUM
PriorAll gamma 3.3 0.3
Seed 8383
"

outCor1 = gsub("BURNINNUM",burnin,outCor1)
outCor1 = gsub("ITERATIONNUM",iterations,outCor1)

# (No fossilisation, because we don't have a reconstruction)
outCor1 = paste(outCor1,"\nRun",sep="")

# Write commands
cat(outCor1,file="GTree_Anc_VOWELSRAW_continuous_makeModels.txt")

# Second model to estimate data
outCor2 =
  "4
2
LoadModels ../../results/BayesTraitsOutput/AncestralStateReconstruction/Anc_VOWELSRAW_continuous_models.bin
LogFile ../../results/BayesTraitsOutput/AncestralStateReconstruction/Anc_VOWELSRAW_continuous_est
Burnin BURNINNUM
Iterations ITERATIONNUM
PriorAll gamma 3.3 0.3
Seed 373
"

outCor2 = gsub("BURNINNUM",burnin,outCor2)
outCor2 = gsub("ITERATIONNUM",iterations,outCor2)

# Add tags to track prob at each node
outTag = ""
for(nodeId in rownames(nodeData(tree.chosen))){
  anc = names(descendants(tree.chosen,as.numeric(nodeId),type="tips"))
  tag = paste0("N",nodeId)
  nodeName = paste0("A",nodeId)
  outTag = paste(outTag,
                 paste("AddTag",tag,paste(anc,collapse=" "),collapse=" "),
                 sep="\n")
  outTag = paste(outTag,
                 paste("AddMRCA",nodeName,tag,collapse=" "),
                 sep="\n")
  
}


outCor2 = paste(outCor2,outTag,"Run","\n",sep="\n")


cat(outCor2,file="GTree_Anc_VOWELSRAW_continuous_estimateAnc.txt")

# SH file (same tree as vowels)
runFileA = "
/Applications/BayesTraitsV4.0.0-OSX/BayesTraitsV4 ../../data/processed/GTree_phylo4g_combined_trimmed_VOWELS.nex ../../data/processed/AncestralStateReconstruction/GTree_Data_VOWELSRAW_continuous.txt < GTree_Anc_VOWELSRAW_continuous_makeModels.txt
/Applications/BayesTraitsV4.0.0-OSX/BayesTraitsV4 ../../data/processed/GTree_phylo4g_combined_trimmed_VOWELS.nex ../../data/processed/AncestralStateReconstruction/GTree_Data_VOWELSRAW_continuous.txt < GTree_Anc_VOWELSRAW_continuous_estimateAnc.txt"

cat(paste(runFileA,"\n",sep="\n"),file="p6_runBayesTraits_Anc_VOWELSRAW_continuous.sh")


#####################
# Raw Consonants

# First command file to create model (continuous MCMC)
outCor1 = 
  "4
2
SaveModels ../../results/BayesTraitsOutput/AncestralStateReconstruction/Anc_CONSRAW_continuous_models.bin
LogFile ../../results/BayesTraitsOutput/AncestralStateReconstruction/Anc_CONSRAW_continuous_make
Burnin BURNINNUM
Iterations ITERATIONNUM
PriorAll gamma 6.85 0.286
Seed 284
"

outCor1 = gsub("BURNINNUM",burnin,outCor1)
outCor1 = gsub("ITERATIONNUM",iterations,outCor1)

# (No fossilisation, because we don't have a reconstruction)
outCor1 = paste(outCor1,"\nRun",sep="")

# Write commands
cat(outCor1,file="GTree_Anc_CONSRAW_continuous_makeModels.txt")

# Second model to estimate data
outCor2 =
  "4
2
LoadModels ../../results/BayesTraitsOutput/AncestralStateReconstruction/Anc_CONSRAW_continuous_models.bin
LogFile ../../results/BayesTraitsOutput/AncestralStateReconstruction/Anc_CONSRAW_continuous_est
Burnin BURNINNUM
Iterations ITERATIONNUM
PriorAll gamma 6.85 0.286
Seed 1093
"

outCor2 = gsub("BURNINNUM",burnin,outCor2)
outCor2 = gsub("ITERATIONNUM",iterations,outCor2)

# Add tags to track prob at each node
outTag = ""
for(nodeId in rownames(nodeData(tree.chosen))){
  anc = names(descendants(tree.chosen,as.numeric(nodeId),type="tips"))
  tag = paste0("N",nodeId)
  nodeName = paste0("A",nodeId)
  outTag = paste(outTag,
                 paste("AddTag",tag,paste(anc,collapse=" "),collapse=" "),
                 sep="\n")
  outTag = paste(outTag,
                 paste("AddMRCA",nodeName,tag,collapse=" "),
                 sep="\n")
  
}


outCor2 = paste(outCor2,outTag,"Run","\n",sep="\n")


cat(outCor2,file="GTree_Anc_CONSRAW_continuous_estimateAnc.txt")

# SH file
runFileA = "
/Applications/BayesTraitsV4.0.0-OSX/BayesTraitsV4 ../../data/processed/GTree_phylo4g_combined_trimmed_VOWELS.nex ../../data/processed/AncestralStateReconstruction/GTree_Data_CONSRAW_continuous.txt < GTree_Anc_CONSRAW_continuous_makeModels.txt
/Applications/BayesTraitsV4.0.0-OSX/BayesTraitsV4 ../../data/processed/GTree_phylo4g_combined_trimmed_CONS.nex ../../data/processed/AncestralStateReconstruction/GTree_Data_CONSRAW_continuous.txt < GTree_Anc_CONSRAW_continuous_estimateAnc.txt"

cat(paste(runFileA,"\n",sep="\n"),file="p6_runBayesTraits_Anc_CONSRAW_continuous.sh")


#######################
# ASJP ratio

load("../../data/processed/GTree_phylo4g_withClimateSimData.RDat")
tree.chosen = tree4d

tips.to.remove  = tree.chosen@label[is.na(tree.chosen@data[1:length(tree.chosen@label),]$ASJPVowelRatio)]
tree.chosen = prune(tree.chosen,tips.exclude = tips.to.remove, trim.internal=T)

write.nexus(as(tree.chosen,'phylo'),file = "../../data/processed/GTree_phylo4g_combined_trimmed_ASJP.nex")

write.table(tipData(tree.chosen[,c("ASJPVowelRatio")]),
            sep=" ",col.names=F, quote=F,
            file="../../data/processed/AncestralStateReconstruction/GTree_Data_ASJP_continuous.txt")

# First command file to create model (continuous MCMC)
outCor1 = 
  "4
2
SaveModels ../../results/BayesTraitsOutput/AncestralStateReconstruction/Anc_ASJP_continuous_models.bin
LogFile ../../results/BayesTraitsOutput/AncestralStateReconstruction/Anc_ASJP_continuous_make
Burnin BURNINNUM
Iterations ITERATIONNUM
PriorAll normal 0.4659 0.0556
Seed 1000
"

outCor1 = gsub("BURNINNUM",burnin,outCor1)
outCor1 = gsub("ITERATIONNUM",iterations,outCor1)

# (No fossilisation, because we don't have a reconstruction)
outCor1 = paste(outCor1,"\nRun",sep="")

# Write commands
cat(outCor1,file="GTree_Anc_ASJP_continuous_makeModels.txt")

# Second model to estimate data
outCor2 =
  "4
2
LoadModels ../../results/BayesTraitsOutput/AncestralStateReconstruction/Anc_ASJP_continuous_models.bin
LogFile ../../results/BayesTraitsOutput/AncestralStateReconstruction/Anc_ASJP_continuous_est
Burnin BURNINNUM
Iterations ITERATIONNUM
PriorAll normal 0.4659 0.0556
Seed 3783
"

outCor2 = gsub("BURNINNUM",burnin,outCor2)
outCor2 = gsub("ITERATIONNUM",iterations,outCor2)

# Add tags to track prob at each node
outTag = ""
for(nodeId in rownames(nodeData(tree.chosen))){
  anc = names(descendants(tree.chosen,as.numeric(nodeId),type="tips"))
  tag = paste0("N",nodeId)
  nodeName = paste0("A",nodeId)
  outTag = paste(outTag,
                 paste("AddTag",tag,paste(anc,collapse=" "),collapse=" "),
                 sep="\n")
  outTag = paste(outTag,
                 paste("AddMRCA",nodeName,tag,collapse=" "),
                 sep="\n")
  
}


outCor2 = paste(outCor2,outTag,"Run","\n",sep="\n")


cat(outCor2,file="GTree_Anc_ASJP_continuous_estimateAnc.txt")

# SH file
runFileA = "
/Applications/BayesTraitsV4.0.0-OSX/BayesTraitsV4 ../../data/processed/GTree_phylo4g_combined_trimmed_ASJP.nex ../../data/processed/AncestralStateReconstruction/GTree_Data_ASJP_continuous.txt < GTree_Anc_ASJP_continuous_makeModels.txt
/Applications/BayesTraitsV4.0.0-OSX/BayesTraitsV4 ../../data/processed/GTree_phylo4g_combined_trimmed_ASJP.nex ../../data/processed/AncestralStateReconstruction/GTree_Data_ASJP_continuous.txt < GTree_Anc_ASJP_continuous_estimateAnc.txt"

cat(paste(runFileA,"\n",sep="\n"),file="p6_runBayesTraits_Anc_ASJP_continuous.sh")
