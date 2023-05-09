setwd("~/OneDrive - Cardiff University/Research/MPI/ClimateAndLanguage/Grollemund/analysis/IndependentContrasts/")

source("../bayestraitr_helper.R")
library(diagram)
library(plotrix)
library(phylobase)
#BiocManager::install("ggtree")
library("ggtree")
library(ggimage)
library(phangorn)


d = bt_read.log("../../results/BayesTraitsOutput/ContinuousA/ContA_Tones_SpecHSimBi.Log.txt")

plot(d$`R Trait 1 2`)
hist(d$`R Trait 1 2`,xlim=c(0,0.1))
mean(d$`R Trait 1 2`)

g = ggplot(d,aes(y=`R Trait 1 2`,x="A")) +
  geom_violin() +
  coord_cartesian(ylim=c(-0.05,0.1)) +
  geom_hline(yintercept = 0) +
  ylab("r")
g

pdf("../../results/BayesTraitsOutput/ContinuousA/Continuous_ModelA_R.pdf",width=2,height=3)
g
dev.off()

getLogMarginalLikelihoodFromStonesFile = function(f){
  a = readLines(f)
  return(as.numeric(strsplit(tail(a,1),"\t")[[1]][2]))
}

logMarginalLikelihood_a = getLogMarginalLikelihoodFromStonesFile("../../results/BayesTraitsOutput/ContinuousA/ContA_Tones_SpecHSimBi.Stones.txt")
logMarginalLikelihood_b = getLogMarginalLikelihoodFromStonesFile("../../results/BayesTraitsOutput/ContinuousA/ContA_Tones_SpecHSimBi_TestCorrel.Stones.txt")

LogBF = 2*(logMarginalLikelihood_a - logMarginalLikelihood_b)
LogBF


### 
# No fossilised Nodes:
dNoF = bt_read.log("../../results/BayesTraitsOutput/ContinuousA/ContA_Tones_SpecHSimBi_NoNodeFossils.Log.txt")

hist(dNoF$`R Trait 1 2`,xlim=c(0,0.1))
mean(dNoF$`R Trait 1 2`)

logMarginalLikelihood_aNoF = getLogMarginalLikelihoodFromStonesFile("../../results/BayesTraitsOutput/ContinuousA/ContA_Tones_SpecHSimBi_NoNodeFossils.Stones.txt")
logMarginalLikelihood_bNoF = getLogMarginalLikelihoodFromStonesFile("../../results/BayesTraitsOutput/ContinuousA/ContA_Tones_SpecHSimBi_NoNodeFossils_TestCorrel.Stones.txt")

LogBFNoF = 2*(logMarginalLikelihood_aNoF - logMarginalLikelihood_bNoF)
LogBFNoF
