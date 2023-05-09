setwd("~/Documents/MPI/ClimateAndLanguage/Grollemund/analysis/IndependentContrasts/")

source("../bayestraitr_helper.R")
library(diagram)
library(plotrix)
library(phylobase)
#BiocManager::install("ggtree")
library("ggtree")
library(ggimage)
library(phangorn)

getLogMarginalLikelihoodFromStonesFile = function(f){
  a = readLines(f)
  return(as.numeric(strsplit(tail(a,1),"\t")[[1]][2]))
}

logMarginalLikelihood_a = getLogMarginalLikelihoodFromStonesFile("../../results/BayesTraitsOutput/IndependentContrasts/IC_Tones_SpecHSimBi.Stones.txt")
logMarginalLikelihood_b = getLogMarginalLikelihoodFromStonesFile("../../results/BayesTraitsOutput/IndependentContrasts/IC_Tones_SpecHSimBi_TestCorrel.Log.txt")

LogBF = 2*(logMarginalLikelihood_a - logMarginalLikelihood_b)
LogBF
