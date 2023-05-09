# TODO: is this using the most up to date files?

library(ape)
library(phytools)
library(caper)
library(phytools)
library(nlme)
library(geiger)

setwd("~/Documents/MPI/ClimateAndLanguage/Grollemund/analysis")


tree = read.nexus("../data/processed/tree_combined.csv")
a = read.csv("../data/processed/tones_combined.csv", row.names = 1)
a = a[!is.na(a$Tones),]
a = a[!is.na(a$glottolog_id),]

ASE <- ace(x = a$Tones, phy = tree, method = "ML")

plot.phylo(tree, label.offset = 0.01)
tiplabels(a$Tones, width = 0.5)
nodelabels(round(ASE$ace, 0), width = 0.5)


tones = a$Tones
names(tones) = a$iso
contMap(tree,tones, outline = F)


########
library(geiger)
mod_l_1 <- fitContinuous(phy = tree, dat = tones, model = "BM")
mod_l_est <- fitContinuous(phy = tree, dat = tones, model = "lambda")

mod_l_est$opt$lambda

LR <- 2 * (logLik(mod_l_est) - logLik(mod_l_1))
LR 
P <- pchisq(LR, df = mod_l_est$opt$k - mod_l_1$opt$k, lower.tail = FALSE)
P

#########

bm<-corBrownian(phy=tree)
bm.gls<-gls(Tones~specH.mean,
            correlation=bm,data=a)
summary(bm.gls)




########

# Different rate models
makeQ.linear = function(n){
  X = matrix(0, nrow = n, ncol=n)
  for(i in 1:n){
    if(i<n){
      X[i,i+1] = i
    }
    if(i>1){
      X[i,i-1] = i-1
    }
  }
  return(X)
}

makeQ.linear.ARD = function(n){
  X = matrix(0, nrow = n, ncol=n)
  varNum = 1
  for(i in 1:n){
    if(i>1){
      X[i,i-1] = varNum
      varNum = varNum+ 1
    }
    if(i<n){
      X[i,i+1] = varNum
      varNum = varNum+ 1
    }
  }
  return(X)
}

categories = unique(a$Tones)

linear.Q = makeQ.linear(length(categories))
# single-rate model
iQ<-matrix(c(0,1,2,0,3,0,0,2,4,0,0,1,0,4,3,0),4,4,byrow=TRUE)

# differet rates for complex -> simple in humid and dry conditions
dQ =matrix(c(0,1,2,0,3,0,0,5,4,0,0,1,0,4,3,0),4,4,byrow=TRUE) 

a = a[a$iso6393 %in% tree$tip.label,]
a = a[!duplicated(a$iso6393),]
tones = as.factor(a$Tones)
names(tones) = a$iso

fit.linear <- fitDiscrete(tree,tones,model = linear.Q)

plot(fit.linear, show.zeros=F)


fit.linear.kappa <- fitDiscrete(tree,tones,model = linear.Q, transform='kappa')

plot(fit.linear.kappa, show.zeros=F)

# 3 categories

tones.3 = cut(a$Tones,c(-1,1,2,3,7))
names(tones.3) = a$iso

fit.linear.3 = fitDiscrete(tree, tones.3, model = makeQ.linear.ARD(length(unique(tones.3))))
plot(fit.linear.3, show.zeros=F)

fit.3 = fitDiscrete(tree,tones.3,model="ARD")
plot(fit.3, show.zeros=F)

# ARD

fit.linear.ARD = fitDiscrete(tree, tones, model = makeQ.linear.ARD(length(unique(tones))))

plot(fit.linear.ARD)

