library(phylobase)
library(phytools)
library(ggtree)

try(setwd("~/OneDrive - Cardiff University/Research/MPI/ClimateAndLanguage/Grollemund_public/processing/"))

load("../data/processed/GTree_phylo4g_withClimateSimData.RDat")
tree.chosen = tree4d

# Prune tree
# This breaks if you give it too many things to trim at once, so do it the slow way:
tips.to.remove  = tree4d@label[is.na(tree4d@data[1:length(tree4d@label),]$Tones)]
tree.chosen = prune(tree.chosen,tips.exclude = tips.to.remove, trim.internal=T)

simTone = function (tree, anc = 0, sdX = 0.001, ngen = 1000, 
                    effectSize = 1, newVarName="sim", ...) 
{
#  if (!inherits(tree, "phylo")) 
#    stop("tree should be an object of class \"phylo\".")
  return.tree <- TRUE
  trO <- reorder(tree)
  tr <- as(trO,"phylo")
  H <- nodeHeights(tr)
  tr$edge.length <- round(tr$edge.length/max(H) * ngen)
  tr <- di2multi(tr)
  h <- nodeHeights(tr)[, 1]
  
  bmSim <- function(start, n, humidityDiff){
    #xh = (parentHumidity - 0.0148)/0.006 # (between -1 and 1 roughly)
    xh = humidityDiff / 0.009# (between -1 and 1 roughly)
    ret = cumsum(c(start, rnorm(n, mean = xh*effectSize, sd = sdX)))
    ret[ret<0] = 0
    ret[ret>12] = 12
    return(ret)
  }
  X <- T <- list()
  N <- length(tr$tip)
  trO@data[,paste0("Tones.",newVarName)] = NA
  for (i in 1:nrow(tr$edge)) {
    if (tr$edge[i, 1] == (N + 1)) 
      X[[i]] <- bmSim(anc, tr$edge.length[i], 0)
    else {
      parent <- match(tr$edge[i, 1], tr$edge[, 2])
      parentHumidity = trO@data[parent,]$specH.mean.sim.gam
      childHumidity = trO@data[i,]$specH.mean.sim.gam
      simulatedValues <- bmSim(X[[parent]][length(X[[parent]])], 
                      tr$edge.length[i],
                      childHumidity - parentHumidity)
      X[[i]] <- simulatedValues
      finalValue = simulatedValues[length(simulatedValues)]
      trO@data[i,paste0("Tones.",newVarName)] = finalValue
    }
    T[[i]] <- h[i] + 0:tr$edge.length[i]
  }
#  if (type == "BM") 
#    cols <- phytools:::bm(X, T, ...)
#  else if (type == "threshold") 
#    cols <- th(X, T, ...)
#  if (return.tree) {
#    if (type == "BM") 
#      tr$X <- X
#  }
  Y <- matrix(NA, nrow(tr$edge), 2)
  Y[, 1] <- sapply(X, function(x) x[1])
  Y[, 2] <- sapply(X, function(x) x[length(x)])
  x <- Y[tr$edge[, 2] %in% 1:N, 2]
  names(x) <- tr$tip.label[tr$edge[tr$edge[, 2] %in% 1:N, 2]]
  a <- c(anc, Y[tr$edge[, 2] %in% (N + 2:tr$Nnode), 2])
  names(a) <- c(N + 1, tr$edge[tr$edge[, 2] %in% (N + 2:tr$Nnode), 
                               2])
  xx <- c(x[tr$tip], a[as.character(N + 1:tr$Nnode)])
  if (!return.tree) 
    return(xx)
  else return(list(x = xx, tree = tr, treeO = trO))
}

set.seed(45678)
st1 = simTone(tree.chosen,anc=2,effectSize=1/10)
st2 = simTone(tree.chosen,anc=2,effectSize=2/10)
st10 = simTone(tree.chosen,anc=2,effectSize=1/2)

tree4d = st1$treeO

tree4d@data[,"Tones.sim2"] = st2$treeO@data$Tones.sim
tree4d@data[,"Tones.sim10"] = st10$treeO@data$Tones.sim

#pdf("~/Downloads/tt.pdf", height=100,width=7)
  ggtree(tree4d) + geom_label(aes(fill=Tones.sim,label=""))
#dev.off()
ggtree(tree4d) + geom_label(aes(fill=Tones.sim2,label=""))
ggtree(tree4d) + geom_label(aes(fill=Tones.sim10,label=""))

ggtree(tree4d) + geom_label(aes(fill=specH.mean.sim.gam,label=""))


save(tree4d,file="../data/processed/GTree_phylo4g_withClimateSimData_withToneSim.RDat")
save(tree4d,file="../../Grollemund/data/GTree_phylo4g_withClimateSimData_withToneSim.RDat")

