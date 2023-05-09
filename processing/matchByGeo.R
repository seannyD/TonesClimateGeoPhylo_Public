library(ape)
library(fields)

setwd("~/Documents/MPI/ClimateAndLanguage/Grollemund/processing/")

t = read.nexus("../data/trees/BP424_M1P_100_cv2_relaxed.trees")

iso = read.csv("../data/trees/bantu_iso_glottolog_links_28Junr2016.csv", stringsAsFactors = F)

#isob = read.csv("../data/trees/taxa.csv", stringsAsFactors = F)
#isob = isob[match(t$tip.label, isob$taxon),]

# From original
tree.iso.codes = iso[match(t$tip.label,iso$Name_on_tree_tip2),]$iso6393

tree.glotto.codes = iso[match(t$tip.label,iso$Name_on_tree_tip2),]$glottolog_id
length(unique(tree.iso.codes))
length(unique(tree.glotto.codes))


tdata = read.delim("../data/trees/data.txt", stringsAsFactors = F)
tdata = tdata[tdata$Node.Group!='Int',]

tdata$iso = iso[match(tdata$Node.Name, iso$Name_on_tree_tip2),]$iso6393

tdata = tdata[is.na(tdata$iso) | tdata$iso=='',]

g = read.csv("../data/glottolog-languoid/languoid.csv", stringsAsFactors = F)
g = g[g$latitude<40 & g$longitude>-35 & g$longitude <59,]
plot(g$longitude, g$latitude)
g = g[!is.na(g$id),]
rownames(g) = g$id

dist = rdist.earth(cbind(g$longitude, g$latitude), cbind(tdata$Node.Long, tdata$Node.Lat), miles = F)


matches = which(dist< 5, arr.ind = T)

g[matches[1,1],]
tdata[matches[1,2],]
# A15A_Mbuu, mboc1235, mbo
# Already in data


g[matches[2,1],]
tdata[matches[2,2],]
# "A54_Tibea", tibe1274, ngy
# Already in data


g[matches[3,1],]
tdata[matches[3,2],]
# A62B_Mmala, elip1238, ekm


matches2 = which(dist>=5 & dist<10, arr.ind = T)

g[matches2[4,1],]
tdata[matches2[4,2],]
