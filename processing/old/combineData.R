library(ape)

# TODO.  The file phoible-aggregated.tsv actually has all the data we need, no need to use the whole database.

setwd("~/Documents/MPI/ClimateAndLanguage/Grollemund/processing/")

t = read.nexus("../data/trees/BP424_M1P_100_cv2_relaxed.trees")

iso = read.csv("../data/trees/bantu_iso_glottolog_links_28Junr2016.csv", stringsAsFactors = F)

# From original
tree.iso.codes = iso[match(t$tip.label,iso$Name_on_tree_tip2),]$iso6393

tree.glotto.codes = iso[match(t$tip.label,iso$Name_on_tree_tip2),]$glottolog_id
length(unique(tree.iso.codes))
length(unique(tree.glotto.codes))

# From dplace (identical)
#dplace = read.csv("../data/trees/taxa.csv", stringsAsFactors = F)
#tree.glotto.codes.dp = dplace[match(t$tip.label,dplace$taxon),]$glottocode
#length(unique(tree.glotto.codes.dp))




a = read.csv("../data/EBR_2016/ANU_numTones_SpecificHumidity_GlottoFams.csv", stringsAsFactors = F)

#TODO: some iso codes pick out dialects (khg-4)
a$iso = sapply(a$iso, function(X){strsplit(X,"-")[[1]][1]})
#a$iso[a$iso=='xxx'] = NA

a = a[a$iso!="",]
#a = a[a$iso!="xxx",]
a = a[a$iso %in% tree.iso.codes,]

a$GTREE.name = iso[match(a$iso,iso$iso6393),]$Name_on_tree_tip2

#t$tip.label[!tree.iso.codes %in% a$iso]

sum(tree.iso.codes %in% a$iso)

t = drop.tip(t,t$tip.label[!tree.iso.codes %in% a$iso])
t = drop.tip(t, "")

tree.iso.codes = iso[match(t$tip.label,iso$Name_on_tree_tip2),]$iso6393

sample.1 = function(X){
  if(length(X)==1){
    return(X)
  } else{
    sample(X,size=1)
  }
}

# The tree has multiple observations per iso code
# TODO: pick by geographic proximity?

set.seed(123)
random.langs = as.numeric(tapply(1:length(tree.iso.codes),as.character(tree.iso.codes),sample.1))
random.langs.iso = tree.iso.codes[random.langs]
random.langs.name = t$tip.label[random.langs]

toRemove = (1:length(t$tip.label))[-random.langs]
t = drop.tip(t,toRemove)

anu.tree = t

t$tip.label = iso[match(t$tip.label,iso$Name_on_tree_tip2),]$iso6393

# There are multiple sources in the anu data, too
a = a[a$iso %in% t$tip.label,]
set.seed(456)
random.a = tapply(1:nrow(a), a$iso, sample.1)
a = a[random.a,]
rownames(a) = a$iso
# sort in correct order
a = a[t$tip.label,]

a$Language.in.tree = anu.tree$tip.label

write.csv(a,"../data/processed/tones.csv")
write.nexus(t, file = "../data/processed/tree.csv")


##############################
# From phoible

library(reshape)

load("../data/phonology/phoible-by-phoneme.rdata")


# How well do different sources agree?
# (actually, Sources don't overlap much)
sx = levels(final.data$Source)
sx = sx[!sx %in% c("ea",'upsid','saphon','ra')]
mx = matrix(NA, nrow=length(sx),ncol=length(sx))
rownames(mx) = sx
colnames(mx) = sx
for(i in sx){
  for(j in sx){
    ax = unique(final.data[final.data$Source==i,]$LanguageCode)
    bx = unique(final.data[final.data$Source==j,]$LanguageCode)
    mx[i,j] = length(intersect(ax,bx))
  }
}
mx

# ph and gm intersect very little
phoible.ph = final.data[final.data$Source=='gm',]
phoible.gm = final.data[final.data$Source=='aa',]

dx.ph = melt(phoible.ph[,c("LanguageCode","SpecificDialect",'InventoryID','tone')], id=c("LanguageCode","SpecificDialect",'InventoryID'))
dx.ph = cast(dx.ph, InventoryID+LanguageCode~variable, function(X){sum(X=='+')})
tones.ph = tapply(dx.ph$tone, dx.ph$LanguageCode, max)

dx.gm = melt(phoible.gm[,c("LanguageCode","SpecificDialect",'InventoryID','tone')], id=c("LanguageCode","SpecificDialect",'InventoryID'))
dx.gm = cast(dx.gm, InventoryID+LanguageCode~variable, function(X){sum(X=='+')})
tones.gm = tapply(dx.gm$tone, dx.gm$LanguageCode, max)

table(tones.gm,tones.ph[names(tones.gm)])

#####

p = final.data
p = p[!is.na(p$LanguageCode),]
p = p[p$LanguageCode!='0',]
p[is.na(p$Trump),]$Trump = FALSE

# Some sources have no info on tone segments
p = p[!p$Source %in% c("ea",'upsid','saphon','ra'),]

# Some languages have more than one source:
table(tapply(p$Source,p$LanguageCode, function(X){(length(unique(X)))}))

# Phoible includes a 'trump' variable that picks out the best source:
p = p[p$Trump,]

# Still not perfect, so:
# Pick the largest inventory?
bestSource = tapply(p$InventoryID,p$LanguageCode, function(X){
  x = table(X)
  names(x)[which(x==max(x))]
})

keep = apply(p[,c("InventoryID","LanguageCode")], 1, function(X){
  bestSource[as.character(X[2])]==X[1]
})

p = p[keep,]
#check again
table(tapply(p$InventoryID,p$LanguageCode, function(X){(length(unique(X)))}))

tonePhonemes.all = tapply(p$tone,p$LanguageCode,function(X){sum(X=="+")})

# Load ANU data
aP = read.csv("../data/EBR_2016/ANU_numTones_SpecificHumidity_GlottoFams.csv", stringsAsFactors = F)

anu.tones = tapply(aP$Number.of.tones, aP$iso, max)

sum(unique(tree.iso.codes) %in% names(tonePhonemes.all))
sum(unique(tree.iso.codes) %in% names(anu.tones))
sum(unique(tree.iso.codes) %in% setdiff(names(anu.tones),names(tonePhonemes.all)))
sum(unique(tree.iso.codes) %in% setdiff(names(tonePhonemes.all),names(anu.tones)))

intx = intersect(names(anu.tones), names(tonePhonemes.all))
anu.tones2 = anu.tones[intx]
tonePhonemes.all2 = tonePhonemes.all[intx]
familyX = tapply(aP$Family, aP$iso,head,n=1)[intx]

plot(anu.tones2,tonePhonemes.all2)
hist(anu.tones2)
hist(tonePhonemes.all2)

ck = cohen.kappa(cbind(anu.tones2,tonePhonemes.all2))
ck


library(lme4)
library(lattice)

dx = data.frame(anu.tones2,tonePhonemes.all2,familyX)

m0 = glmer(anu.tones2~tonePhonemes.all2 + (1+tonePhonemes.all2|familyX), family=poisson, data=dx)

dotplot(ranef(m0))
summary(m0)

predict(m0, newdata=data.frame(tonePhonemes.all2=0:10,familyX="Niger-Kongo"))

plot(anu.tones2[familyX=="Niger-Kongo"], tonePhonemes.all2[familyX=="Niger-Kongo"], col=rgb(0,0,0,0.1),pch=16)





##########
# Limit dataset to languages in the tree
tP = read.nexus("../data/trees/BP424_M1P_100_cv2_relaxed.trees")

iso = read.csv("../data/trees/bantu_iso_glottolog_links_28Junr2016.csv", stringsAsFactors = F)
tree.iso.codes = iso[match(tP$tip.label,iso$Name_on_tree_tip2),]$iso6393

length(unique(tree.iso.codes))
# Get rid of languages without iso codes
tP = drop.tip(tP, tP$tip.label[tree.iso.codes == ""])

# Recalculate tip iso codes
tree.iso.codes = iso[match(tP$tip.label,iso$Name_on_tree_tip2),]$iso6393

x = intersect(tree.iso.codes, p$LanguageCode)


# Count tone phonemes
tonePhonemes = tapply(p[p$LanguageCode %in% x,]$tone,p[p$LanguageCode %in% x,]$LanguageCode,function(X){sum(X=="+")})

table(tonePhonemes)

#TODO: 1 language only has one tone category?

sum(unique(tree.iso.codes) %in% unique(names(tonePhonemes)))

tP2 = drop.tip(tP, tP$tip.label[!tree.iso.codes %in% names(tonePhonemes)])

phoible.tree = tP2

tP2$tip.label = iso[match(tP2$tip.label,iso$Name_on_tree_tip2),]$iso6393


set.seed(123)
random.langs = as.numeric(tapply(1:length(tP2$tip.label),as.character(tP2$tip.label),sample.1))

toRemove = (1:length(tP2$tip.label))[-random.langs]
phoible.tree = drop.tip(phoible.tree,toRemove)
tP2 = drop.tip(tP2,toRemove)


tonePhonemes = tonePhonemes[names(tonePhonemes) %in% tP2$tip.label]
tonePhonemes = tonePhonemes[tP2$tip.label]

out = data.frame(Language.in.tree = phoible.tree$tip.label ,iso=names(tonePhonemes), Number.of.tones = tonePhonemes)

write.csv(out,"../data/processed/tones_phoible.csv")
write.nexus(tP2, file = "../data/processed/tree_phoible.csv")


#### Combined

tC = read.nexus("../data/trees/BP424_M1P_100_cv2_relaxed.trees")

commonNames = tC$tip.label[tC$tip.label %in% c(phoible.tree$tip.label,anu.tree$tip.label)]

other = tC$tip.label[!tC$tip.label %in% commonNames]

tC = drop.tip(tC,other)

ax = a[!a$Language.in.tree %in% phoible.tree$tip.label,]
#ax = ax[!ax$iso %in% tP2$tip.label,]

out2 = rbind(out, ax[,c("Language.in.tree",'iso',"Number.of.tones")])
out2$Language.in.tree =  as.character(out2$Language.in.tree)
out2$iso =  as.character(out2$iso)
out2$Source = c(rep("Phoible",nrow(out)), rep("ANU",nrow(out2)-nrow(out)))

dup = out2$Language.in.tree[duplicated(out2$iso)]
tC = drop.tip(tC,dup)
out2 = out2[!duplicated(out2$iso),]

tx = out2[match(tC$tip.label, out2$Language.in.tree),]$iso
tx2 = iso[match(tC$tip.label, iso$Name_on_tree_tip2),]$iso6393
tC$tip.label = tx

# Sort by tip label order
rownames(out2) = out2$iso
out2 = out2[tC$tip.label,]

write.csv(out2,"../data/processed/tones_combined.csv")
write.nexus(tC, file = "../data/processed/tree_combined.csv")
