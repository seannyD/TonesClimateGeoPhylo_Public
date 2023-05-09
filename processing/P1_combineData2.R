library(dplyr)
try(setwd("~/OneDrive - Cardiff University/Research/MPI/ClimateAndLanguage/Grollemund/processing/"))

# Data from the Grollemund et al. tree
#tree.data = read.csv("../data/trees/bantu_iso_glottolog_links_28Junr2016.csv", stringsAsFactors = F)
# Updated taxa in Feb 2020 that maps more glottolog ids
tree.data = read.csv("../data/trees/taxa.csv",stringsAsFactors = F)
tree.data[tree.data$glottocode=="",]$glottocode = NA

# Load Phoible data
#library(dplyr)
phoible = read.csv("../data/phonology/phoible.csv",stringsAsFactors = F,encoding = "utf-8", fileEncoding = "utf-8")

# While we're here, let's work out the prior for vowel ratio
numVowels = tapply(phoible$SegmentClass ,
                   paste(phoible$Glottocode,
                         phoible$InventoryID,sep="#"),
                   function(X){sum(X=="vowel")})
numConsonants = tapply(phoible$SegmentClass ,
                   paste(phoible$Glottocode,
                         phoible$InventoryID,sep="#"),
                   function(X){sum(X=="consonant")})
vowelRatio = numVowels/(numVowels+numConsonants)
# A gamma distribution seems sensible.
library(fitdistrplus)
fitdist(as.vector(vowelRatio), "beta")
fitdist(as.vector(vowelRatio), "gamma")
hist(vowelRatio,freq = F,breaks = 40)
p = seq(0,1, length=100)
lines(p, dbeta(p, 4.1, 9.6))
lines(p, dgamma(p, 5.8, 19.7),col=2)

fitdist(as.vector(numVowels),"gamma")
hist(numVowels,freq = F,breaks = 20)
p = seq(0,50, length=100)
lines(p, dgamma(p, 3.3, 0.3),col=2)

fitdist(as.vector(numConsonants),"gamma")
hist(numConsonants,freq = F,breaks = 20)
p = seq(0,120, length=100)
lines(p, dgamma(p, 6.85, 0.286),col=2)

phoible = phoible[phoible$Glottocode %in% tree.data$glottocode[!is.na(tree.data$glottocode)],]

# Filter sources so that each language has one source, with the following priority
sourcePriority = c("ph","aa","upsid","spa", # Main sources
                   "saphon", "ra","er","ea")# these others don't appear in our data
for(gl in unique(phoible$Glottocode)){
  sources = unique(phoible[phoible$Glottocode==gl,]$Source)
  if(length(sources)>1){
    mx = match(sources,sourcePriority)
    trumpSource = sources[which(mx==min(mx))]
    phoible = phoible[(phoible$Glottocode != gl) | (phoible$Glottocode==gl & phoible$Source==trumpSource),]
  }
}

# There are still multiple sources for some languages
numSources = tapply(phoible$InventoryID,phoible$Glottocode,function(X){length(unique(X))})
multiSources = numSources[numSources>1]
# Check some examples manually:
dx = phoible[phoible$Glottocode=="mboc1235",]
table(dx$SegmentClass,dx$InventoryID)

# For tone, the decision about which to keep only affects the estimate for one language. By coin flip, I chose InventoryID 1258.
# Variation in vowel/consonant ratio is very small.
phoible[phoible$InventoryID!=1444,]

# Pick the first inventory listed for the rest
iidsToKeep = sapply(names(multiSources),function(X){phoible[phoible$Glottocode==X,]$InventoryID[1]})
for(gx in names(multiSources)){
  phoible = phoible[phoible$Glottocode!=gx | (phoible$Glottocode==gx & phoible$InventoryID==iidsToKeep[gx]),]
}

# Old aggregated data
#p2langs = read.delim("../data/phonology/phoible-aggregated.tsv", sep='\t', stringsAsFactors = F, fileEncoding = 'utf-8')

pget = function(X){
  tapply(phoible[,X],phoible$Glottocode,head,n=1)
}

p2 = data.frame(
  "Language" = pget("LanguageName"),
  "InventoryID"= pget("InventoryID"),
  'iso'= pget("ISO6393"),
  "Latitude"= NA,
  "Longitude"= NA,
  'Tones'=  tapply(phoible$SegmentClass,phoible$Glottocode,function(X){sum(X=="tone")}),
  'Phonemes'= tapply(phoible$SegmentClass,phoible$Glottocode,length),
  'Consonants' = tapply(phoible$SegmentClass,phoible$Glottocode,function(X){sum(X=="consonant")}),
  "Vowels" = tapply(phoible$SegmentClass,phoible$Glottocode,function(X){sum(X=="vowel")}),
  "Source" = paste0("Phoible (",pget("Source"),")"),
  "Glottocode" = pget("Glottocode")
)

# Take languages from phoible, and top up from ANU

a = read.csv("../data/phonology/ANU_numTones_SpecificHumidity_GlottoFams_utf8.csv", stringsAsFactors = F, fileEncoding = 'utf-8')
# Glottolog to match by ISO
g = read.csv("../data/glottolog/languages-and-dialects-geo.csv", stringsAsFactors = F)
g = g[g$isocodes!="",]

a2 = a[,c("Language","OID_","iso","Latitude",'Longitude','Number.of.tones')]
a2$Phonemes = NA
a2$Consonants = NA
a2$Vowels = NA
a2$Source = "ANU"
a2$Glottocode = g[match(a2$iso,g$isocodes),]$glottocode
names(a2) = c("Language","InventoryID","iso","Latitude",'Longitude','Tones','Phonemes','Consonants',"Vowels","Source","Glottocode")

combined = rbind(p2[!is.na(p2$Glottocode),],a2[!a2$Glottocode %in% p2$Glottocode & (!is.na(a2$Glottocode)),])

# Match language data onto tree data
tree.data[,c("Language","Tones","Phonemes","Consonants","Vowels",'Lang.Source','InventoryID')] = combined[match(tree.data$glottocode, combined$Glottocode), c("Language","Tones","Phonemes","Consonants","Vowels",'Source','InventoryID')]


# The number of tones was found for 193 languages. For the remaining 232 languages, we attempted to find phonological references listed in Glottolog. We found additional data on the number of tones in 19 languages (see SI) to bring the total to 212.

############
# Random sample to code by hand
tree.data[is.na(tree.data$Tones),c("taxon",'glottocode')]
newG = unique(p2langs$Glottocode)
tree.data[is.na(tree.data$Tones),]$glottcode %in% newG

set.seed(2389)
x = sample(tree.data[is.na(tree.data$Tones),]$glottocode)
x

# Load data found by hand:
extra = read.delim("../data/phonology/extraToneData.tab",stringsAsFactors = F,sep='\t')

tree.data[match(extra$glotto,tree.data$glottocode),]$Tones = extra$tones
tree.data[match(extra$glotto,tree.data$glottocode),]$Consonants = extra$consonants
tree.data[match(extra$glotto,tree.data$glottocode),]$Vowels = extra$vowels
tree.data[match(extra$glotto,tree.data$glottocode),]$Lang.Source = "Extra"

# Rename fields to keep compatibility with other scripts
names(tree.data)[names(tree.data)=="glottocode"] = "glottolog_id"
names(tree.data)[names(tree.data)=="taxon"] = "Name_on_tree_tip2"

#############
## ASJP data

asjp = read.csv("../data/asjp/forms.csv",stringsAsFactors = F,encoding = "UTF-8",fileEncoding = "UTF-8")
asjpLangs = read.csv('../data/asjp/languages.csv',stringsAsFactors = F)

asjp$Glottocode = asjpLangs[match(asjp$Language_ID,asjpLangs$ID),]$Glottocode
asjp$LangName = asjpLangs[match(asjp$Language_ID,asjpLangs$ID),]$Name
treeGlotto = tree.data$glottolog_id
treeGlotto = treeGlotto[!is.na(treeGlotto)]
#asjp = asjp[asjp$Glottocode %in% treeGlotto,]

# Language names are unique, while there may be multiple sources 
#  per glottocode
vowels = c("a"  , "ã"  , "ɐ" ,  "ɐ̃", "e",   "ẽ" ,  "ə" ,  "ə̃", "i",   "ĩ",
"o" ,  "õ" , "u"  , "ũ")
vowelRatiosByLangNames = 
  tapply(asjp$Segments, asjp$LangName, function(X){
  segs = unlist(strsplit(X," "))
  vowels = sum(segs %in% vowels)
  return(vowels/length(segs))
})
vrGlotto = asjp[match(names(vowelRatiosByLangNames),asjp$LangName),]$Glottocode
meanVowelRatiosByGlottoCodes = tapply(vowelRatiosByLangNames,vrGlotto,mean)

# Look at vowel ratio distirbution: looks very normal.
hist(meanVowelRatiosByGlottoCodes,probability = T,breaks=40)
mean(meanVowelRatiosByGlottoCodes)
sd(meanVowelRatiosByGlottoCodes)
lines(seq(0,1,length.out=100), dnorm(seq(0,1,length.out=100), mean = mean(meanVowelRatiosByGlottoCodes), 
             sd = sd(meanVowelRatiosByGlottoCodes)))


tree.data$ASJPVowelRatio = meanVowelRatiosByGlottoCodes[tree.data$glottolog_id]

write.csv(tree.data,"../data/processed/tones_combined.csv",row.names = F)
