setwd("~/Documents/MPI/ClimateAndLanguage/Grollemund/processing/")
a = read.csv("~/Desktop/Stuff/Everett/ANU_numTones_SpecificHumidity_GlottoFams_utf8.csv", stringsAsFactors = F, fileEncoding = 'utf-8')


a2 = a[,c("Language","OID_","iso","Latitude",'Longitude')]
a2$Source = "ANU"
a2$Glottocode = NA
names(a2) = c("Language","InventoryID","iso","Latitude",'Longitude',"Source","Glottocode")

p2langs = read.delim("../data/phonology/phoible-aggregated.tsv", sep='\t', stringsAsFactors = F, fileEncoding = 'utf-8')

p2 = p2langs[,c("LanguageName","InventoryID",'LanguageCode',"Latitude","Longitude")]
p2$Source = "Phoible"
p2$Glottocode = p2langs$Glottocode
names(p2) = c("Language","InventoryID","iso","Latitude",'Longitude',"Source","Glottocode")

out = rbind(p2,a2)

out = out[!is.na(out$Latitude),]
out = out[out$Latitude!="NULL",]
out = out[out$Latitude!="",]

write.csv(out, "../data/climate/Valdes/PresentLanguages_LatLong.csv", fileEncoding = 'utf-8')



tdata = read.delim("../data/trees/data.txt", stringsAsFactors = F)
