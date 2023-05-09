try(setwd("~/OneDrive - Cardiff University/Research/MPI/ClimateAndLanguage/Grollemund/processing/"))

library(chron)
library(RColorBrewer)
library(lattice)
library(ncdf4)
library(akima)
library(phylobase)
library(splines)
library(spacetime)
library(fields)
library(mgcv)

climateReferenceYear = 1850
# Valdes: Basically it should be thought of as 1850 but there is considerable debate about what is pre-industrial (some people saying it is 1750). For some strange model history reasons, the orbit is never changed so the orbit is for 1950 but the CO2/CH4/N2O values are for 1850. The differences between 1750/1800/1850 are quite small. The time attribute is largely meaningless and is related to how long the model has been run for.

linguisticReferenceYear = 2000


getSimulationData = function(filename, age, extract.var = "q_mm_1_5m", plotData=F){
  # Load simulation data. Filenames indicate age:
  # teiia1 - 0k (i.e. the pre-industrial climate)
  # teiib1 - 1 kaBP 
  # teiic1 - 2 kaBP
  # etc
  # The files containing the specific humidity are called, for example fo 6kaBP - teiig1, 

  d = nc_open(filename)
  #attributes(d$var)$names
  
  qcat = ncatt_get(d, extract.var)
  qcat$date
  qvar = ncvar_get(d, extract.var)
  
  # Check no missing data
  if(min(qvar) < qcat$valid_min){
    warning(paste0("Missing data in ",filename))
  }
  #attributes(d$dim)$names

  nc_lat <- ncvar_get( d, "latitude")
  nc_lon <- ncvar_get( d, "longitude")
  
  #These dimensions match those of the chl-a data which is arranged in a nc_lon(rows) by nc_lat(cols) matrix. This also matches what we saw in the text file.
  
  dimnames(qvar) <- list(lon=nc_lon, lat=nc_lat)
  
  if(plotData){
    plot(c(0,360),c(-90,90),type='n')
    
    colx = heat.colors(20)
    
    for(i in 1:length(nc_lon)){
      for(j in 1:length(nc_lat)){
        points(nc_lon[i],nc_lat[j], col=colx[1+(qvar[i,j]/max(qvar))*(length(colx)-1)],pch=15)
      }
    }
  }
  
  dx = data.frame(
    x = as.numeric(rep(rownames(qvar),ncol(qvar))),
    y = as.numeric(rep(colnames(qvar),each=nrow(qvar))),
    z = as.vector(qvar),
    year = climateReferenceYear - age
  )
  
  return(dx)
  
}

# Process first file
simData = getSimulationData("../data/climate/Valdes/teii/teiiaa.pdclann.nc",0)
# Process other files
simData = rbind(simData, getSimulationData("../data/climate/Valdes/teii/teiiba.pdclann.nc",1000))
simData = rbind(simData, getSimulationData("../data/climate/Valdes/teii/teiica.pdclann.nc",2000))
simData = rbind(simData, getSimulationData("../data/climate/Valdes/teii/teiida.pdclann.nc",3000))
simData = rbind(simData, getSimulationData("../data/climate/Valdes/teii/teiiea.pdclann.nc",4000))
simData = rbind(simData, getSimulationData("../data/climate/Valdes/teii/teiifa.pdclann.nc",5000))
simData = rbind(simData, getSimulationData("../data/climate/Valdes/teii/teiiga.pdclann.nc",6000))
simData = rbind(simData, getSimulationData("../data/climate/Valdes/teii/teiiha.pdclann.nc",7000))

# Demo of spline fitting for a particular location through time
selx = simData$x==18.75 & simData$y==-2.5
fit = interpSpline(z~year, simData[selx,])
pdf("../results/graphs/SplineInterpolation.pdf",width=5,height=4)
par(mar=c(4,5,1,1))
plot( predict( fit, seq(min(simData$year),max(simData$year) )), type = "l",xlab="Time",ylab="",col=2,las=1)
title(ylab="Humidity",line=4)
points(simData[selx,]$z~simData[selx,]$year,pch=19)
dev.off()
# 
# simData = simData[order(simData$year),]
# 
# sp = SpatialPoints(simData[!duplicated(simData[,c("x",'y')]),c("x",'y')],proj4string = CRS("+proj=longlat"))
# time = as.POSIXct("0000-01-01") + (unique(simData$year)*365.25*24*60*60)
# stfdf = STFDF(sp, time, simData)
# 
# empVgm <- variogramST(z~1, stfdf)
# 
# linStAni <- estiStAni(empVgm, c(50000,200000))
# 
# # rescale empVgm and linStAni to km for estimation
# empVgm$dist  <- empVgm$dist/1000
# empVgm$avgDist  <- empVgm$avgDist/1000
# empVgm$spacelag <- empVgm$spacelag/1000
# 
# linStAni <- linStAni/1000
# 
# # separable
# 
# separableModel <- vgmST(NA,"separable",NA,NA)
# fit.StVariogram(empVgm,separableModel)
# 
# separableModel <- vgmST("separable", 
#                         space=vgm(0.9,"Exp", 200, 0.1),
#                         time =vgm(0.9,"Sph", 3.5, 0.1),
#                         sill=120)
# fitSepModel <- fit.StVariogram(empVgm, separableModel, fit.method = 7, 
#                                stAni = linStAni, method = "L-BFGS-B", 
#                                control = list(parscale=c(100,1,10,1,100)),
#                                lower = c(10,0,.1,0,0.1), 
#                                upper = c(2000,1,12,1,200))


# Load phylogeny to get time/location for nodes
# (loads `tree4d`)
load("../data/processed/GTree_phylo4g.Rdat")
tree4d@data$year = linguisticReferenceYear - tree4d@data$Node.Age


selx = simData$x>=0 & simData$x <=52.5 & simData$y >= -35.0 & simData$y <= 17.5

#########
#  GAM  #
#########
# Use a gam to predict humidity given year, lat, long and interactions
k.x = length(unique(simData[selx,]$x))
k.y = length(unique(simData[selx,]$x))
b2 <- gam(z ~ s(year,k=8,fx = T) + te(x, y, year,fx = T) + s(x,y,fx = T) +
            s(x,by=year,fx = T,k=k.x) + s(y,by=year,fx = T,k=k.y),
          data=simData[selx,])
summary(b2)
plot(simData[selx,]$z,predict(b2))

tree4d@data$specH.mean.sim.gam = predict(b2,
      newdata = data.frame(
        x = tree4d@data$Longitude,
        y = tree4d@data$Latitude,
        year = tree4d@data$year))

hist(tree4d@data$specH.mean.sim.gam)

########################
# Linear interpolation #
########################
# Given the geographic sample points,
#  interpolate the humidity for the target year for the sample points,
#  then use this surface to interpolate the target geographic position.

x.vals = as.character(sort(unique(simData[selx,]$x)))
y.vals = as.character(sort(unique(simData[selx,]$y),decreasing =T))
yearVals = as.character(sort(unique(tree4d@data$year)))

m = array(dim=c(length(x.vals),length(y.vals),length(yearVals)),
          dimnames = list(x.vals,y.vals,yearVals))
for(targetYear in yearVals){
  # For each location ...
  for(x in x.vals){
    for(y in y.vals){
      # ... interpolate values and work out z for each target year
      dx = simData[simData$x==as.numeric(x) & simData$y==as.numeric(y),]
      fit = interpSpline(z~year, dx)
      #plot(fit)
      #points(dx$year,dx$z)
      interpolated.z = predict( fit, as.numeric(yearVals))$y
      m[x,y,] = interpolated.z
    }
  }
}

getBilinearInterpolation = function(tree.x,tree.y,tree.targetYear){
  simHForTargetYear = m[,,as.character(tree.targetYear)]
  
  dx = data.frame(
    x = as.numeric(rep(rownames(simHForTargetYear),ncol(simHForTargetYear))),
    y = as.numeric(rep(colnames(simHForTargetYear),each=nrow(simHForTargetYear))),
    z = as.vector(simHForTargetYear))
  
  bl.est = interp.surface(list(x=as.numeric(rownames(simHForTargetYear)),
                               y=as.numeric(colnames(simHForTargetYear)),
                               z=simHForTargetYear),
                          loc = matrix(c(tree.x,tree.y),ncol=2))
  return(bl.est)
}


tree4d@data$specH.mean.sim.bilinear = NA
for(i in 1:nrow(tree4d@data)){
  tree.x = tree4d@data$Longitude[i]
  tree.y = tree4d@data$Latitude[i]
  tree.targetYear = tree4d@data$year[i]
  
  tree4d@data[i,]$specH.mean.sim.bilinear = getBilinearInterpolation(tree.x,tree.y,tree.targetYear)
}

# Bilinear
cor(tree4d@data$specH.mean.sim.bilinear,tree4d@data$specH.mean,use = "complete.obs")
# GAM
cor(tree4d@data$specH.mean.sim.gam,tree4d@data$specH.mean,use = "complete.obs")

save(tree4d,file="../data/processed/GTree_phylo4g_withClimateSimData.RDat")
