rm(list=ls())
require(sf)
require(raster)
require(MASS)
require(dismo)
require(rJava)
require(virtualspecies)
require(dplyr)

#To deal with java issues
.jinit()

wd <- "~/Desktop/GVTP/"
#wd <- "/home1/alsimons/GVTP"

setwd(wd)
set.seed(1)

#Generate project study area if it has not already been generated.
if(length(list.files(wd,"GVTPStudyArea.shp"))<1){
  #Read in US ecoregions and CA counties.
  Ecoregions <- sf::st_read("S_USA.EcomapSections.shp")
  Counties <- sf::st_read("CA_Counties_TIGER2016.shp")
  
  #Select Santa Barbra and San Luis Obispo counties and merge them.
  CountySubset <- Counties[Counties$NAME %in% c("Santa Barbara","San Luis Obispo"),]
  CountySubset <- sf::st_union(CountySubset)
  #Reproject to EPSG:4326
  CountySubset <- sf::st_transform(CountySubset,4326)
  #Clip to remove Channel islands
  mainland <- sf::st_bbox(c(xmin = -121.4, xmax = -119.4, ymax = 35.8, ymin = 34.3), crs = st_crs(4326)) %>% st_as_sfc()
  CountySubset <- sf::st_intersection(CountySubset,mainland)
  
  #Selection ecoregions 261A and 261B and merge them.
  EcoregionSubset <- Ecoregions[Ecoregions$MAP_UNIT_S %in% c("261A","261B"),]
  EcoregionSubset <- sf::st_union(EcoregionSubset)
  #Reproject to EPSG:4326
  EcoregionSubset <- sf::st_transform(EcoregionSubset,4326)
  
  #Find the overlap between counties and ecoregions, this is the study area.  Reproject to EPSG:2229
  StudyArea <- sf::st_intersection(CountySubset,EcoregionSubset)
  StudyArea <- sf::st_transform(StudyArea,2229)
  sf::st_write(StudyArea,"GVTPStudyArea.shp")
}

StudyArea <- sf::st_read("GVTPStudyArea.shp")

#Determine the environmental variables, with a low degree of collinearity, to use
#in this analysis.
if(length(list.files(pattern="EnvFiltered.txt",full.names=TRUE)) != 1){
  #Read in environmental layers for project.
  EnvironmentalLayers <- raster::stack(Sys.glob("GVTP*.tif"))
  #names(EnvironmentalLayers) <- gsub("GVTP","",names(EnvironmentalLayers))
  #Remove environmental layers with a high degree of multicollinearity.
  #Using r=0.7 as the cutoff.
  env.filtered <- removeCollinearity(EnvironmentalLayers,nb.points = 10000,sample.points = T,select.variables = T,multicollinearity.cutoff = 0.7,plot=T)
  env.filtered <- as.data.frame(env.filtered)
  colnames(env.filtered) <- c("Variable")
  #Save filtered list of environmental variables.
  write.table(env.filtered,"EnvFiltered.txt",quote=FALSE,sep="\t",row.names = FALSE) 
}
#Get environmental rasters which have low collinearity.
env.filtered <- read.table("EnvFiltered.txt", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote="", encoding = "UTF-8")
#Filter rasters to only keep those with low collinearity.
env.names <- gsub(".tif","",Sys.glob("GVTP*.tif"))
env.names <- paste(env.names[env.names %in% env.filtered$Variable],".tif",sep="")
#Stack environmental layers
EnvironmentalLayers <- stack(c(env.names))

#Randomly generate specific gaviota tarplant observation points.
#Use the number of observations recorded per patch, and randomly generate those points within them.
GVTPObservationAreas <- sf::st_read("GaviotaTarplant_CNDDB2021.shp")
GVTPObservations <- sf::st_sfc()
for(patch in unique(GVTPObservationAreas$MAPNDX)){
  GVTPObservationArea <- GVTPObservationAreas[GVTPObservationAreas$MAPNDX==patch,]
  PatchObservations <- sf::st_sample(GVTPObservationArea,GVTPObservationArea$OCCNUMBER)
  GVTPObservations <- c(GVTPObservations,PatchObservations)
}
GVTPObservations <- as_Spatial(GVTPObservations)

#Create a bias file for illustrating the observation density bias in species observations.
#This is for downstream use in Maxent modeling for generating background points.
if(length(list.files(wd,"Bias.tif"))<1){
  occurrence <- as.data.frame(GVTPObservations@coords)
  colnames(occurrence) <- c("longitude","latitude")
  dens.ras <- raster(kde2d(occurrence$longitude, occurrence$latitude, n=c(nrow(EnvironmentalLayers[[1]]), ncol(EnvironmentalLayers[[1]])),lims=c(sf::st_bbox(StudyArea)[1],sf::st_bbox(StudyArea)[2]),range(sf::st_bbox(StudyArea)[3],sf::st_bbox(StudyArea)[4])))
  dens.ras <- raster(kde2d(occurrence$longitude, occurrence$latitude, n=c(nrow(EnvironmentalLayers[[1]]), ncol(EnvironmentalLayers[[1]])),lims=c(range(sf::st_bbox(StudyArea)[1],sf::st_bbox(StudyArea)[3]),range(sf::st_bbox(StudyArea)[2],sf::st_bbox(StudyArea)[4]))))
  writeRaster(dens.ras, "tmp.tif",overwrite=T)
  #Clip and align bias raster to the LA area and projection.
  #Get gdalwarp path in Unix by typing 'which gdalwarp' and copying the path.
  gdalwarpPath <- "/opt/homebrew/bin/gdalwarp"
  gdalwarpCommand <- "-srcnodata -3.39999995214436425e+38 -dstnodata -9999 -of GTiff -co COMPRESS=LZW -s_srs EPSG:2229 -cutline GVTPStudyArea.shp -crop_to_cutline -co BIGTIFF=YES -tap -tr 30 30 tmp.tif Bias.tif"
  system2(gdalwarpPath,gdalwarpCommand)
}
dens.ras <- raster("Bias.tif")

if(length(list.files(wd,"EnvironmentalObservations.txt"))<1){
  #Extract environmental values at the observation locations.
  EnvironmentalObservations <- cbind(GVTPObservations@coords,raster::extract(EnvironmentalLayers,GVTPObservations))
  #Remove missing data
  EnvironmentalObservations <- as.data.frame(EnvironmentalObservations[complete.cases(EnvironmentalObservations),])
  #Save output.
  write.table(EnvironmentalObservations,"EnvironmentalObservations.txt",quote=FALSE,sep="\t",row.names = FALSE)
}
EnvironmentalObservations <- read.table("EnvironmentalObservations.txt", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote="", encoding = "UTF-8")

if(length(list.files(wd,"EnvironmentalBackground.txt"))<1){
  #Create a background point set.
  #Create background points using bias background layer.
  bg <- as.data.frame(xyFromCell(dens.ras, sample(which(!is.na(values(dens.ras))), 10*nrow(EnvironmentalObservations), prob=values(dens.ras)[!is.na(values(dens.ras))])))
  colnames(bg) <- c("lon","lat")
  #Extract environmental values at the observation locations.
  EnvironmentalBackground <- cbind(bg,raster::extract(EnvironmentalLayers,bg))
  #Remove missing data
  EnvironmentalBackground <- as.data.frame(EnvironmentalBackground[complete.cases(EnvironmentalBackground),])
  #Save output.
  write.table(EnvironmentalBackground,"EnvironmentalBackground.txt",quote=FALSE,sep="\t",row.names = FALSE)
}
EnvironmentalBackground <- read.table("EnvironmentalBackground.txt", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote="", encoding = "UTF-8")

#Define categorical environmental variables.
env.factors <- c("GVTPGeology","GVTPLandform","GVTPLandCover","GVTPSoilTexture")

#Generate and evaluate MaxEnt species distribution models using randomly subsampled presence/background data.
if(length(list.files(wd,"GVTPModelEvaluationSummary.txt"))<1){
  XMEvaluations <- data.frame()
  for(i in 1:100){
    #Randomly subsample observation and background points.
    sampleNum <- floor(0.8*nrow(EnvironmentalObservations))
    ObservationsSubsample <- EnvironmentalObservations[sample(nrow(EnvironmentalObservations),sampleNum),]
    ObservationsSubsample$pa <- 1
    BackgroundSubsample <- EnvironmentalBackground[sample(nrow(EnvironmentalBackground),10*sampleNum),]
    BackgroundSubsample$pa <- 0
    #Create a merged data input and set certain columns to factor.
    ModelInput <- rbind(ObservationsSubsample,BackgroundSubsample)
    ModelInput[env.factors] <- lapply(ModelInput[env.factors], factor)
    #Run Maxent model
    xm <- dismo::maxent(x=ModelInput[,grep("GVTP",colnames(ModelInput))],p=ModelInput$pa,factors=env.factors)
    #Get relative importance of environmental layers.
    XMImportance <- as.data.frame(plot(xm))
    XMImportance <- as.data.frame(t(XMImportance))
    #Evaluate maxent model.
    exm <- suppressWarnings(evaluate(p=ModelInput[ModelInput$pa==1,grep("GVTP",colnames(ModelInput))],a=ModelInput[ModelInput$pa==0,grep("GVTP",colnames(ModelInput))],model=xm))
    tmp <- data.frame(matrix(nrow=1,ncol=4))
    colnames(tmp) <- c("AUC","TSS","r","p")
    tmp$AUC <- exm@auc
    #Calculate the true skill statistic
    a <- mean(exm@TPR,na.rm=T)
    b <- mean(exm@FPR,na.rm=T)
    c <- mean(exm@FNR,na.rm=T)
    d <- mean(exm@TNR,na.rm=T)
    H <- a/(a+c)
    F <- b/(b+d)
    tmp$TSS <- H+F-1
    tmp$r <- exm@cor
    tmp$p <- exm@pcor
    tmp <- cbind(tmp,XMImportance)
    XMEvaluations <- rbind(XMEvaluations,tmp)
    print(paste("auc:",tmp$AUC,"TSS:",tmp$TSS,"cor:",tmp$r,"p:",tmp$p)) 
  }
  #Summarize model evaluations.
  XMEvaluationsSD <- XMEvaluations %>% dplyr::summarise_if(is.numeric, sd)
  rownames(XMEvaluationsSD) <- c("Standard Deviation (% Importance)")
  XMEvaluationsMean <- XMEvaluations %>% dplyr::summarise_if(is.numeric, mean)
  rownames(XMEvaluationsMean) <- c("Mean (% Importance")
  XMEvaluationsSummary <- rbind(XMEvaluationsMean,XMEvaluationsSD)
  colnames(XMEvaluationsSummary) <- gsub(x = colnames(XMEvaluationsSummary), pattern = "\\GVTP", replacement = "")
  XMEvaluationsSummary <- as.data.frame(t(XMEvaluationsSummary))
  write.table(XMEvaluationsSummary,"GVTPModelEvaluationSummary.txt",quote=FALSE,sep="\t",row.names = TRUE)
}
XMEvaluationsSummary <- read.table("GVTPModelEvaluationSummary.txt", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote="", encoding = "UTF-8")

#Plot percent relative importance
require(ggplot2)
tmp <- XMEvaluationsSummary
colnames(tmp) <- c("Mean","SD")
tmp$Variable <- as.factor(rownames(tmp))
tmp <- tmp[!(tmp$Variable %in% c("AUC","TSS","r","p")),]
tmp <- dplyr::arrange(tmp,-Mean)
ggplot(tmp,aes(x=Variable,y=Mean))+
  scale_x_discrete(limits = tmp$Variable)+
  geom_bar(position=position_dodge(),stat="identity",color="gray55")+
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD, color="red"),width=.2,position=position_dodge(.9))+
  ggtitle("Percent relative importance of environmental variables\n(Mean and standard deviation, 100 MaxEnt models)")+
  xlab("Environmental variables")+ylab("Percent relative importances")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),strip.background=element_blank(),strip.text.x=element_blank(),legend.position = "none")

#Run Maxent model
EnvironmentalObservations$pa <- 1
EnvironmentalBackground$pa <- 0
ModelInput <- rbind(EnvironmentalObservations,EnvironmentalBackground)
xm <- dismo::maxent(x=ModelInput[,grep("GVTP",colnames(ModelInput))],p=ModelInput$pa,factors=env.factors)
#Evaluate maxent model.
exm <- suppressWarnings(dismo::evaluate(EnvironmentalObservations[,c("lon","lat")],EnvironmentalBackground[,c("lon","lat")],xm,EnvironmentalLayers))
#Evaluate probability threshold of species detection
bc.threshold <- dismo::threshold(x = exm, stat = "spec_sens")
#Generate prediction maps for the Gaviota tarplant, both presence/absence and by probability of observations.
r <- dismo::predict(object=xm,x=EnvironmentalLayers,progress='text')
#r <- dismo::predict(object=xm,x=test,progress='text')
#Convert prediction probability raster to a presence/absence prediction.
rPA <- r > bc.threshold
writeRaster(r,"GVTP_ObservationProbability.tif",overwrite=T)
writeRaster(rPA,"GVTP_PresenceAbsence.tif",overwrite=T)

#Map Gaviota tarplant species distribution model.
require(ggmap)
require(maptools)
require(maps)
require(mapdata)
require(webshot)
require(leaflet)
require(viridis)
require(ggrastr)
MapDF <- as.data.frame(raster::raster("GVTP_PresenceAbsenceEPSG4326.tif"),xy=T, na.rm=T)
MapDF <- as.data.frame(raster::rasterToPoints(raster::raster("GVTP_PresenceAbsenceEPSG4326.tif")))
colnames(MapDF) <-c ("Presence","longitude","latitude")
MapCoordinates <- MapDF[seq(1,nrow(MapDF),floor(nrow(MapDF)/100000)),]
ggmap(get_stamenmap(bbox = c(left=min(MapCoordinates$longitude)-buffer,bottom=min(MapCoordinates$latitude)-buffer,right=max(MapCoordinates$longitude)+buffer,top=max(MapCoordinates$latitude)+buffer),zoom=8,maptype="terrain-background"))+
  geom_point(data=MapCoordinates,aes(x=longitude,y=latitude,color=Presence))+
  labs(title="Gaviota tarplant suitable habitat",x="longitude",y="latitude")+
  theme(legend.position="none")+
  coord_sf(crs = st_crs(4326))

#Map of Gaviota tarplant observations areas.
GVTPObservationAreas <- sf::st_read("GaviotaTarplant_CNDDB2021.shp")
bbox <- sf::st_bbox(sf::st_bbox(sf::st_transform(GVTPObservationAreas,crs=st_crs(4326))))
bbox <- c(left=bbox[1],bottom=bbox[2],right=bbox[3],top=bbox[4])
buffer <- 0.1
ggmap(get_stamenmap(bbox = c(left=as.numeric(bbox[1])-buffer,bottom=as.numeric(bbox[2])-buffer,right=as.numeric(bbox[3])+buffer,top=as.numeric(bbox[4])+buffer),zoom=12,maptype="terrain-background"))+
  geom_sf(data=GVTPObservationAreas[,"CNAME"],aes(fill=factor(CNAME)),inherit.aes=F)+
  scale_fill_discrete(na.value="white",na.translate=F)+
  labs(title="Gaviota tarplant observation areas",x="longitude",y="latitude",fill='OBJECTID')+
  theme(legend.position="none")+
  coord_sf(crs = st_crs(4326))


#Generate partial response plots of importance variables in Gaviota tarplant SDM.
require(zoo)
require(ggplot2)
#Partial response plot of geology
test <- as.data.frame(response(x=xm,var="GVTPGeology"))
test <- as.data.frame(approx(x=test$V1,y=test$p,xout=1:68))
test <- dplyr::mutate(test,y=na.spline(y))
GeologyCategories <- read.table("GeologyMetadata.tsv", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote="", encoding = "UTF-8")
test <- dplyr::left_join(test,GeologyCategories,by=c("x"="GeologicalCategoryNumber"))
test <- dplyr::arrange(test,y) %>%
  mutate(name=factor(GeologicalCategory,levels=GeologicalCategory))
ggplot(test,aes(x=name,y=y))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  xlab("Geological categories")+ylab("Probability of suitable habitat being observed")+
  ggtitle("Probability of suitable habitat for Gaviota tarplant\nbeing observed by geological category")+geom_point()
#Partial response of interannual cloud cover.
test <- as.data.frame(response(x=xm,var="GVTPInterAnnualCloudCover"))
ggplot(test,aes(x=V1,y=p))+geom_line()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  xlab("Interannual cloud cover\nMean of the 12 monthly standard deviations on % cloud cover")+ylab("Probability of suitable habitat being observed")+
  ggtitle("Probability of suitable habitat for Gaviota tarplant\nbeing observed by interannual cloud cover")
#Partial response plot of bioclim3
test <- as.data.frame(response(x=xm,var="GVTPBioclim3"))
ggplot(test,aes(x=V1,y=p))+geom_line()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  xlab("Isothermality (%)\nMean Diurnal Temperature Range / Mean Annual Temperature Range")+ylab("Probability of suitable habitat being observed")+
  ggtitle("Probability of suitable habitat for Gaviota tarplant\nbeing observed by isothermality")
#Partial response plot of bioclim3
test <- as.data.frame(response(x=xm,var="GVTPBioclim12"))
ggplot(test,aes(x=V1,y=p))+geom_line()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  xlab("Annual precipitation (mm water)")+ylab("Probability of suitable habitat being observed")+
  ggtitle("Probability of suitable habitat for Gaviota tarplant\nbeing observed by annual rainfall") 
#Partial response plot of landform
test <- as.data.frame(response(x=xm,var="GVTPLandform"))
test <- as.data.frame(approx(x=test$V1,y=test$p,xout=seq(1000,9000,1000)))
test <- dplyr::mutate(test,y=na.spline(y))
LandformCategories <- read.table("Landforms.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote="", encoding = "UTF-8")
test <- dplyr::left_join(test,LandformCategories,by=c("x"="Value"))
test <- dplyr::arrange(test,y) %>%
  mutate(name=factor(Landform,levels=Landform))
ggplot(test,aes(x=name,y=y))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  xlab("Landform categories")+ylab("Probability of suitable habitat being observed")+
  ggtitle("Probability of suitable habitat for Gaviota tarplant\nbeing observed by landform category")+geom_point()
