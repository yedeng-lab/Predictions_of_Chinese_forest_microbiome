# R code for "Warming-driven migration of core microbiota projects soil property changes at continental scale"
# Written by Kai Feng, please contact by kaifeng_st@rcees.ac.cn
# Corresponding author: Prof. Ye Deng, yedeng@rcees.ac.cn
# Research Center for Eco-Environmental Sciences, Chinese Academy of Sciences

library(randomForest)
library(parallel)
library(ppcor)
library(pheatmap)
library(raster)
library(sp)
library(ggplot2)
library(maptools)
library(gstat)
library(Cubist)

#---------------------------------Instruction------------------------------------------------------------------#
# 1. Prepare directories with following:                                                                       #   
# --------RF     # Store the results from random forest                                                        # 
# --------map    # Store the maps based on projection                                                          # 
# --------map/raster # Store the raster files of the projected map                                             # 
# 2. Please prepare the OTU table, taxonomy information file and the environmental variables table for samples # 
# 3. Please prepare the extracted covariates table for Chinese gridded map based on current climate model and  # 
#    future climate change models, e.g. BCC-CSM1-1                                                             # 
# 4. Please change the parameters or values according to own case                                              # 
#--------------------------------------------------------------------------------------------------------------# 


# Import the three table of OTU, taxonomy and environmental variables
tax.orig <- read.table("Classifier_of_16srrna_1654_uparse_5reads.txt",sep = "\t",header = T,row.names = 1)
otu.orig <- read.table("resample_normalized_uparse_table_1654_5reads_10000.txt",sep="\t",header=T,row.names=1)
env.orig <- read.table("Meta_info_1411_4RF_ecological_envs.txt",sep="\t",header=T,row.names=1)

# remove the chloroplast and mitochondria
sort(unique(tax.orig$Phylum)) 
ex <- rownames(tax.orig)[grep("Chloroplast",tax.orig$Phylum)]

otu.norm <- otu.orig[-match(ex,rownames(otu.orig)),]
otu.norm <- otu.norm[,match(rownames(env.orig),colnames(otu.norm))]
otu.norm <- otu.norm[-c(which(rowSums(otu.norm)==0)),]

tax <- tax.orig[match(rownames(otu.norm),rownames(tax.orig)),]
otu <- otu.norm
env <- env.orig

##---------a. Cubist model construction for core phylotypes based on ecological cluster ------------------------
## find core phylotypes
# select 1411 samples, 10% top dominant OTUs, >40% occurred in 1411 samples 564 samples
otu1 <- otu[,rownames(env)]
otu.core1 <- otu1[names(sort(rowSums(otu1),decreasing = T))[1:ceiling(0.1*nrow(otu1))],]
sum(otu.core1)/sum(otu)  # 0.8794186
otu.core2 <- otu.core1[names(which(sort(rowSums(otu.core1>0),decreasing = T)>floor(0.4*ncol(otu.core1)))),]
sum(otu.core2)/sum(otu) # 0.3970043

## Random forest model to select OTUs associated with environmental preference
if(!dir.exists("RF")){ dir.create("RF") } # store all results from random forest
# Covariates Combination #1
sl <- c("pH","SOC","AMT","MDR","TemperatureSeasonality","MTWM","MTCM","AMP",
        "PWM","PDM","PrecipitationSeasonality","CEC","NDVI","NPP","im", "Aridity",
        "NPP2010_CN","NDVI2015_CN","Clay_CN","Sand_CN","Silt_CN")
env.st <- env[,sl]
env.sd <- decostand(env.st,method = "standardize")
library(randomForest)
library(parallel)
report <- c()
for(i in 1:nrow(otu.core2)){
  dat <- cbind(env.sd,t(otu.core2)[,i])
  t.rf <- randomForest(x=dat[,-ncol(dat)],y=dat[,ncol(dat)],importance=TRUE, proximity=TRUE, ntree=5000, norm.votes = FALSE)
  imp <- t(importance(t.rf)[,1])
  explained <- mean(t.rf$rsq *100)
  best_predictor <- names(sort(importance(t.rf)[,1],decreasing = T)[1])
  res <- cbind(imp,explained,best_predictor)
  report <- rbind(report,res)
}
colnames(report) <- c(sl,"explained","best_predictor")
rownames(report) <- rownames(otu.core2)
parallel::stopCluster(cl)
colnames(report)[which(colnames(report)=="im")] <- "IM"

## Calculate semi-partial spearman correlations
library(ppcor)
# Filtering OTUs above the 30% explanation by random forest results
otu.rf <- rownames(report)[which(report[,length(sl)+1]>30)]
otu.rf.tab <- otu.core2[otu.rf,]
spc.env <- function(j,e,r,sl){
  library(ppcor)
  res <- rbind(lapply(as.list(1:length(sl)),function(i,x1,y) {
    i <- as.numeric(i)
    x.inv <- try(ppcor::spcor.test(y,x1[,i],x1[,-i],method = "spearman"),silent = TRUE)
    if ('try-error' %in% class(x.inv)) {
      return(NULL)
    }else{
      if(x.inv$p.value <= 0.001){
        x.inv$estimate
      }
    }
  },x1=e,y=r[,j]))
  res <- as.data.frame(res)
  return(res)
}
cl <- parallel::makeCluster(2,type = "PSOCK")
spc.report <- parLapply(cl,1:nrow(otu.rf.tab),spc.env,e=env.sd,r=t(otu.rf.tab),sl=sl)
names(spc.report) <- rownames(otu.rf.tab)
parallel::stopCluster(cl)
spcor.res <- data.frame()
for (i in 1:length(spc.report)) {
  spcor.res <- rbind(spcor.res,spc.report[[i]])
}
spcor.res[spcor.res=="NULL"] <- as.numeric(0)
spcor.res <- matrix(as.numeric(unlist(spcor.res)),nrow=nrow(otu.rf.tab))
colnames(spcor.res) <- sl
rownames(spcor.res) <- names(spc.report)
spcor.res[spcor.res=="0"] <-  NA

## Clustering the semi-partial spearman correlation coefficient
library(pheatmap)
spcor.tab <- spcor.res
spcor.tab <- as.matrix(spcor.tab)
spcor.plot <- matrix(as.numeric(as.vector(spcor.tab)),byrow=F,nrow = nrow(spcor.tab))
colnames(spcor.plot) <- colnames(spcor.tab)
rownames(spcor.plot) <- rownames(spcor.tab)
spcor.plot[is.na(spcor.plot)] = 0
row.clus <- hclust(dist(spcor.plot,method="maximum"),method="ward.D2")
Cluster <- cutree(row.clus,k=10)
jpeg("RF/semi-partial correlation 0.001 maximum distance 10 Cluster ward.D2 hclust 764 OTUs.jpg",width = 7,height = 10,units = "in",res = 600)
annotation_row = data.frame(Cluster)
pheatmap(spcor.tab,color= colorRampPalette(c("blue","white", "red"), space = "rgb")(15),
         cluster_rows = row.clus, cluster_cols =T,scale = "none",
         fontsize_row=12,  cutree_rows = length(unique(Cluster)),
         fontsize_col=12, show_rownames = F,main = "",border_color = NA,
         treeheight_row = 50,treeheight_col = 30,
         legend = T,height = NA,width = NA,legend_labels =seq(-0.4,0.45,length.out = 40) ,na_col = "white",
         annotation_row =  annotation_row,annotation_colors = list(rainbow(nk)),annotation_names_row = T)
dev.off()
write.csv(cbind(as.matrix(spcor.tab),Cluster),"RF/764 OTUs 10 cluster groups maximum distance ward.D2 hclust with 0.001 semi-correlation.csv")

# Classifying taxa into ecological clusters.
# two columns in the table, one is the cluster number and another is the ecological cluster names
# e.g. Low pH, High pH, Low SOC etc.
# this process should be done byself
eco.tab <- read.table("RF/764 OTUs 10 cluster groups maximum distance ward.D2 hclust with 0.001 semi-correlation ecological culster preference.txt",sep="\t",header = T,row.names = 1)
eco.grp <- list()
for (i in 1:length(unique(eco.tab[,2]))) {
  eco.grp[[i]] <- rownames(eco.tab)[which(eco.tab[,2]==unique(eco.tab[,2])[i])]
}
names(eco.grp) <- unique(eco.tab[,2])
eco <- names(eco.grp) 
eco <- eco[-c(which(eco==""))]

## use covariates combination #2 to build cubist model
library(raster)
library(sp)
library(ggplot2)
library(maptools)
library(gstat)
if(!dir.exists("map")){ dir.create("map") }
if(!dir.exists("map/raster")){ dir.create("map/raster") }
# load covariates combination #3 based on Chinese gridded map with 0.1 degree resolution, extracted by ArcGis software
CN.grid <- read.csv("RF/china_grid_01_degree_extract_CN_4predict.csv",header = T,row.names = 3)

# load Chinese forest distribution map and boundary
cn.forest <- raster("map/CNforestRas1/CNforestRas1.tif")
cn.boundary <- readShapePoly("map/chinaprovinceborderdata/bou2_4p.shp")

##------------b. Current distribution map for each ecological cluster over forested region -----------------------------
reg.method <- "cubist"
standard = "Relative"
for (sl in eco) {
  otu.sub <- otu.core2[eco.grp[[sl]],]
  
  # select covariates combination #2 to build prediction model
  pre.dat <- cbind(env.st[,c("pH","SOC","AMT","MDR","TemperatureSeasonality","MTWM","MTCM","AMP","PWM","PDM",
                             "PrecipitationSeasonality","CEC","NDVI","NPP","Clay_CN","Sand_CN")],colSums(otu.sub)/colSums(otu)*100)
  colnames(pre.dat)[ncol(pre.dat)] <- "abundance"
  
  # select covariates combination #3 to predict ecological abundance across Chinese forest map, extracted by ArcGis software
  dat4pre <- CN.grid[,c("LONG","LAT","SoilpH","SoilSOC","AMT","MDR","TemperatureSeasonality","MTWM","MTCM","AMP","PWM","PDM",
                        "PrecipitationSeasonality","CEC","NDVI","NPP","CLAY","SAND")]
  dat4pre[dat4pre==0] <- NA
  colnames(dat4pre)[3:ncol(dat4pre)] <- colnames(pre.dat)[1:(ncol(pre.dat)-1)]  
  
  # using cubuist model to predict abundance
  library(Cubist)
  cub.mod <- cubist(x=pre.dat[,-ncol(pre.dat)],y=pre.dat[,ncol(pre.dat)],committees = 100)
  sink(paste("map/raster/",reg.method," predicted ",standard," abundance of ",sl," summary.txt",sep=""))
  print(summary(cub.mod))
  sink()
  pre.abund<-cbind(dat4pre[which(rowSums(is.na(dat4pre))==0),],
                   predict(cub.mod,dat4pre[which(rowSums(is.na(dat4pre))==0),],neighbors = 0))
  colnames(pre.abund)[ncol(dat4pre)+1] <- paste("stardard ",sl," abundance",sep="")
  data.obs <- pre.abund[,c(1,2,ncol(pre.abund))]
  colnames(data.obs) <- c("X","Y","VALUE")
  raster_to_dataframe<- function(inRaster, name){
    rp<- rasterToPoints(inRaster)
    df = data.frame(rp)
    colnames(df) = c("X", "Y", "VALUE")
    df$NAME<- name
    return(df)
  }
  
  # Kriging and raster export
  v<- variogram(object = VALUE~1,data = data.obs,locations =~X+Y)
  v.fit=fit.variogram(v,vgm(model="Exp",psill=1,range=max(v$dist),kappa=30),fit.sills=50)
  kingnewdata<- raster_to_dataframe(cn.forest, "Interpolate location")
  aa1 <- krige(formula=VALUE~1,locations=~X+Y,model=vgm(v.fit$psill,"Exp",v.fit$range,v.fit$kappa),data=data.obs,newdata= kingnewdata, nmax=12, nmin=10)
  normal.map <- rasterFromXYZ(aa1[,c(1,2,3)])
  projection(normal.map) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  writeRaster(normal.map,paste("map/raster/",reg.method," predicted ",standard," abundance of ",sl," for 1411 samples.tif",sep=""),overwrite=TRUE)
}

##------------c. Future distribution map for each ecological cluster over forested regions -------------------------------------
reg.method <- "cubist"
standard = "Relative"
cc.model <- "bcc_csm1"
# four different scenarios for BCC-CSM1-1 model
map.dir <- c("bcc_csm1_1_rcp2_6_2050s_bio_30s_r1i1p1_no_tile_asc","bcc_csm1_1_rcp2_6_2080s_bio_30s_r1i1p1_no_tile_asc",
             "bcc_csm1_1_rcp8_5_2050s_bio_30s_r1i1p1_no_tile_asc","bcc_csm1_1_rcp8_5_2080s_bio_30s_r1i1p1_no_tile_asc")
for (i in map.dir) {
  if(length(grep("2050",i))>0 && length(grep("rcp2_6",i))>0){
    sw.dir <- paste(cc.model,"rcp2_6_2050",sep="")
    if(!dir.exists(sw.dir)){ dir.create(sw.dir) }
    year <- "2050"
    rcp <- "RCP2.6"
  }else if(length(grep("2080",i))>0 && length(grep("rcp2_6",i))>0){
    sw.dir <- paste(cc.model,"rcp2_6_2080",sep="")
    if(!dir.exists(sw.dir)){ dir.create(sw.dir) }
    year <- "2080"
    rcp <- "RCP2.6"
  }else if(length(grep("2050",i))>0 && length(grep("rcp8_5",i))>0){
    sw.dir <- paste(cc.model,"rcp8_5_2050",sep="")
    if(!dir.exists(sw.dir)){ dir.create(sw.dir) }
    year <- "2050"
    rcp <- "RCP8.5"
  }else if(length(grep("2080",i))>0 && length(grep("rcp8_5",i))>0){
    sw.dir <- paste(cc.model,"rcp8_5_2080",sep="")
    if(!dir.exists(sw.dir)){ dir.create(sw.dir) }
    year <- "2080"
    rcp <- "RCP8.5"
  }
  # load covariates combination #4 based on gridded Chinese map, extracted by ArcGis software
  CN.grid.future <- read.table(paste(i,"_cn_grid.txt",sep=""),sep="\t",header=T,row.names=1)
  
  # predict abundance for each ecological clusters under each climate change condition
  for (sl in eco) {
    otu.sub <- otu.core2[eco.grp[[sl]],]
    pre.dat <- cbind(env.st[,c("pH","SOC","AMT","MDR","TemperatureSeasonality","MTWM","MTCM","AMP","PWM","PDM",
                               "PrecipitationSeasonality","CEC","NDVI","NPP","Clay_CN","Sand_CN")],colSums(otu.sub)/colSums(otu)*100)
    colnames(pre.dat)[ncol(pre.dat)] <- "abundance"
    dat4pre <- CN.grid.future[,c("long","lat","SoilpH","SoilSOC","AMT","MDR","TemperatureSeasonality","MTWM","MTCM","AMP","PWM","PDM",
                                 "PrecipitationSeasonality","CEC","NDVI","NPP","CLAY","SAND")]
    dat4pre[dat4pre==0] <- NA
    colnames(dat4pre)[3:ncol(dat4pre)] <- colnames(pre.dat)[1:(ncol(pre.dat)-1)]
    
    library(Cubist)
    cub.mod <- cubist(x=pre.dat[,-ncol(pre.dat)],y=pre.dat[,ncol(pre.dat)],committees = 100)
    sink(paste(sw.dir,"/",reg.method," predicted ",year,rcp," ",standard," abundance of ",sl," summary.txt",sep=""))
    print(summary(cub.mod))
    sink()
    pre.abund<-cbind(dat4pre[which(rowSums(is.na(dat4pre))==0),],
                     predict(cub.mod,dat4pre[which(rowSums(is.na(dat4pre))==0),-c(1,2)],neighbors = 0))
    colnames(pre.abund)[ncol(dat4pre)+1] <- paste("stardard ",sl," abundance",sep="")
    write.csv(pre.abund,paste(sw.dir,"/",reg.method," predicted ",year,rcp," ",standard," abundance of ",sl," for 1411 samples.csv",sep=""))
    data.obs <- pre.abund[,c(1,2,ncol(pre.abund))]
    colnames(data.obs) <- c("X","Y","VALUE")
    v<- variogram(object = VALUE~1,data = data.obs,locations =~X+Y)
    v.fit=fit.variogram(v,vgm(model="Sph",psill=1,range=max(v$dist),kappa=30),fit.sills=50)
    kingnewdata<- raster_to_dataframe(cn.forest, "Interpolate location")
    aa1 <- krige(formula=VALUE~1,locations=~X+Y,model=vgm(v.fit$psill,"Exp",v.fit$range,v.fit$kappa),data=data.obs,newdata= kingnewdata, nmax=12, nmin=10)
    map.pre <- rasterFromXYZ(aa1[,c(1,2,3)])
    projection(map.pre) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
    writeRaster(map.pre,paste(sw.dir,"/",reg.method," predicted ",year,rcp," ",standard," abundance of ",sl," for 1411 samples.tif",sep=""),overwrite=TRUE)
    
    # make difference for predicted future data and current data
    map.now <- raster(paste("map/raster/",reg.method," predicted ",standard," abundance of ",sl," for 1411 samples.tif",sep=""))
    map.diff <- overlay(map.pre, map.now, fun=function(x,y){return(x-y)})
    projection(map.diff) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
    writeRaster(map.diff,paste(sw.dir,"/",reg.method," predicted ",year,rcp," ",standard," difference abundance of ",sl," for 1411 samples.tif",sep=""),overwrite=TRUE)
  }
}

