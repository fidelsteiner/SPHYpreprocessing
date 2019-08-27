################################################################################
# SPHY_glacmap - create IDmaps and glac_IDs for the SPHY model
# 
# SPHY_glacmap.R
#
# ReadMe:
# Uses RGI6.0 glacier inventory + associated debris cover and ice thickness products to create raster product for a distributed hydrological model as well as the 
# glac_table.csv necessary for the SPHY3.0 model
#
# RGI6.0 (http://www.glims.org/RGI/rgi60_dl.html) and the two datasets supplied there (Scherler2018, Farinotti2019) are neccessary to run the code and need to be stored locally
#
# Additionally a clone.map file marking the domain and resolution as well as a dem.map file with the DEM needs to be available
#
# Outputs are saved as .asc files that can be converted in PCRaster to .map files by the following command: asc2map --clone \clone.map -S -a \glacfrac.asc \glacfrac.map
# The files glacfrac.asc / glacfrac_db.asc / glacfrac_ic.asc / glacID.asc / icedepth.asc / modID.asc need to be converted
# Created:          2019/08/15
# Latest Revision:  2019/08/22
#
# Jakob F Steiner| PhD candidate | Faculty of Geosciences | Universiteit Utrecht | Princetonlaan 8a, 3584 CB Utrecht 
# Vening Meinesz building, room 4.30 | P.O. Box 80.115, 3508 TC Utrecht | j.f.steiner@uu.nl | www.uu.nl/staff/jfsteiner | www.mountainhydrology.org 
################################################################################
# clear entire workspace (excl. packages)
rm(list = ls())
gc()

# define &-sign for pasting string-elements
'&' <- function(...) UseMethod('&')
'&.default' <- .Primitive('&')
'&.character' <- function(...) paste(...,sep='')

# install necessary packages if not available yet via install.packages()
library(pacman)
p_load(rgdal,rgeos,maptools,raster)

########################## 
# PATHS FOR GLACIER RAW DATA (adapt and make sure all data is available)
path_maps <- 'F:\\PhD\\Research\\SPHY\\trisulitestSPHY\\Code\\SPHY-master\\input'             # Folder of SPHY input data
path_output <- 'F:\\PhD\\Research\\SPHY\\trisulitestSPHY\\Code\\SPHY-master\\input'             # Folder for all outputs
path_RGI <- 'F:\\PhD\\GeoSpatialData\\RGI60_Asia'                                             #folder with (RGI 6.0)
path_debris <- 'F:\\PhD\\GeoSpatialData\\DCG_Scherler\\S2_2015-2017_NDSI\\S2_2015-2017_NDSI'  # folder debris data (Scherler 2018)
path_thick <- 'F:\\PhD\\GeoSpatialData\\IceThickness_Farinotti\\RGI60-15\\RGI60-15'           # folder with ice thickness data (Farinotti 2019)

RGI_filename <- '15_rgi60_SouthAsiaEast.shp'
RGI_debris_filename <- '15_rgi60_SouthAsiaEast_S2_DC_2015_2017_NDSI.shp'

##########################

##########################
# Load SPHY data
projec<-'+proj=utm +zone=45N +datum=WGS84'
domain <- raster(paste(path_maps,'\\clone.map',sep=''))   # Load the domain
projection(domain) <- projec

dem <- raster(paste(path_maps,'\\dem.map',sep=''))        # Load the DEM
projection(dem) <- projec

domain_deg <- projectRaster(domain,crs = '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')

modID_raster <- domain
modID_raster[seq(1,dim(domain)[1]*dim(domain)[2],1)] <- seq(1,dim(domain)[1]*dim(domain)[2],1)
##########################

##########################
# Load RGI data
ogrInfo(path_RGI&'\\'&RGI_filename)
RGI60_15<-readOGR(dsn=path_RGI&'\\'&RGI_filename)
ogrInfo(path_debris&'\\'&RGI_debris_filename)
RGI60_15_debris<-readOGR(dsn=path_debris&'\\'&RGI_debris_filename)
projection(RGI60_15_debris) <- projection(RGI60_15)

RGI60_15 <- spTransform(RGI60_15, projec)
RGI60_15_debris <- spTransform(RGI60_15_debris, projec)
  
# Restrict the datasets to the domain
sub_15 <- subset(RGI60_15, CenLon >= extent(domain_deg)[1] & CenLon <= extent(domain_deg)[2] & CenLat >= extent(domain_deg)[3] & CenLat <= extent(domain_deg)[4])
sub_15_debris <- subset(RGI60_15_debris, CenLon >= extent(domain_deg)[1] & CenLon <= extent(domain_deg)[2] & CenLat >= extent(domain_deg)[3] & CenLat <= extent(domain_deg)[4])
  
# Find ice thickness for all glaciers and combine into 1 raster
mergeRaster <- domain * 0
  for(i in 1:length(sub_15$RGIId)){
    thick_ <- raster(path_thick&'\\'&sub_15$RGIId[i]&'_thickness.tif')
    thick_ <- aggregate(thick_,fact=res(domain)[1]/res(thick_)[1])
    thick_ <- projectRaster(thick_,crs = projection(domain))
    if(is.null(intersect(extent(thick_), extent(domain)))){}
    else{
    thick_ <- resample(thick_,domain)
    thick_[is.na(thick_)] <- 0
    mergeRaster <- mergeRaster + thick_
    }
  }
  
# discard ice thickness < 2m
mergeRaster[mergeRaster < 2] <- 0
########################## 
 
##########################
# Rasterize all datasets and save
r_glaciermask_cover <- rasterize(sub_15, domain,getCover=TRUE) # relative cover of glacier over pixel
r_glaciermask <- rasterize(sub_15, domain,mask=F) # Individual glacier numbers (only pixels that lie on the polygon)
r_glaciermask[is.na(r_glaciermask)] <- 0
  
r_glaciermask_debris_cover <- rasterize(sub_15_debris, domain,getCover=TRUE) # relative cover of debris over pixel
r_glaciermask_debris <- rasterize(sub_15_debris, domain,mask=F) # Individual glacier numbers (only pixels that lie on the polygon)
r_glaciermask_debris[is.na(r_glaciermask_debris)] <- 0
  
# Save all rasters (GTiff is an option but ascii files are needed for PCRaster/SPHY)
#writeRaster(r_glaciermask_cover, filename=file.path(path_maps, "glacfrac"), format="GTiff",overwrite=T)
#writeRaster(r_glaciermask, filename=file.path(path_maps, "glacID"), format="GTiff",overwrite=T) 
writeRaster(r_glaciermask_cover, file.path(path_maps, "glacfrac.asc"), format="ascii")
writeRaster(r_glaciermask, file.path(path_maps, "glacID.asc"), format="ascii")

writeRaster(r_glaciermask_cover - r_glaciermask_debris_cover, file.path(path_maps, "glacfrac_ci.asc"), format="ascii")
writeRaster(r_glaciermask_debris, file.path(path_maps, "glacID_debris.asc"), format="ascii")

#writeRaster(r_glaciermask_debris_cover, filename=file.path(path_maps, "debrisfrac"), format="GTiff",overwrite=T)
writeRaster(r_glaciermask_debris_cover, file.path(path_maps, "glacfrac_db.asc"), format="ascii")


#writeRaster(mergeRaster, filename=file.path(path_maps, "icedepth"), format="GTiff",overwrite=T)
writeRaster(mergeRaster, file.path(path_maps, "icedepth.asc"), format="ascii")

writeRaster(modID_raster, file.path(path_maps, "modID.asc"), format="ascii")
##########################  
 
########################## 
# Make glac_table for SPHY input
  U_ID <- seq(1,length(r_glaciermask[which(r_glaciermask[]>0)]),1)
  MOD_ID <- which(r_glaciermask[]>0)
  GLAC_ID <- r_glaciermask[which(r_glaciermask[]>0)]
  MOD_H <- dem[which(r_glaciermask[]>0)]
  
  # If SPHY Model domain is coarser than 1km, Model elevations and glacier elevations can be different, otherwise it is assumed that they are the same
  if(res(domain)[1]>1000){
  GLAC_H <- levels(r_glaciermask)[[1]]$Zmed[GLAC_ID]
  } else {
    GLAC_H <- MOD_H  
  }
  
  FRAC_DEBRIS <- r_glaciermask_debris_cover[which(r_glaciermask[]>0)]
  DEBRIS <- FRAC_DEBRIS
  DEBRIS[DEBRIS>0] <- 1
  
  FRAC_GLAC <- r_glaciermask_cover[which(r_glaciermask[]>0)]
  
  ICE_DEPTH <- mergeRaster[which(r_glaciermask[]>0)]
 
  
# Save Table for SPHY
glac_table <- cbind(U_ID,MOD_ID,GLAC_ID,MOD_H,GLAC_H,FRAC_DEBRIS,DEBRIS,FRAC_GLAC,ICE_DEPTH)
glac_table <- glac_table[-which(is.na(MOD_H)),] # remove all rows where no DEM is available to save space
colnames(glac_table) <- c('U_ID','MOD_ID','GLAC_ID','MOD_H','GLAC_H','FRAC_DEBRIS','DEBRIS','FRAC_GLAC','ICE_DEPTH')
write.table(glac_table, file = path_output&'\\glac_table.csv', append = FALSE,col.names = T,row.names = F,sep=',')
##########################

##########################
# Visualize domain data
  png(file=path_output&'\\Domain_Glaciers.png', res = 160,width=dim(domain)[1]*2,height=dim(domain)[2]*2)
  par(mar=c(4,5,2,2),cex.lab=1.5,cex.axis=1.5)
  #par(mfrow=c(1))
  layout(matrix(c(1), nrow = 1, ncol = 1, byrow = FALSE))
  plot(dem,xlab="Easting [m]",ylab="Northing [m]",legend.args=list(text='Elevation [m asl]', side=2, font=1, line=0.3, cex=1.5))
  plot(sub_15,add=T)
  plot(sub_15_debris,add=T,col='red')
  legend('bottomright',bty="n",c('glaciers','debris cover'),pch = 1,col = c('black','red'))
  #plot(r_glaciermask_cover - r_glaciermask_debris_cover,xlab="Easting [m]",ylab="Northing [m]",legend.args='')
  dev.off()
##########################  