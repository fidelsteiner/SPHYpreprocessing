################################################################################
#SPHY_glacmap - create IDmaps and glac_IDs for the SPHY model
# 
# SPHY_glacmap.R
#
# ReadMe:
# Uses RGI6.0 glacier inventory + associated debris cover and ice thickness products to create raster product for a distributed hydrological model as well as the 
# glac_table.csv necessary for the SPHY3.0 model
#
# Created:          2018/02/05
# Latest Revision:  2017/02/05
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

library(pacman)
p_load(rgdal,rgeos,maptools,raster)

# PATHS FOR GLACIER RAW DATA
path_maps <- 'F:\\PhD\\Research\\SPHY\\trisulitestSPHY\\Code\\SPHY-master\\input'             # Folder of SPHY input data
path_output <- 'F:\\PhD\\Research\\SPHY\\trisulitestSPHY\\Code\\SPHY-master\\input'             # Folder for all outputs
path_RGI <- 'F:\\PhD\\GeoSpatialData\\RGI60_Asia'                                             #(RGI 6.0)
path_debris <- 'F:\\PhD\\GeoSpatialData\\DCG_Scherler\\S2_2015-2017_NDSI\\S2_2015-2017_NDSI'  # debris data (Scherler 2018)
path_thick <- 'F:\\PhD\\GeoSpatialData\\IceThickness_Farinotti\\RGI60-15\\RGI60-15'           # ice thickness data (Farinotti 2019)

# Make mod_id and glac_id maps

projec<-'+proj=utm +zone=45N +datum=WGS84'
domain <- raster(paste(path_maps,'\\clone.map',sep=''))
projection(domain) <- projec

dem <- raster(paste(path_maps,'\\dem.map',sep=''))
projection(dem) <- projec

domain_deg <- projectRaster(domain,crs = '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')

ogrInfo(path_RGI&'\\14_rgi60_SouthAsiaWest.shp')
RGI60_14<-readOGR(dsn=path_RGI&'\\14_rgi60_SouthAsiaWest.shp')
ogrInfo(path_RGI&'\\15_rgi60_SouthAsiaEast.shp')
RGI60_15<-readOGR(dsn=path_RGI&'\\15_rgi60_SouthAsiaEast.shp')
ogrInfo(path_debris&'\\15_rgi60_SouthAsiaEast_S2_DC_2015_2017_NDSI.shp')
RGI60_15_debris<-readOGR(dsn=path_debris&'\\15_rgi60_SouthAsiaEast_S2_DC_2015_2017_NDSI.shp')

# Read in Surface Features with different melt properties

  projection(RGI60_15_debris) <- projection(RGI60_15)
  RGI60_14 <- spTransform(RGI60_14, projec)
  RGI60_15 <- spTransform(RGI60_15, projec)
  RGI60_15_debris <- spTransform(RGI60_15_debris, projec)
  
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
  
  mergeRaster[mergeRaster < 2] <- 0
  
  r_glaciermask_cover <- rasterize(sub_15, domain,getCover=TRUE) # relative cover of cliff over pixel
  r_glaciermask <- rasterize(sub_15, domain,mask=F) # Individual cliff numbers (only pixels that lie on the polygon)
  r_glaciermask[is.na(r_glaciermask)] <- 0
  
  r_glaciermask_debris_cover <- rasterize(sub_15_debris, domain,getCover=TRUE) # relative cover of cliff over pixel
  r_glaciermask_debris <- rasterize(sub_15_debris, domain,mask=F) # Individual cliff numbers (only pixels that lie on the polygon)
  r_glaciermask_debris[is.na(r_glaciermask_debris)] <- 0
  
  
  writeRaster(r_glaciermask_cover, filename=file.path(path_maps, "glacfrac_1"), format="GTiff",overwrite=T)
  writeRaster(r_glaciermask, filename=file.path(path_maps, "glacID"), format="GTiff",overwrite=T) 
  
  writeRaster(r_glaciermask_debris_cover, filename=file.path(path_maps, "debrisfrac"), format="GTiff",overwrite=T)
  writeRaster(r_glaciermask_debris, filename=file.path(path_maps, "glacID_debris"), format="GTiff",overwrite=T) 

  writeRaster(mergeRaster, filename=file.path(path_maps, "icedepth"), format="GTiff",overwrite=T)
  
  
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
  
  glac_table <- cbind(U_ID,MOD_ID,GLAC_ID,MOD_H,GLAC_H,FRAC_DEBRIS,DEBRIS,FRAC_GLAC,ICE_DEPTH)
  glac_table <- glac_table[-which(is.na(MOD_H)),]
  colnames(glac_table) <- c('U_ID','MOD_ID','GLAC_ID','MOD_H','GLAC_H','FRAC_DEBRIS','DEBRIS','FRAC_GLAC','ICE_DEPTH')
  write.table(glac_table, file = path_output&'\\glac_table.csv', append = FALSE,col.names = T,row.names = F,sep=',')
  
  png(file=path_output&'\\Domain_Glaciers.png', res = 160,width=dim(domain)[1]*2,height=dim(domain)[2]*2)
  par(mar=c(4,5,2,2),cex.lab=1.5,cex.axis=1.5)
  #par(mfrow=c(3,1))
  layout(matrix(c(1), nrow = 1, ncol = 1, byrow = FALSE))
  plot(dem,xlab="Easting [m]",ylab="Northing [m]",legend.args=list(text='Elevation [m asl]', side=2, font=1, line=0.3, cex=1.5))
  plot(sub_15,add=T)
  plot(sub_15_debris,add=T,col='red')
  legend('bottomright',bty="n",c('glaciers','debris cover'),pch = 1,col = c('black','red'))
  dev.off()