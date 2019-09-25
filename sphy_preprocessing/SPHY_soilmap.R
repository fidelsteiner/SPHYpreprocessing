################################################################################
# SPHY_soilmap - preprocess all necessary soil related data for the SPHY mode
# 
# SPHY_soilmap.R
#
# ReadMe: 
# Necessary data for the SPHY model are (a) soil depth, (b) field capacity of the root zone, (c) wilting point (pF3), (d) permanent wilting point (pF 4.2), (e) saturated water content and (f) saturated hydrualic conductivity.
# All variables if available for the top soil (0-30 cm), as well as the subsoil (<30 cm)
#
# Soil depth is freely available from https://soilgrids.org or can be computed from epirical relations with slope and curvature, all provided here
#
# All soil properties are available from the HiHydroSoil product, which is based on data from soilgrids (https://www.futurewater.eu/2015/07/soil-hydraulic-properties/)
#
# The code below uses this data for any location according to the provided DEM and saves it as an ascii file that needs to be converted into a .map file using the following command in python
#
# asc2map --clone clone.map -S -a -m -9999 ---.asc ---.map
#
# Output: 
# (1) soil thickness (3 versions, based on soilgrids, a slope and a curvature model)
# (2) saturated water content in top and sub soil (wc_top_sat.asc, wc_sub_sat.asc)
# (3) hydraulic conductivity top and sub soil (k_top_sat.asc, k_sub_sat.asc)
# (4) field capacity of the top and subsoil (rootfield.asc, sub_field.asc)
# (5) wilting and permanent wilting point for the rootzone (root_wilt.asc and root_dry.asc)
#
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
p_load(rgdal,rgeos,maptools,raster,rasterVis)
library(spatialEco)
##########################
# SPECIFY FILENAMES AND DESIRED PROJECTION
##########################
RGI_filename <- '15_rgi60_SouthAsiaEast.shp' 

# Define Projection
projec_lambers <- "+proj=laea +lat_0=28 +lon_0=95 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
projec_utm <- '+proj=utm +zone=45N +datum=WGS84'

projec <- projec_utm

# PATHS FOR GLACIER RAW DATA (adapt and make sure all data is available)
path_maps <- 'F:\\PhD\\Research\\SPHY\\SPHYLangtang\\Code\\SPHY-master\\input'    # Folder of SPHY input data, used as output folder for all data
path_output <- 'F:\\PhD\\Research\\SPHY\\SPHYLangtang\\figures'                   # Folder for all figures
path_soildata <- 'F:\\PhD\\Research\\SPHY\\data\\soil'                            # Folder for raw data
path_RGI <- 'F:\\PhD\\GeoSpatialData\\RGI60_Asia'                                 # Folder for glacier outlines

file_soil <- 'dominantsoil_npl.shp'                                               # soil data from nepal survey (optional)
file_rock <- 'dominantrock_npl.shp'                                               # rock data from nepal survey (optional)
file_1km <- 'BDRICM_M_250m.tif'                                                   # soil thickness data
 
# if available, provide .shp for catchment outline for visualisation, set 'cSHP' to 1
cSHP <- 1
file_catchmentoutline <- 'F:\\PhD\\Research\\SPHY\\CatchmentMapping\\OutlineLangtang\\catch_proj_1.shp' # Catchment Outline (shp file)

# if available, provide raster file with data where no soil available, set 'bare' to 1
bare <- 1
file_barerock <- 'F:\\PhD\\Research\\SPHY\\data\\barerock_langtang100m.tif'

##########################
# Load SPHY data
##########################
domain <- raster(paste(path_maps,'\\clone.map',sep=''))   # Load the domain
projection(domain) <- projec

dem <- raster(paste(path_maps,'\\dem.map',sep=''))        # Load the DEM
projection(dem) <- projec
contourDEM <- rasterToContour(dem,levels = pretty(range(dem[], na.rm = TRUE), 10))

domain_deg <- projectRaster(domain,crs = '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')

modID_raster <- domain
modID_raster[seq(1,dim(domain)[1]*dim(domain)[2],1)] <- seq(1,dim(domain)[1]*dim(domain)[2],1)

if(bare ==1){
  bareRaster <- raster(file_barerock)
  projection(bareRaster) <- projec
  bareRaster <- resample(bareRaster,domain)
  bareRaster <- mask(bareRaster,dem)
}

##########################
# Load RGI and other catchment data
##########################
ogrInfo(path_RGI&'\\'&RGI_filename)
RGI60_15<-readOGR(dsn=path_RGI&'\\'&RGI_filename)
RGI60_15 <- spTransform(RGI60_15, projec)

# Restrict the datasets to the domain
sub_15 <- subset(RGI60_15, CenLon >= extent(domain_deg)[1] & CenLon <= extent(domain_deg)[2] & CenLat >= extent(domain_deg)[3] & CenLat <= extent(domain_deg)[4])

r_glaciermask <- rasterize(sub_15, domain,mask=F) # Individual glacier numbers (only pixels that lie on the polygon)
r_glaciermask[is.na(r_glaciermask)] <- 0
r_glaciermask <- mask(r_glaciermask,dem)

# Load Catchment outline data if available
if(cSHP == 1) {
  ogrInfo(file_catchmentoutline)
  catch<-readOGR(dsn=file_catchmentoutline)
  catch <- spTransform(catch, projec)
}

# Load soil data from Nepal survey if provided
#ogrInfo(path_soildata&'\\'&file_soil)
#soil_map<-readOGR(dsn=path_soildata&'\\'&file_soil)
#soil_map <- spTransform(soil_map, projec)
#r_soilmask <- rasterize(soil_map, domain,mask=F,field='LCODE_GEN')
#r_soilmask <- mask(r_soilmask,dem)

#r_soilmask_fac<-as.factor(r_soilmask)
#tar<-levels(r_soilmask_fac)[[1]]
#tar[["soil"]]<-c("6","7","8","9","18")
#levels(r_soilmask_fac)<-tar

##########################
# Load HiHydroSoilmaps and process to domain
##########################
WCSat_top <- raster('F:\\PhD\\Research\\SPHY\\data\\soil\\wcsat_topsoil_gapfilled_HMA.tif') 
WCSat_top <- projectRaster(WCSat_top,crs=projec,method="bilinear")
WCSat_top <- resample(WCSat_top,dem)        # Saturated Water Content top soil (0 - 0.3 m)
WCSat_top <- mask(WCSat_top,dem) / 10000       # [m3/m3]
WCSat_top[!is.na(bareRaster)] <- 0
WCSat_top[r_glaciermask > 0 ] <- 0

WCSat_sub <- raster('F:\\PhD\\Research\\SPHY\\data\\soil\\wcsat_subsoil_gapfilled_HMA.tif') 
WCSat_sub <- projectRaster(WCSat_sub,crs=projec,method="bilinear")
WCSat_sub <- resample(WCSat_sub,dem)        # Saturated Water Content bottom soil (<0.3 m)
WCSat_sub <- mask(WCSat_sub,dem) / 10000       # [m3/m3]
WCSat_sub[!is.na(bareRaster)] <- 0
WCSat_sub[r_glaciermask > 0 ] <- 0

ksat_top <- raster('F:\\PhD\\Research\\SPHY\\data\\soil\\ksat_topsoil_gapfilled_HMA.tif') 
ksat_top <- projectRaster(ksat_top,crs=projec,method="bilinear")
ksat_top <- resample(ksat_top,dem)        # hydraulic conductivity top soil (0 - 0.3 m)
ksat_top <- mask(ksat_top,dem)  / 10000 / 100 # [m /d]
ksat_top[!is.na(bareRaster)] <- 0
ksat_top[r_glaciermask > 0 ] <- 0

ksat_sub <- raster('F:\\PhD\\Research\\SPHY\\data\\soil\\ksat_subsoil_gapfilled_HMA.tif') 
ksat_sub <- projectRaster(ksat_sub,crs=projec,method="bilinear")
ksat_sub <- resample(ksat_sub,dem)        # hydraulic conductivity bottom soil (<0.3 m)
ksat_sub <- mask(ksat_sub,dem)  / 10000 / 100 # [m /d]
ksat_sub[!is.na(bareRaster)] <- 0
ksat_sub[r_glaciermask > 0 ] <- 0

fieldCap_top <- raster('F:\\PhD\\Research\\SPHY\\data\\soil\\wcpf2_topsoil_gapfilled_HMA.tif') 
fieldCap_top <- projectRaster(fieldCap_top,crs=projec,method="bilinear")
fieldCap_top <- resample(fieldCap_top,dem)        # field capacity top soil (0 - 0.3 m)
fieldCap_top <- mask(fieldCap_top,dem)  / 10000  # [m3 /m3]
fieldCap_top[!is.na(bareRaster)] <- 0
fieldCap_top[r_glaciermask > 0 ] <- 0

fieldCap_sub <- raster('F:\\PhD\\Research\\SPHY\\data\\soil\\wcpf2_subsoil_gapfilled_HMA.tif') 
fieldCap_sub <- projectRaster(fieldCap_sub,crs=projec,method="bilinear")
fieldCap_sub <- resample(fieldCap_sub,dem)        # field capacity bottom soil (<0.3 m)
fieldCap_sub <- mask(fieldCap_sub,dem)  / 10000  # [m3 / m3]
fieldCap_sub[!is.na(bareRaster)] <- 0
fieldCap_sub[r_glaciermask > 0 ] <- 0

wilting_top <- raster('F:\\PhD\\Research\\SPHY\\data\\soil\\wcpf3_topsoil_gapfilled_HMA.tif') 
wilting_top <- projectRaster(wilting_top,crs=projec,method="bilinear")
wilting_top <- resample(wilting_top,dem)        # wilting point top soil (0 - 0.3 m)
wilting_top <- mask(wilting_top,dem)  / 10000  # [m3 /3]
wilting_top[!is.na(bareRaster)] <- 0
wilting_top[r_glaciermask > 0 ] <- 0

wilting_sub <- raster('F:\\PhD\\Research\\SPHY\\data\\soil\\wcpf3_subsoil_gapfilled_HMA.tif') 
wilting_sub <- projectRaster(wilting_sub,crs=projec,method="bilinear")
wilting_sub <- resample(wilting_sub,dem)        # wilting point bottom soil (<0.3 m)
wilting_sub <- mask(wilting_sub,dem)  / 10000  # [m3 / m3]
wilting_sub[!is.na(bareRaster)] <- 0
wilting_sub[r_glaciermask > 0 ] <- 0

pwilting_top <- raster('F:\\PhD\\Research\\SPHY\\data\\soil\\wcpf4.2_topsoil_gapfilled_HMA.tif') 
pwilting_top <- projectRaster(pwilting_top,crs=projec,method="bilinear")
pwilting_top <- resample(pwilting_top,dem)        # wilting point top soil (0 - 0.3 m)
pwilting_top <- mask(pwilting_top,dem)  / 10000  # [m3 /3]
pwilting_top[!is.na(bareRaster)] <- 0
pwilting_top[r_glaciermask > 0 ] <- 0

pwilting_sub <- raster('F:\\PhD\\Research\\SPHY\\data\\soil\\wcpf4.2_subsoil_gapfilled_HMA.tif') 
pwilting_sub <- projectRaster(pwilting_sub,crs=projec,method="bilinear")
pwilting_sub <- resample(pwilting_sub,dem)        # wilting point bottom soil (<0.3 m)
pwilting_sub <- mask(pwilting_sub,dem)  / 10000  # [m3 / m3]
pwilting_sub[!is.na(bareRaster)] <- 0
pwilting_sub[r_glaciermask > 0 ] <- 0

##########################
# Load soil depth and create additional soil depth maps
##########################
depth1kmdata <- raster(path_soildata&'\\'&file_1km)
depth1kmdata <- projectRaster(depth1kmdata,crs=projec,method="bilinear")

depth1kmdata <- resample(depth1kmdata,dem) / 100
depth1kmdata <- mask(depth1kmdata,dem)
#depth1kmdata[is.na(depth1kmdata)] <- 0
depth1kmdata[!is.na(bareRaster)] <- 0
depth1kmdata[r_glaciermask > 0 ] <- 0

slopeMap <- terrain(dem, opt='slope', unit='degrees', neighbors=8)
slopeMap <- mask(slopeMap,dem)
curvatureMap <- curvature(dem,type = c("profile"))
curvatureMap <- mask(curvatureMap,dem)

# Pelletier/Heimsath slope thickness models
h_0 <- 0.5                      # characteristic soil depth (standard literature value)
P_0 <- 1 / 1000                 # denudation (m / yr) (see for example https://www.earth-surf-dynam-discuss.net/esurf-2019-7/)
U <- 0.3 / 1000                 # orogenic uplift (m / yr)
kappa <- 0.5468                 # Hillslope diffusivity [L^2/T] (Culling 1963)

# based on slope, Heimsath, 1997 (https://www.nature.com/articles/41056) and Pelletier 2009 (https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2008WR007319)  
H <- h_0 * cos(slopeMap*pi/180) * log(P_0 / U  / cos(slopeMap*pi/180))
H[!is.na(bareRaster)] <- 0
H[r_glaciermask > 0 ] <- 0

# based on curvature, Pelletier, 2016 - https://agupubs.onlinelibrary.wiley.com/doi/10.1002/2015MS000526
curvPel <- -P_0 / kappa / curvatureMap
curvPel[curvPel<1] <- 1
H_n <- h_0 * log(curvPel)
H_n[is.infinite(H_n)]<- 0
H_n <- focal(H_n, w=matrix(1/9,nrow=3,ncol=3))
H_n[!is.na(bareRaster)] <- 0
H_n[r_glaciermask > 0 ] <- 0

##########################
# Save final data in ascii format
##########################
writeRaster(depth1kmdata, file.path(path_maps, "H_SoilGrid.asc"), format="ascii",overwrite=T)
writeRaster(H, file.path(path_maps, "H_slope.asc"), format="ascii",overwrite=T)
writeRaster(H_n, file.path(path_maps, "H_curvature.asc"), format="ascii",overwrite=T)

writeRaster(WCSat_top, file.path(path_maps, "wc_top_sat.asc"), format="ascii",overwrite=T)
writeRaster(WCSat_sub, file.path(path_maps, "wc_sub_sat.asc"), format="ascii",overwrite=T)
writeRaster(fieldCap_top, file.path(path_maps, "rootfield.asc"), format="ascii",overwrite=T)
writeRaster(pwilting_top, file.path(path_maps, "root_dry.asc"), format="ascii",overwrite=T)
writeRaster(wilting_top, file.path(path_maps, "root_wilt.asc"), format="ascii",overwrite=T)
writeRaster(ksat_top, file.path(path_maps, "k_top_sat.asc"), format="ascii",overwrite=T)
writeRaster(fieldCap_sub, file.path(path_maps, "sub_field.asc"), format="ascii",overwrite=T)
writeRaster(ksat_sub, file.path(path_maps, "k_sub_sat.asc"), format="ascii",overwrite=T)

##########################
# Visualize domain data
##########################

png(file=path_output&'\\soildepth.png', res = 300,width=dim(dem)[2]*9,height=dim(dem)[1]*6)
par(mar=c(4,5,2,1),cex.lab=1.5,cex.axis=1.5)
par(mfrow=c(2,4))
layout(matrix(c(1,2,3,,4,5,6,), nrow = 2, ncol = 4, bycol = T))
plot(depth1kmdata,xlab="Easting [m]",ylab="Northing [m]",zlim=c(0,2),legend=F)
plot(sub_15,add=T,border='blue',lwd=0.5)
plot(contourDEM,add=T,lwd=0.5)
if(cSHP == 1) {
  plot(catch,add=T)
} 
plot(H,xlab="Easting [m]",ylab="",zlim=c(0,2),legend=F)
plot(sub_15,add=T,border='blue',lwd=0.5)
plot(contourDEM,add=T,lwd=0.5)
if(cSHP == 1) {
  plot(catch,add=T)
} 
plot(H_n,xlab="Easting [m]",ylab="",zlim=c(0,2),legend.args=list(text='Soil Depth [m]', side=2, font=1, line=-4, cex=1))
plot(sub_15,add=T,border='blue',lwd=0.5)
plot(contourDEM,add=T,lwd=0.5)
if(cSHP == 1) {
  plot(catch,add=T)
} 
plot.new()
plot(depth1kmdata[],dem[],xlab="Depth [m]",ylab="Elevation [m]")
plot(H[],slopeMap[],xlab="Depth [m]",ylab="slope [deg]")
plot(H_n[],curvatureMap[],xlab="Depth [m]",ylab="curvature [-]")
dev.off()

#png(file=path_output&'\\soiltype.png', res = 160,width=dim(domain)[2]*2,height=dim(domain)[1]*2)
#par(mar=c(4,5,2,2),cex.lab=1,cex.axis=1)
##par(mfrow=c(1))
#layout(matrix(c(1), nrow = 1, ncol = 1, byrow = FALSE))
#levelplot(r_soilmask_fac,xlab="Easting [m]",ylab="Northing [m]",col.regions=brewer.pal(5,"Set3")) + layer(sp.polygons(catch))
#if(cSHP == 1) {
#  plot(catch,add=T)
#} 
#dev.off()

png(file=path_output&'\\SoilProps.png', res = 300,width=dim(dem)[2]*9,height=dim(dem)[1]*3)
par(mar=c(4,5,2,5),cex.lab=1.5,cex.axis=1.5)
par(mfrow=c(1,4))
layout(matrix(c(1,2,3,4), nrow = 1, ncol = 4, bycol = T))
plot(ksat_top,zlim=c(0,1),xlab="Easting [m]",ylab="Northing [m]",legend.args=list(text='k [m / d]', side=2, font=1, line=-5, cex=1))
plot(sub_15,add=T,border='blue',lwd=0.5)
plot(contourDEM,add=T,lwd=0.5)
if(cSHP == 1) {
  plot(catch,add=T)
} 

plot(fieldCap_top,zlim=c(0.2,0.4),xlab="Easting [m]",ylab="",legend.args=list(text='field capacity [m3 / m3]', side=2, font=1, line=-5, cex=1))
plot(sub_15,add=T,border='blue',lwd=0.5)
plot(contourDEM,add=T,lwd=0.5)
if(cSHP == 1) {
  plot(catch,add=T)
} 

plot(WCSat_top,zlim=c(0.5,0.6),xlab="Easting [m]",ylab="",legend.args=list(text='saturated water content [m3 / m3]', side=2, font=1, line=-5, cex=1))
plot(sub_15,add=T,border='blue',lwd=0.5)
plot(contourDEM,add=T,lwd=0.5)
if(cSHP == 1) {
  plot(catch,add=T)
} 
dev.off()