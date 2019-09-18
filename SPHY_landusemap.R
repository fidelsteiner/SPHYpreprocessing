################################################################################
# SPHY_landusemap - create create landuse map for SPHY model from CCI global landuse dataset
# 
# SPHY_landusemap.R
#
# ReadMe: 
# 
# Uses .tiff files from https://www.esa-landcover-cci.org/, updates it with latest glacier outline data and produces raster file that can 
# then be converted into crop factors
# Uses additional input for bare rock surfaces if available

# Output:
# ascii file  - needs to be converted into .map files by 
# asc2map --clone clone.map -S -a -m -9999 ---.asc ---.map
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

##########################
# SPECIFY FILENAMES AND DESIRED PROJECTION
RGI_filename <- '15_rgi60_SouthAsiaEast.shp' 

# Define Projection
projec_lambers <- "+proj=laea +lat_0=28 +lon_0=95 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
projec_utm <- '+proj=utm +zone=45N +datum=WGS84'

projec <- projec_utm

# PATHS FOR GLACIER RAW DATA (adapt and make sure all data is available)
path_maps <- 'F:\\PhD\\Research\\SPHY\\SPHYLangtang\\Code\\SPHY-master\\input'             # Folder of SPHY input data
path_output <- 'F:\\PhD\\Research\\SPHY\\SPHYLangtang\\figures'             # Folder for all figures
path_ccidata <- 'F:\\PhD\\Research\\SPHY\\data\\landuse\\cci-landcover-2015'
filename <- 'ESACCI-LC-L4-LCCS-Map-300m-P1Y-2015-v2.0.7.tif'                # CCI file name
path_RGI <- 'F:\\PhD\\GeoSpatialData\\RGI60_Asia' 

# if available, provide .shp for catchment outline for visualisation, set 'cSHP' to 1
cSHP <- 1
file_catchmentoutline <- 'F:\\PhD\\Research\\SPHY\\CatchmentMapping\\OutlineLangtang\\catch_proj_1.shp' # Catchment Outline (shp file)

# if available, provide raster file with data where no soil available, set 'bare' to 1
bare <- 1
file_barerock <- 'F:\\PhD\\Research\\SPHY\\data\\barerock_langtang100m.tif'
##########################

##########################
# Load SPHY data
domain <- raster(paste(path_maps,'\\clone.map',sep=''))   # Load the domain
projection(domain) <- projec

dem <- raster(paste(path_maps,'\\dem.map',sep=''))        # Load the DEM
projection(dem) <- projec

domain_deg <- projectRaster(domain,crs = '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')

modID_raster <- domain
modID_raster[seq(1,dim(domain)[1]*dim(domain)[2],1)] <- seq(1,dim(domain)[1]*dim(domain)[2],1)

if(bare ==1){
  bareRaster <- raster(file_barerock)
  projection(bareRaster) <- projec
  bareRaster <- resample(bareRaster,domain)
}
##########################
# Load RGI data
ogrInfo(path_RGI&'\\'&RGI_filename)
RGI60_15<-readOGR(dsn=path_RGI&'\\'&RGI_filename)
RGI60_15 <- spTransform(RGI60_15, projec)

# Restrict the datasets to the domain
sub_15 <- subset(RGI60_15, CenLon >= extent(domain_deg)[1] & CenLon <= extent(domain_deg)[2] & CenLat >= extent(domain_deg)[3] & CenLat <= extent(domain_deg)[4])

# Load Catchment outline data if available
if(cSHP == 1) {
ogrInfo(file_catchmentoutline)
catch<-readOGR(dsn=file_catchmentoutline)
catch <- spTransform(catch, projec)
}
##########################
# Load global CCI data

cci_global <- raster(paste(path_ccidata,'\\',filename,sep=''))
cci_local <- crop(cci_global,domain_deg)

cci_local <- projectRaster(cci_local,crs=projec , method = "ngb")
cci_local <- resample(cci_local,domain)

r_glaciermask <- rasterize(sub_15, domain,mask=F) # Individual glacier numbers (only pixels that lie on the polygon)
r_glaciermask[is.na(r_glaciermask)] <- 0

cci_local <- resample(cci_local,domain)

cci_local[which(cci_local[] < 10)] <- 0
cci_local[which(cci_local[] >=10 & cci_local[]<=40)] <- 1   # cropland
cci_local[which(cci_local[] >40 & cci_local[]<120)] <- 2    # trees
cci_local[which(cci_local[] >=120 & cci_local[]<130)] <- 3  # shrubland
cci_local[which(cci_local[] >=130 & cci_local[]<140)] <- 4  # grassland
cci_local[which(cci_local[] >=140 & cci_local[]<150)] <- 5  # lichen/moss
cci_local[which(cci_local[] >=150 & cci_local[]<160)] <- 6  # sparse vegetation
cci_local[which(cci_local[] >=160 & cci_local[]<190)] <- 7  # flooded land
cci_local[which(cci_local[] >=190 & cci_local[]<200 )] <- 8 # urban land
cci_local[which(cci_local[] >=200  & cci_local[]<210)] <- 9 # bare land
cci_local[which(cci_local[] >=210 & cci_local[]<220)] <- 10 # water body
cci_local[which(cci_local[] >=220)] <- NA                   # old snow/ice (superimposed by actual RGI)
if(bare == 1){
cci_local[bareRaster>0] <- 9 # Convert all pixels identified as bare rock as bare land
}
cci_local[r_glaciermask>0] <- 11 # Convert all pixels covered by glacier mask to 220 (ID for ice)

cci_local <- cci_local + 1

cci_local_fac<-as.factor(cci_local)
tar<-levels(cci_local_fac)[[1]]
tar[["landcover"]]<-c("no data / 0","cropland / 1", "trees / 1","shrub / 1.05","grass / 0.9 ","lichen/moss / 1","sparse veg / 1","flooded / 1","urban / 0","bare / 0","water / 1.05","ice / 0")
levels(cci_local_fac)<-tar

writeRaster(cci_local, file.path(path_maps, "luse.asc"), format="ascii",overwrite=T)

##########################
# Visualize domain data

  png(file=path_output&'\\landuse.png', res = 160,width=dim(domain)[2]*2,height=dim(domain)[1]*2)
  par(mar=c(4,5,2,2),cex.lab=1,cex.axis=1)
  #par(mfrow=c(1))
  layout(matrix(c(1), nrow = 1, ncol = 1, byrow = FALSE))
  if(cSHP == 1) {
    levelplot(cci_local_fac, col.regions=brewer.pal(12,"Set3"),xlab="Easting [m]",ylab="Northing [m]") + layer(sp.polygons(catch))
  } else{
    levelplot(cci_local_fac, col.regions=brewer.pal(12,"Set3"),xlab="Easting [m]",ylab="Northing [m]")
  }
  dev.off()
##########################  