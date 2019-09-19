################################################################################
# SPHY_snowmap - create snowmaps for validation of the SPHY-snow model
# 
# SPHY_snowmap.R
#
# ReadMe:
# 
# Uses .tiff files from MODIS / Landsat or / Sentinel snow cover maps
# To be extended for snow depth and SWE

# Output:
# - Snowcover: ascii files with 0 = No Snow, 1 = Snow and NA = NO DATA/CLOUDS/UNKNOWN - needs to be converted into .map files by 
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
p_load(rgdal,rgeos,maptools,raster)

########################## 
# Define Paths and Projection
##########################
projec_lambers <- "+proj=laea +lat_0=28 +lon_0=95 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
projec_utm <- '+proj=utm +zone=45N +datum=WGS84'

projec <- projec_utm

# PATHS  (adapt and make sure all data is available)
path_maps <- 'F:\\PhD\\Research\\SPHY\\SPHYLangtang\\Code\\SPHY-master\\input'        # Folder of SPHY input data (used as output folder)
path_output <- 'F:\\PhD\\Research\\SPHY\\SPHYLangtang\\figures'                       # Folder for all figures
path_snow_output <- 'F:\\PhD\\Research\\SPHY\\SPHYLangtang\\validation\\snow'
path_ccidata <- 'F:\\PhD\\Research\\SPHY\\data\\landuse\\cci-landcover-2015'          # Folder for raw landcover data
filename <- 'ESACCI-LC-L4-LCCS-Map-300m-P1Y-2015-v2.0.7.tif'                          # CCI file name
path_RGI <- 'F:\\PhD\\GeoSpatialData\\RGI60_Asia'                                     # Folder for RGI glacier outlines
RGI_filename <- '15_rgi60_SouthAsiaEast.shp'                                          # RGI filename

path_snowdata <- 'F:\\PhD\\Research\\SPHY\\data\\Trisuli\\snow\\MODIS_maximum_snow_extent'           # folder with snow data, needs to be a list of folders that have as a name the numerical time and contain a tiff

# Specify the satelite data used (currently supported Modis 8day composite for snow cover extent)
satellite <- 'modis8day'
##########################
# Load SPHY data
##########################
domain <- raster(paste(path_maps,'\\clone.map',sep=''))   # Load the domain
projection(domain) <- projec

dem <- raster(paste(path_maps,'\\dem.map',sep=''))        # Load the DEM
projection(dem) <- projec

domain_deg <- projectRaster(domain,crs = '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')

modID_raster <- domain
modID_raster[seq(1,dim(domain)[1]*dim(domain)[2],1)] <- seq(1,dim(domain)[1]*dim(domain)[2],1)

##########################
# Load Snow Data for each time step
##########################

listFil <- list.files(path = path_snowdata, pattern = NULL, all.files = T,
           full.names = T, recursive = F,
           ignore.case = FALSE, include.dirs = F)

listDate <- list.files(path = path_snowdata, pattern = NULL, all.files = F,
                      full.names = F, recursive = F,
                      ignore.case = FALSE, include.dirs = F)
           
for (i in 1:length(listDate)){

  date <- as.numeric(listDate[i])
  getTiff <- list.files(path = paste(path_snowdata,'\\',listDate[i],sep=''), pattern = "\\.tif$",full.names = T)
  snowmap <- raster(getTiff)        # Load the extent map
  projection(snowmap) <- projec
  snowmap_resampled <- resample(snowmap,dem)
  snowmap_resampled <- mask(snowmap_resampled,dem)
  
  switch(satellite,
         'modis8day' = {
           snowmap_resampled2 <- snowmap_resampled
           snowmap_resampled[snowmap_resampled>=190] <- 1 # Snow
           snowmap_resampled[snowmap_resampled>=22 & snowmap_resampled<=27] <- 0  # no snow
           snowmap_resampled[snowmap_resampled>0 & snowmap_resampled<1] <- NA
           snowmap_resampled[snowmap_resampled>1] <- NA

           # Create new folder
           ifelse(!dir.exists(file.path(path_output,paste('\\',listDate[i],sep=''))), dir.create(file.path(path_snow_output,paste('\\',listDate[i],sep=''))), FALSE)
           
           # Save all rasters (GTiff is an option but ascii files are needed for PCRaster/SPHY)
           writeRaster(snowmap_resampled, file.path(file.path(path_snow_output,paste('\\',listDate[i],sep='')), "snowcover.asc"), format="ascii",overwrite=TRUE)
         })
}

##########################
# Visualize domain data
##########################
pal <- colorRampPalette(c("red","blue"))
  png(file=path_output&'\\Domain_SnowCover.png', res = 160,width=dim(dem)[2]*2,height=dim(dem)[1]*2)
  par(mar=c(4,5,2,2),cex.lab=1.5,cex.axis=1.5)
  #par(mfrow=c(1))
  layout(matrix(c(1), nrow = 1, ncol = 1, byrow = FALSE))
  plot(dem,xlab="Easting [m]",ylab="Northing [m]",legend.args=list(text='Elevation [m asl]', side=2, font=1, line=0.3, cex=1.5))
  plot(snowmap_resampled,add=T,col=pal(3),legend=FALSE)
  legend('bottomright',bty="n",c('no snow','snow'),pch = 1,col = c('red','blue'))
  dev.off()
##########################  