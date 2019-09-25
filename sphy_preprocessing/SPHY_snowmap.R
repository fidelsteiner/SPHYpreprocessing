################################################################################
# SPHY_snowmap - create snowmaps for validation of the SPHY-snow model
# 
# SPHY_snowmap.R
#
# ReadMe:
# 
# (1) Uses files from MODIS snow cover maps
#     (a) TIFF Snowmaps from Modis
#     (b) HDF Snowmaps (direct output from datadownload; python code provided separately)
# To be extended for snow depth and SWE
# (2) Generates output_snow.asc file with locations where SnowDepth and SWE need to be reported 

# Output:
# - Snowcover: ascii or tiff files with 0 = No Snow, 1 = Snow and NA = NO DATA/CLOUDS/UNKNOWN - ascii needs to be converted into .map files by 
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
library(ncdf4)
library(gdalUtils)

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

formatMODIS <- 'HDF'
path_snowdata <- 'F:\\PhD\\Research\\SPHY\\data\\langtang\\Snow\\MODIS10A2\\HDF'       # folder with snow data (for .tiff needs to be a list of folders that have as a name the numerical time and contain a tiff)

# individual locations for snow measurements - still hardwired
path_snowloc <- 'F:\\PhD\\Research\\SPHY\\CatchmentMapping'                           # Location of file with coordinates of snow measurements 
file_snowloc <- 'StationLocsMap.csv'                                                  # File with locations, including Lat/Lon column as 2nd and 3d column

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

switch(formatMODIS,
       'TIFF' = {
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
           snowmap_resampled[snowmap_resampled>=190] <- 1 # Snow
           snowmap_resampled[snowmap_resampled>=22 & snowmap_resampled<=27] <- 0  # no snow
           snowmap_resampled[snowmap_resampled>0 & snowmap_resampled<1] <- NA
           snowmap_resampled[snowmap_resampled>1] <- NA

           # Create new folder
           ifelse(!dir.exists(file.path(path_snow_output,paste('\\',listDate[i],sep=''))), dir.create(file.path(path_snow_output,paste('\\',listDate[i],sep=''))), FALSE)
           
           # Save all rasters (GTiff is an option but ascii files are needed for PCRaster/SPHY)
           writeRaster(snowmap_resampled, file.path(file.path(path_snow_output,paste('\\',listDate[i],sep='')), "snowcover.asc"), format="ascii",overwrite=TRUE)
         })
}
       },
'HDF' = {
  listFil <- list.files(path = path_snowdata, pattern = "\\.hdf$", all.files = F,
                        full.names = T, recursive = F,
                        ignore.case = FALSE, include.dirs = F)
  
  listDate <- list.files(path = path_snowdata, pattern = "\\.hdf$", all.files = F,
                         full.names = F, recursive = F,
                         ignore.case = FALSE, include.dirs = F)
  
  doyDate <- substr(listDate, 10, 16) # Date in YYYYDDD format
  
  for (i in 1:length(listDate)){
    
    Actdate <- as.numeric(as.Date(as.numeric(substr(doyDate[i],5,7)), origin = paste(substr(doyDate[i],1,4),"-01-01",sep='')))
    sds <- get_subdatasets(listFil[i])
    gdal_translate(sds[1], dst_dataset = "temp_snowmaxextent.tif")
    snowmap <- raster("temp_snowmaxextent.tif")        # Load the extent map
    snowmap_resampled <- projectRaster(snowmap,dem)
    snowmap_resampled <- mask(snowmap_resampled,dem)
    
    switch(satellite,
           'modis8day' = {
             snowmap_resampled[snowmap_resampled>250] <- NA # detector saturated
             snowmap_resampled[snowmap_resampled>=190] <- 1 # Snow
             snowmap_resampled[snowmap_resampled>=22 & snowmap_resampled<=27] <- 0  # no snow
             snowmap_resampled[snowmap_resampled>0 & snowmap_resampled<1] <- NA
             snowmap_resampled[snowmap_resampled>1] <- NA
             
             # Create new folder
             ifelse(!dir.exists(file.path(path_snow_output,'MODIS')), dir.create(file.path(path_snow_output,'MODIS')), FALSE)
             
             # Save all rasters (GTiff is an option but ascii files are needed for PCRaster/SPHY)
             writeRaster(snowmap_resampled, file.path(file.path(path_snow_output,'MODIS'), paste(Actdate,"max8dayMODIS.tiff",sep="")), format="GTiff",overwrite=TRUE)
           })
  }
}
)
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
# Get Locations for Snow Measurements and return as .asc File
##########################  
SnowLocs <- read.csv(path_snowloc&'\\'&file_snowloc)

SnowDepth <- SnowLocs[c(2:13),2:3]

valLocs <- as.data.frame(seq(1,12,1)) # individual ID for each location
colnames(valLocs) <- 'siteID'
spSnowMeas <- SpatialPointsDataFrame(as.data.frame(cbind(SnowDepth[,2],SnowDepth[,1])), valLocs, coords.nrs = numeric(0),proj4string = CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'),bbox = NULL)
spSnowMeas <- spTransform(spSnowMeas, projec)

snowLocRaster <- rasterize(spSnowMeas,dem,'siteID')
writeRaster(snowLocRaster, file.path(path_maps, "output_snow.asc"), format="ascii",overwrite=T)