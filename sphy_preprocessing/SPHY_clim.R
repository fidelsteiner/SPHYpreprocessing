################################################################################
# SPHY_discharge - prepare hydrology data needed for SPHY clibration/validation
# 
# SPHY_discharge.R
#
# ReadMe: 
# 
# (1) Creates outlet.asc file to define where discharge is computed

# Output:
# ascii file for landuse  - needs to be converted into .map files by 
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
##########################
# Define Projection
projec_lambers <- "+proj=laea +lat_0=28 +lon_0=95 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
projec_utm <- '+proj=utm +zone=45N +datum=WGS84'

projec <- projec_utm

# PATHS  (adapt and make sure all data is available)
path_maps <- 'F:\\PhD\\Research\\SPHY\\SPHYLangtang\\Code\\SPHY-master\\input'        # Folder of SPHY input data (used as output folder)
path_output <- 'F:\\PhD\\Research\\SPHY\\SPHYLangtang\\figures'                       # Folder for all figures
path_ccidata <- 'F:\\PhD\\Research\\SPHY\\data\\landuse\\cci-landcover-2015'          # Folder for raw landcover data
filename <- 'ESACCI-LC-L4-LCCS-Map-300m-P1Y-2015-v2.0.7.tif'                          # CCI file name
path_RGI <- 'F:\\PhD\\GeoSpatialData\\RGI60_Asia'                                     # Folder for RGI glacier outlines
RGI_filename <- '15_rgi60_SouthAsiaEast.shp'                                          # RGI filename
path_watershed <- 'F:\\PhD\\Research\\SPHY\\SPHYLangtang\\RawData\\WaterSheds'        # Location for all watersheds

outlet1 <- 'SyafruBesi.shp'
outlet2 <- 'Kyanjing_Meas.shp'


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
dem_deg <- projectRaster(dem,crs = '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')

lat_dem <- as.data.frame(coordinates(dem_deg)[!is.na(dem_deg[]),2])
lon_dem <- as.data.frame(coordinates(dem_deg)[!is.na(dem_deg[]),1])
rad_table <- cbind(seq(1,length(lat_dem[,1]),1),lat_dem,lon_dem)
colnames(rad_table) <- c('label','lat','lon')

write.table(rad_table,'F:\\PhD\\Research\\SPHY\\SPHYLangtang\\RawData\\radIDs.csv',sep = ",",row.names = F)