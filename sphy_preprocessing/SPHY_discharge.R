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
library(RColorBrewer)

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
path_debris <- 'F:\\PhD\\GeoSpatialData\\DCG_Scherler\\S2_2015-2017_NDSI\\S2_2015-2017_NDSI'  # folder debris data (Scherler 2018)
path_thick <- 'F:\\PhD\\GeoSpatialData\\IceThickness_Farinotti\\RGI60-15\\RGI60-15'           # folder with ice thickness data (Farinotti 2019)
RGI_debris_filename <- '15_rgi60_SouthAsiaEast_S2_DC_2015_2017_NDSI.shp'


path_watershed <- 'F:\\PhD\\Research\\SPHY\\SPHYLangtang\\RawData\\WaterSheds'        # Location for all watersheds

outlets <- c('SyafruBesi','Kyanjing_Meas','Langshisha_Meas','Lirung_Meas','Lama','LangtangVillage','LangtangGlacier','Shalbachum','Numthang','Yala','Kimoshung','Ganja','Ganja2')
watersheds <- c('ws_syafru.tif','ws_kyan.tif','ws_las.tif','ws_lir.tif','ws_lam.tif','ws_vil.tif','ws_lan.tif','ws_shal.tif','ws_num.tif','ws_yala.tif','ws_kimoshun.tif','ws_gan1.tif','ws_gan2.tif')
outletType <- c(1,2,2,2,3,3,4,4,5,6,6,6,6) # 1 - main outlet, 2 - measurements, 3 - ungauged catchments, 4 - ungauged subcatchments with debris covered glacier 5 - unguaged subcatchment with nearly no glaciers, 6 - unguaged catchment with clean ice cover


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
writeOGR(contourDEM, dsn = path_watershed&'\\contoursDomain.shp', layer='contourDEM',driver = "ESRI Shapefile")


domain_deg <- projectRaster(domain,crs = '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')
dem_deg <- projectRaster(dem,crs = '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')

dem_lat <- domain
dem_lat[seq(1,dim(domain)[1]*dim(domain)[2],1)] <- coordinates(domain_deg)[seq(1,dim(domain)[1]*dim(domain)[2],1),2]

modID_raster <- domain
modID_raster[seq(1,dim(domain)[1]*dim(domain)[2],1)] <- seq(1,dim(domain)[1]*dim(domain)[2],1)

##########################
# Load RGI and thickness data
##########################
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

if(bare ==1){
  bareRaster <- raster(file_barerock)
  projection(bareRaster) <- projec
  bareRaster <- resample(bareRaster,dem)
  bareRaster <- mask(bareRaster,dem)
}

##########################
# Load locations for outlets and catchment data
##########################
allOutlets <- SpatialPoints(data.frame(x = 0, y = 0))[-1,]
for(subc in 1:length(outlets)){
  ogrInfo(path_watershed&'\\'&outlets[subc]&'.shp')
  out1 <- readOGR(dsn=path_watershed&'\\'&outlets[subc])
  out1 <- spTransform(out1, projec)
  out1$Id <- subc
  allOutlets <- allOutlets + out1
}

outletRaster <- rasterize(allOutlets,dem,'Id') # Combine all outlet locations
outletRaster[is.na(outletRaster[])] <- 0

writeRaster(outletRaster, file.path(path_maps, "outlet.tif"), format="GTiff",overwrite=T)

##########################
# Get Properties of subcatchments
##########################

r_glaciermask_cover <- rasterize(sub_15, dem,getCover=TRUE) # relative cover of glacier over pixel
r_glaciermask_cover <- mask(r_glaciermask_cover,dem)
r_glaciermask_debris_cover <- rasterize(sub_15_debris, dem,getCover=TRUE) # relative cover of debris over pixel
r_glaciermask_debris_cover <- mask(r_glaciermask_debris_cover,dem)
bareRaster[r_glaciermask_cover>0] <- NA

for(subc in 1:length(outlets)){
op <- raster(path_watershed&'\\'&watersheds[subc]) + 1
projection(op) <- projec

elev_sub <- mask(dem,op)
elev_sub_his <- hist(elev_sub,breaks=seq(1000,7000,by=500))

outline_sub <- buffer(rasterToPolygons(elev_sub, dissolve=TRUE), width = 1, dissolve = T)
outline_sub <- as(outline_sub, "SpatialPolygonsDataFrame")

# Make new subfolder for subcatchment
newPath <- file.path(path_watershed&'\\'&outlets[subc])
dir.create(newPath)
writeOGR(outline_sub, dsn = newPath&'\\outline_'&outlets[subc]&'.shp', layer='outline_sub',driver = "ESRI Shapefile")
writeRaster(elev_sub, newPath&'\\'&outlets[subc]&'_Dem.tif', format="GTiff",overwrite=T)

glac_sub <- mask(r_glaciermask_cover,elev_sub)
glac_sub[glac_sub>0] <- 1
glac_sub[glac_sub<1] <- NA
glac_sub_his <- hist(glac_sub*elev_sub,breaks=seq(1000,7000,by=500))

deb_sub <- mask(r_glaciermask_debris_cover,elev_sub)
deb_sub[deb_sub>0] <- 1
deb_sub[deb_sub<1] <- NA
deb_sub_his <- hist(deb_sub*elev_sub,breaks=seq(1000,7000,by=500))

bare_sub <- mask(bareRaster,elev_sub)
bare_sub[bare_sub>0] <- 1
bare_sub[bare_sub<1] <- NA
bare_sub_his <- hist(bare_sub*elev_sub,breaks=seq(1000,7000,by=500))

elevHist <- elev_sub_his$counts/100 - glac_sub_his$counts/100 - bare_sub_his$counts/100 # Area per elevation band in km2 not covered by glaciers or rock
glacHist <- glac_sub_his$counts/100 - deb_sub_his$counts/100                            # Area only covered by clean ice
relGlac <- floor(sum(glacHist,na.rm=T) /sum(elev_sub_his$counts/100,na.rm=T) * 100)
debHist <- deb_sub_his$counts/100
relDeb <- floor(sum(debHist,na.rm=T) /sum(elev_sub_his$counts/100,na.rm=T) * 100)
bareHist <- bare_sub_his$counts/100
relBare <- floor(sum(bareHist,na.rm=T) /sum(elev_sub_his$counts/100,na.rm=T) * 100)
relFree <- floor(sum(elevHist,na.rm=T)/sum(elev_sub_his$counts/100,na.rm=T) * 100)

png(file=path_watershed&'\\ElevBins'&outlets[subc]&'.png', res = 300,width=1000,height=3000)
par(mar=c(4,5,2,2),cex.lab=1.5,cex.axis=1.5)
#par(mfrow=c(1))
xlimMax <- ceiling(max(rbind(elevHist+glacHist+debHist+bareHist)))
barplot(rbind(bareHist,debHist,glacHist,elevHist), xaxt="n",
        names.arg = elev_sub_his$mids,
        las = 1,
        xlab = expression('km' ^ 2),
        horiz = TRUE,
        border = "grey",
        xlim = c(0, xlimMax),
        xaxp = c(0, xlimMax, 2),col=c(brewer.pal(1, 'Greys')[c(3,2)],brewer.pal(1, 'BuGn')[1:2]))
if(xlimMax>=10){
abline(v=seq(0, xlimMax, by=floor(xlimMax/10)*10/5), col='white', lwd=2)
axis(side=1, at=seq(0, xlimMax, by=floor(xlimMax/10)*10/5), labels = seq(0, xlimMax, by=floor(xlimMax/10)*10/5))
}else{
  axis(side=1, at=seq(0, xlimMax, by=1), labels = seq(0, xlimMax, by=1))
}
grid(NULL,NULL)
legend('bottomright',bty="n",c(parse(text=paste("Area:~",sum(elev_sub_his$counts/100), "*km^2")),paste('bare rock ',relBare,'%'),paste('debris cover ',relDeb,'%'),paste('clean ice ',relGlac,'%'),paste('vegetated ',relFree,'%')),pch = c(NA,1,1,1,1),col = c(NA,brewer.pal(3, 'Greys')[c(3,2)],brewer.pal(3, 'BuGn')[1:2]))
dev.off()

}

##########################
# Get Discharge data
##########################
