# LaikipiaSpatialAnalysis_Test
Analysis of vegetation index change images exported from Google Earth Engine.

# Load required packages
library(raster) # rasters/grids
library(rgdal) # readORG?
library(sp) # spatial data classes
library(rasterVis) #raster visualisation
#library(maptools)
library(rgeos) #for gIntersection
#library(dismo) #google maps
library(units) #set units
library(smoothr) #for removing small polygons
library(lwgeom) #units packge doesn't seem to work without it?
#library(foreach) #faster for loops?
library(ggplot2)
library(plyr) #for ddply
library(splitstackshape) #for stratified sampling
#library(olsrr) #ols_test_normality()

setwd("~/Desktop/LaikipiaSpatialAnalysisTest") #set working directory
#options(stringsAsFactors = FALSE) #prevent strings from being set as factors

#KenyaMap <- gmap("Kenya", type = "satellite", exp = 1)
#plot(KenyaMap)

# Define spatial projection
#crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84") #geographical, datum WGS84
#projection(x) <- crs.geo # 'x' is the raster

# Load Laikipia vegetation cover map (Taiti 1992)
vegCoverMap <- readOGR(dsn=".", layer="vegetation_taiti_1992_EPSG4326")

# Load Laikipia land management unit boundaries map with customised 'MGMT' column
boundaries <- readOGR(dsn=".", layer="Property_boundary_2015_EPSG4326_valid")

# Load 1984-2018 VI change slope image
viSlope <- raster('./NDVI_pcChange_2000-2018_Laikipia_syntheticLandsat-7_0pcWaterMask.tif')

vegCoverTypes <- c("Acacia drep. Bushland",
                   #"Acacia Seyal Bushland", #absent in WO
                   "Arid Zone Acacia Bushland",
                   #"Dwarf Bush Grassland", #absent in WO
                   "Grassland",
                   #"Leafy Bushland", #absent in 
                   #"Leafy Bushland and Thicket", #absent in
                   "Open Acacia brevispica thicket")
boundaries.df <- as(boundaries, "data.frame") #convert to dataframe for unique()
mgmtTypes <- na.omit(unique(boundaries.df$"MGMT")) #remove NAs

final <- data.frame() #empty dataframe to house results of for loop
# For loop to create all vegetation cover by management type combinations
for (mgmtType in mgmtTypes) {
  for (vegCoverType in vegCoverTypes) {
    vegSubset <- subset(vegCoverMap, DESCRIPTIO==vegCoverType)
    mgmtSubset <- subset(boundaries, MGMT==mgmtType)
    mgmtDissolved <- gUnaryUnion(mgmtSubset) #dissolve internal lines within polygons
    intersection <- gIntersection(mgmtDissolved, vegSubset)
    #drop polygons smaller than 1km2
    smoothed <- if (!is.null(intersection)) drop_crumbs(intersection, threshold = set_units(1, km^2)) else NULL
    extracted <- if (!is.null(smoothed)) extract(viSlope, smoothed) else NA #extract slope values
    tab <- table(extracted, exclude = NULL) #keep NAs to ensure consistent row numbers for cbind
    df <- as.data.frame(tab) #prepare for cbind()
    result <- cbind(vegCoverType=vegCoverType, mgmtType=mgmtType, df)
    final <- rbind(final, result)
  }
}

# Convert from factor to numeric
final$extracted.num <- as.numeric(levels(final$extracted))[final$extracted] 
# Rename for clarity in plot
final$vegCoverType <- revalue(final$vegCoverType, 
                              c("Acacia drep. Bushland" = "AdB",
                                "Acacia Seyal Bushland" = "AsB",
                                "Arid Zone Acacia Bushland" = "AZAB",
                                "Dwarf Bush Grassland" = "DBG",
                                "Grassland" = "G",
                                "Leafy Bushland" = "LB",
                                "Leafy Bushland and Thicket" = "LBT",
                                "Open Acacia brevispica thicket" = "OAbT"))


globalM <- median(final$extracted.num, na.rm=TRUE) #compute mean/median for whole dataset
final$anomaly <- final$extracted.num-globalM #compute anomaly
final$anomalyProp <- (final$extracted.num-globalM)/globalM #convert to % anomaly

# Subsample 10% using stratified random sampling to account for spatial autocorrelation
sampled <- stratified(final, c("mgmtType", "vegCoverType"), 0.1) #0.1=10%

ggplot(data = sampled, aes(x=vegCoverType, y=extracted.num, fill=vegCoverType))+
  geom_boxplot(outlier.size = 0.1)+
  facet_grid(~mgmtType)+
  geom_hline(aes(yintercept=globalM), linetype="dashed")+
  ylab("2000-2018 NDVI change of 10% sample pixels (%)")+
  xlab(NULL)+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"), name="vegetation\ncover type")+ #colour blind-friendly palette
  theme(axis.text.x = element_text(angle = 90))

# Compute summary statistics
#stats <- ddply(.data = final,
               .variables = c("mgmtType", "vegCoverType"),
               summarize,
               M = mean(extracted.num, na.rm=T), SD = sd(extracted.num, na.rm=T), SEM = SD/(sqrt(length(extracted.num))))
 
#stats$normAnomaly <- (stats$M-globalM)/stats$SD #compute normalized anomaly 

# Plot VI change by management & vegetation type
ggplot(data=stats, aes(x=vegCoverType, y=normAnomaly, fill=mgmtType))+
  geom_bar(stat="identity", position=position_dodge())+
  #geom_errorbar(aes(ymin=norm-SEM, ymax=norm+SEM), width=0.2, position = position_dodge(0.9))+
  ylab("Normalised mean 1984-2018\nVI change slope anomaly")+
  xlab("Vegetation cover type")+
  scale_fill_manual(values=c("#009E73", "#e79f00", "#CC79A7", "#0072B2"), name="managemet\ntype") #colour blind-friendly palette

# Test for statistical differences between vegCover-mgmt combinations
aov.fit <- aov(extracted.num~mgmtType*vegCoverType, data = sampled) #'final'/'sampled' log-transform y?
summary(aov.fit)
print(model.tables(aov.fit, "means"), digits = 3) # report means & # of pixels per cell
TukeyHSD(aov.fit)
plot(aov.fit)
hist(resid(aov.fit))

# Export plot
#png("plot1.png")
#plot1 <- ggplot(data = final, aes(x=vegCoverType, y=extracted, fill=mgmtType))+
#geom_boxplot()
#print(plot1)
#dev.off()
