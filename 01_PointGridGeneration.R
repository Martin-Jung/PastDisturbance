## Package loading ##
library(dplyr)
library(reshape2)
library(stringr)
library(lubridate)
library(tidyr)
# Over
#library(roquefort)
library(ggplot2)
# Spatial files
library(sp)
library(gstat)
library(rgdal)
library(rgeos)
# Custom packages
source("000_HelperFunction.R")

#### Prepare PREDICTS site data ####
# Load
# Limit estimates to only those in published version
#r <- readRDS("../../Data/diversity-2016-02-03-03-37-46.rds") 
#r <- r[which(r$Source_ID %in% unique(r.pub$Source_ID)),];rm(r.pub)

## THIS POINTS TO THE PREDICTS FULL PUBLIC DATABSE
# To be downloaded!
r <- readRDS("../../../Data/PREDICTS_v1/database.rds")
# Convert all SSBS names to correct encoding
r$SSBS <- iconv(r$SSBS,"UTF-8","latin1")
#r <- DropInvalidMetricsAndMethods(r)
r <- CorrectSamplingEffort(r)

sites <- SiteMetrics(diversity=r,
                     extra.cols=c("SSB","SSBS","Longitude","Latitude","Sample_start_earliest","Sample_end_latest","Sample_midpoint","Sample_date_resolution",
                                  "Ecoregion","Biome","Country","UN_subregion","Site_name","Order","Family",
                                  "Sampling_method","Study_common_taxon","Max_linear_extent","Coordinates_precision_metres"
                     ))

# Calculate pairwise dissimilarity
ss <- CompDissim2(r, metric="SorVeg",binary=T)
saveRDS(ss,"resSaves/Out_MatricesSor.rds")

# Calculate PIE
r$Is_abundance <- "Abundance" == r$Diversity_metric_type
site.abundance <- tapply(r$Is_abundance, r$SSS, unique)
sites$PIE <- SiteHurlbertsPie(r,site.abundance)
rm(r)
# ---- # 
# Set up
sites$startyear <- year(ymd(sites$Sample_start_earliest))
sites$endyear <- year(ymd(sites$Sample_end_latest))

# Encode taxonomic group category
sites$TGrouping <- as.character(sites$Study_common_taxon)
sites$TGrouping[grep("Hymenoptera",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Insecta",x = sites$TGrouping,ignore.case = T)]<- "Invertebrates"
sites$TGrouping[grep("Chordata",x = sites$TGrouping,ignore.case = T)] <- "Other"
sites$TGrouping[grep("Animalia",x = sites$TGrouping,ignore.case = T)] <- "Other"
sites$TGrouping[grep("Formicidae",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Scarabaeidae",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Tracheophyta",x = sites$TGrouping,ignore.case = T)] <- "Plantae"
sites$TGrouping[grep("Strigiformes",x = sites$TGrouping,ignore.case = T)] <- "Aves"
sites$TGrouping[grep("Isoptera",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Coleoptera",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Anogeissus",x = sites$TGrouping,ignore.case = T)] <- "Plantae"
sites$TGrouping[grep("Poaceae",x = sites$TGrouping,ignore.case = T)] <- "Plantae"
sites$TGrouping[grep("Colubridae",x = sites$TGrouping,ignore.case = T)] <- "Reptilia"
sites$TGrouping[grep("Chiroptera",x = sites$TGrouping,ignore.case = T)] <- "Mammalia"
sites$TGrouping[grep("Ascomycota",x = sites$TGrouping,ignore.case = T)] <- "Fungi"
sites$TGrouping[grep("Lepidoptera",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Bryophyta",x = sites$TGrouping,ignore.case = T)] <- "Plantae"
sites$TGrouping[grep("Sarcoptiformes",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[which(str_length(sites$TGrouping)==0)] <- "Other"
sites$TGrouping[grep("Bombus",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Apidae",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Arthropoda",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Drosophilidae",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Colletes floralis",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Phasianidae",x = sites$TGrouping,ignore.case = T)] <- "Aves"
sites$TGrouping[grep("Lophophorus impejanus",x = sites$TGrouping,ignore.case = T)] <- "Aves"
sites$TGrouping[grep("Gastropoda",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Araneae",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Arachnida",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Clitellata",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Anura",x = sites$TGrouping,ignore.case = T)] <- "Amphibia"
sites$TGrouping[grep("Carabidae",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Hemiptera",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Isopoda",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Collembola",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Agaricomycetes",x = sites$TGrouping,ignore.case = T)] <- "Fungi"
sites$TGrouping[grep("Curculionidae",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Pongo pygmaeus",x = sites$TGrouping,ignore.case = T)] <- "Mammalia"
sites$TGrouping[grep("Squamata",x = sites$TGrouping,ignore.case = T)] <- "Reptilia"
sites$TGrouping[grep("Culicidae",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Phyllostomidae",x = sites$TGrouping,ignore.case = T)] <- "Mammalia"
sites$TGrouping[grep("Maerua subcordata",x = sites$TGrouping,ignore.case = T)] <- "Plantae"
sites$TGrouping[grep("Oryctolagus cuniculus",x = sites$TGrouping,ignore.case = T)] <- "Mammalia"
sites$TGrouping[grep("Arecaceae",x = sites$TGrouping,ignore.case = T)] <- "Plantae"
sites$TGrouping[grep("Pteropus tonganus",x = sites$TGrouping,ignore.case = T)] <- "Mammalia"
sites$TGrouping[grep("Nymphalidae",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Diptera",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Staphylinidae",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Opiliones",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Orthoptera",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Swietenia macrophylla",x = sites$TGrouping,ignore.case = T)] <- "Plantae"
sites$TGrouping[grep("Aenictus",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Dorylus",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Vespertilionidae",x = sites$TGrouping,ignore.case = T)] <- "Mammalia"
sites$TGrouping[grep("Primates",x = sites$TGrouping,ignore.case = T)] <- "Mammalia"
sites$TGrouping[grep("Panthera pardus",x = sites$TGrouping,ignore.case = T)] <- "Mammalia"
sites$TGrouping[grep("Odocoileus virginianus",x = sites$TGrouping,ignore.case = T)] <- "Mammalia"
sites$TGrouping[grep("Cephalophus",x = sites$TGrouping,ignore.case = T)] <- "Mammalia"
sites$TGrouping[grep("Geometridae",x = sites$TGrouping,ignore.case = T)] <- "Mammalia"
sites$TGrouping[grep("Rodentia",x = sites$TGrouping,ignore.case = T)] <- "Mammalia"
sites$TGrouping[grep("Magnoliopsida",x = sites$TGrouping,ignore.case = T)] <- "Plantae"
sites$TGrouping[grep("Sciomyzidae",x = sites$TGrouping,ignore.case = T)] <- "Plantae"
sites$TGrouping[grep("Liolaemus",x = sites$TGrouping,ignore.case = T)] <- "Reptilia"
sites$TGrouping[grep("Dolichopus",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Muridae",x = sites$TGrouping,ignore.case = T)] <- "Mammalia"
sites$TGrouping[grep("Soricidae",x = sites$TGrouping,ignore.case = T)] <- "Mammalia"
sites$TGrouping[grep("Lumbricidae",x = sites$TGrouping,ignore.case = T)] <- "Mammalia"
sites$TGrouping[grep("Lecanoromycetes",x = sites$TGrouping,ignore.case = T)] <- "Fungi"
sites$TGrouping[grep("Clethrionomys gapperi",x = sites$TGrouping,ignore.case = T)] <- "Mammalia"
sites$TGrouping[grep("Passeriformes",x = sites$TGrouping,ignore.case = T)] <- "Aves"
sites$TGrouping[grep("Dipteryx oleifera",x = sites$TGrouping,ignore.case = T)] <- "Plantae"
sites$TGrouping[grep("Nematoda",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Diprotodontia",x = sites$TGrouping,ignore.case = T)] <- "Mammalia"
sites$TGrouping[grep("Glomeromycetes",x = sites$TGrouping,ignore.case = T)] <- "Fungi"
sites$TGrouping[grep("Strabomantidae",x = sites$TGrouping,ignore.case = T)] <- "Amphibia"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "CC1_2008__Schon 1")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "CC1_2008__Schon 2")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "DG1_2012__Ge 1")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "DL1_2009__Woinarski 2")] <- "Reptilia"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "DL1_2012__Dominguez 1")] <- "Plantae"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "HP1_2010__Bicknell 1")] <- "Aves"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "HP1_2010__Bicknell 2")] <- "Aves"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "MJ1_2009__Lehouck 2")] <- "Sturnidae"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "SC1_2014__Kurz 1")] <- "Amphibia"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "SE2_2010__Gaigher 1")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "SE2_2012__Craig 1")] <- "Amphibia"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "SE2_2014a_Craig 1")] <- "Amphibia"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "SE2_2014b_Craig 1")] <- "Amphibia"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "SH1_2011__Todd 1")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "SH1_2013__Peri 1")] <- "Plantae"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "SH1_2014__Walker 3")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "VB1_2012__Carpenter 1")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "VB1_2012__Carpenter 2")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "VB1_2012__Carpenter 6")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "VB1_2015__Mumme 1")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "DI1_2012__Muchane 1")] <- "Invertebrates"
sites$TGrouping[grep("Sturnidae",x = sites$TGrouping,ignore.case = T)] <- "Aves"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "AR1_2008__Basset 1")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "CC1_2005__Barratt 1")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "CC1_2005__Barratt 2")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "YP1_2012__Sung 1")] <- "Reptilia"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "CC1_2006__Norton 1")] <- "Plantae"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "VK1_2007__StLaurent 1")] <- "Mammalia"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "VB1_2013__Burton 3")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "VB1_2013__Burton 2")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "VB1_2012a_Carpenter 1")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "VB1_2009__Boutin 3")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "VB1_2009__Boutin 2")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "VB1_2009__Boutin 1")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "VB1_2012__LeightonGoodall 1")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "VB1_2008__Smith 1")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "VB1_2008a_Smith 1")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "VB1_2006__Smith 2")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "VB1_2006__Smith 1")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "VB1_2005__Eggleton 2")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "VB1_2005__Eggleton 1")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "SE2_2010__McCarthy 1")] <- "Mammalia"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "SC1_2009__Christensen 1")] <- "Plantae"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "TN1_2008__Ngai 1")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "TN1_2007__Gardner 1")] <- "Amphibia"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "SC1_2006__UrbinaCardona 1")] <- "Reptilia"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "SC1_2005__Richardson 1")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "MJ1_2009__Lehouck 1")] <- "Aves"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "CC1_2010__Schon 1")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "CC1_2010__Schon 2")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "SE2_2009__Craig 1")] <- "Amphibia"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "SE1_2012__Lopez 1")] <- "Fungi"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "MH1_2010__CATIE 5")] <- "Reptilia"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "MH1_2010__CATIE 4")] <- "Reptilia"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "MH1_2010__CATIE 3")] <- "Amphibia"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "MH1_2010__CATIE 2")] <- "Reptilia"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "MH1_2010__CATIE 1")] <- "Amphibia"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "MG1_2011__Schon 1")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "MG1_2008__Buscardo 1")] <- "Plantae"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "LK1_2010__Endo 1")] <- "Mammalia"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "LK1_2009__Hayward 1")] <- "Mammalia"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "JD1_2004__Alcala 1")] <- "Reptilia"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "HZ1_2012__Kutt 1")] <- "Amphibia"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "HW1_2011__Robinson 1")] <- "Fungi"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "GP1_2009__Vasconcelos 2")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "HB1_2009__Parry 1")] <- "Mammalia"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "GP1_2007__Kutt 1")] <- "Reptilia"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "DL1_2013__deThoisy 1")] <- "Mammalia"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "DL1_2008__MacipRios 1")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "CC1_2010__Schon 2")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "DB1_2010__Garden 1")] <- "Mammalia"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "DI1_2008__Noeske 1")] <- "Plantae"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "DL1_2008__MacipRios 1")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "AD1_2002__Vazquez 1")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "HB1_2002__Woinarski 1")] <- "Aves"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "SH1_2002__Bonham 1")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "HZ1_2006__Peres 1")] <- "Mammalia"

# Check
stopifnot( length(sites$SS[which(sites$TGrouping=="Other")]) == 0 )
# # Check number of studies

length(which(is.na(sites$Max_linear_extent)) ) / nrow(sites)

d <- sites
temp <- d %>% dplyr::select(Max_linear_extent,Sampling_method,Study_common_taxon) %>% 
  dplyr::group_by(Study_common_taxon,Sampling_method) %>% 
  dplyr::summarise(MLE_avg = mean(Max_linear_extent,na.rm=T))
d <- dplyr::left_join(d,temp) # Join back
# Insert where empty
d$Max_linear_extent <- ifelse(is.na(d$Max_linear_extent),d$MLE_avg,d$Max_linear_extent)
d$MLE_avg <- NULL # Kickout the previously calc. average

# Then check again this time with Higher Taxa at the remaining empty values
temp <- d %>% dplyr::select(Max_linear_extent,Sampling_method,TGrouping) %>% 
  group_by(TGrouping,Sampling_method) %>% 
  dplyr::summarise(MLE_avg = mean(Max_linear_extent,na.rm=T))
d <- left_join(d,temp)
d$Max_linear_extent <- ifelse(is.na(d$Max_linear_extent),d$MLE_avg,d$Max_linear_extent)
d$MLE_avg <- NULL # Kickout the previously calc. average

# Last try
temp <- d %>% dplyr::select(Max_linear_extent,Sampling_method,Study_common_taxon) %>% 
  group_by(Study_common_taxon) %>% 
  dplyr::summarise(MLE_avg = mean(Max_linear_extent,na.rm=T))
d <- left_join(d,temp)
d$Max_linear_extent <- ifelse(is.na(d$Max_linear_extent),d$MLE_avg,d$Max_linear_extent)
d$MLE_avg <- NULL # Kickout the previously calc. average
rm(temp)
# How many left
print(paste("How many remaining without MaxLinExtent:", length(which(is.na(d$Max_linear_extent)))))

# # # # All MLI filledin # # # #
cat( ((length(which(d$Max_linear_extent<3000)) / nrow(d) ) ) *100,"% smaller than 3000m" )
sites <- d
rm(d)
myLog("Preperation of PREDICTS done!")
myLog("-----------------------------")
saveRDS(sites,"PREDICTS_allsites.rds")

#### Generating site-based point grids ####
# This 
source('https://raw.githubusercontent.com/Martin-Jung/Icarus/master/R/miscellaneous/latlong2UTMzone.R')
library(rgeos)
sites <- readRDS("PREDICTS_allsites.rds")
myLog("Creating site based point grid...")
sp = subset(sites,select = c("SS","SSBS","Longitude","Latitude","Max_linear_extent"))

# First buffer radially
sp.long <- list()
for(site in unique(sp$SSBS)){
  print(site)
  sub.site <- subset(sp,SSBS == site)
  # Get UTM zone projection
  zone = CRS(paste0("+proj=utm +zone=",latlong2UTMzone(lon = sub.site$Longitude,lat = sub.site$Latitude )," +datum=WGS84 +units=m +no_defs"))
  
  # Make spatial file
  coordinates(sub.site) <- ~Longitude+Latitude
  proj4string(sub.site) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  
  # Transform to a metric meter based projection
  sub.site <- spTransform(sub.site,CRSobj = zone )
  
  buf <- gBuffer(sub.site,width = sub.site$Max_linear_extent / 2,quadsegs = 50)  
  # Backtransform
  buf <- spTransform(buf,CRSobj = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") )
  buf@polygons[[1]]@ID <- site
  sp.long[[site]] <- buf
}
stopifnot(nrow(sp) == length(sp.long))
joined <- do.call(rbind,sp.long)
jdata = SpatialPolygonsDataFrame(Sr = joined, data = as.data.frame(sp),match.ID = FALSE)
writeOGR(jdata,"PREDICTS_RadialBuffers.kml","PREDICTS_RadialBuffers",driver = "KML",overwrite_layer = T)
# Also write jdata as shapefile while we at it.
writeOGR(jdata,"PREDICTS_RadialBuffers.shp","PREDICTS_RadialBuffers",driver = "ESRI Shapefile",overwrite_layer = T)
rm(joined,jdata,sp.long);gc(verbose = T)

## According to http://landsat.usgs.gov/band_designations_landsat_satellites.php ##
## with few exceptions all Landsat 4-8 bands are at 30m resolution ## 
cs = 30 # Resolution
co = 101 # Cut off
sp.long = list()

# Make a site based point identifier #
# The encoding will be increasing string with p%i where %i = i++
sp$SSBSp <- NA
# If extent is smaller than demandend extent, take given coordinate only
#sp$SSBSp[which(sp$Max_linear_extent <= cs )] <- paste0("p",1)

# The offset at both ends of a rectangular buffer increases by two at both ends
# If multiplied with the cellsize that gives me 
x = seq(from=1,to = co,by = 2)
# Remove studies greater than x times 30 of the cutoff
sp <- subset(sp,Max_linear_extent < (max(x) * 30))
# For each size
for(ext in x ){
  # Get all those sites that fall within the given zone
  sub <- subset(sp,Max_linear_extent <= (ext*cs) & Max_linear_extent > ((ext-2)*cs) )
  # Now for each site SSBS
  for(site in unique(sub$SSBS)){
    myLog("Processing ",site, " - ",ext)
    sub.site <- subset(sub,SSBS == site)
    # Get UTM zone projection
    zone = CRS(paste0("+proj=utm +zone=",latlong2UTMzone(lon = sub.site$Longitude,lat = sub.site$Latitude )," +datum=WGS84 +units=m +no_defs"))
    
    # Make spatial file
    coordinates(sub.site) <- ~Longitude+Latitude
    proj4string(sub.site) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
    
    # Transform to a metric meter based projection
    sub.site <- spTransform(sub.site,CRSobj = zone )
    
    #define the plot boundaries based upon the plot radius. 
    #NOTE: this assumes that plots are oriented North and are not rotated. 
    yPlus <- sub.site@coords[2]+ (ext*cs)/2
    xPlus <- sub.site@coords[1]+ (ext*cs)/2
    yMinus <- sub.site@coords[2]- (ext*cs)/2
    xMinus <- sub.site@coords[1]- (ext*cs)/2
    
    # Calculate a square
    square = cbind(xMinus,yPlus, xPlus,yPlus, xPlus,yMinus, xMinus,yMinus,xMinus,yPlus,xMinus,yPlus) 
    
    # Create spatial polygons out of square
    polys <- SpatialPolygons(mapply(function(poly, id) {
      xy <- matrix(poly, ncol=2, byrow=TRUE)
      Polygons(list(Polygon(xy)), ID=sub.site$SSBS)
    }, split(square, row(square)), sub.site$SSBS),proj4string=zone) 
    
      # Get data to polygon and insert into list
      o <- spTransform(polys,CRSobj = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
      od <- as.data.frame(sub.site) %>% mutate(zone=latlong2UTMzone(lon = coordinates(o)[1],lat = coordinates(o)[2] ) ) %>% 
        rename(x = Longitude, y = Latitude) %>% 
        # Correct coordinates to long,lat
        mutate(Longitude = coordinates(o)[1], Latitude = coordinates(o)[2])
      row.names(od) <- od$SSBS
      sp.long[[sub.site$SSBS]] <- SpatialPolygonsDataFrame(o,data = od)
      rm(o,od,polys)
  }
}
# Control check
stopifnot(length(sp.long)==nrow(sp))
print("Done!");rm(x)

# Now join back the multiple polygons
joined = SpatialPolygons(lapply(sp.long, function(x){x@polygons[[1]]}),proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
jdata = SpatialPolygonsDataFrame(Sr=joined, data=do.call(rbind,(lapply(sp.long,as.data.frame))),match.ID = FALSE)

# Finally write result to KML
library(rgdal)
writeOGR(jdata,"PREDICTS_ExtentPolygons.kml","PREDICTS_ExtentPolygons",driver = "KML",overwrite_layer = T)
# Also write jdata as shapefile while we at it.
writeOGR(jdata,"PREDICTS_ExtentPolygons.shp","PREDICTS_ExtentPolygons",driver = "ESRI Shapefile",overwrite_layer = T)
print("Finished exporting! Full Done!")

### Export again ###
library(rgdal)
sp <- readOGR(".","PREDICTS_ExtentPolygons")

library(plotKML)
plotKML(sp, filename = " /PREDICTS_ExtentPolygons.kml")

# FROM ABOVE
# Save center coordiates rather than extents
library(plotKML)
library(rgdal)
coordinates(sp) <- ~Longitude+Latitude
proj4string(sp) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
writeOGR(sp," PREDICTS_CenterPoints.kml","PREDICTS_CenterPoints",driver = "KML",overwrite_layer = T)

plotKML(sp, filename = " PREDICTS_CenterPoints.kml")