# Package loading
par.ori <- par(no.readonly = T)
library(dplyr)
library(reshape2)
library(stringr)
library(lubridate)
library(tidyr)
library(jsonlite)
library(zoo)
library(xts)
library(bfast)
library(imputeTS)
library(sf)
library(progress)

source("000_HelperFunction.R")
myLog <- function(...) {
  cat(paste0("[AssembleTimeSeries] ", Sys.time(), " | ", ..., "\n"))
}

# The PREDICTS sites 
sites <- readRDS("PREDICTS_allsites.rds")

# This is in the input directory where the downloaded geojson files
# from Google Earth Engine have to be placed
input <- "LandsatPREDICTS/"

#### Assemble time series ####
# Get list of input files
ff <- sort( list.files(input,"*.json",full.names=T) );stopifnot(length(ff)>0)
type = "EVI2FINAL"
bands = "mean" # EVI2 aggregate
#bands = c("red","green","blue","nir","swir1","swir2","EVI2","NBR","NDMI","NDVI")
product = "LEDAPS"# "#LEDAPS" #"TOA"
bounds = T # Filter out extreme values
outlier = T # extreme outliers
fillgap = F #Save without gapfilling 
maxgap = 5

# Remove those files as they are too big
#ff <- ff[-which(str_detect(ff,'11_227__EVI_LEDAPS_MEAN_30_POLY.geojson'))]
#ff <- ff[-which(str_detect(ff,'11_233__EVI_LEDAPS_MEAN_30_POLY.geojson'))]

res <- list()
# Loop through
for(entry in ff){
  # Need to skip the empty ones
  a <- try( jsonlite::flatten( jsonlite::fromJSON(entry)$features ), silent = T )
  if(class(a)=="try-error"){ print(paste("Couldn't read",basename(entry))); next()}
  if( nrow(a) < 12) { next()}
  # read in

  for(site in unique(a$properties.SSBS)){
    sub <-  subset(a,properties.SSBS == site)
    if( nrow(sub) < 12) { next()}

    myLog("Processing ",unique(sub$properties.SSBS))
    subs <- sites %>% dplyr::select(SS,SSBS,Sample_start_earliest) %>% dplyr::filter(SSBS==site)
    if(nrow(subs)==0){myLog(site,"   NOT FOUND!   "); next()}
    # Get Julian date from ID
    id <- str_split(sub$id,"_",simplify = F)
    sub$date <- unlist(
      lapply(id, function(x){
        # Get the longest bit but without the last one
        d = which.max(str_length(x[1:(length(x)-1)]))
        return( x[[d]] ) # Longest bit is the date
      } )
    )
    # Format the date
    #sub$date <- as.POSIXct( strptime( sub$date, "%Y%j",tz = "GMT" ) ) # Y - j (pre-collection 1)
    sub$date <- as.POSIXct( strptime( sub$date, "%Y%m%d",tz = "GMT" ) ) # Y - m - d
    
    # Security check for date recognition
    if(any(is.na(sub$date))) stop("Something has gone wrong with the date recognition")
    
    # Get the satellite
    sub$Sat <- unlist(
      lapply(id, function(x){
        if(length(x)==6){
          return( x[3] )
        } else if(length(x) == 5){ return( x[2] ) } else { return(x[4])}
      })
    )

    sub$geometry <- NULL;sub$type <- NULL
    names(sub) <- str_replace_all(names(sub),"properties.","")
    
    sub <- merge.data.frame(x = sub,y = subs,by=c("SSBS") )

    # Only consider sites before sampling start
    sub <- sub[which( ymd(sub$date) <= ymd(unique(sub$Sample_start_earliest)) ),]
    
    if(!assertthat::has_name(sub,'mean')){ next()}
    # Build time series
    zz <- suppressWarnings( zooreg(data = sub[,bands], order.by = sub$date ) )
    
    if(length(time(zz)) < 12) {next()} # Useless time series
    # Filter out extreme values
    if(bounds){ # Filter by bounds to remove sensor errors
      zz[union(which(zz < -1),which(zz > 1))] = NA
    }
    
    # Filter out extreme values if existant
    # Detect and remove outliers based on MAD defined threshold
    # http://eurekastatistics.com/using-the-median-absolute-deviation-to-find-outliers/
    if(outlier){
      x = abs(zz - median(zz,na.rm = T)) / mad(zz,na.rm = T)
      tr = ifelse(quantile(x,.99,na.rm=T)>2,quantile(x,.99,na.rm=T),NA) # Determine threshold
      if(!is.na(tr)){
        zz[ which(x > tr )] <- NA
      }
    }

    # Aggregate to monthly data using a max value composite
    xx = suppressWarnings( aggregate(zz,as.yearmon, function(x) max(x,na.rm = T)) )
    xx[which(is.infinite(xx))] <- NA
    
    # Theoretical possible measurements
    # March 1, 1984 - Landsat 5 (landsat 4 = July 16, 1982 )
    x.na <- zooreg(data=NA,order.by =  as.yearmon(1982 + seq(0, (2014-1982)*12)/12),frequency = 12 )
    if(length(bands)== 1) xx <- merge(xx,x.na)[,1] else xx <- merge(xx,x.na)[,bands]
    
    xx.ori = xx
    
    if( fillgap ) {
      # Next fill all gaps using a kalman smoother of an ARIMA state-space model
      sub1_i <- try(
        na.kalman(as.ts(xx),model = "StructTS",smooth = T),
        silent=T)
      
      if(class(sub1_i)=="try-error") next() # If this fails that very likely the site is not suitable
      
      # The time series of NDVI was smoothed using the Savitzky–Golay filter with a length of 5 to reduce noise 
      # but still expose abrupt change events that might occur in the series (Jönsson & Eklundh, 2004)
      # Note this queries the signal package sgolayfilter function, which has other default parameters
      # that might defer from the default as set in the TIMESAT software
      sub1_i <- signal::sgolayfilt(sub1_i, n = 5)
      
      # Overwrite previous values
      coredata(xx) <- sub1_i
      
      # Reset larger gaps
      s.na <- na.approx(xx.ori,maxgap=maxgap)
      sub1_i2 <- merge(xx,s.na) # prev. cbind
      # Equalize gaps where continious gap length too long
      sub1_i2[which(is.na(apply(sub1_i2,1,mean))),1] <- NA
      xx <- sub1_i2[,1];rm(sub1_i2,s.na)
      
      # Furthermore reset permanent gaps in the time series to original values
      pg <- greenbrown::IsPermanentGap(as.ts(xx),min.gapfrac = 0.5,lower = NULL)
      xx[which(pg==T)] <- xx.ori[which(pg==T)]
      
      # Finally skip if more than 50% of possible measurements are missing
      gapt <- imputeTS::statsNA(xx,printOnly = F)
      # If consequetive NA gap is longer than half of the monitoring period
      if( ((gapt$naGapLongest / gapt$lengthTimeSeries) >=.5) ) { next() }
      rm(gapt)
    }
    
    # Get Samp.ling
    startd <- subset(sites,SSBS == site)
    # Finally create a window of only those values relevant that occured in the past
    out <- window(xx,end = as.yearmon(startd$Sample_start_earliest) )
    
    # Lastly trim
    out <- na.trim(out,sides = "both",is.na = "all")
    
    if(length(which(!is.na(out))) < 12){ next()}
    
    # Save in list
    res[[as.character(unique(sub$SS))]][[site]] <- out
    rm(out)
  }
  rm(a)
}

if( fillgap) {
  #saveRDS(res,paste0("resSaves/LS_",type,"_",product,"_GAPF_PolMean.rds")) 
  print('Not done in the paper!')
  } else{
  saveRDS(res,paste0("resSaves/LS_",type,"_",product,"_NONGAPF_PolMean.rds")) 
}

stop("DONE!")
