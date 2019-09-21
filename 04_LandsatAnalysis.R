# Package loading
par.ori <- par(no.readonly = T)
library(dplyr)
library(reshape2)
library(stringr)
library(lubridate)
library(tidyr)
library(TSdist)
library(zoo)
library(xts)
library(bfast)
library(strucchange)
library(stlplus)
library(broom)
library(imputeTS)
source("000_HelperFunction.R")
par.ori <- par(no.readonly = T)

# This is the downloaded and prepared list of EVI values per site
spast <- readRDS("resSaves/LS_EVI2FINAL_LEDAPS_NONGAPF_PolMean.rds") # Load data

# Previously prepared PREDICTS sites
sites <- readRDS("PREDICTS_allsites.rds")
sites <- subset(sites, SS %in% names(spast))


#### Apply BFAST on Landsat time series ####
stopifnot(exists("spast"))
# Parameters
band = "EVI2" 
season = "harmonic"
maxgap = 5 #
mvcpt = F
gapfill = F # Gap filling or not?
bfastF = F # Should TS be gapfilled and a normal bfast model be fitted ?
gapFTtest = T # Should a paired t-test on gap-filled data be computed
outlier = F # Should outliers be removed before LS processing?
bound = T # Should time series be bounded to 0-1 ?
shorten = T # Apply criteria for time series shortening in case of extensive gaps
reg = "lm"

# Gap fill
py_short <- list()
for(id in names(spast)){
  py = spast[[id]]
  # Also trim the time series to remove non-existing early data
  py <- lapply(py, function(x) lapply( x, na.trim  ) )
  
  # Finally remove studies from spast that have only one site
  if( length( names(py) ) < 2 ) {next()}
  
  for(SSBS in names(py) ){
    myLog(SSBS)
    x = py[[SSBS]]
    x = do.call("merge",x) # Merge all bands
    if(is.null(x)) { next()}
    if(length(x)<12) { next()}
    
    if(gapfill) {
      # Trim 
      x = na.trim(x)
      x.ori = x

      # Next fill all gaps using a kalman smoother of an ARIMA state-space model
      sub1_i <- try(
        na.kalman(coredata(x),model = "StructTS",smooth = TRUE,type="BSM"),
        silent=T)
      
      if(class(sub1_i)=="try-error") next() # If this fails that very likely the site is not suitable
      
      # The time series of NDVI was smoothed using the Savitzky–Golay filter with a length of 5 to reduce noise 
      # but still expose abrupt change events that might occur in the series (Jönsson & Eklundh, 2004)
      # Note this queries the signal package sgolayfilter function, which has other default parameters
      # that might defer from the default as set in the TIMESAT software
      sub1_i <- signal::sgolayfilt(sub1_i, n = 5)
      
      # Overwrite previous values
      coredata(x) <- sub1_i
      
      # Reset larger gaps
      s.na <- na.approx(x.ori,maxgap=maxgap)
      sub1_i2 <- merge(x,s.na) # prev. cbind
      # Equalize gaps where continious gap length too long
      sub1_i2[which(is.na(apply(sub1_i2,1,mean))),1] <- NA
      subi_1 <- sub1_i2[,1];rm(sub1_i2,s.na)
      
      py_short[[SSBS]] <- subi_1
      
    } else {
      py_short[[SSBS]] <- x
    }
  }
  rm(py)
}
rm(spast)
stop("Gap-filled!!!")

# Now calculate for all studies 
# if a breakpoint occurs / how many
# the biggest and last magnitude change
# and since the last disturbance:
# average values, trend direction before and after the tend
res <- data.frame() # Save output2
res2 <- data.frame()
# Counter
cou = length(names(py_short))

library(parallel)
library(doParallel)
cl <- parallel::makeCluster(detectCores()-2)
doParallel::registerDoParallel(cl)

res <- foreach(id =  names(py_short),
        .combine = rbind,
        .multicombine = FALSE,
        .errorhandling = 'stop',
        .packages = c("dplyr",'zoo','imputeTS','bfast','xts','lubridate'),
        .verbose = FALSE
        ) %dopar% {
  #  print(id) # id = "HP1_2005__Kumar 1  28" # id = "AD1_2001__Liow 1  2"
  myLog("Processing ",id," | ", cou," remaining...");cou = cou-1
  xo = py_short[[id]]
  if(is.null(xo)){ return( data.frame() ) }
  if(length(xo)==0){  return( data.frame() ) }

  # For each band
  for(b in band){
    x = xo[,b]
    
    if(outlier){
      zz = x
      x = abs(zz - median(zz,na.rm = T)) / mad(zz,na.rm = T)
      tr = ifelse(quantile(x,.95,na.rm=T)>2,quantile(x,.95,na.rm=T),NA) # Determine threshold
      if(!is.na(tr)){
        zz[ which(x > tr )] <- NA
      }
      x = na.trim(zz)
    }
    if(bound){ x[x < 0 | x > 1] <- NA; x = na.trim(x) }
    ## --- ##
    # Shorten if long extensive gaps
    if(shorten){
      # If tere is a 5 year gap
        if( ((imputeTS::statsNA(x,printOnly = F)$naGapLongest)/12) >5 ){
        xx = na.trim(window(x,start = 1999)) # Clip since Landsat 7
        # Did it help?
        #((imputeTS::statsNA(xx,printOnly = F)$naGapLongest)/12)
        x = xx
        if( length(x) < 24 ) return( data.frame() ) # Shorter than 2 year time series
        lg = ((imputeTS::statsNA(x,printOnly = F)$naGapLongest)/12)
        if( is.na(lg) | lg > 5)  return( data.frame() )
      }
      print(x)
    } 
    
    # Best largest break based on deseasonalized trend and /or trend + harmonic
    # First trend and harmonic model
    rdist = 12/(length(x)) # calculate h relative to total sample length
    bf1 = try(bfast01(as.ts(x),formula = response ~  trend + harmon, test = c("OLS-MOSUM","BIC"),order = 2,aggregate=any, bandwidth = rdist ),silent=T) # In trend + harmon
    if(class(bf1)!= "try-error") { l = list(bf1); model.deseasonalized = 0} else {
      # Otherwise deseasonalize the time series and run normal model
      x_g <- try( stlplus::stlplus(x,s.window = "periodic")$data,silent = T )
      if(class(x_g)!="try-error") {
        # Can a seasonal term be determined
        x2 <- x;coredata(x2) <- x_g$trend+x_g$remainder
        bf2 = try(bfast01(as.ts(x2),formula = response ~  trend, test = c("OLS-MOSUM","BIC"),aggregate=any,bandwidth =rdist),silent=T) # In trend
        l = list(bf2)
        model.deseasonalized = 1
        rm(x_g,x2)
      } else {
        # T-test results empty
        tt <- data.frame(tt.estimate = NA,tt.statistic = NA, tt.parameter = NA, tt.pvalue = NA )
        tt2 <- data.frame(tt.estimate.gf = NA,tt.statistic.gf = NA,tt.parameter.gf = NA, tt.pvalue.gf = NA)
        # BF failed
        bfc = data.frame(flag_type=NA,pct_segment1=NA,pct_segment2=NA,flag_significance=NA)
        model.deseasonalized = 0
        return(
                     data.frame(SSBS = id,band = b,
                                PastPeriod=NA,PastPeriodNA = NA,Break=0,LongestGap = NA, # General info
                                befpoint = NA,bf.AIC = NA,
                                largest_year = NA,largest_mag = NA , largest_mag_prop = NA,
                                tt,tt2,
                                na.1ybefore = NA,na.1yafter = NA,
                                largest_trendbef =  NA, largest_trendaft = NA,rec_time = NA,
                                AvgOverall=NA,AvgBefore=NA,AvgAfter=NA,
                                la_nrbf=NA,la_sdMag=NA,la_avgMag=NA,la_largest_mag=NA, la_largest_year = NA, 
                                bfc,reg,
                                test="Largest Breakpoint", model = NA,model.deseasonalized  )
        )
        rm(tt,tt2,bfc)
        next() # Continue
      }
    }
    
    # Loop through
    for(bf in l){
      if(class(bf)!="try-error"){
        
        if(bf$breaks > 0){
          # Time of largest change
          la_year <- bf$data$time[bf$breakpoints]
          # Also write magnitude of biggest change - # As predicted change in trends
          #la_mag <- as.numeric( predict(bf$model[[2]])[bf$breakpoints+1] - predict(bf$model[[2]])[bf$breakpoints] )
          la_mag <- as.numeric( (predict(lm(response~time,data=bf$data[which(bf$data$segment==2),]))[1]) - (predict(lm(response~time,data=bf$data[which(bf$data$segment==1),]))[bf$breakpoints]) )
          
          # Magnitude loss/gain in proportion to relative mean before
          #la_mag_prop <- (bf$Mags[which.max(abs(bf$Mags[,3])),2] - bf$Mags[which.max(abs(bf$Mags[,3])),1]) /  bf$Mags[which.max(abs(bf$Mags[,3])),1]
          la_mag_prop <-  as.numeric(la_mag /  abs( predict(lm(response~time,data=bf$data[which(bf$data$segment==1),]))[bf$breakpoints] ) )
          
          befpoint = predict(lm(response~time,data=bf$data[which(bf$data$segment==1),]))[bf$breakpoints]
          # Linear Trend before largest breakpoint
          #y = bf$output[[1]]$Tt[1:bf$Time-1]
          #la_tbef <- try( coef(lm(y~time(y)))[2],silent = T)
          #la_tbef <- ifelse(class(la_tbef)=="try-error", NA,la_tbef) # in case the regression doesn't fit
          la_tbef = coef(bf)[1,2]
          # Trend after largest breakpoint
          #y = bf$output[[1]]$Tt[(bf$Time+1):length(bf$output[[1]]$Tt)]
          #la_taft <- try( coef(lm(y~time(y)))[2], silent = T )
          #la_taft <- ifelse(class(la_taft)=="try-error", NA,la_taft) # in case the regression doesn't fit
          la_taft = coef(bf)[2,2]
          
          # In the trend after break: Assess recovery. When (if at all) does it reach pre-disturbance levels again ?
          rec = (predict(lm(response~time,data=bf$data[which(bf$data$segment==1),]))[bf$confint[1]] - predict(lm(response~time,data=bf$data[which(bf$data$segment==2),] )) )
          if(la_mag < 0){ # Negative
            rec = bf$data$time[which(bf$data$segment==2)][which(rec < 0)[1]] # First date when pre-dist is reached again
          } else{
            rec = bf$data$time[which(bf$data$segment==2)][which(rec > 0)[1]] # First date when pre-dist is reached again
          }
          rec_time = ifelse(length(rec) == 0,NA,rec - la_year)
          
          # Ratio of year before to year at time of sampling
          # Follow https://www.sciencedirect.com/science/article/pii/S0034425717305308#bb0305
          # Ratio between before and after magnitude
          #mean(x[(length(x)-11):length(x)],na.rm=T) / mean(w.earlier,na.rm=T) 
          
          # Missing observations before and after disturbance
          w.earlier = window(x,start = (as.yearmon(la_year)+1/12)-1,end = as.yearmon(la_year)) # One year including break date
          w.later = window(x,start = (as.yearmon(la_year)+ 1/12),end = (as.yearmon(la_year)) +1 )
          na.1ybefore = length(which(is.na(w.earlier)))
          na.1yafter = length(which(is.na(w.later)))
          
          # Overall, what is the minimum of the deseasonalized trend
          #MinimumTrend = min(stl(x,s.window = "periodic")$time.series[,2])
          # Or maximum
          #MaximumTrend = max(stl(x,s.window = "periodic")$time.series[,2])
          # Average of the time series
          AvgOverall = mean(x,na.rm=T)
          # Use from earlier
          AvgBefore = mean( coredata(w.earlier),na.rm=T ) 
          AvgAfter = mean( coredata(w.later),na.rm = T)

          # Normal T-test
          # Enough observations ?
          if( sum(!is.na(coredata(w.earlier))) < 2 | sum(!is.na(coredata(w.later))) < 2 ) {
            tt <- data.frame(tt.estimate = NA, tt.statistic = NA, tt.parameter = NA, tt.pvalue = NA )
          } else {
            k = inner_join(
              data.frame(m = month(w.earlier), x1 = coredata(w.earlier)),
              data.frame(m = month(w.later), x2 = coredata(w.later)),
              by = "m"
            ) %>% subset(.,complete.cases(.))
            # Construct a paired two-sample t.test
            try(ptt <- t.test(k$x1,k$x2,paired = TRUE),silent=T)
            if(!exists("ptt")) {tt <- data.frame(tt.estimate = NA, tt.statistic = NA, tt.parameter = NA, tt.pvalue = NA )} else
              if(class(ptt) == "try-error"){
                tt <- data.frame(tt.estimate = NA, tt.statistic = NA, tt.parameter = NA, tt.pvalue = NA )
              } else {
                tt <- ptt %>% broom::glance() %>% dplyr::select(estimate,statistic,parameter,p.value) %>%
                  dplyr::rename(tt.estimate = estimate, tt.statistic = statistic, tt.parameter = parameter, tt.pvalue = p.value )
                rm(ptt)
              }
            rm(k)
          }
          
          # t-test using interpolated data
          if(gapFTtest){
            xx = na.kalman(x,model = "StructTS",smooth = FALSE)
            w.earlier2 = window(xx,start = (as.yearmon(la_year)+1/12)-1,end = as.yearmon(la_year)) # One year including break date
            w.later2 = window(xx,start = (as.yearmon(la_year)+ 1/12),end = (as.yearmon(la_year)) +1 )
            
            # Enough observations ?
            if( sum(!is.na(coredata(w.earlier2))) < 2 | sum(!is.na(coredata(w.later2))) < 2 ) {
              tt2 <- data.frame(tt.estimate = NA, tt.statistic = NA, tt.parameter = NA, tt.pvalue = NA )
            } else {
              k = inner_join(
                data.frame(m = month(w.earlier2), x1 = coredata(w.earlier2)),
                data.frame(m = month(w.later2), x2 = coredata(w.later2)),
                by = "m"
              ) %>% subset(.,complete.cases(.))
              # Construct a paired two-sample t.test
              try(ptt <- t.test(k$x1,k$x2,paired = TRUE),silent=T)  
              if(!exists("ptt")) {tt2 <- data.frame(tt.estimate = NA, tt.statistic = NA, tt.parameter = NA, tt.pvalue = NA )} else 
                if(class(ptt) == "try-error"){
                  tt2 <- data.frame(tt.estimate = NA, tt.statistic = NA, tt.parameter = NA, tt.pvalue = NA )
                } else {
                  tt2 <- ptt %>% broom::glance() %>% dplyr::select(estimate,statistic,parameter,p.value) %>% 
                    dplyr::rename(tt.estimate = estimate, tt.statistic = statistic, tt.parameter = parameter, tt.pvalue = p.value )
                  rm(ptt)
                }
              rm(k)
            }
            tt2 <- tt2 %>% rename(tt.estimate.gf = tt.estimate, tt.statistic.gf = tt.statistic,tt.parameter.gf = tt.parameter, tt.pvalue.gf = tt.pvalue)
          } else {
            tt2 <- data.frame(tt.estimate = NA, tt.statistic = NA, tt.parameter = NA, tt.pvalue = NA )%>% rename(tt.estimate.gf = tt.estimate, tt.statistic.gf = tt.statistic,tt.parameter.gf = tt.parameter, tt.pvalue.gf = tt.pvalue)
          }
          
          # ALso do the normal bfast model
          if( bfastF ){
            xx = na.kalman(x,model = "StructTS",smooth = FALSE)
            # Generally number of breakpoints
            bf3 <-  bfast(as.ts(xx),h = rdist,season = "harmonic",type = "OLS-MOSUM",max.iter = 10)
            la_nrbf = if(any(bf3$Mags==0)){ NA } else { nrow(bf3$Mags) }
            la_sdMag = if(any(bf3$Mags==0)){ NA } else { sd(bf3$Mags[,3]) } 
            la_avgMag = if(any(bf3$Mags==0)){ NA } else { mean(bf3$Mags[,3]) }
            la_largest_mag = if(any(bf3$Mags==0)){ NA } else { bf3$Magnitude }
            la_largest_year = if(any(bf3$Mags==0)){ NA } else { time(xx)[bf3$Time]  }
          } else{
            la_nrbf = NA;la_sdMag = NA;la_avgMag=NA; la_largest_mag=NA; la_largest_year = NA;  
          }
          
          # Classify the change type
          # pct_segment = percent change per unit time (0-100)
          bfc = bfast01classify(bf,alpha = 0.05) %>% dplyr::select(flag_type,pct_segment1,pct_segment2,flag_significance)
          
          # Save all
          return(
                       data.frame(SSBS = id,band = b, Break=1, # General info
                                  PastPeriod=length(x),PastPeriodNA = length(which(!is.na(x))),LongestGap = imputeTS::statsNA(x,printOnly = F)$naGapLongest,
                                  befpoint,bf.AIC = AIC(bf),
                                  largest_year = la_year,largest_mag = la_mag , largest_mag_prop = la_mag_prop,
                                  na.1ybefore,na.1yafter,
                                  tt,tt2,
                                  largest_trendbef = la_tbef, largest_trendaft = la_taft,rec_time,
                                  AvgOverall,AvgBefore,AvgAfter,
                                  la_nrbf,la_sdMag,la_avgMag,la_largest_mag,la_largest_year,
                                  bfc,reg,
                                  test="Largest Breakpoint", model = ifelse("harmon" %in% all.vars(bf$formula),"harmon","trend"),model.deseasonalized)
                       )
          rm(tt)
        } else {
          # No changepoint detected
          # Just get the overall trend and segment
          trend = as.numeric( coef(bf)["trend"] )
          AvgOverall = mean(bf$data$response,na.rm=T)
          bfc = bfast01classify(bf,alpha = 0.05) %>% dplyr::select(flag_type,pct_segment1,pct_segment2,flag_significance)
          # T-test results empty
          tt <- data.frame(tt.estimate = NA,tt.statistic = NA, tt.parameter = NA, tt.pvalue = NA )
          tt2 <- data.frame(tt.estimate.gf = NA,tt.statistic.gf = NA,tt.parameter.gf = NA, tt.pvalue.gf = NA)
          return(
                       data.frame(SSBS = id,band = b,
                                  PastPeriod=length(x),PastPeriodNA = length(which(!is.na(x))),Break=0,LongestGap = imputeTS::statsNA(x,printOnly = F)$naGapLongest, # General info
                                  befpoint = NA,bf.AIC = AIC(bf),
                                  largest_year = NA,largest_mag = NA , largest_mag_prop = NA,
                                  na.1ybefore = NA,na.1yafter = NA,
                                  tt,tt2,
                                  largest_trendbef = trend, largest_trendaft = NA,rec_time = NA,
                                  AvgOverall,AvgBefore=NA,AvgAfter=NA,
                                  la_nrbf=NA,la_sdMag=NA,la_avgMag=NA,la_largest_mag=NA, la_largest_year = NA, 
                                  bfc,reg,
                                  test="Largest Breakpoint", model = ifelse("harmon" %in% all.vars(bf$formula),"harmon","trend"),model.deseasonalized )
                       )
          rm(tt)
        }
        
      } else {
        # T-test results empty
        tt <- data.frame(tt.estimate = NA, tt.statistic = NA, tt.parameter = NA, tt.pvalue = NA )
        tt2 <- data.frame(tt.estimate.gf = NA,tt.statistic.gf = NA,tt.parameter.gf = NA, tt.pvalue.gf = NA)
        # BF failed
        bfc = data.frame(flag_type=NA,pct_segment1=NA,pct_segment2=NA,flag_significance=NA)
        return( 
                     data.frame(SSBS = id,band = b,
                                PastPeriod=NA,PastPeriodNA = NA,Break=0,LongestGap = NA, # General info
                                befpoint = NA,bf.AIC = NA,
                                largest_year = NA,largest_mag = NA , largest_mag_prop = NA,
                                tt,tt2,
                                na.1ybefore = NA,na.1yafter = NA,
                                largest_trendbef =  NA, largest_trendaft = NA,rec_time = NA,
                                AvgOverall=NA,AvgBefore=NA,AvgAfter=NA,
                                la_nrbf=NA,la_sdMag=NA,la_avgMag=NA,la_largest_mag=NA, la_largest_year = NA, 
                                bfc,reg,
                                test="Largest Breakpoint", model = NA,model.deseasonalized  )
              )
        rm(tt)
        } 
      }
  }
  
    rm(bf,bf1,bf2,l)
    
} # End for foreach loop

stopCluster(cl)
stop("DONE!!!")
res$bestmodel <- 1
saveRDS(res,paste0("resSaves/LS_EVI2_BFast",ifelse(outlier,"outlierREM",""),"_",reg,"_","BOUNDS_LEDAPS_SEASON.rds")) # <- File for all

# ----------------------- #
# Go through all and save the best model!
fname = paste0("resSaves/LS_EVI2_BFast",ifelse(outlier,"outlierREM",""),"_",reg,"_","BOUNDS_LEDAPS_SEASON.rds")
# Load in files
res_ls_bfast <- readRDS(fname) 
res_ls_bfast <- subset(res_ls_bfast,band=="EVI2") # %>% dplyr::filter(model == "trend") 
res_ls_bfast$bestmodel <- 0
res_ls_bfast2  <- data.frame()
for(i in unique(res_ls_bfast$SSBS)){
  myLog(i)
  ss <- subset(res_ls_bfast,SSBS == i)
  if(nrow(ss)>2){ stop()}
  ss$bestmodel[ which.min((ss$bf.AIC))] <- 1
  res_ls_bfast2 <- rbind(res_ls_bfast2, ss )
}
stopifnot( nrow(res_ls_bfast) == nrow(res_ls_bfast2) )
saveRDS(res_ls_bfast2,fname) 
