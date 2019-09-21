### Helper functions ####
library(dplyr)


# Open png after saving it
ggpreview <- function (..., device = "png", dpi = 400) {
  fname <- tempfile(fileext = paste0(".", device))
  ggplot2::ggsave(filename = fname,device = device,dpi = dpi,  ...)
  system2("open", fname)
  invisible(NULL) 
}

myLog <- function(...) {
  cat(paste0("[Processing] ", Sys.time(), " | ", ..., "\n"))
}

# Mode
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


# Standard error of the mean
sem <- function (x) 
{
  sqrt(var(x, na.rm = TRUE)/length(na.omit(x)))
}

vif.mer <- function (fit) 
{
  require(lme4)
  v <- vcov(fit)
  nam <- names(fixef(fit))
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)]
  }
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v
}

# Get standardized coefficients for lmer
lm.beta.lmer <- function(mod) {
  b <- fixef(mod)[-1]
  sd.x <- apply(getME(mod,"X")[,-1],2,sd)
  sd.y <- sd(getME(mod,"y"))
  b*sd.x/sd.y
}

# Describe function
describe <- function(x) {Hmisc::describe(x)} # Convenience describe function
# Not in function
'%!in%' <- function(x,y)!('%in%'(x,y))

# ---------------------------------------

# Function to isolate stable periods from time series - cpt.varmean
stabPer <- function(x,cp){
  cp.l <- cpts(cp) # Location of changepoints
  # Return a list of subsets
  o <- list(x[1:cp.l[1]]) # The first one
  cp.l <- c(cp.l, length(x)) # Remove the first and add the last
  for(i in 1:(length(cp.l)-1)){ 
    o[[i+1]] <- x[cp.l[i]:cp.l[i+1]] 
  }
  return(o)
}

new.freq <- function(d, freq = 46) {
  y <- as.Date(cut(range(d), "years")) + c(0, 366)
  yd <- seq(y[1], y[2], "year")
  yy <- as.numeric(format(yd, "%Y"))
  floor(freq * approx(yd, yy, xout = d)$y) / freq
}

# Load in al MODIS BRDF Bands
readInFormat <- function(x,idv = "SSBS"){
  require(data.table)
  val = as.data.frame(data.table::fread(x,showProgress = T))[,-1] %>% dplyr::select(-.geo,-latitude, -longitude) %>% 
    reshape2::melt(., id.vars = idv) %>%
    distinct() %>% mutate(variable = as.character(variable)) %>% 
    mutate(year = as.numeric( str_sub(variable,13,16)), # Get year out of file name
           month = as.numeric( str_sub(variable,18,19) ), # get month
           day = as.numeric( str_sub(variable,21,22) )) %>%  # get day
    # Make a date column
    mutate(date = ymd(paste(year,month,day,sep="-")))
  return(val)
}

# Function to convert from MODIS QA scores - BRDF Albedo band scores
# Inspired by MODISTools -> Tuck et al.
BRDFconvertQAScores <- function(qsc,band=NA,QualityThreshold=3){
  #MCD43A4 = c(0,4294967294,4294967295), # BRDF albedo band quality, taken from MCD43A2, for reflectance data
  QualityScores <- qsc
  if(max(QualityScores,na.rm = T) == 0) num.binary.digits <- 1
  if(max(QualityScores,na.rm = T) != 0) num.binary.digits <- floor(log(max(QualityScores,na.rm = T), base = 2)) + 1
  
  binary.set<- matrix(nrow = length(QualityScores), ncol = num.binary.digits)
  for(n in 1:num.binary.digits){
    binary.set[ ,(num.binary.digits - n) + 1] <- QualityScores %% 2
    QualityScores <- QualityScores %/% 2
  }
  # Construct quality binary score
  quality.binary <- apply(binary.set, 1, function(x) paste(x, collapse = ""))
  rm(binary.set)
  
  if(!is.na(band)){
    band.num <- as.numeric(substr(band, nchar(band), nchar(band)))
    
    # Select respective subset of quality score
    qa.binary <- substr(quality.binary, (nchar(quality.binary) - (((band.num - 1) * 2) + 2)),
                        (nchar(quality.binary) - ((band.num - 1) * 2)))
    
    # Make a result vector
    qa.int <- numeric(length(qa.binary))
    qa.int[qa.binary == "000"] <- 0 # best quality, full inversion (WoDs, RMSE majority good)
    qa.int[qa.binary == "001"] <- 1 # good quality, full inversion
    qa.int[qa.binary == "010"] <- 2 # Magnitude inversion (numobs >=7)
    qa.int[qa.binary == "011"] <- 3 # Magnitude inversion (numobs >=3&<7)
    qa.int[qa.binary == "100"] <- 4 # Fill value
    qa.int[qa.binary == "ANA" | is.na(qa.binary) ] <- NA # NA
    
    # Finally replace everything above the treshold with NA
    qa.int[qa.int > QualityThreshold] <- NA
    # And return
    return(qa.int)
  } else {
   # If band score is NA, it is assumed that we intend to calcualte the solar zenith angle
    
    # Solar zenith information is stored between 8-14
    qa.binary <- substr(quality.binary,8,14)
    qa.binary[qa.binary == "ANANANA"] <- NA
    
    x = strtoi(qa.binary, base = 2)
    # Which values have to high zenith values -> Encode with 1 otherwise 0
    x = ifelse(x>QualityThreshold,1,0)
  }
}

# Backtransformation of asin(sqrt)
bt.asin <- function(x) {
  z <- sin(x)^2
  return(z)
}

#### Max CCF ####
# Returns the maximal value of the crosscorrelation function,
# thus indicating the time when there is a certain lag
aMaxCCF <- function(a,b,lt = 5)
{
  d <- ccf(a, b, plot = FALSE, lag.max = length(a)-lt)
  cor = d$acf[,,1]
  abscor = abs(d$acf[,,1])
  lag = d$lag[,,1]
  res = data.frame(cor,lag)
  absres = data.frame(abscor,lag)
  absres_max = res[which.max(absres$abscor),]
  return(absres_max)
}

#### LME4 functions ####
# Non-Convergence
checkConv <- function(mod){
  gg <- mod@optinfo$derivs$grad
  hh <- mod@optinfo$derivs$Hessian
  vv <- sqrt(diag(solve(hh/2)))
  cat("Should be smaller or around ~ e-05 --> ",mean(abs(gg*vv)),"\n")
  #  the numbers printed here should be very small (examples that the developers deemed acceptable were around e-05.
  relgrad <- with(mod@optinfo$derivs,solve(Hessian,gradient))
  cat("Should be smaller than 0.001 --> ",max(abs(relgrad)),"\n") # this number should ideally be < 0.001.        
}

# Convergence boolean test
checkConvTest <- function(mod){
  gg <- mod@optinfo$derivs$grad
  hh <- mod@optinfo$derivs$Hessian
  # In case system is exactly singular 
  try( vv <- sqrt(diag(solve(hh/2))), silent=T )
  if(exists("vv")){
    if( any(is.nan(vv))) { return(FALSE)}
    # Otherwise test for degenerated Hessians
    try( relgrad <- with(mod@optinfo$derivs,solve(Hessian,gradient)),silent=T )
    if(!exists("relgrad")) {return(FALSE)}
    if( all( mean(abs(gg*vv)) < 0.00001,max(abs(relgrad)) < 0.001 ) ){
      return(TRUE)
    } else {
      return(FALSE)
    }
  } else {return(FALSE)} 
}

## Check for overdispersion
overdisp_fun <- function(model) {
  ## number of variance parameters in 
  ##   an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
  rdf <- nrow(model.frame(model))-model.df
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}


#### Compdis update #####

# Requires a full predicts subset
CompDissim2 <- function (data, metric,binary=F)
{
  require(vegan)
  require(data.table)
  data <- as.data.table(data) # For speed improvement
  if (metric == "SorAbd") {
    data <- data[data$Diversity_metric_type == "Abundance", 
                 ]
  }
  if (metric == "SorCorr") {
    data <- data[((data$Diversity_metric == "abundance") & 
                    (data$Diversity_metric_unit == "individuals")), ]
  }
  if (metric == "BCVeg" & binary == F) {
    data <- data[data$Diversity_metric_type == "Abundance",]
  }
  if (metric == "BrayCurt" & binary == F) {
    data <- data[data$Diversity_metric_type == "Abundance",]
  }
  
  data <- subset(data, select = c("SS", "SSBS", "Measurement", 
                                  "Taxon_name_entered"))
  data <- na.omit(data)
  if (metric == "SorCorr") {
    study.all.int.meas <- tapply(data$Measurement, data$SS, 
                                 function(m) all(floor(m) == m))
    int.meas <- study.all.int.meas[match(data$SS, names(study.all.int.meas))]
    data <- data[int.meas, ]
  }
  
  # Results
  result <- list()
  for (st in unique(data$SS)) {
    cat(paste("Processing ", st, "\n", sep = ""))
    sub.data <- data[data$SS == st, ]
    if (metric != "SorCorr") {
      sub.data <- sub.data[sub.data$Measurement > 0, ]
    }
    if(metric == "SorVeg"){
      if(length(unique(sub.data$SSBS)) < 2) next()
      m <- reshape2::acast(data=sub.data,SSBS~Taxon_name_entered,value.var = "Measurement",fill = 0) # Fill with zero assuming absence
      sites.matrix <- as.matrix( suppressWarnings( vegan::betadiver(m,method = "sor",binary=binary) ) )

      #sites.matrix[upper.tri(sites.matrix)] <- -9999
      diag(sites.matrix) <- NA  # Set diagonal to NA
      #mat_a <- reshape2::melt(sites.matrix) %>% rename(SSBS_x = Var1, SSBS_y = Var2 ) %>% #Melt and rename
      #  mutate(SSBS_x = as.character(SSBS_x), SSBS_y = as.character(SSBS_y)) # format again
      #mat_a <- mat_a[-which(mat_a$value==-9999),] # remove triangle and missing data
      #mat_a$SS <- st # Finally append the study ID
      
      result[[st]] <- sites.matrix#mat_a
      
    } else if(metric == "BCVeg"){
      if(binary==T) stop("Binary should be FALSE for BcVEG")
      if(length(unique(sub.data$SSBS)) < 2) next()
      m <- reshape2::acast(data=sub.data,SSBS~Taxon_name_entered,value.var = "Measurement",fill = 0) # Fill with zero assuming absence
      sites.matrix <- as.matrix( suppressWarnings( vegdist(m,method="bray",binary=binary,na.rm = T) ) )
      
      #sites.matrix[upper.tri(sites.matrix)] <- -9999
      diag(sites.matrix) <- NA  # Set diagonal to NA
      #mat_a <- reshape2::melt(sites.matrix) %>% rename(SSBS_x = Var1, SSBS_y = Var2 ) %>% #Melt and rename
      #  mutate(SSBS_x = as.character(SSBS_x), SSBS_y = as.character(SSBS_y)) # format again
      #mat_a <- mat_a[-which(mat_a$value==-9999),] # remove triangle and missing data
      #mat_a$SS <- st # Finally append the study ID
      
      result[[st]] <- sites.matrix#mat_a
      
    } else {
      sites.matrix <- matrix(nrow = length(unique(sub.data$SSBS)), 
                             ncol = length(unique(sub.data$SSBS)))
      i1 <- 1
      for (s1 in unique(sub.data$SSBS)) {
        i2 <- 1
        for (s2 in unique(sub.data$SSBS)) {
          if (metric == "Sor") {
            # sorrensen
            u <- length(union(sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                            s1], sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                                               s2]))
            i <- length(intersect(sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                s1], sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                                                   s2]))
            sor <- (2 * i)/((2 * i) + (u - i))
          }
          else if (metric == "SorAbd") {
            # abundance corrected sorrensen
            u <- sum(sub.data$Measurement[(sub.data$SSBS == 
                                             s1) & (sub.data$Taxon_name_entered %in% union(sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                                                                         s1], sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                                                                                                            s2]))])
            v <- sum(sub.data$Measurement[(sub.data$SSBS == 
                                             s2) & (sub.data$Taxon_name_entered %in% union(sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                                                                         s1], sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                                                                                                            s2]))])
            sor <- (2 * u * v)/(u + v)
          }
          else if (metric == "SorBM") {
            # abundance corrected sorrensen
            # Sum of estimates of species recorded at both sites
            A <- sum(sub.data$Measurement[union(sub.data$SSBS == s1,sub.data$SSBS == s2) &
                                            (sub.data$Taxon_name_entered %in% union(sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                                                                         s1], sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                                                                                                            s2]))])
            
            B <- sum(sub.data$Measurement[(sub.data$SSBS == s1) & 
                                            (sub.data$Taxon_name_entered %in% union(sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                                                                  s1], sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                                                                                                     s2]))])
            C <- sum(sub.data$Measurement[(sub.data$SSBS == s2) & 
                                            (sub.data$Taxon_name_entered %in% union(sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                                                                  s1], sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                                                                                                     s2]))])
            sor <- (2 * A)/(2*A + B + C)
          }
          else if (metric == "BrayCurt"){
           # Bray curtis similarity 
            u <- sub.data$Measurement[(sub.data$SSBS == 
                                             s1) & (sub.data$Taxon_name_entered %in% union(sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                                                                         s1], sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                                                                                                            s2]))]
            v <- sub.data$Measurement[(sub.data$SSBS == 
                                             s2) & (sub.data$Taxon_name_entered %in% union(sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                                                                         s1], sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                                                                                                            s2]))]
            sor <- ( sum( abs(u - v) )  / ( sum(u) + sum(v) ) ) # Similarity 1- following Faith
            # Tested and completly identical to vegan vegdist
            
          }
          else if (metric == "SorCorr") {
            # Sampling corrected sorrensen
            n <- sum(sub.data$Measurement[sub.data$SSBS == 
                                            s1])
            m <- sum(sub.data$Measurement[sub.data$SSBS == 
                                            s2])
            if ((n > 0) & (m > 0)) {
              xi <- sub.data$Measurement[sub.data$SSBS == 
                                           s1][(match(union(sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                                          s1], sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                                                                             s2]), sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                                                                                                                 s1]))]
              yi <- sub.data$Measurement[sub.data$SSBS == 
                                           s2][(match(union(sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                                          s1], sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                                                                             s2]), sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                                                                                                                 s2]))]
              xi[is.na(xi)] <- 0
              yi[is.na(yi)] <- 0
              f1. <- length(which((xi == 1) & (yi > 0)))
              f2. <- max(1, length(which((xi == 2) & (yi > 
                                                        0))))
              f.1 <- length(which((xi > 0) & (yi == 1)))
              f.2 <- max(1, length(which((xi > 0) & (yi == 
                                                       2))))
              p1 <- sum(xi[yi > 0]/n)
              p2 <- ((m - 1)/m) * (f.1/(2 * f.2))
              p3 <- sum(xi[yi == 1]/n)
              u <- min(1, p1 + p2 * p3)
              q1 <- sum(yi[xi > 0]/m)
              q2 <- ((n - 1)/n) * (f1./(2 * f2.))
              q3 <- sum(yi[xi == 1]/m)
              v <- min(1, q1 + q2 * q3)
              if ((u > 0) & (v > 0)) {
                sor <- (2 * u * v)/(u + v)
              }
              else {
                sor <- 0
              }
            }
            else {
              sor <- 0
            }
          }
          else {
            stop("Error: specfied dissimilarity metric is not supported")
          }
          if (s1 != s2) 
            sites.matrix[i1, i2] <- sor
          i2 <- i2 + 1
        }
        i1 <- i1 + 1
      }
      #sites.matrix[upper.tri(sites.matrix)] <- -9999
      #diag(sites.matrix) <- -9999  # Remove duplicates and set diagonal to -9999
      rownames(sites.matrix) <- unique(sub.data$SSBS); colnames(sites.matrix) <- unique(sub.data$SSBS)
      #mat_a <- reshape2::melt(sites.matrix) %>% rename(SSBS_x = Var1, SSBS_y = Var2 ) %>% #Melt and rename
      #  mutate(SSBS_x = as.character(SSBS_x), SSBS_y = as.character(SSBS_y)) # format again
      #mat_a <- mat_a[-which(mat_a$value==-9999),] # remove triangle and missing data
      #mat_a$SS <- st # Finally append the study ID
      
      
      result[[st]] <- sites.matrix
    } 
  }
  return(result)
}


#### Monthly functions ####

# File monthlyfunction.R
# Part of the hydroTSM R package, http://www.rforge.net/hydroTSM/ ; 
#                                 http://cran.r-project.org/web/packages/hydroTSM/
# Copyright 2008-2013 Mauricio Zambrano-Bigiarini
# Distributed under GPL 2 or later

################################################################################
# monthlyfunction: Generic function for applying any R function to             #
#                  ALL the values in 'x' belonging to a given month            #
################################################################################
# Author : Mauricio Zambrano-Bigiarini                                         #
################################################################################
# Started: May 15th, 2009;                                                     #
# Updates: 31-Aug-2009 ; 25-Jul-2011 ; 08-Aug-2011                             #
#          08-Apr-2013                                                         #
################################################################################
# 'x   '    : variable of type 'zoo' or 'data.frame', with daily or monthly frequency
# 'FUN'      : Function that will be applied to ALL the values in 'x' belonging to each one of the 12 months of the year
# 'na.rm'    : Logical. Should missing values be removed?
#              TRUE : the monthly values  are computed considering only those values in 'x' different from NA
#              FALSE: if there is AT LEAST one NA within a month, the FUN and monthly values are NA
monthlyfunction <- function(x, ...) UseMethod("monthlyfunction")

monthlyfunction.default <- function(x, FUN, na.rm=TRUE,...) {
  
  # Checking that 'x' is a zoo object
  if ( !is.zoo(x) ) stop("Invalid argument: 'class(x)' must be in c('zoo', 'xts')")
  
  monthlyfunction.zoo(x=x, FUN=FUN, na.rm=na.rm, ...)
  
} # 'monthlyfunction.default' end


################################################################################
# Author : Mauricio Zambrano-Bigiarini                                         #
################################################################################
# Started: 25-Jul-2011                                                         #
# Updates: 08-Aug-2011                                                         #
#          05-Jun-2012                                                         #
#          08-Apr-2013                                                         #
#          28-Nov-2015                                                         #
################################################################################
monthlyfunction.zoo <- function(x, FUN, na.rm=TRUE,...) {
  
  # Checking that the user provied a valid argument for 'FUN'
  if (missing(FUN)) stop("Missing argument: 'FUN' must be provided")
  
  # Checking the user provide a valid value for 'x'
  if (sfreq(x) %in% c("quarterly", "annual"))
    stop("Invalid argument: 'x' is not a sub-daily, daily, weekly or monthly ts. 'x' is a ", sfreq(x), " ts" )
  
  # Monthly index for 'x'
  dates  <- time(x)
  m      <- as.numeric(format( dates, "%m" ))
  months <- factor( month.abb[m], levels=unique(month.abb[m]) )
  
  # 'as.numeric' is necessary for being able to change the names to the output
  totals <- aggregate(x, by= months, FUN=FUN, na.rm= na.rm, ... ) 
  
  # Replacing the NaNs by 'NA.
  # NaN's are obtained when using the FUN=mean with complete NA values
  nan.index          <- which(is.nan(totals))
  if ( length(nan.index) > 0 )  totals[ nan.index] <- NA
  
  # Replacing all the Inf and -Inf by NA's
  # min(NA:NA, na.rm=TRUE) == Inf  ; max(NA:NA, na.rm=TRUE) == -Inf
  inf.index <- which(is.infinite(totals))
  if ( length(inf.index) > 0 ) totals[inf.index] <- NA
  
  # Giving meaningful names to the output
  if ( (is.matrix(x)) | (is.data.frame(x)) ) {
    totals <- t(totals) # For having the months' names as column names
    colnames(totals) <- levels(months)
  } #IF end
  
  return(totals)
  
} # 'monthlyfunction.zoo' end



################################################################################
# Author : Mauricio Zambrano-Bigiarini                                         #
################################################################################
# Started: 25-Jul-2011                                                         #
# Updates: 08-Aug-2011                                                         #
#          29-May-2013 ; 03-Jun-2013                                           #
################################################################################
# 'dates'   : "numeric", "factor", "Date" indicating how to obtain the
#             dates for correponding to the 'sname' station
#             If 'dates' is a number, it indicates the index of the column in
#                'x' that stores the dates
#             If 'dates' is a factor, it have to be converted into 'Date' class,
#                using the date format  specified by 'date.fmt'
#             If 'dates' is already of Date class, the following line verifies that
#                the number of days in 'dates' be equal to the number of element in the
#                time series corresponding to the 'st.name' station
# 'date.fmt': format in which the dates are stored in 'dates'.
#             ONLY required when class(dates)=="factor" or "numeric"
# 'out.type': string that define the desired type of output. Possible values are
#             -) "data.frame": a data.frame, with 12 columns representing the months,
#                              and as many rows as stations are included in 'x'
#             -) "db"        : a data.frame, with 4 colums will be produced.
#                              The first column stores the ID of the station
#                              The second column stores the Year,
#                              The third column stores the ID of the station,
#                              The fourth column contains the monthly value corresponding to the year specified in the second column
# 'verbose'      : logical; if TRUE, progress messages are printed
monthlyfunction.data.frame <- function(x, FUN, na.rm=TRUE,
                                       dates=1, date.fmt="%Y-%m-%d",
                                       out.type="data.frame",
                                       verbose=TRUE,...) {
  
  # Checking that the user provied a valid argument for 'out.type'
  if (is.na(match( out.type, c("data.frame", "db") ) ) )
    stop("Invalid argument: 'out.type' must be in c('data.frame', 'db'")
  
  # Checking that the user provied a valid argument for 'FUN'
  if (missing(FUN)) stop("Missing argument: 'FUN' must be provided")
  
  # Checking that the user provied a valid argument for 'dates'
  if (is.na(match(class(dates), c("numeric", "factor", "Date"))))
    stop("Invalid argument: 'dates' must be of class 'numeric', 'factor', 'Date'")
  
  # If 'dates' is a number, it indicates the index of the column of 'x' that stores the dates
  # The column with dates is then substracted form 'x' for easening the further computations
  if ( class(dates) == "numeric" ) {
    tmp   <- dates
    dates <- as.Date(x[, dates], format= date.fmt)
    x     <- x[-tmp]
  }  # IF end
  
  # If 'dates' is a factor, it have to be converted into 'Date' class,
  # using the date format  specified by 'date.fmt'
  if ( class(dates) == "factor" ) dates <- as.Date(dates, format= date.fmt)
  
  # If 'dates' is already of Date class, the following line verifies that
  # the number of days in 'dates' be equal to the number of element in the
  # time series corresponding to the 'st.name' station
  if ( ( class(dates) == "Date") & (length(dates) != nrow(x) ) )
    stop("Invalid argument: 'length(dates)' must be equal to 'nrow(x)'")
  
  # Transforming 'x' into a zoo object
  x <- zoo(x, dates)
  
  ##############################################################################
  if (out.type == "data.frame") {
    
    monthlyfunction.zoo(x=x, FUN=FUN, na.rm=na.rm, ...)
    
  } else if (out.type == "db") {
    
    # Amount of stations in 'x'
    nstations <- ncol(x)
    
    # ID of all the stations in 'x'
    snames <- colnames(x)
    
    if (is.null(snames)) snames <- paste("V", 1:nstations, sep="")
    
    # Computing the Starting and Ending Year of the analysis
    Starting.Year <- as.numeric(format(start(x), "%Y"))
    Ending.Year   <- as.numeric(format(end(x), "%Y"))
    
    # Amount of Years belonging to the desired period
    nyears <- Ending.Year - Starting.Year + 1
    
    # Computing the numeric index of the resulting months
    month.index <- unique(as.numeric(format( time(x), "%m" )))
    
    # Amount of different months belonging to the desired period
    nmonths <- length(month.index)
    
    # Total amount of months belonging to the desired period
    totalmonths <- nmonths*nyears
    
    # Creating a vector with the names of the field that will be used for storing the results
    field.names <- c("StationID", "Year", "Month", "Value" )
    
    # Creating the data.frame that will store the computed averages for each subcatchment
    z <- as.data.frame(matrix(data = NA, nrow = totalmonths*nstations, ncol = 4, byrow = TRUE, dimnames = NULL) )
    
    for (j in 1:nstations) {
      
      if (verbose) message( "[ Station: ", format(snames[j], width=10, justify="left"),
                            " : ", format(j, width=3, justify="left"), "/",
                            nstations, " => ",
                            format(round(100*j/nstations,2), width=6, justify="left"),
                            "% ]" )
      
      # Computing the annual values
      tmp <- monthlyfunction.default(x= x[,j], FUN=FUN, na.rm=na.rm, ...)
      
      # Putting the annual/monthly values in the output data.frame
      # The first column of 'x' corresponds to the Year
      row.ini <- (j-1)*totalmonths + 1
      row.fin <-  j*totalmonths
      
      z[row.ini:row.fin, 1] <- snames[j] # it is automatically repeted 'totalmonths' times
      z[row.ini:row.fin, 2] <- rep(Starting.Year:Ending.Year, each=nmonths)
      z[row.ini:row.fin, 3] <- month.abb[month.index]
      z[row.ini:row.fin, 4] <- tmp
      
    } # FOR end
    
    colnames(z) <- field.names
    
    return( z )
    
  } # ELSE end
  
} #'monthlyfunction.data.frame' END


################################################################################
# Author : Mauricio Zambrano-Bigiarini                                         #
################################################################################
# Started: 25-Jul-2011                                                         #
# Updates: 08-Aug-2011                                                         #
#          29-May-2013                                                         #
################################################################################
monthlyfunction.matrix <- function(x, FUN, na.rm=TRUE,
                                   dates=1, date.fmt="%Y-%m-%d",
                                   out.type="data.frame",
                                   verbose=TRUE,...) {
  x <- as.data.frame(x)
  #NextMethod("monthlyfunction")
  monthlyfunction.data.frame(x=x, FUN=FUN, na.rm=na.rm,
                             dates=dates, date.fmt=date.fmt,
                             out.type=out.type,
                             verbose=verbose,...)
  
} # 'monthlyfunction.matrix' END  

# File sfreq.R
# Part of the hydroTSM R package, http://www.rforge.net/hydroTSM/ ; 
#                                 http://cran.r-project.org/web/packages/hydroTSM/
# Copyright 2009-2013 Mauricio Zambrano-Bigiarini
# Distributed under GPL 2 or later

################################################################################
# sfreq: Sampling frequency of a ts/zoo object                                 #
################################################################################
# This function generates a table indicating the number of days                #
# with information (<>NA's) within a data.frame                                #
################################################################################
# Author : Mauricio Zambrano-Bigiarini                                         #
################################################################################
# Started: 13-May-2009                                                         #
# Updates: Mar 2009                                                            #
#          Nov 2010                                                            #
#          Apr 2011 ; 09-Aug-2011                                              #
#          18-Oct-2012                                                         #
#          29-May-2013                                                         #
################################################################################
sfreq <- function(x, min.year=1800) {
  
  # Checking that 'class(x)'
  valid.class <- c("xts", "zoo")    
  if (length(which(!is.na(match(class(x), valid.class )))) <= 0) 
    stop("Invalid argument: 'x' must be in c('xts', 'zoo')" )
  
  out <- periodicity(x)$scale # xts::periodicity
  
  if (out == "yearly") out <- "annual"
  
  return(out)
  
} # 'sfreq' END


### Weekly function ###
weeklyfunction <- function(x, FUN, na.rm=TRUE,...) {
  
  # Checking that the user provied a valid argument for 'FUN'
  if (missing(FUN)) stop("Missing argument: 'FUN' must be provided")
  
  # Checking the user provide a valid value for 'x'
  if (sfreq(x) %in% c("quarterly", "annual"))
    stop("Invalid argument: 'x' is not a sub-daily, daily, weekly or monthly ts. 'x' is a ", sfreq(x), " ts" )
  
  # Monthly index for 'x'
  dates  <- time(x)
  m      <- as.numeric(format( dates, "%W" ))

  # 'as.numeric' is necessary for being able to change the names to the output
  totals <- aggregate(x, by= m, FUN=FUN, na.rm= na.rm, ... ) 
  
  # Replacing the NaNs by 'NA.
  # NaN's are obtained when using the FUN=mean with complete NA values
  nan.index          <- which(is.nan(totals))
  if ( length(nan.index) > 0 )  totals[ nan.index] <- NA
  
  # Replacing all the Inf and -Inf by NA's
  # min(NA:NA, na.rm=TRUE) == Inf  ; max(NA:NA, na.rm=TRUE) == -Inf
  inf.index <- which(is.infinite(totals))
  if ( length(inf.index) > 0 ) totals[inf.index] <- NA
  
  return(totals)
  
}

simpson <- function(y, a, b, n = 100) {
  # numerical integral of y from a to b
  # using Simpson's rule with n subdivisions
  #
  # y is a function of a single variable
  # we assume a < b and n is a positive even integer
  
  n <- max(c(2*(n %/% 2), 4),na.rm=T)
  h <- (b-a)/n
  x.vec1 <- seq(a+h, b-h, by = 2*h)
  x.vec2 <- seq(a+2*h, b-2*h, by = 2*h)
  f.vec1 <- y[x.vec1]
  f.vec2 <- y[x.vec2]
  S <- h/3*(y[a] + y[b] + 4*sum(f.vec1,na.rm = T) + 2*sum(f.vec2,na.rm = T))
  return(S)
}
eudis = function(x, y) { sqrt( sum( (x-y)^2 ) ) } # define Euclidean distance

trapezoid = function(x, y) 
{ # computes the integral of y with respect to x using trapezoidal integration. 
  idx = 2:length(x)
  return (as.double( (x[idx] - x[idx-1]) %*% (y[idx] + y[idx-1])) / 2)
}

# Jet Colors
jet.colors <-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

myvis.gam <- function (x, view = NULL, cond = list(), n.grid = 30, too.far = 0, 
                       col = NA, color = "heat", contour.col = NULL, se = -1, type = "link", 
                       plot.type = "persp", zlim = NULL, nCol = 50, ...) 
{
  fac.seq <- function(fac, n.grid) {
    fn <- length(levels(fac))
    gn <- n.grid
    if (fn > gn) 
      mf <- factor(levels(fac))[1:gn]
    else {
      ln <- floor(gn/fn)
      mf <- rep(levels(fac)[fn], gn)
      mf[1:(ln * fn)] <- rep(levels(fac), rep(ln, fn))
      mf <- factor(mf, levels = levels(fac))
    }
    mf
  }
  dnm <- names(list(...))
  v.names <- names(x$var.summary)
  if (is.null(view)) {
    k <- 0
    view <- rep("", 2)
    for (i in 1:length(v.names)) {
      ok <- TRUE
      if (is.matrix(x$var.summary[[i]])) 
        ok <- FALSE
      else if (is.factor(x$var.summary[[i]])) {
        if (length(levels(x$var.summary[[i]])) <= 1) 
          ok <- FALSE
      }
      else {
        if (length(unique(x$var.summary[[i]])) == 1) 
          ok <- FALSE
      }
      if (ok) {
        k <- k + 1
        view[k] <- v.names[i]
      }
      if (k == 2) 
        break
    }
    if (k < 2) 
      stop("Model does not seem to have enough terms to do anything useful")
  }
  else {
    if (sum(view %in% v.names) != 2) 
      stop(gettextf("view variables must be one of %s", 
                    paste(v.names, collapse = ", ")))
    for (i in 1:2) if (!inherits(x$var.summary[[view[i]]], 
                                 c("numeric", "factor"))) 
      stop("Don't know what to do with parametric terms that are not simple numeric or factor variables")
  }
  ok <- TRUE
  for (i in 1:2) if (is.factor(x$var.summary[[view[i]]])) {
    if (length(levels(x$var.summary[[view[i]]])) <= 1) 
      ok <- FALSE
  }
  else {
    if (length(unique(x$var.summary[[view[i]]])) <= 1) 
      ok <- FALSE
  }
  if (!ok) 
    stop(gettextf("View variables must contain more than one value. view = c(%s,%s).", 
                  view[1], view[2]))
  if (is.factor(x$var.summary[[view[1]]])) 
    m1 <- fac.seq(x$var.summary[[view[1]]], n.grid)
  else {
    r1 <- range(x$var.summary[[view[1]]])
    m1 <- seq(r1[1], r1[2], length = n.grid)
  }
  if (is.factor(x$var.summary[[view[2]]])) 
    m2 <- fac.seq(x$var.summary[[view[2]]], n.grid)
  else {
    r2 <- range(x$var.summary[[view[2]]])
    m2 <- seq(r2[1], r2[2], length = n.grid)
  }
  v1 <- rep(m1, n.grid)
  v2 <- rep(m2, rep(n.grid, n.grid))
  newd <- data.frame(matrix(0, n.grid * n.grid, 0))
  for (i in 1:length(x$var.summary)) {
    ma <- cond[[v.names[i]]]
    if (is.null(ma)) {
      ma <- x$var.summary[[i]]
      if (is.numeric(ma)) 
        ma <- ma[2]
    }
    if (is.matrix(x$var.summary[[i]])) 
      newd[[i]] <- matrix(ma, n.grid * n.grid, ncol(x$var.summary[[i]]), 
                          byrow = TRUE)
    else newd[[i]] <- rep(ma, n.grid * n.grid)
  }
  names(newd) <- v.names
  newd[[view[1]]] <- v1
  newd[[view[2]]] <- v2
  if (type == "link") 
    zlab <- paste("linear predictor")
  else if (type == "response") 
    zlab <- type
  else stop("type must be \"link\" or \"response\"")
  fv <- predict.gam(x, newdata = newd, se.fit = TRUE, type = type)
  z <- fv$fit
  if (too.far > 0) {
    ex.tf <- exclude.too.far(v1, v2, x$model[, view[1]], 
                             x$model[, view[2]], dist = too.far)
    fv$se.fit[ex.tf] <- fv$fit[ex.tf] <- NA
  }
  if (is.factor(m1)) {
    m1 <- as.numeric(m1)
    m1 <- seq(min(m1) - 0.5, max(m1) + 0.5, length = n.grid)
  }
  if (is.factor(m2)) {
    m2 <- as.numeric(m2)
    m2 <- seq(min(m1) - 0.5, max(m2) + 0.5, length = n.grid)
  }
  if (se <= 0) {
    old.warn <- options(warn = -1)
    av <- matrix(c(0.5, 0.5, rep(0, n.grid - 1)), n.grid, 
                 n.grid - 1)
    options(old.warn)
    max.z <- max(z, na.rm = TRUE)
    z[is.na(z)] <- max.z * 10000
    z <- matrix(z, n.grid, n.grid)
    surf.col <- t(av) %*% z %*% av
    surf.col[surf.col > max.z * 2] <- NA
    if (!is.null(zlim)) {
      if (length(zlim) != 2 || zlim[1] >= zlim[2]) 
        stop("Something wrong with zlim")
      min.z <- zlim[1]
      max.z <- zlim[2]
    }
    else {
      min.z <- min(fv$fit, na.rm = TRUE)
      max.z <- max(fv$fit, na.rm = TRUE)
    }
    surf.col <- surf.col - min.z
    surf.col <- surf.col/(max.z - min.z)
    surf.col <- round(surf.col * nCol)
    con.col <- 1
    if (color == "heat") {
      pal <- heat.colors(nCol)
      con.col <- 3
    }
    else if (color == "topo") {
      pal <- topo.colors(nCol)
      con.col <- 2
    }
    else if (color == "cm") {
      pal <- cm.colors(nCol)
      con.col <- 1
    }
    else if (color == "terrain") {
      pal <- terrain.colors(nCol)
      con.col <- 2
    }
    else if (color == "gray" || color == "bw") {
      pal <- gray(seq(0.1, 0.9, length = nCol))
      con.col <- 1
    }
    ### customized here
    else if (color == 'jet') {
      pal <- jet.colors(nCol)
      con.col = 1
    }
    else if (color == "viridis"){
      require(viridis)
      pal <- viridis::viridis(nCol)
      con.col = 1
    }
    else if (color == "plasma"){
      require(viridis)
      pal <- viridis::plasma(nCol)
      con.col = 1
    }
    ####
    else stop("color scheme not recognised")
    if (is.null(contour.col)) 
      contour.col <- con.col
    surf.col[surf.col < 1] <- 1
    surf.col[surf.col > nCol] <- nCol
    if (is.na(col)) 
      col <- pal[as.array(surf.col)]
    z <- matrix(fv$fit, n.grid, n.grid)
    if (plot.type == "contour") {
      stub <- paste(ifelse("xlab" %in% dnm, "", ",xlab=view[1]"), 
                    ifelse("ylab" %in% dnm, "", ",ylab=view[2]"), 
                    ifelse("main" %in% dnm, "", ",main=zlab"), ",...)", 
                    sep = "")
      if (color != "bw") {
        txt <- paste("image(m1,m2,z,col=pal,zlim=c(min.z,max.z)", 
                     stub, sep = "")
        eval(parse(text = txt))
        txt <- paste("contour(m1,m2,z,col=contour.col,zlim=c(min.z,max.z)", 
                     ifelse("add" %in% dnm, "", ",add=TRUE"), ",...)", 
                     sep = "")
        eval(parse(text = txt))
      }
      else {
        txt <- paste("contour(m1,m2,z,col=1,zlim=c(min.z,max.z)", 
                     stub, sep = "")
        eval(parse(text = txt))
      }
    }
    else {
      stub <- paste(ifelse("xlab" %in% dnm, "", ",xlab=view[1]"), 
                    ifelse("ylab" %in% dnm, "", ",ylab=view[2]"), 
                    ifelse("zlab" %in% dnm, "", ",zlab=zlab"), ",...)", 
                    sep = "")
      if (color == "bw") {
        op <- par(bg = "white")
        txt <- paste("persp(m1,m2,z,col=\"white\",zlim=c(min.z,max.z) ", 
                     stub, sep = "")
        eval(parse(text = txt))
        par(op)
      }
      else {
        txt <- paste("persp(m1,m2,z,col=col,zlim=c(min.z,max.z)", 
                     stub, sep = "")
        eval(parse(text = txt))
      }
    }
  }
  else {
    if (color == "bw" || color == "gray") {
      subs <- paste("grey are +/-", se, "s.e.")
      lo.col <- "gray"
      hi.col <- "gray"
    }
    else {
      subs <- paste("red/green are +/-", se, "s.e.")
      lo.col <- "green"
      hi.col <- "red"
    }
    if (!is.null(zlim)) {
      if (length(zlim) != 2 || zlim[1] >= zlim[2]) 
        stop("Something wrong with zlim")
      min.z <- zlim[1]
      max.z <- zlim[2]
    }
    else {
      max.z <- max(fv$fit + fv$se.fit * se, na.rm = TRUE)
      min.z <- min(fv$fit - fv$se.fit * se, na.rm = TRUE)
      zlim <- c(min.z, max.z)
    }
    z <- fv$fit - fv$se.fit * se
    z <- matrix(z, n.grid, n.grid)
    if (plot.type == "contour") 
      warning("sorry no option for contouring with errors: try plot.gam")
    stub <- paste(ifelse("xlab" %in% dnm, "", ",xlab=view[1]"), 
                  ifelse("ylab" %in% dnm, "", ",ylab=view[2]"), ifelse("zlab" %in% 
                                                                         dnm, "", ",zlab=zlab"), ifelse("sub" %in% dnm, 
                                                                                                        "", ",sub=subs"), ",...)", sep = "")
    txt <- paste("persp(m1,m2,z,col=col,zlim=zlim", ifelse("border" %in% 
                                                             dnm, "", ",border=lo.col"), stub, sep = "")
    eval(parse(text = txt))
    par(new = TRUE)
    z <- fv$fit
    z <- matrix(z, n.grid, n.grid)
    txt <- paste("persp(m1,m2,z,col=col,zlim=zlim", ifelse("border" %in% 
                                                             dnm, "", ",border=\"black\""), stub, sep = "")
    eval(parse(text = txt))
    par(new = TRUE)
    z <- fv$fit + se * fv$se.fit
    z <- matrix(z, n.grid, n.grid)
    txt <- paste("persp(m1,m2,z,col=col,zlim=zlim", ifelse("border" %in% 
                                                             dnm, "", ",border=hi.col"), stub, sep = "")
    eval(parse(text = txt)) 
    }
}


#### Subsetting functions - pairwise ####
# Subsetting function
pairwiseRandomCellPick <- function(sout,cluster=T, pr=T){
  if(cluster) {
    library(doParallel)
    cl <- makeCluster(getOption("cl.cores", 2))
    registerDoParallel(cl)
  }
  ### Picks a random combination of cells for biodiversity values
  res <- data.frame()
  # Do a progress bar
  if(pr) pb <- txtProgressBar(min = 0, max = length(unique(unique(sout$SS))), style = 3)
  
  for(study in unique(sout$SS)){
    sub <- subset(sout,SS==study)
    if(nrow(sub) < 2) next() # To the rare case of only a single study site
    
    # Kickout those where cells.x and cell.y are identical
    sub <- sub[which(sub$cells.x!=sub$cells.y),]
    if(nrow(sub)==0) next()
    # Get unique pairwise combinations
    sub <- transform(sub, cells.u = paste(pmin(cells.x,cells.y), pmax(cells.x,cells.y), sep="_")) 
    
    if(cluster){
      df2 <- parLapply(cl,split(sub, sub$cells.u),
                       function(subdf) subdf[sample(1:nrow(subdf), 1),]
      )
    } else {
      df2 <- lapply(split(sub, sub$cells.u),
                       function(subdf) subdf[sample(1:nrow(subdf), 1),]
      )
    }
    
    res <- rbind(res, do.call('rbind', df2) )
    if(pr) setTxtProgressBar(pb, getTxtProgressBar(pb)+1)
  }
  if(pr) close(pb) # close progressbar
  if(cluster) stopCluster(cl)
  return(res)
}


# Get a Single permutation round of the subdiagonal 
pairwiseSubDiagonalPermutation <- function(sout,pr=T) {
  res <- data.frame()
  # Do a progress bar
  if(pr) pb <- txtProgressBar(min = 0, max = length(unique(unique(sout$SS))), style = 3)
  
  # Take one value per cell
  for(study in unique(sout$SS)){
    sub <- subset(sout,SS==study)

    # Kickout those where cells.x and cell.y are identical
    sub <- sub[which(sub$cells.x!=sub$cells.y),]
    if(nrow(sub)==0) next()
    # Get unique pairwise combinations
    sub <- transform(sub, cells.u = paste(pmin(cells.x,cells.y), pmax(cells.x,cells.y), sep="_")) 
    
    # Sample from combinations
    df2 <- do.call('rbind', lapply(split(sub, sub$cells.u),
                  function(subdf) subdf[sample(1:nrow(subdf), 1),]) )
    
    # if(nrow(df2)>2) {
    #   df2$SSBS_c <- paste0(df2$SSBS_x,"_",df2$SSBS_y) # variable to merge on
    #   # Build pairwise matrix
    #   m <- reshape2::acast(data=df2,SSBS_y~SSBS_x,value.vdar = "value",fill = NA) # Fill with zero assuming absence
    #   # Permute / rotate the matrix
    #   #p <- sample.int(dim(m)[1])
    #   #m <- m[p, p]
    #   # take diagonal. A bit fishy to get the row.names
    #   diag(m) <- -9999 # Take diagonal
    #   m <- reshape2::melt(m)
    #   m <- m[which(m$value == -9999),];m$SSBS_c <- paste0(m$Var1,"_",m$Var2);m$Var1 <- NULL;m$Var2 <- NULL;m$value <- NULL
    #   # Merge back and remove merging columns
    #   df2 <- merge.data.frame(m,df2,by="SSBS_c")
    #   df2$SSBS_c <- NULL; df2$cells.u <- NULL
    # }
    # df2$cells.u <- NULL
    res <- rbind(res,df2)
    if(pr) setTxtProgressBar(pb, getTxtProgressBar(pb)+1)
    
  }
  if(pr) close(pb) # close progressbar
  return(res)
}

singleRandomCellPick <- function(sout){
  ### Picks a random combination of cells for biodiversity values
  res <- data.frame()
  # Do a progress bar
  pb <- txtProgressBar(min = 0, max = length(unique(unique(sout$SS))), style = 3)
  
  for(study in unique(sout$SS)){
    sub <- subset(sout,SS==study)
    # For each cell individual cell pick a random entry
    df2 <- lapply(split(sub, sub$cells),
                  function(subdf) subdf[sample(1:nrow(subdf), 1),]
    )
    res <- rbind(res, do.call('rbind', df2) )
    setTxtProgressBar(pb, getTxtProgressBar(pb)+1)
  }
  close(pb) # close progressbar
  return(res)
}

#' Theme to make creating journal ready figures easier
#'
#'
#' @param base_size Font size
#' @param base_family Font used for everything, unless axis fonts etc. are defined
#' @param plot_grid = Do you want a grid? (all = major/minor, none = no grid, major = major grid)
#' @param lines_lwd = width of lines
#' @param title_size = font size of title
#' @param legend_size = font size of legend
#' @param legend_title = Do you want the legend title? (TRUE/FALSE)
#' @param bg_col = background colour
#' @param title_font = font use for title
#' @param base_col  = font colour
#' @param horz_grid = Do you want a horizontal grid?
#' @param bord_size = width of a rectangular border
#' @param alpha_leg = opacity of the legend. 0 = totally transparent
#' @param strip_bg = colour background for facets
#' @param grid_thick = A multiplier to apply to the grid lines.
#' @param grid_type = Grid type. Default is a solid line
#' @param ticks_xy = Do you want ticks on the x or y axis? "x" = x-axis only, "y" = y-axis only, "xy" = both axes.
#' @param grid_cols = Colour of the grid. 2 element vector. First element is major grid colour. If only one element, the first will be used for minor grid.
#' @return ggplot2 theme
#'
#' @export


theme_pub <- function(base_size = 11,
                      base_family = "Gill Sans MT",
                      plot_grid = "none",
                      lines_lwd = 0.50,
                      bg_col = "white",
                      axis_lines = TRUE,
                      horz_grid = ifelse(plot_grid == "none", FALSE, TRUE),
                      vert_grid = ifelse(plot_grid == "none", FALSE, TRUE),
                      ticks_type = "outer",
                      ticks_xy  = "xy",
                      alpha_leg = 0.1,
                      bord_size = lines_lwd,
                      legend_bg = "white",
                      strip_bg = "white",
                      title_size = base_size * 1.2,
                      legend_size = base_size * 0.9,
                      grid_thick = 1,
                      grid_type = "solid",
                      grid_cols = c("grey10", "grey30"),
                      legend_title = TRUE,
                      plot_box = TRUE) {
  minor_grid = ifelse(plot_grid == "all", TRUE, FALSE)
  
  if (legend_title)
    title_element <- element_text()
  else
    title_element <- element_blank()
  if (axis_lines)
    axis_line =  element_line(size = ifelse(axis_lines, grid::unit(lines_lwd, "mm"), 0),
                              color = "black")
  else
    axis_line <- element_blank()
  
  if (plot_box)
    panel_border <-
      ggplot2::element_rect(colour = "black",
                            fill = NA,
                            size = bord_size)
  else
    panel_border <- element_blank()
  
  
  theme_bw() +
    ggplot2::theme(
      text = ggplot2::element_text(family = base_family, size = base_size),
      axis.line =  axis_line,
      axis.ticks.length = grid::unit(ifelse(ticks_type == "outer", 0.15,-0.15), "cm"),
      axis.ticks.x =  element_line(
        size = ifelse(
          stringr::str_detect(ticks_xy, "x"),
          grid::unit(lines_lwd, "cm"),
          0
        ),
        color = "black"
      ),
      axis.ticks.y =  element_line(
        size = ifelse(
          stringr::str_detect(ticks_xy, "y"),
          grid::unit(lines_lwd, "cm") ,
          0
        ),
        color = "black"
      ),
      axis.text.x = ggplot2::element_text(
        size = base_size,
        family = base_family,
        margin = ggplot2::margin(ifelse(ticks_type == "inner", 11, 5), 5, 10, 5, "pt")
      ),
      axis.text.y = ggplot2::element_text(
        size = base_size,
        family = base_family,
        margin = ggplot2::margin(5, ifelse(ticks_type == "inner", 11, 5), 10, 5, "pt")
      ),
      axis.title.y = ggplot2::element_text(
        size =  base_size,
        vjust = 1.5,
        family = base_family
      ),
      axis.title.x = ggplot2::element_text(
        size = base_size,
        vjust = -.5,
        family = base_family
      ),
      panel.background = ggplot2::element_rect(fill = bg_col),
      plot.background = ggplot2::element_rect(fill = bg_col),
      panel.border = panel_border,
      panel.grid.major.x = ggplot2::element_line(
        linetype = grid_type,
        colour = ifelse(vert_grid, grid_cols[1], bg_col),
        size = ifelse(vert_grid, 0.25 * grid_thick, 0)
      ),
      panel.grid.minor.x = ggplot2::element_line(
        linetype = grid_type,
        colour = ifelse(vert_grid, ifelse(minor_grid, grid_cols[2 - (length(grid_cols) == 1)], bg_col), bg_col),
        size = ifelse(vert_grid, 0.15 * grid_thick, 0)
      ),
      panel.grid.major.y = ggplot2::element_line(
        linetype = grid_type,
        colour = ifelse(horz_grid, grid_cols[1], bg_col),
        size = ifelse(horz_grid, 0.25 * grid_thick, 0)
      ),
      panel.grid.minor.y = ggplot2::element_line(
        linetype = grid_type,
        colour = ifelse(horz_grid, ifelse(minor_grid, grid_cols[2 - (length(grid_cols) == 1)], bg_col), bg_col),
        size = ifelse(horz_grid, 0.15 * grid_thick, 0)
      ),
      plot.title = ggplot2::element_text(
        face = "bold",
        vjust = 2,
        size = title_size,
        family = base_family
      ),
      legend.background = ggplot2::element_rect(fill = scales::alpha(legend_bg, alpha_leg)),
      legend.key = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = legend_size, family = base_family),
      legend.title = title_element,
      strip.background =  ggplot2::element_rect(colour = strip_bg, fill = strip_bg),
      strip.text.x = ggplot2::element_text(size = base_size + 1),
      strip.text.y = ggplot2::element_text(size = base_size + 1),
      plot.caption = element_text(hjust = 0)
    )
}


#### ---- GRID arrange #####
grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
  require(grid);require(gridExtra)
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  
  grid.newpage()
  grid.draw(combined)
  
  # return gtable invisibly
  invisible(combined)
  
}


#### GGPLOT 2 helper functions ####
# Stat_n
stat_n_text <-
  function (mapping = NULL, data = NULL, geom = ifelse(text.box, 
                                                       "label", "text"), position = "identity", na.rm = FALSE, show.legend = NA, 
            inherit.aes = TRUE, y.pos = NULL, y.expand.factor = 0.1, 
            text.box = FALSE, alpha = 1, angle = 0, color = "black", 
            family = "", fontface = "plain", hjust = 0.5, label.padding = ggplot2::unit(0.25, 
                                                                                        "lines"), label.r = ggplot2::unit(0.15, "lines"), label.size = 0.25, 
            lineheight = 1.2, size = 4, vjust = 0.5, ...) 
  {
    geom <- match.arg(geom, c("label", "text"))
    params <- list(y.pos = y.pos, y.expand.factor = y.expand.factor, 
                   alpha = alpha, angle = angle, color = color, family = family, 
                   fontface = fontface, hjust = hjust, lineheight = lineheight, 
                   size = size, vjust = vjust)
    if (geom == "label") {
      params <- c(params, list(label.padding = label.padding, 
                               label.r = label.r, label.size = label.size, na.rm = na.rm, 
                               ...))
    }
    else {
      params <- c(params, na.rm = na.rm, ...)
    }
    ggplot2::layer(stat = StatNText, data = data, mapping = mapping, 
                   geom = geom, position = position, show.legend = show.legend, 
                   inherit.aes = inherit.aes, params = params)
  }

#### Function to make plot transparent ####
# Transparent
transparent <- theme(plot.background = element_rect(fill = "transparent",colour = NA),panel.background = element_rect(fill = "transparent",colour = NA),
                     axis.text = element_text(colour="white",size = 18),axis.text.x =element_text(colour="white"),
                     axis.text.y =element_text(colour="white"),axis.ticks = element_line(colour = "white"),
                     axis.line = element_line(colour="white"),axis.title = element_text(colour="white"),
                     legend.key = element_rect(fill=NA,colour = NA),legend.background = element_rect(fill = NA,colour=NA),legend.text = element_text(colour="white"),
                     title = element_text(colour="white",size=18),strip.text = element_text(colour="white",size=18),
                     strip.background = element_rect(fill=NA),
                     panel.grid.major = element_blank(),panel.grid.minor = element_blank()
) 


#### Do a qucik MapView lookup ####
spatialLookup <- function(spp){
  paste0("Coordinates (Lon|Lat): ",spp$Longitude,", ",spp$Latitude)
  library(mapview)
  library(sp)
  coordinates(spp) <- ~Longitude + Latitude
  proj4string(spp) <- CRS("+proj=longlat +datum=WGS84")
  mapView(spp) 
}

#### Filter out permanent water ####
# Returns site ID and if site falls into water (MODE) or not
waterFilter <- function(type="mode"){
  library(jsonlite);require(dplyr)
  if(type == "mode"){
    wa <- fromJSON("resSaves/PREDICTS_Permanentwater_mode.geojson.json",flatten=T)$features
    wa <- wa %>% dplyr::select(properties.SS,properties.SSBS,properties.mode) %>% 
      dplyr::rename(SS = properties.SS, SSBS = properties.SSBS, WaterMode = properties.mode)
    return(wa)
  } else {
    print('Not implemented!')
  }
}


#==============================================================================
#### 'Probability of Interspecific Encounter (PIE) ####
# PIE <- function(x){
#   pie_part1 <- sum(x) / (sum(x) - 1)
#   pie_part2 <- (x / sum(x)) ^ 2
#   pie_part3 <- ( 1 - sum(pie_part2) )
#   pie_part1 * pie_part3
# }
# A= c(rep(1,4),rep(2,4), rep(3,4),rep(4,4),rep(5,4))
# B= c(rep(1,16),2,3,4,5)
# A = table(A)
# B = table(B)
# PIE(B)
SiteHurlbertsPie <- function(diversity, site.abundance){
  yarg:::.Log("Computing Hurlberts PIE")
  sd <- rep(NA, length(site.abundance))
  sd[site.abundance] <- tapply(diversity$Measurement[diversity$Is_abundance], 
                               droplevels(diversity$SSS[diversity$Is_abundance]), function(m) {
                                 if (any(m > 2)) {
                                   #Hurlbert's (1971) Probability of Interspecific Encounter (PIE):
                                   #'PIE = (N/(N-1)) * (1 - sum(p^2))
                                   #'where:
                                   #'p = the propotion a specific taxon composes of a sample
                                   #'N = the total abundance of organisms
                                   #'
                                   #'This measurement is equivalent to Simpson's Index but includes a correction
                                   #'factor base on the total abundance of organisms.
                                   
                                   pie_part1 <- sum(m) / (sum(m) - 1)
                                   pie_part2 <- (m / sum(m)) ^ 2
                                   pie_part3 <- ( 1 - sum(pie_part2) )
                                   
                                   return(pie_part1 * pie_part3)
                                   #   1/(sum((m/sum(m))^2))
                                 } else NA
                               })
  return(sd)
}

# Correct Sampling effort figure
CorrectSamplingEffort <- function(diversity){
  ## Namespace - From yarg package
  missing <- is.na(diversity$Sampling_effort)
  myLog("Correcting ", sum(missing), " missing sampling effort values\n")
  diversity$Sampling_effort[missing] <- 1
  myLog("Rescaling sampling effort\n")
  diversity$Sampling_effort <- do.call("c", tapply(diversity$Sampling_effort, 
                                                   diversity$SS, function(se) return(se/max(se))))
  sensitive <- diversity$Diversity_metric_is_effort_sensitive
  myLog("Correcting ", sum(sensitive), " values for sensitivity to sampling", 
       "effort\n")
  diversity$Measurement[sensitive] <- diversity[sensitive, 
                                                "Measurement"]/diversity[sensitive, "Sampling_effort"]
  return(diversity)
}

# SiteMetrics
source('https://raw.githubusercontent.com/timnewbold/predicts-demo/master/predictsFunctions/R/SiteMetrics.R')
