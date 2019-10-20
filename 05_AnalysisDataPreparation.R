library(tidyverse)
library(xts)
source('000_HelperFunction.R')

# Load prepared sites data
sites <- readRDS("resSaves/PREDICTS_allsites.rds")
# Correct for certain issues and remove problematic studies
# A wrong classified study
sites$TGrouping[which(sites$SS == "DL1_2008__MacipRios 1" )] <- "Reptilia"
sites <- subset(sites, !SS %in% c("DL1_2012__CalvinoCancela 1","DL1_2012__CalvinoCancela 2") ) # Remove that one problematic study
sites <- subset(sites, Max_linear_extent <= 3000) #Remove sites with  very large extents

# Metrics
sites$logabund <- log10(sites$Total_abundance+1)
sites$logsimpson <- log10(sites$Simpson_diversity)
sites$Richness_rarefied <- round(sites$Richness_rarefied,0) # Round estimates as they are based on 1000 random drawn samples
sites$LogPopDen <- log1p(sites$PopDen)
sites$asPIE <- asin(sqrt(sites$PIE))

# Refactorize land use
sites$LCLU <- as.character(sites$Predominant_habitat)
sites$LCLU[which(sites$LCLU %in% c("Mature secondary vegetation","Intermediate secondary vegetation","Young secondary vegetation","Secondary vegetation (indeterminate age)"))] <- "Secondary vegetation"
sites <- sites[which(sites$LCLU != "Cannot decide"),]
# Reclass
sites$LCLU <- factor(sites$LCLU,levels = c("Primary forest","Primary non-forest","Secondary vegetation","Plantation forest","Cropland","Pasture","Urban"))
sites$LCLU2 <- sites$LCLU
# Replace estimates that are based on abundance values with NA if only occurence has been estimated
# Do this to ensure that abundance-based metrics are not calculated on occurence data per accident
# They should all be NA already, but just to be sure
sites[which(sites$Diversity_metric_type %in% c("Occurrence","Species richness")),c("logabund","asPIE","PIE","Richness_rarefied")] <- NA

## Assign taxonomic grouping 
# Classify based on order ("flying" to "non-flying"). This is usually reflected by differing sampling methodologies ( pitfall traps vers )
sites$TGrouping[which(sites$Order %in%
                        c("Lepidoptera","Diptera","Hymenoptera","Odonata","Thysanoptera" ) )] <- "Invertebrates.flying"
# Get those where study common taxon  and look up manually
sites$TGrouping[which(sites$SS == "AD1_2009__Farwig 1")] <- "Invertebrates.flying"
sites$TGrouping[which(sites$TGrouping == "Invertebrates")] <- "Invertebrates.ground"

# Reorder taxonomic grouping
sites$TGrouping <- factor(sites$TGrouping, levels = c(
  "Plantae", "Fungi","Invertebrates.ground","Invertebrates.flying", "Amphibia", "Reptilia","Aves", "Mammalia"
),labels = c(
  "Plants","Fungi","Invertebrates.ground","Invertebrates.flying","Amphibians","Reptiles","Birds","Mammals"
)
)
# Check if Taxonomic grouping is unique within studies
stopifnot( any( rowSums( table(sites$SS,sites$TGrouping)>0 ) !=2 ) )

# ----------------------------------------------------- #
# Load BFAST results
res_ls_bfast <- readRDS("resSaves/LS_EVI2_BFast_lm_BOUNDS_LEDAPS_SEASON.rds") %>%  # Load the bounded LEDAPS processed time series
  subset(.,band=="EVI2") %>% dplyr::filter(bestmodel == 1) 

# Get the study id's
res_ls_bfast$SS <- sites$SS[match(res_ls_bfast$SSBS,sites$SSBS)]
# Remove those which HDVe to large sampling extent and where removed earlier
res_ls_bfast <- res_ls_bfast[which(!is.na(res_ls_bfast$SS)),]
# Merge back into sites and clean up
out <- merge.data.frame(sites,res_ls_bfast, by = c("SS","SSBS"),all = TRUE) %>% mutate(Sat = "Landsat")
# Filter out those studies that had no detected break among all their sites
tSS <- out %>% dplyr::group_by(SS) %>% dplyr::summarise( BPS = length(which(Break==1)) ) %>% 
  dplyr::filter(BPS > 0 ) %>% dplyr::select(SS) %>% distinct()
out1 <- subset(out,SS %in% tSS$SS)
out <- out1
rm(out1,res_ls_bfast)

# --- Set ----- #           
# Get only those studies with at least one breakpoint
out <- rename(out,Break_binom = Break)
out$Break_binom[which(is.na(out$Break_binom))] <- "0"# Assume NA values to HDVe no break ?
out$Break_binom <- factor(out$Break_binom)
out$Break_binom <-relevel(out$Break_binom,ref = "0")
# Furthmore classify those with Breakpoints into Pos|Neg breaks
out$Break_direction <- ifelse(out$largest_mag>0,"P","N")
out$Break_direction[which(is.na(out$Break_direction))] <- "S"
out$Break_direction <- factor(out$Break_direction,levels = c("S","N","P"))
out$Break_direction <- relevel(out$Break_direction,ref="S")
# Also calculate the time since disturbance
out$LargeTimeAgo <- as.numeric(as.yearmon(out$Sample_start_earliest)) - out$largest_year

# Has the site recovered from a past disturbance (in trend)
out$Recovered <- factor(ifelse(is.na(out$rec_time),0,1),levels=c(0,1))
out$Recovered[which(out$Break_binom==0)] <- NA # Overwrite non-disturbed sites with NA

# Make sure that reference estimates come from definetly checked time series (so sites with an estimate)!
out <- out[which(!is.na(out$largest_trendbef)),] # Use only those as spatial comparison were we HDVe spec. tested for stable time series

# Check permanent water sites and filter them out.
out <- waterFilter(type = "mode") %>% dplyr::right_join(.,out,by=c("SS","SSBS"))
out <- out %>% dplyr::filter(WaterMode==0)

# Filter out those that are followed / superseded by a full year of missing data
out$na.1ybefore[which(is.na(out$na.1ybefore))] <- 9999 # Placeholder value
out$na.1yafter[which(is.na(out$na.1yafter))] <- 9999 # Placeholder value

# Assign land use
out$LCLU <- factor(out$LCLU,levels = levels(sites$LCLU),labels = c("PV","PNV","SV","PF","CL","PA","UR"))
# New classification
out$LCLU <- as.character(out$LCLU)
out$LCLU[which( out$LCLU %in% c("PV","PNV") )] <- "PV"
out$LCLU[which( out$LCLU %in% c("PF","CL","PA","UR") )] <- "HDV"
out$LCLU <- factor(out$LCLU,levels = c("PV","SV","HDV"),labels = c("PV","SV","HDV"))

# If largest_map_prop overshoots (regression forced lower than zero), set to zero implying 100 % loss
# This has no effect on the grouping per se
out[which(out$largest_mag_prop < -0.99),"largest_mag_prop"] <- -1

# Do a binning for magnitude and time span passed
tmc <- c(min(out$largest_mag_prop,na.rm=T)-0.001,-0.5,-0.25,0,0.25,0.5,max(out$largest_mag_prop,na.rm=T))
out$BinMagn <- cut(out$largest_mag_prop,breaks = tmc,labels = c("< -50%", "-50% <> -25%","-25% <> 0%","0% <> 25%","25% <> 50%","> 50%"))
out$BinMagn <- as.character(out$BinMagn)
out$BinMagn[which(is.na(out$BinMagn))] <- "ND" # Stable/non-disturbed reference
out$BinMagn <- factor(out$BinMagn,levels = c("ND","< -50%", "-50% <> -25%","-25% <> 0%","0% <> 25%","25% <> 50%","> 50%")) # Backtransform to factor
table(out$BinMagn,out$LCLU)
rm(tmc)

# Coarse binner based on terciles
tmc = c( quantile(out$largest_mag_prop[which(out$largest_mag_prop < 0)],c(0,.33,.66)), 0,quantile(out$largest_mag_prop[which(out$largest_mag_prop > 0)],c(.33,.66)),max(out$largest_mag_prop,na.rm = T) );tmc[1] <- tmc[1] - 0.001
#tmc <- c(min(out$largest_mag_prop,na.rm=T)-0.001,-0.25,0,0.25,max(out$largest_mag_prop,na.rm=T))
out$BinMagn2 <- cut(out$largest_mag_prop,breaks = tmc,labels = c("---","--","-","+","++","+++"))
#out$BinMagn <- cut(out$largest_mag_prop,breaks = tmc,labels = c("< -33%", "-33% <> 0%","0% <> 33%","> 33%"))
out$BinMagn2 <- as.character(out$BinMagn2)
out$BinMagn2[which(is.na(out$BinMagn2))] <- "ND" # Stable/non-disturbed reference
out$BinMagn2 <- factor(out$BinMagn2,levels = c("ND","---","--","-","+","++","+++") ) # Relevel
table(out$BinMagn2,out$LCLU)

# Cuts for time
tgc <- c(0,5,10,30) # Time cuts
out$BinTime <- cut(out$LargeTimeAgo,breaks = tgc,labels = c(">0-5y","5-10y",">10y"))
out$BinTime<- as.character(out$BinTime)
out$BinTime[which(is.na(out$BinTime))] <- "ND" # Stable/non-disturbed reference
out$BinTime <- factor(out$BinTime,levels = c("ND",">0-5y","5-10y",">10y")) # Backtransform to factor

# Bin Trend
before = out$largest_trendbef*12
after = out$largest_trendaft*12
trendchange <- after-before

tgt <- c(min(trendchange,na.rm = T)-0.0001,-0.05,-0.01,0,0.01,0.05,max(trendchange,na.rm = T))
out$BinTrend <- cut(trendchange,breaks = tgt,labels = c("Largest NTC","Large NTC","Small NTC","Small PTC","Large PTC","Largest PTC"))
out$BinTrend <- as.character(out$BinTrend)
out$BinTrend[which(is.na(out$BinTrend))] <- "ND" # Stable/non-disturbed reference
out$BinTrend <- factor(out$BinTrend,levels = c("ND","Largest NTC","Large NTC","Small NTC","Small PTC","Large PTC","Largest PTC")) # Backtransform to factor
table(out$BinTrend,useNA = "ifany")
# ---- #
# Only the columns necessary
out <- out %>% dplyr::select(Source_ID,SS,SSS,SSB,SSBS,LCLU,Longitude,Latitude,
                             Sample_start_earliest,Sample_midpoint,Sample_end_latest,Diversity_metric_type,startyear,endyear,
                             Total_abundance, logabund, Species_richness,Richness_rarefied, PIE,asPIE,
                             TGrouping,Koeppen,
                             largest_year,largest_mag_prop,largest_trendbef, largest_trendaft,
                             Break_binom,Break_direction, LargeTimeAgo,Recovered, BinMagn, BinTime, BinTrend
                             )

# Save output for analysis
write_rds(out,'resSaves/PREDICTS_prepared_data.rds')
