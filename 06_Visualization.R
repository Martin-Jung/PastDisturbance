#### Load packages ####
# --------------------------------------------------------- #
library(dplyr)
library(forcats) # Factor reorganization
library(stringr)
library(tidyr);library(tidyverse)
library(lubridate)
library(zoo)
library(lme4)
library(MuMIn)
library(effects)
library(coda)
library(pbkrtest)
library(arm)
library(ggplot2);theme_set(theme_bw()) # simpler theme for the plots
library(extrafont) # font_import() # Run first to import all fonts on system
loadfonts(device="win",quiet = T) # Load fonts
library(scales)
library(ggthemes);library(ggsci)
library(viridisLite);library(viridis)
library(grid);library(gridExtra)
source("000_HelperFunction.R");par.ori <- par(no.readonly = T)
redBlue <- c("#d73027","#fc8d59","#fee090","black","#e0f3f8","#91bfdb","#4575b4")
# --------------------------------------------------------- #

# Load in prepared data 
out <- read_rds('resSaves/PREDICTS_prepared_data.rds')

# Overview
o <- out
# Some summary statistics
myLog("Number of studies: ", length(unique(o$SS)) )
# How much is negative vs positive change?
table(o$Break_binom)
table(o$Break_direction)[2:3]/(sum(table(o$Break_direction)[2:3]) )

# Disturbance magnitude
o %>% dplyr::filter(Break_direction == "N") %>% dplyr::summarise(m = median(largest_mag_prop*100,na.rm=T),s = mad(largest_mag_prop*100,na.rm = T))
o %>% dplyr::filter(Break_direction == "P") %>% dplyr::summarise(m = median(largest_mag_prop*100,na.rm=T),s = mad(largest_mag_prop*100,na.rm = T))

# Trend difference
before = out$largest_trendbef*12
after = out$largest_trendaft*12
o$trendchange <- (after-before)
o$Trendchange_direction = ifelse(o$trendchange < 0,"N","P");o$Trendchange_direction[which(is.na(o$trendchange))] <- "S"
o$Trendchange_direction <- factor(o$Trendchange_direction,c("S","N","P"))
table(o$Trendchange_direction)[2:3]/(sum(table(o$Trendchange_direction)[2:3]) )
o %>% dplyr::filter(Trendchange_direction == "N") %>% dplyr::summarise(m = median(trendchange,na.rm=T),s = mad(trendchange,na.rm = T))
o %>% dplyr::filter(Trendchange_direction == "P") %>% dplyr::summarise(m = median(trendchange,na.rm=T),s = mad(trendchange,na.rm = T))

mean(out$LargeTimeAgo,na.rm=T);sd(out$LargeTimeAgo,na.rm=T)

# --------------------------------------------------#
# The code below recreates the main figures of the manuscript
# Figures might slightly differ visually as they are edited posthoc

#### Overall tests and results ####
options(na.action = "na.exclude")

o <- out

# Test for best random intercept
# Find SR random
f1 <- glmer(Species_richness ~ Break_binom + (1|SS),data=o,family = "poisson")
f2 <- glmer(Species_richness ~ Break_binom + (1|SS)+(1|LCLU),data=o,family = "poisson")
f3 <- glmer(Species_richness ~ Break_binom + (1|SS)+(1+Koeppen|LCLU),data=o,family = "poisson")
f4 <- glmer(Species_richness ~ Break_binom + (1|SS)+(Break_binom|LCLU),data=o,family = "poisson")
f5 <- glmer(Species_richness ~ Break_binom + (1|SS)+(1|Koeppen/LCLU),data=o,family = "poisson")
f6 <- glmer(Species_richness ~ Break_binom + (1|SS)+(1|LCLU) + (1|SSB),data=o,family = "poisson")
f7 <- glmer(Species_richness ~ Break_binom + (1|SS)+(1|Koeppen) + (1|LCLU),data=o,family = "poisson")
f8 <- glmer(Species_richness ~ Break_binom + (Break_binom|SS)+(1|LCLU)+(1|SSB),data=o,family = "poisson")

AICcmodavg::aictab(cand.set = list(f1,f2,f3,f4,f5,f6,f7,f8),modnames = c("(1|SS)","(1|SS)(1|LCLU)","(1|SS)(1Koeppen|LCLU)","(1|SS)(Break_binom|LCLU)",
                                                                         "(1|SS)(1|Biome/LCLU)","(1|SS)(1|LCLU)(1|SSB)","(1|SS)(1|Biome)(1|LCLU)","(Break_binom|SS)(1|LCLU)(1|SSB)"),sort = T)
# (Break_binom|SS)+(1|LCLU)+(1|SSB), then (1|SS)(1|LCLU)(1|SSB) <- Best random 

# For abundance
f1 <- glmer(logabund ~ 1 + (1|SS),data=out,family = "gaussian",REML=F)
f2 <- glmer(logabund ~ 1 + (1|SS)+(1|LCLU),data=out,family = "gaussian",REML=F)
f3 <- glmer(logabund ~ 1 + (1|SS)+(1+Koeppen|LCLU),data=out,family = "gaussian",REML=F)
f4 <- glmer(logabund ~ 1 + (1|SS)+(Break_binom|LCLU),data=out,family = "gaussian",REML=F)
f5 <- glmer(logabund ~ 1 + (1|SS)+(1|Koeppen/LCLU),data=out,family = "gaussian",REML=F)
f6 <- glmer(logabund ~ 1 + (1|SS)+(1|LCLU) + (1|SSB),data=out,family = "gaussian",REML=F)
f7 <- glmer(logabund ~ 1 + (1|SS)+(1|Koeppen) + (1|LCLU),data=out,family = "gaussian",REML=F)
f8 <- glmer(logabund ~ 1 + (Break_binom|SS)+(1|LCLU)+(1|SSB),data=out,family = "gaussian",REML=F)
f9 <- glmer(logabund ~ 1 + (1|SS)+(Break_binom|LCLU)+(1|SSB),data=out,family = "gaussian",REML=F)

AICcmodavg::aictab(cand.set = list(f1,f2,f3,f4,f5,f6,f7,f8,f9),modnames = c("(1|SS)","(1|SS)(1|LCLU)","(1|SS)(1Koeppen|LCLU)","(1|SS)(Break_binom|LCLU)","(1|SS)(1|Biome/LCLU)","(1|SS)(1|LCLU)(1|SSB)","(1|SS)(1|Biome)(1|LCLU)","(Break_binom|SS)(1|LCLU)(1|SSB)","(1|SS)(Break_binom|LCLU)(1|SSB)"),sort = T)

# ------------------------------------------ #
# Do the models for overall disturbance magnitude bins
# ------------------------------------------ #
library(glmmTMB);library(sjPlot)
# First for sites without any break
options(na.action = "na.exclude")

o <- out
before = out$largest_trendbef*12
after = out$largest_trendaft*12
o$trendchange <- (after-before)
o$Trendchange_direction = ifelse(o$trendchange < 0,"N","P");o$Trendchange_direction[which(is.na(o$trendchange))] <- "ND"
o$Trendchange_direction <- factor(o$Trendchange_direction,c("ND","N","P"))
o$Break_binom <- factor(o$Break_binom,labels = c("ND","D"))

# --- #
# Overall effect
sfit.overall <- glmer(Species_richness ~ Break_binom + (Break_binom|SS)+(1|LCLU)+(1|SSB) +(1|SSBS),data=o,family = "poisson");checkConv(sfit.overall)
sfit.null <- glmer(Species_richness ~ 1 + (Break_binom|SS)+(1|LCLU)+(1|SSB) +(1|SSBS),data=o,family = "poisson")
anova(sfit.overall,sfit.null);summary(sfit.overall)
# LA 
afit.overall <- glmer(logabund ~ Break_binom + (Break_binom|SS)+(1|LCLU) + (1|SSB),data=o %>% dplyr::filter(Diversity_metric_type == "Abundance"),family = "gaussian");checkConv(afit.overall)
afit.null <- glmer(logabund ~ 1 + (1|SS)+(1|LCLU)+(1|SSB),data=o %>% dplyr::filter(Diversity_metric_type == "Abundance"),family = "gaussian")
anova(afit.overall,afit.null)
summary(afit.overall)
# PIE
piefit.overall <- glmer(asPIE ~ Break_binom + (Break_binom|SS)+(1|LCLU) + (1|SSB),data=o %>% dplyr::filter(Diversity_metric_type == "Abundance"),family = "gaussian");checkConv(piefit.overall)
pfit.null <- glmer(asPIE ~ 1 + (1|SS)+(1|LCLU)+(1|SSB),data=o %>% dplyr::filter(Diversity_metric_type == "Abundance"),family = "gaussian")
summary(piefit.overall)
anova(piefit.overall,pfit.null)

## Difference between Models ? ##
sfit1 <- glmer(Species_richness ~ BinMagn + (Break_binom|SS)+(1|LCLU)+(1|SSB) +(1|SSBS),data=o,family = "poisson");checkConv(sfit1)
sfit2 <- glmer(Species_richness ~ BinTrend + (Break_binom|SS)+(1|LCLU)+(1|SSB) +(1|SSBS),data=o,family = "poisson");checkConv(sfit2)
sfit0 <- glmer(Species_richness ~ Break_binom + (Break_binom|SS)+(1|LCLU)+(1|SSB) +(1|SSBS),data=o,family = "poisson");checkConv(sfit0)
diff( AIC(sfit1,sfit2)$AIC )
anova(sfit0,sfit1)
anova(sfit0,sfit2)
cor.test(as.vector(fixef(sfit1)[-1]) , as.vector(fixef(sfit2)[-1]),method = "pear")

# LA
afit1 <- glmer(logabund ~ BinMagn + (Break_binom|SS)+(1|LCLU)+(1|SSB),data=o %>% dplyr::filter(Diversity_metric_type == "Abundance"),family = "gaussian")
afit2 <- glmer(logabund ~ BinTrend + (Break_binom|SS)+(1|LCLU)+(1|SSB),data=o %>% dplyr::filter(Diversity_metric_type == "Abundance"),family = "gaussian")
afit0 <- glmer(logabund ~ Break_binom + (1|SS)+(1|LCLU)+(1|SSB),data=o %>% dplyr::filter(Diversity_metric_type == "Abundance"),family = "gaussian")
afit3 <- glmer(logabund ~ BinMagn + BinTrend + (1|SS)+(1|LCLU)+(1|SSB),data=o %>% dplyr::filter(Diversity_metric_type == "Abundance"),family = "gaussian")
diff( AIC(afit1,afit2)$AIC )
anova(afit0,afit2)
anova(afit1,afit2)
r.squaredGLMM(afit0);sjstats::icc(afit0)
cor.test(as.vector(fixef(afit1)[-1]) , as.vector(fixef(afit2)[-1]),method = "pear")

# PIE
pfit1 <- glmer(asPIE ~ BinMagn + (Break_binom|SS)+(1|LCLU)+(1|SSB),data=o %>% dplyr::filter(Diversity_metric_type == "Abundance"),family = "gaussian")
pfit2 <- glmer(asPIE ~ BinTrend + (Break_binom|SS)+(1|LCLU)+(1|SSB),data=o %>% dplyr::filter(Diversity_metric_type == "Abundance"),family = "gaussian")
pfit0 <- glmer(asPIE ~ Break_binom + (1|SS)+(1|LCLU)+(1|SSB),data=o %>% dplyr::filter(Diversity_metric_type == "Abundance"),family = "gaussian")
diff( AIC(pfit1,pfit2)$AIC )
anova(pfit0,pfit1,pfit2)
anova(pfit0,pfit2)
cor.test(as.vector(fixef(pfit1)[-1]) , as.vector(fixef(pfit2)[-1]),method = "pear")

## For time passed ##
sfit1 <- glmer(Species_richness ~ interaction(Break_direction,BinTime) + (1|SS)+(1|LCLU)+(1|SSB) +(1|SSBS),data=o,family = "poisson");checkConv(sfit1)
sfit2 <- glmer(Species_richness ~ interaction(Trendchange_direction,BinTime) + (1|SS)+(1|LCLU)+(1|SSB) +(1|SSBS),data=o,family = "poisson");checkConv(sfit2)
sfit0 <- glmer(Species_richness ~ BinTime + (1|SS)+(1|LCLU)+(1|SSB) +(1|SSBS),data=o,family = "poisson");checkConv(sfit0)
diff( AIC(sfit1,sfit2)$AIC )
anova(sfit0,sfit1)
anova(sfit0,sfit2)
cor.test(as.vector(fixef(sfit1)[-1]) , as.vector(fixef(sfit2)[-1]),method = "pear")

# LA
afit1 <- glmer(logabund ~ interaction(Break_direction,BinTime) + (1|SS)+(1|LCLU)+(1|SSB),data=o %>% dplyr::filter(Diversity_metric_type == "Abundance"),family = "gaussian")
afit2 <- glmer(logabund ~ interaction(Trendchange_direction,BinTime) + (1|SS)+(1|LCLU)+(1|SSB),data=o %>% dplyr::filter(Diversity_metric_type == "Abundance"),family = "gaussian")
afit0 <- glmer(logabund ~ BinTime + (1|SS)+(1|LCLU)+(1|SSB),data=o %>% dplyr::filter(Diversity_metric_type == "Abundance"),family = "gaussian")
diff( AIC(afit1,afit2)$AIC )
anova(afit0,afit1)
anova(afit1,afit2)
cor.test(as.vector(fixef(afit1)[-1]) , as.vector(fixef(afit2)[-1]),method = "pear")

pfit1 <- glmer(asPIE ~ interaction(Break_direction,BinTime) + (Break_binom|SS)+(1|LCLU)+(1|SSB),data=o %>% dplyr::filter(Diversity_metric_type == "Abundance"),family = "gaussian")
pfit2 <- glmer(asPIE ~ interaction(Trendchange_direction,BinTime) + (Break_binom|SS)+(1|LCLU)+(1|SSB),data=o %>% dplyr::filter(Diversity_metric_type == "Abundance"),family = "gaussian")
pfit0 <- glmer(asPIE ~ BinTime + (Break_binom|SS)+(1|LCLU)+(1|SSB),data=o %>% dplyr::filter(Diversity_metric_type == "Abundance"),family = "gaussian")
diff( AIC(pfit1,pfit2)$AIC )
anova(pfit0,pfit1)
anova(pfit0,pfit2)
cor.test(as.vector(fixef(pfit1)[-1]) , as.vector(fixef(pfit2)[-1]),method = "pear")

# ------------------------------------------ #
#### Figure 2 ####
# ------------------------------------------ #
library(glmmTMB);library(sjPlot)
# First for sites without any break
options(na.action = "na.exclude")
o <- out

sfit.overall <- glmer(Species_richness ~ BinMagn + (Break_binom|SS)+(1|LCLU)+(1|SSB) +(1|SSBS),data=o,family = "poisson");checkConv(sfit.overall)
sfit.null <- glmer(Species_richness ~ 1 + (1|SS)+(1|LCLU)+(1|SSB) + (1|SSBS),data=o,family = "poisson");checkConv(sfit.overall)
anova(sfit.overall,sfit.null)
overdisp_fun(sfit.overall) # Overdispersion, thus site observation level random effect
#sjPlot::plot_model(sfit.overall,type = 'est',show.p = T,show.values = T)
sfit.overall.x <- format.results(sfit.overall,"BinMagn",maxlevels = 6) 

# --- #
afit.overall <- glmer(logabund ~ BinMagn + (Break_binom|SS)+(1|LCLU) + (1|SSB),data=o %>% dplyr::filter(Diversity_metric_type == "Abundance"),family = "gaussian");checkConv(afit.overall)
afit.null <- glmer(logabund ~ 1 + (1|SS)+(1|LCLU)+(1|SSB),data=o %>% dplyr::filter(Diversity_metric_type == "Abundance"),family = "gaussian")
anova(afit.overall,afit.null)
#sjPlot::plot_model(afit.overall,type = 'est',show.p = T,show.values = T)
afit.overall.x <- format.results(afit.overall,"BinMagn",maxlevels = 6) 

# --- #
piefit.overall <- glmer(asPIE ~ BinMagn + (Break_binom|SS)+(1|LCLU) + (1|SSB),data=o %>% dplyr::filter(Diversity_metric_type == "Abundance"),family = "gaussian");checkConv(piefit.overall)
pfit.null <- glmer(asPIE ~ 1 + (1|SS)+(1|LCLU)+(1|SSB),data=o %>% dplyr::filter(Diversity_metric_type == "Abundance"),family = "gaussian")
anova(piefit.overall,pfit.null)
#sjPlot::plot_model(piefit.overall,type = 'est',show.p = T,show.values = T)
piefit.overall.x <- format.results(piefit.overall,"BinMagn",maxlevels = 6) 

# Combine for plot 
sfit.overall.x$wrap.facet <- NULL
full1 <- rbind(sfit.overall.x %>% mutate(type="Species richness"),
              afit.overall.x %>% mutate(type = "Total abundance"),
              piefit.overall.x %>% mutate(type= "Probability of interspecific encounter") )
# ---- #
### ############################################################### ###
# Next for trend
sfit.overall <- glmer(Species_richness ~ BinTrend + (Break_binom|SS)+(1|LCLU)+(1|SSB) +(1|SSBS),data=o,family = "poisson");checkConv(sfit.overall)
sfit.null <- glmer(Species_richness ~ 1 + (1|SS)+(1|LCLU)+(1|SSB) +(1|SSBS),data=o,family = "poisson");checkConv(sfit.overall)
anova(sfit.overall,sfit.null)
overdisp_fun(sfit.overall) # Overdispersion, thus site observation level random effect
#sjPlot::plot_model(sfit.overall,type = 'est',show.p = T,show.values = T)
sfit.overall.x <- format.results(sfit.overall,"BinTrend",maxlevels = 6) 
# --- #
afit.overall <- glmer(logabund ~ BinTrend + (Break_binom|SS)+(1|LCLU) + (1|SSB),data=o %>% dplyr::filter(Diversity_metric_type == "Abundance"),family = "gaussian");checkConv(afit.overall)
afit.null <- glmer(logabund ~ 1 + (1|SS)+(1|LCLU)+(1|SSB),data=o %>% dplyr::filter(Diversity_metric_type == "Abundance"),family = "gaussian")
anova(afit.overall,afit.null)
#sjPlot::sjp.glmer(afit.overall,type = "fe",prnt.plot=T,fade.ns=T,show.intercept = F)
#sjPlot::plot_model(afit.overall,type = 'est',show.p = T,show.values = T)
afit.overall.x <- format.results(afit.overall,"BinTrend",maxlevels = 6) 
# --- #
piefit.overall <- glmer(asPIE ~ BinTrend + (Break_binom|SS)+(1|LCLU) + (1|SSB),data=o %>% dplyr::filter(Diversity_metric_type == "Abundance"),family = "gaussian");checkConv(piefit.overall)
pfit.null <- glmer(asPIE ~ 1 + (1|SS)+(1|LCLU)+(1|SSB),data=o %>% dplyr::filter(Diversity_metric_type == "Abundance"),family = "gaussian")
anova(piefit.overall,pfit.null)
#sjPlot::plot_model(piefit.overall,type = 'est',show.p = T,show.values = T)
piefit.overall.x <- format.results(piefit.overall,"BinTrend",maxlevels = 6) 

# Plot 
sfit.overall.x$wrap.facet <- NULL
full2 <- rbind(sfit.overall.x %>% mutate(type = "Species richness"),
               afit.overall.x %>% mutate(type = "Total abundance"),
               piefit.overall.x %>% mutate(type = "Probability of interspecific encounter") )

# ------------------ #
# Now combine both
full = rbind(
  full1 %>% mutate(model = "Magnitude") %>% 
      mutate(term = factor( str_split(term,"BinMagn",simplify = T)[,2] ) ) %>% 
      mutate(term = fct_collapse(term,
                           "ND" = "",
                          "---" = c("< -50%","Largest NTC"),
                          "--" = c("-50% <> -25%","Large NTC"),
                          "-" = c("-25% <> 0%","Small NTC"),
                          "+" = c("0% <> 25%","Small PTC"),
                          "++" = c("25% <> 50%","Large PTC"),
                          "+++" = c("> 50%","Large PTC")
                )
      ),
  full2 %>% mutate(model = "Trend") %>% 
    mutate(term = factor( str_split(term,"BinTrend",simplify = T)[,2] ) ) %>% 
    mutate(term = fct_collapse(term,
                               "ND" = "",
                               "---" = c("< -50%","Largest NTC"),
                               "--" = c("-50% <> -25%","Large NTC"),
                               "-" = c("-25% <> 0%","Small NTC"),
                               "+" = c("0% <> 25%","Small PTC"),
                               "++" = c("25% <> 50%","Large PTC"),
                               "+++" = c("> 50%","Largest PTC")
    )
    )
)

# Other labels
full$type <- factor(full$type,levels=c("Species richness","Total abundance","Probability of interspecific encounter"),
                    labels = c("Species richness","Total abundance","Probability of\n interspecific encounter") 
#                    labels = c("Species richness","Total abundance","Rarefied richness","Probability\n of interspecific encounter") 
                    )
full$term <- factor(full$term,levels=c("---","--","-","ND","+","++","+++"),labels = c("---","--","-","ND","+","++","+++") ) # Backtransform to factor
full$direction <- ifelse(str_detect(as.character(full$term),"-"),"N","P");full$direction[which(full$term=="ND")] <- "ND"
full$direction <- factor(full$direction,levels = c("ND","N","P"))
# Format p-string value, remove dots and whitespaces
full$p.stars[which(full$term==0)] <- ""
redBlue <- c("#d73027","#fc8d59","#fee090","black","#e0f3f8","#91bfdb","#4575b4")
# Correct for axis
full[,c("estimate","se.low","se.high")] <- (full[,c("estimate","se.low","se.high")] -1 ) *100

# Remove duplicate intercept
full <- full[-which(full$term=="ND" & full$model=="Trend"),]
# Correct intercept symbol
full$model[which(full$term=="ND")] <- "Intercept"
# Rename
full$term <- fct_recode(full$term, "UC" = "ND")

full.magnitude <- full # Security copy

g1 <- full %>%  dplyr::filter(type != c("Rarefied richness") ) %>% 
  ggplot(.,aes(x=term,y=estimate, ymin=se.low, ymax=se.high,group=model,shape=model,color=term) ) + 
  theme_tufte(base_size = 22,base_family="Arial",ticks = T) +
  theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank() ) +
  # Point range
  geom_hline(yintercept = 0,color="grey") +
  #geom_hline(yintercept = 0,linetype="dotted", color = "white") + # For transparent
  #geom_point(fatten = 3,size = 4, stroke = 6,position = position_dodge(.5)) + geom_linerange(fatten = 3,size = 4, position = position_dodge(.5)) +
  geom_pointrange(fatten = 3,size=2.5,position = position_dodge(.6)) + 
  scale_y_continuous(breaks=pretty_breaks(5), limits = c(-28, 10)) + 
  scale_color_manual(values = redBlue) + 
  scale_shape_manual(values = c(20,15,18)) +
  facet_wrap(~type,scales = "fixed",as.table = T,ncol = 1,strip.position = "left") + 
  theme(strip.text = element_text(size = 18), strip.background = element_blank(), strip.placement = "outside") +#, ,panel.spacing.y = unit(-.25, "lines"),panel.spacing.x = unit(-.5, "lines")) +
  # Add text
  #geom_text(aes(label=nSSBS, y=max(estimate) + 0.035*max(estimate)), colour="grey20", size=3,angle=90) +
  # Add significance symbols
  geom_text(aes(x = term, y=se.low-3.5, label = p.stars ),position = position_dodge(.6), size = 6,colour = "black") +
#  geom_text(aes(x = term, y=se.low-3.5, label = p.stars ),position = position_dodge(.5), size = 6,colour = "white") + # For transparent
  labs(x="Shift in magnitude or trend",y="Difference (% \u00B1 1 SE)") + 
#  annotate("text",x = .9,y = 10,label = c("a","b","c"),size=8,fontface = "bold", family = "Gill Sans MT") +   # Facet Label
  guides(color = "none", shape = "none") +  
  theme(axis.text.x = element_text(size = 22,face = 'bold'),axis.ticks.x = element_blank()) +
  theme(panel.spacing.x =  unit(1, "lines")) +
  # Add Y -axis again
  theme(axis.line.y = element_line(size = .5)) 
g1
ggsave("Figure1_a.png",plot = g1,width = 7,height = 9,dpi = 400)


# --------------------------------------------- #
# For Time bin
o <- out
sfit.overall <- glmer(Species_richness ~ BinTime + (Break_binom|SS)+(1|LCLU)+(1|SSB) +(1|SSBS),data=o,family = "poisson");checkConv(sfit.overall)
sfit.null <- glmer(Species_richness ~ 1 + (1|SS)+(1|LCLU)+(1|SSB) +(1|SSBS),data=o,family = "poisson");checkConv(sfit.overall)
anova(sfit.overall,sfit.null)
format.results(sfit.overall,"BinTime",maxlevels = 3) 

afit.overall <- glmer(logabund ~ BinTime + (Break_binom|SS)+(1|LCLU) + (1|SSB),data=o %>% dplyr::filter(Diversity_metric_type == "Abundance"),family = "gaussian");checkConv(afit.overall)
afit.null <- glmer(logabund ~ 1 + (1|SS)+(1|LCLU)+(1|SSB),data=o %>% dplyr::filter(Diversity_metric_type == "Abundance"),family = "gaussian")
format.results(afit.overall,"BinTime",maxlevels = 3)[,'p.value'][2]
anova(afit.overall,afit.null);format.results(afit.overall,"BinTime",maxlevels = 3) 

piefit.overall <- glmer(asPIE ~ BinTime + (Break_binom|SS)+(1|LCLU) + (1|SSB),data=o %>% dplyr::filter(Diversity_metric_type == "Abundance"),family = "gaussian");checkConv(piefit.overall)
pfit.null <- glmer(asPIE ~ 1 + (1|SS)+(1|LCLU)+(1|SSB),data=o %>% dplyr::filter(Diversity_metric_type == "Abundance"),family = "gaussian")
anova(piefit.overall,pfit.null);format.results(piefit.overall,"BinTime",maxlevels = 3) 
format.results(piefit.overall,"BinTime",maxlevels = 3)[,'p.value'][2]

# ---- #
o <- out
# Filter only to extremes for both 
before = out$largest_trendbef*12
after = out$largest_trendaft*12
o$trendchange <- (after-before)
o$Trendchange_direction = ifelse(o$trendchange < 0,"N","P");o$Trendchange_direction[which(is.na(o$trendchange))] <- "S"
o$Trendchange_direction <- factor(o$Trendchange_direction,c("S","N","P"))
o$Inter <- interaction(o$Break_direction,o$Trendchange_direction)
# Make new factor level for Break_direction and time
o <- o %>% mutate(NewComb = str_c(Break_direction,"_",BinTime) ) %>% mutate(NewComb = fct_relevel(NewComb,"S_ND"))

# Make two separate fits
sfit.overall <- glmer(Species_richness ~ NewComb + (Break_binom|SS)+(1|LCLU) + (1|SSB)+(1|SSBS),data=o,family = "poisson")
sfit.null <- glmer(Species_richness ~ 1 + (1|SS)+(1|LCLU) + (1|SSB)+(1|SSBS),data=o,family = "poisson")
anova(sfit.overall,sfit.null)
#plot_model(sfit.overall, type = 'est',show.values = T)
sfit.overall.x <- format.results(sfit.overall,"NewComb",maxlevels = 6) # BinMagn = 6 levels
sfit.overall.x$wrap.facet <- NULL

# --- #
afit.overall <- glmer(logabund ~ NewComb + (Break_binom|SS) +(1|LCLU) + (1|SSB),data=o %>% dplyr::filter(Diversity_metric_type == "Abundance"),family = "gaussian")
afit.null <- glmer(logabund ~ 1 + (1|SS) +(1|LCLU) + (1|SSB),data=o %>% dplyr::filter(Diversity_metric_type == "Abundance"),family = "gaussian")
anova(afit.overall,afit.null)
checkConv(afit.overall);#;sjPlot::sjp.glmer(afit.overall,type="re.qq")
#sjPlot::sjp.glmer(afit.overall,type = "fe",prnt.plot=T,fade.ns=T,show.intercept = F)
afit.overall.x <- format.results(afit.overall,"NewComb",maxlevels = 6) 

piefit.overall <- glmer(asPIE ~ NewComb + (Break_binom|SS)+(1|LCLU) + (1|SSB),data=o %>% dplyr::filter(Diversity_metric_type == "Abundance"),family = "gaussian");checkConv(piefit.overall)
piefit.null <- glmer(asPIE ~ 1 + (1|SS)+(1|LCLU) + (1|SSB),data=o %>% dplyr::filter(Diversity_metric_type == "Abundance"),family = "gaussian")
anova(piefit.overall,piefit.null)
#sjPlot::sjp.glmer(piefit.overall,type = "fe",prnt.plot=T,fade.ns=T,show.intercept = F)
piefit.overall.x <- format.results(piefit.overall,"NewComb",maxlevels = 6) 

# ---------------------------------------------- #
full1 <- rbind(sfit.overall.x %>% mutate(type="Species richness"),
              afit.overall.x %>% mutate(type = "Total abundance"),
              piefit.overall.x %>% mutate(type= "Probability of interspecific encounter") )
# ---------------------------------------------- #
## Now do for trend changes  ##
o <- o %>% mutate(NewComb = str_c(Trendchange_direction,"_",BinTime) ) %>% mutate(NewComb = fct_relevel(NewComb,"S_ND"))

# Make two separate fits
sfit.overall <- glmer(Species_richness ~ NewComb + (1|SS)+(1|LCLU) + (1|SSB)+(1|SSBS),data=o,family = "poisson")
sfit.null <- glmer(Species_richness ~ 1 + (1|SS)+(1|LCLU) + (1|SSB)+(1|SSBS),data=o,family = "poisson")
anova(sfit.overall,sfit.null)
overdisp_fun(sfit.overall);checkConv(sfit.overall) #;sjPlot::sjp.glmer(sfit.overall,type="re.qq")
#plot_model(sfit.overall, type = 'est',show.values = T)
sfit.overall.x <- format.results(sfit.overall,"NewComb",maxlevels = 6) # BinMagn = 6 levels
sfit.overall.x$wrap.facet <- NULL

# --- #
afit.overall <- glmer(logabund ~ NewComb + (Break_binom|SS) +(1|LCLU) + (1|SSB),data=o %>% dplyr::filter(Diversity_metric_type == "Abundance"),family = "gaussian")
afit.null <- glmer(logabund ~ 1 + (1|SS) +(1|LCLU) + (1|SSB),data=o %>% dplyr::filter(Diversity_metric_type == "Abundance"),family = "gaussian")
anova(afit.overall,afit.null)
checkConv(afit.overall);#;sjPlot::sjp.glmer(afit.overall,type="re.qq")
#sjPlot::sjp.glmer(afit.overall,type = "fe",prnt.plot=T,fade.ns=T,show.intercept = F)
afit.overall.x <- format.results(afit.overall,"NewComb",maxlevels = 6) 

piefit.overall <- glmer(asPIE ~ NewComb + (Break_binom|SS)+(1|LCLU) + (1|SSB),data=o %>% dplyr::filter(Diversity_metric_type == "Abundance"),family = "gaussian");checkConv(piefit.overall)
piefit.null <- glmer(asPIE ~ 1 + (1|SS)+(1|LCLU) + (1|SSB),data=o %>% dplyr::filter(Diversity_metric_type == "Abundance"),family = "gaussian")
anova(piefit.overall,piefit.null)
#sjPlot::sjp.glmer(piefit.overall,type = "fe",prnt.plot=T,fade.ns=T,show.intercept = F)
piefit.overall.x <- format.results(piefit.overall,"NewComb",maxlevels = 6) 

# ---------------------------------------------- #
full2 <- rbind(sfit.overall.x %>% mutate(type="Species richness"),
               afit.overall.x %>% mutate(type = "Total abundance"),
               piefit.overall.x %>% mutate(type= "Probability of interspecific encounter") )

## NOW COMBINE ##
full <- rbind(full1 %>% mutate(Model = "Magnitude"),full2 %>% mutate(Model = "Trend"))

full$term <- str_split(full$term,"NewComb",simplify = T)[,2]
full$term[which(full$term=="")] <- "S_S"
full <- separate(full,term,c("direction","Comparison"),"_")

full$type <- factor(full$type,levels=c("Species richness","Total abundance","Probability of interspecific encounter"),
                    labels = c("Species richness","Total abundance","Probability\n of interspecific encounter") 
                    )
full$Comparison <- factor(full$Comparison,levels = c("S",">0-5y","5-10y",">10y"),labels = c("ND","\U2264 5","5-10",">10") ) # Backtransform to factor
full$direction <- factor(full$direction,levels = c("S","N","P"),labels = c("ND","N","P"))

# Correct for axis
full[,c("estimate","se.low","se.high")] <- (full[,c("estimate","se.low","se.high")] -1 )*100
#full$p.string2 <- gsub('[[:digit:]]+', '', full$p.string); full$p.string2 <- gsub("\\.","",full$p.string2); full$p.string2 <- gsub("^[[:space:]]+|[[:space:]]+$", "", full$p.string2)
full$p.stars[which(full$direction==0)] <- ""

full$Type_term <- paste(full$Model,full$direction) # New Type term
# Remove trend intercept
full <- full[-which(full$Comparison=="ND" & full$Model=="Trend"),]

full$Model[which(full$Comparison=="ND")] <- "ND"
full$Model <- factor(full$Model,levels = c("ND","Magnitude","Trend"))
full$Type_term <- factor(full$Type_term,levels = c("Magnitude ND","Magnitude N","Trend N","Magnitude P","Trend P"))

full$NewX <- interaction(full$direction,full$Type_term)
# Recode
full$Model <- fct_recode(full$Model,UC = 'ND')

# Backup copy
full.time <- full

g2 <- full  %>% 
  dplyr::filter(type != "Rarefied richness", Comparison != 'ND'
                ) %>% 
  ggplot(.,aes(x=Comparison,y=estimate, ymin=se.low, ymax=se.high,group=Type_term,shape=Model,color=direction
               ) ) + 
  theme_tufte(base_size = 22,base_family="Arial",ticks = T) +
  theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank() ) +
  # Point range
#  geom_pointrange(fatten = 3,size=1.5,position = position_dodge(.75)) +
#  geom_point(x=0.5,y=0,size = 6,color="black",shape=20) +
  geom_hline(yintercept = 0,color = "grey") +
  geom_pointrange(fatten = 3,size=2.5,position = position_dodge(.6)) +
  #geom_point(x=0.5,y=0,size = 6,color="black",shape=20) + # FOr transparent view
  coord_cartesian() +
  scale_y_continuous(breaks=pretty_breaks(5), limits = c(-15,6) ) +
  #scale_shape_manual(values = c("\u25CF", "\u25bc","\u25b2", 
  #                              "\u25bc","\u25b2") ) +
  # Shapes
  scale_shape_manual(values = c(15,18),guide = 
                       guide_legend(title = '',reverse = T,ncol = 1,label.position = 'top',override.aes = aes(color = 'grey60'),
                                    label.vjust = .5,label.theme = element_text(size = 24,hjust = .5,angle = 90))) +
  theme(legend.position = 'right',legend.direction = 'vertical') +
  scale_color_manual(values = c(redBlue[1],redBlue[7] ) ) + 
  facet_wrap(~type,scales = "fixed",ncol=1,as.table = T,strip.position = "left") + 
  theme(strip.text = element_text(size = 20), strip.background = element_blank(), strip.placement = "outside") +#, ,panel.spacing.y = unit(-.25, "lines"),panel.spacing.x = unit(-.5, "lines")) +
  # Add text
  #geom_text(aes(label=nSSBS, y=max(estimate) + 0.055*max(estimate)), colour="grey20", size=3,angle=90) +
  # Add significance symbols
  geom_text(aes(x = Comparison, y=se.low-2.5, label = p.stars ),position = position_dodge(.6), size = 6,colour = "black") +
  #  geom_text(aes(x = Comparison, y=se.low-2.5, label = p.stars ),position = position_dodge(.75), size = 6,colour = "white") + # For transparent
  # Add plus and minuses + vlines
  #annotate("text", x = 0.8, y = -14,label = c("","","-"),size=13, fontface = "bold", color = redBlue[1] ) +
  #annotate("text", x = 1.8, y = -14,label = c("","","-"),size=13, fontface = "bold", color = redBlue[1] ) +
  #annotate("text", x = 2.8, y = -14,label = c("","","-"),size=13, fontface = "bold", color = redBlue[1] ) +
  #annotate("text", x = 1.2, y = -14,label = c("","","+"),size=13, fontface = "bold", color = redBlue[7] ) +
  #annotate("text", x = 2.2, y = -14,label = c("","","+"),size=13, fontface = "bold", color = redBlue[7] ) +
  #annotate("text", x = 3.2, y = -14,label = c("","","+"),size=13, fontface = "bold", color = redBlue[7] ) +
  geom_vline(xintercept = c(1,2,3),linetype = "dashed", color = "black") +
  labs(x="Time passed (years)",y="Difference (% \u00B1 1 SE)") + 
  guides(color = "none", shape = 'none') + 
  theme(axis.text.x = element_text(size=20,face = 'bold'),axis.ticks.x = element_blank()) +
  theme(panel.spacing.x =  unit(1, "lines")) + 
  # Add Y -axis again
  theme(axis.line.y = element_line(size = .5)) 
g2
ggsave("Figure2_b.png",plot = g2,width = 7,height = 9,dpi = 400)
# Without strip and y-axis label
ggsave("newFigures/F1Overall_2.png",plot = g2 + labs(y="") + theme(strip.text = element_blank()) + theme(axis.line.y = element_line(size = .5)),
       width = 7,height = 9,dpi = 400)

# Combine both
library(cowplot)
pg = cowplot::plot_grid(g1 + scale_y_continuous(limits = c(-28,7)) ,# + theme(axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0))),
                        g2  + ylab("") + theme(strip.text = element_blank() ) + theme(axis.line.y = element_line(size = .5)) ,
                        align = "h",axis = "l",nrow = 1,labels = "",
                        label_fontfamily = "Arial")
pg
cowplot::ggsave("newFigures/F1Overall_1and2.png",plot=pg,width = 11,height = 9)

# Data for export
Figure2_part1 <- full.magnitude
Figure2_part2 <- full.time

# ------------------------ #
# Diagnostic. Number per Bin and Metric
g.bar <- o %>% 
  dplyr::select(BinMagn,Species_richness,logabund,asPIE) %>% 
  reshape2::melt() %>% subset(.,complete.cases(.)) %>% 
  dplyr::select(-value) %>% # No need for the value as sites are counted
  #dplyr::filter(BinMagn != "ND") %>% 
  mutate(
    variable = factor(variable,levels=c("Species_richness","logabund","Richness_rarefied","asPIE"),labels = c("Species richness","Total abundance","Rarefied richness","Probability\n of interspecific encounter") ),
    BinMagn = factor(BinMagn,levels = c("< -50%","-50% <> -25%","-25% <> 0%","ND","0% <> 25%", "25% <> 50%","> 50%"),labels = c("---","--","-","ND","+","++","+++") ) 
  ) %>% 
  filter(variable != "Rarefied richness") %>% 
  ggplot(.,aes(x=BinMagn,group = BinMagn,fill=BinMagn)) + 
  theme_tufte(base_size = 20,base_family="Arial",ticks = T) +
  geom_bar(position = position_dodge()) +
  scale_fill_manual(values = redBlue) +
  scale_x_discrete(expand = c(0,0)) + theme(axis.text.x = element_text(size = 20)) +
  scale_y_log10(breaks = c(1,10,100,1000),expand = c(0,0)) +
  labs(x = "", y = "Number of sites\n (log10-transformed)") +
  facet_wrap(~variable,as.table = T,ncol=1)+
  guides(fill = "none" ) +
  theme(legend.margin = margin(.5),legend.position = "bottom",legend.spacing.x = unit(5,"in") ) +
  # Add Y -axis again
  theme(axis.line.y = element_line(size = .5))
g.bar
ggsave("newFigures/F1Overall_1_MetricBars.png",plot = g.bar,scale = 1.25,dpi = 400)

# --- # 
# Make a violin plot for time since disturbance for each point
g.time1 <- o %>% 
  dplyr::filter(Break_binom == 1) %>% 
  dplyr::select(BinMagn,LargeTimeAgo) %>% 
  mutate(
    BinMagn = factor(BinMagn,levels = c("< -50%","-50% <> -25%","-25% <> 0%","ND","0% <> 25%", "25% <> 50%","> 50%"),labels = c("---","--","-","ND","+","++","+++") ) 
  ) %>% 
  ggplot(.,aes(x = BinMagn, y= LargeTimeAgo,fill = BinMagn)) +
  theme_tufte(base_size = 20,base_family="Arial",ticks = T) + theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
  geom_violin(position = position_dodge(),scale = "width",color="white",alpha=.6,size=1) +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), geom = "pointrange",size=2, color = "black") +
  scale_fill_manual(values = redBlue[-4],guide = "none") + 
  scale_x_discrete(expand = c(0,0)) + theme(axis.text.x = element_text(size = 20)) +
  # Add text
  geom_text(data= o %>% group_by(BinMagn) %>% summarise(nSSBS = n()) %>% filter(BinMagn != "ND") %>% 
              mutate(BinMagn = factor(BinMagn,levels = c("< -50%","-50% <> -25%","-25% <> 0%","ND","0% <> 25%", "25% <> 50%","> 50%"),labels = c("---","--","-","ND","+","++","+++") ) ) ,
            aes(label=nSSBS, y=26.5), colour="black", size=5) +
  scale_y_continuous(breaks = pretty_breaks(5)) +
  labs(x = "Shift in magnitude", y = "Time passed (years)") +
  theme(legend.margin = margin(.5),legend.position = "bottom",legend.spacing.x = unit(5,"in") ) +
  # Add Y -axis again
  theme(axis.line.y = element_line(size = .5))
g.time1
ggsave("newFigures/F1Overall_1_TimeViolins.png",plot = g.time1,scale = 1.25,dpi = 400)

# Now for trend bins
g.time2 <- o %>% 
  dplyr::filter(Break_binom == 1) %>% 
  dplyr::select(BinTrend,LargeTimeAgo) %>% 
  mutate(
    BinMagn = factor(BinTrend,levels = c("Largest NTC","Large NTC","Small NTC","ND","Small PTC", "Large PTC","Largest PTC"),labels = c("---","--","-","0","+","++","+++") ) 
  ) %>% 
  ggplot(.,aes(x = BinMagn, y= LargeTimeAgo,fill = BinMagn)) +
  theme_tufte(base_size = 20,base_family="Arial",ticks = T) + theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
  geom_violin(position = position_dodge(),scale = "width",color="white",alpha=.6,size=1) +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), geom = "pointrange",size=2, color = "black") +
  scale_fill_manual(values = redBlue[-4],guide = "none") + 
  scale_x_discrete(expand = c(0,0)) + theme(axis.text.x = element_text(size = 20)) +
  # Add text
  geom_text(data= o %>% group_by(BinTrend) %>% summarise(nSSBS = n()) %>% filter(BinTrend != "ND") %>% 
              mutate(BinMagn = factor(BinTrend,levels =c("Largest NTC","Large NTC","Small NTC","ND","Small PTC", "Large PTC","Largest PTC"),labels = c("---","--","-","ND","+","++","+++") ) ) ,
            aes(label=nSSBS, y=26.5), colour="black", size=5) +
  scale_y_continuous(breaks = pretty_breaks(5)) +
  labs(x = "Shift in trend", y = "Time passed (years)") +
  theme(legend.margin = margin(.5),legend.position = "bottom",legend.spacing.x = unit(5,"in") ) +
  # Add Y -axis again
  theme(axis.line.y = element_line(size = .5))
g.time2
ggsave("newFigures/F1Overall_2_TimeViolins.png",plot = g.time2,scale = 1.25,dpi = 400)


# And per Broad-Traxonomic group
g.bar <- o %>% 
  #dplyr::filter(BinMagn != "ND") %>% 
  mutate(
    BinMagn = factor(BinMagn,levels = c("< -50%","-50% <> -25%","-25% <> 0%","ND","0% <> 25%", "25% <> 50%","> 50%"),labels = c("---","--","-","0","+","++","+++") ), 
    BinTrend = factor(BinTrend,levels = c("Largest NTC","Large NTC","Small NTC","ND","Small PTC", "Large PTC","Largest PTC"),labels = c("---","--","-","0","+","++","+++") ) 
  ) %>% 
  dplyr::select(BinMagn,BinTrend,TGrouping) %>% 
  reshape2::melt(id.vars="TGrouping") %>% 
  ggplot(.,aes(x=value,group = variable,fill=variable)) + 
  theme_tufte(base_size = 20,base_family="Gill Sans MT",ticks = T) +
  geom_bar(position = position_dodge()) +
  scale_fill_d3() +
  scale_x_discrete(expand = c(0,0)) + theme(axis.text.x = element_text(size = 20)) +
  scale_y_log10(breaks = c(1,10,100,1000),expand = c(0,0)) +
  facet_wrap(~TGrouping,scales = "free_x",nrow = 3) +
  labs(x = "", y = "Number of sites\n (log10-transformed)") +
  guides(fill = guide_legend(title = "",label.theme = element_text(size=14,angle = 0) ) ) +
  theme(legend.margin = margin(.5),legend.position = "bottom",legend.spacing.x = unit(5,"in") )
g.bar
ggsave("newFigures/F1Overall_1_TGroupingBars.png",plot = g.bar,scale = 1.25,dpi = 400)


# Make violins per time for trend changes
g.time <- o %>% 
  dplyr::filter(Break_binom == 1) %>% 
  dplyr::select(BinTrend,LargeTimeAgo) %>% 
  mutate(
    BinTrend  = factor(BinTrend, levels = c("Largest NTC","Large NTC","Small NTC","ND","Small PTC","Large PTC","Largest PTC"), labels = c("<<<","<<","<","ND",">",">>",">>>"))
  ) %>% 
  ggplot(.,aes(x = BinTrend, y= LargeTimeAgo,fill = BinTrend)) +
  theme_tufte(base_size = 20,base_family="Gill Sans MT",ticks = T) + theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
  geom_violin(draw_quantiles = c(.5),position = position_dodge(),scale = "width",alpha=.6,size=1) +
  scale_fill_manual(values = brownGreen[-4],guide = "none") + 
  scale_x_discrete(expand = c(0,0)) + theme(axis.text.x = element_text(size = 20)) +
  # Add text
  geom_text(data= o %>% group_by(BinTrend) %>% summarise(nSSBS = n()) %>% filter(BinTrend != "ND") %>% 
              mutate(
                BinTrend  = factor(BinTrend, levels = c("Largest NTC","Large NTC","Small NTC","ND","Small PTC","Large PTC","Largest PTC"), labels = c("<<<","<<","<","ND",">",">>",">>>")) ),
            aes(label=nSSBS, y=26.5), colour="grey20", size=5) +
  scale_y_continuous(breaks = pretty_breaks(5)) +
  labs(x = "", y = "Time since disturbance (years)") +
  theme(legend.margin = margin(.5),legend.position = "bottom",legend.spacing.x = unit(5,"in") )
g.time
ggsave("newFigures/F1Overall_3_TimeViolins.png",plot = g.time,scale = 1.25,dpi = 400)


# --------------------------------------------- #
#### Figure 3 ####
library(gridExtra);library(cowplot)
library(ggdendro)
library(cluster)
library(dendextend)

o <- out
# Load pairwise dissimilarity matrices
ss <- readRDS("resSaves/Out_MatricesSor.rds")
ss <- ss[which(names(ss) %in% unique(o$SS)) ]
# Combine all
ol <- data.frame()
for(study in names(ss)){
  print(study)
  sub <- ss[[study]] %>% reshape2::melt() %>% 
    # Filter to only abrupt changes present in comparison
    dplyr::filter(Var1 %in% o$SSBS,Var2  %in% o$SSBS) %>% 
    left_join(., 
              out %>% dplyr::select(SSBS,LCLU,BinMagn) %>% dplyr::rename(Var1 = SSBS),
              by = "Var1") %>% rename(BinMagn_1 = BinMagn) %>% dplyr::select(-Var1) %>% 
    left_join(., 
              out %>% dplyr::select(SSBS,LCLU,BinMagn) %>% dplyr::rename(Var2 = SSBS),
              by = "Var2") %>% rename(BinMagn_2 = BinMagn) %>% dplyr::select(-Var2) %>% 
    # Make a combinated column
    mutate(BinMagnComb = paste0(BinMagn_1,"_",BinMagn_2)) %>% dplyr::select(-BinMagn_1,-BinMagn_2) %>% 
    subset(.,complete.cases(.)) %>% mutate(SS = study)
  # Save
  ol <- rbind(ol, sub)
}
# Now summarize
out_d <- ol %>% mutate(LCLUCc = paste0(LCLU.x,"_",LCLU.y)) %>% 
  dplyr::group_by(BinMagnComb,LCLUCc) %>% 
  dplyr::summarize(avg = mean(value,na.rm=T)) %>% 
  # Only consider within same land use
  dplyr::filter(LCLUCc %in% c("PV_PV","SV_SV","HDV_HDV")) %>% 
  # Group and summarise
  dplyr::group_by(BinMagnComb) %>% 
  dplyr::summarize(avg = mean(avg,na.rm=T) ) %>% 
  # Join sample sizes back in 
  left_join(.,ol %>% mutate(LCLUCc = paste0(LCLU.x,"_",LCLU.y)) %>%
              dplyr::filter(LCLUCc %in% c("PV_PV","SV_SV","HDV_HDV")) %>%
              dplyr::group_by(BinMagnComb) %>%
              dplyr::summarize(N = n_distinct(unique(SS)) )
  ) %>% 
  # Split
  separate(BinMagnComb,c("A","B"),"_") %>% 
  mutate(A = factor(A,levels = c("< -50%","-50% <> -25%","-25% <> 0%","ND", "0% <> 25%","25% <> 50%","> 50%"),labels = c("---","--","-","ND","+","++","+++")  ),
         B = factor(B,levels = c("< -50%","-50% <> -25%","-25% <> 0%","ND", "0% <> 25%","25% <> 50%","> 50%"),labels = c("---","--","-","ND","+","++","+++") ))
stopifnot(max(out_d$N) < n_distinct(ol$SS))  # Security check

# Rescale all values relative tothe comparison between stable sites
for(s1 in unique(out_d$A) ){
  for(s2 in unique(out_d$B) )
    if(s1 != "ND" | s2 != "ND" ){
      out_d$avg[which(out_d$A== s1 & out_d$B == s2)] <- (out_d$avg[which(out_d$A == s1 & out_d$B == s2)] - out_d$avg[which(out_d$A == "ND" & out_d$B == "ND")] )
    }
}
out_d$avg[which(out_d$A=="ND" & out_d$B == "ND")] <- 0

out_d$A <- fct_recode(out_d$A, "UC" = "ND")
out_d$B <- fct_recode(out_d$B, "UC" = "ND")

# Convert to matrix
m <- reshape2::acast(out_d,A~B,value.var = "avg")
# Save for Output
Figure3_part1 <- m

hc <- hclust(dist(m,"man"),method = "complete");plot(hc)
# Reshuffle 
hc <- as.dendrogram(hc)
hc <- rotate(hc, order = c("---","+++","--","++","+","-","UC") )
gtree <- ggdendrogram(hc, rotate = F, size = 4, theme_dendro = FALSE) + 
  theme_tufte(base_size = 16,base_family="Arial",ticks = F) + 
  # Remove labels and ticks - order identical to factor order releveling
  labs(x= "",y="Manhattan distance") +
  scale_y_continuous(breaks = pretty_breaks(5)) +
  # Add Y -axis again
  theme(axis.line.y = element_line(size = .5)) 
gtree  
ggsave("F3_Treeclust.png",plot=gtree,width=5,height=3)

# Construct ggplot manually with sample size inserted
m2 <- m
m2[lower.tri(m2)] <- NA
# Do the same for label
m.lab <- reshape2::acast(out_d,A~B,value.var = "N",fun.aggregate = sum)
m.lab[lower.tri(m.lab)] <- NA

m.lab.col <- ifelse( abs(reshape2::melt(m2)[,"value"]) >= 0.05,"white","black")

gd <- reshape2::melt(m2) %>% subset(.,complete.cases(.)) %>% 
  # Reorder
  mutate(Var2 = factor(Var2,levels = c("+++","++","+","UC","-","--","---") )) %>% 
  # Reorder based on clustering
  ggplot(.,aes(x=Var1,y=Var2,fill=value)) + 
  theme_tufte(base_size = 20,base_family="Arial",ticks = F) + 
  geom_tile() + coord_equal() +# Tiles
  suppressWarnings( geom_text(data=reshape2::melt(m.lab),aes(x=Var1,y=Var2,label=value),
                              inherit.aes = F,size=5,color = m.lab.col,fontface = "bold") ) +
  scale_fill_gradient2(low = "#601200", mid = "white",high = "#001260",
                       na.value = "white",midpoint = 0,breaks = pretty_breaks(5),
                       guide = guide_colorbar(title = "Sørensen similarity index",title.position="top",title.hjust = .5,title.theme = element_text(size=16),
                                              direction = "horizontal",nbin=100,ticks = FALSE, barwidth = 18, barheight = 1)) +
  guides(color="none") + theme(legend.position = "bottom") +
  scale_x_discrete(position = "top") +
  labs(x = "", y = "",title = "") +
  theme(legend.margin=margin(t = -.75, unit='cm')) #  theme(legend.margin=margin(t=0, r=0, b=0, l=0, unit="cm"))
gd
ggsave("F3_grid.png",plot=gd + theme(plot.margin=unit(c(0,0,0,0), "mm")),width=6,height=6)

# --------------------------------------------- #
#                Time Bins                      # 
# --------------------------------------------- #
# Version 2 with time bins and direction
o <- out
ss <- readRDS("resSaves/Out_MatricesSor.rds")
ss <- ss[which(names(ss) %in% unique(o$SS)) ]
# Combine all
ol <- data.frame()
for(study in names(ss)){
  print(study)
  sub <- ss[[study]] %>% reshape2::melt() %>% 
    # Filter to only abrupt changes present in comparison
    dplyr::filter(Var1 %in% o$SSBS,Var2  %in% o$SSBS) %>% 
    left_join(., 
              out %>% dplyr::select(SSBS,LCLU,BinTime,Break_direction) %>% 
                mutate(NewComb = str_c(Break_direction,"_",BinTime) ) %>% mutate(NewComb = fct_relevel(NewComb,"S_ND")) %>% dplyr::select(-BinTime,-Break_direction) %>% 
                # Rename
                dplyr::rename(Var1 = SSBS),
              by = "Var1") %>% rename(NewComb_1 = NewComb) %>% dplyr::select(-Var1) %>% 
    left_join(., 
              out %>% dplyr::select(SSBS,LCLU,BinTime,Break_direction) %>% 
                mutate(NewComb = str_c(Break_direction,"_",BinTime) ) %>% mutate(NewComb = fct_relevel(NewComb,"S_ND")) %>% dplyr::select(-BinTime,-Break_direction) %>% 
                # Rename
                dplyr::rename(Var2 = SSBS),
              by = "Var2") %>% rename(NewComb_2 = NewComb) %>% dplyr::select(-Var2) %>% 
    # Make a combinated column
    mutate(NewCombFull = paste0(NewComb_1,"__",NewComb_2)) %>% dplyr::select(-NewComb_1,-NewComb_2) %>% 
    subset(.,complete.cases(.)) %>% mutate(SS = study)
  # Save
  ol <- rbind(ol, sub)
}

# Now summarize
out_d <- ol %>% mutate(LCLUCc = paste0(LCLU.x,"_",LCLU.y)) %>% 
  dplyr::group_by(NewCombFull,LCLUCc) %>% 
  dplyr::summarize(avg = mean(value,na.rm=T)) %>% 
  # Only consider within same land use
  dplyr::filter(LCLUCc %in% c("PV_PV","SV_SV","HDV_HDV")) %>% 
  # Group and summarise
  dplyr::group_by(NewCombFull) %>% 
  dplyr::summarize(avg = mean(avg,na.rm=T)) %>% # Maximal value to get the number of studies contributing to a combination
  # Join sample sizes back in 
  left_join(.,ol %>% mutate(LCLUCc = paste0(LCLU.x,"_",LCLU.y)) %>%
              dplyr::filter(LCLUCc %in% c("PV_PV","SV_SV","HDV_HDV")) %>%
              dplyr::group_by(NewCombFull) %>%
              dplyr::summarize(N = n_distinct(unique(SS)) )
  ) %>% 
  # Split
  separate(NewCombFull,c("A","B"),"__") %>% 
  mutate(A = factor(A,levels = c("N_>10y","N_5-10y","N_>0-5y","S_ND", "P_>0-5y","P_5-10y","P_>10y"),labels = c(">10 ","5-10 ","0-5 ","ND","0-5","5-10",">10")  ),
         B = factor(B,levels = c("N_>10y","N_5-10y","N_>0-5y","S_ND", "P_>0-5y","P_5-10y","P_>10y"),labels = c(">10 ","5-10 ","0-5 ","ND","0-5","5-10",">10")  ) )
stopifnot(max(out_d$N) < n_distinct(ol$SS))  # Security check

# Rescale all values relative tothe comparison between stable sites
for(s1 in unique(out_d$A) ){
  for(s2 in unique(out_d$B) )
    if(s1 != "ND" | s2 != "ND" ){
      out_d$avg[which(out_d$A== s1 & out_d$B == s2)] <- (out_d$avg[which(out_d$A == s1 & out_d$B == s2)] - out_d$avg[which(out_d$A == "ND" & out_d$B == "ND")] )
    }
}
out_d$avg[which(out_d$A=="ND" & out_d$B == "ND")] <- 0

out_d$A <- fct_recode(out_d$A, "UC" = "ND")
out_d$B <- fct_recode(out_d$B, "UC" = "ND")

# Convert to matrix
m <- reshape2::acast(out_d,A~B,value.var = "avg",fun.aggregate = mean)
# Save for output
Figure3_part2 <- m

m2 <- m
colnames(m2) <- c(">10 ","5-10 ","\U2264 5 ","UC", "\U2264 5", "5-10", ">10")
rownames(m2) <- c(">10 ","5-10 ","\U2264 5 ","UC", "\U2264 5", "5-10", ">10")

hc <- hclust(dist(m2,"man"),method = "complete")
cols <- c(redBlue[1],redBlue[7],redBlue[1],redBlue[7],redBlue[7],redBlue[1],"black")

gtree <- ggdendrogram(hc, rotate = F, size = 4, theme_dendro = FALSE) + 
  theme_tufte(base_size = 16,base_family="Arial",ticks = F) + 
  labs(x= "",y="Manhattan distance") +
  scale_y_continuous(breaks = pretty_breaks(5)) +
  # Custom x-axis labels
  theme(axis.text.x = element_text(colour = cols) ) +
  # Add Y -axis again
  theme(axis.line.y = element_line(size = .5))
gtree  
ggsave("F3_timetree.png",plot=gtree,width=5,height=3)

# Construct ggplot manually with sample size inserted
m2 <- m
m2[lower.tri(m2)] <- NA
# Do the same for label
m.lab <- reshape2::acast(out_d,A~B,value.var = "N")
#m.lab <- m.lab[hc$labels[hc$order],hc$labels[hc$order]] # Reorder
m.lab[lower.tri(m.lab)] <- NA

colnames(m2) <- c(">10 ","5-10 ","\U2264 5 ","UC", "\U2264 5", "5-10", ">10")
rownames(m2) <- c(">10 ","5-10 ","\U2264 5 ","UC", "\U2264 5", "5-10", ">10")

colnames(m.lab) <- c(">10 ","5-10 ","\U2264 5 ","UC", "\U2264 5", "5-10", ">10")
rownames(m.lab) <- c(">10 ","5-10 ","\U2264 5 ","UC", "\U2264 5", "5-10", ">10")
m.lab.col <- ifelse( abs(reshape2::melt(m2)[,"value"]) >= 0.05,"white","black")

cols <- c(redBlue[1],redBlue[1],redBlue[1],"black",redBlue[7],redBlue[7],redBlue[7])

gd <- reshape2::melt(m2) %>% subset(.,complete.cases(.)) %>% 
  # Reorder 
  ggplot(.,aes(x=Var1,y=fct_rev(Var2),fill=value)) + 
  theme_tufte(base_size = 20,base_family="Arial",ticks = F) + 
  geom_tile() + coord_equal() +# Tiles
  suppressWarnings( geom_text(data=reshape2::melt(m.lab),aes(x=Var1,y=Var2,label=value),inherit.aes = F,size=5,color = m.lab.col,fontface = "bold") ) +
  scale_fill_gradient2(low = "#601200", mid = "white",high = "#001260",
                       na.value = "white",midpoint = 0,breaks = pretty_breaks(5),
                       guide = guide_colorbar(title = "Sørensen similarity index",title.position="top",title.hjust = .5,title.theme = element_text(size=16),
                                              direction = "horizontal",nbin=100,ticks = FALSE, barwidth = 18, barheight = 1)) +
  guides(color="none") + theme(legend.position = "bottom") +
  scale_x_discrete(position = "top") +
  labs(x = "", y = "",title = "") +
  # Custom x-axis labels
  theme(axis.text.x = element_text(colour = cols) ) + theme(axis.text.y = element_text(colour = rev(cols)) ) +
  theme(legend.margin=margin(t = -.75, unit='cm'))#  theme(legend.margin=margin(t=0, r=0, b=0, l=0, unit="cm"))
gd
ggsave("F3_SorGrid_TimeDir.png",plot=gd+ theme(plot.margin=unit(c(0,0,0,0), "mm")),width=6,height=6,dpi = 300)

# ------------------------------ #
# Alternative with post-disturbance trend (SI)
o <- out
ss <- readRDS("resSaves/Out_MatricesSor.rds")
ss <- ss[which(names(ss) %in% unique(o$SS)) ]
# Combine all
ol <- data.frame()
for(study in names(ss)){
  print(study)
  sub <- ss[[study]] %>% reshape2::melt() %>% 
    # Filter to only abrupt changes present in comparison
    dplyr::filter(Var1 %in% o$SSBS,Var2  %in% o$SSBS) %>% 
    left_join(., 
              out %>% dplyr::select(SSBS,LCLU,BinTrend) %>% dplyr::rename(Var1 = SSBS),
              by = "Var1") %>% rename(BinTrend_1 = BinTrend) %>% dplyr::select(-Var1) %>% 
    left_join(., 
              out %>% dplyr::select(SSBS,LCLU,BinTrend) %>% dplyr::rename(Var2 = SSBS),
              by = "Var2") %>% rename(BinTrend_2 = BinTrend) %>% dplyr::select(-Var2) %>% 
    # Make a combinated column
    mutate(BinTrendComb = paste0(BinTrend_1,"_",BinTrend_2)) %>% dplyr::select(-BinTrend_1,-BinTrend_2) %>% 
    subset(.,complete.cases(.)) %>% mutate(SS = study)
  # Save
  ol <- rbind(ol, sub)
}

# Now summarize
out_d <- ol %>% mutate(LCLUCc = paste0(LCLU.x,"_",LCLU.y)) %>% 
  dplyr::group_by(BinTrendComb,LCLUCc) %>% 
  dplyr::summarize(avg = mean(value,na.rm=T)) %>% 
  # Only consider within same land use
  dplyr::filter(LCLUCc %in% c("PV_PV","SV_SV","HDV_HDV")) %>% 
  # Group and summarise
  dplyr::group_by(BinTrendComb) %>% 
  dplyr::summarize(avg = mean(avg,na.rm=T) ) %>% 
  # Join sample sizes back in 
  left_join(.,ol %>% mutate(LCLUCc = paste0(LCLU.x,"_",LCLU.y)) %>%
              dplyr::filter(LCLUCc %in% c("PV_PV","SV_SV","HDV_HDV")) %>%
              dplyr::group_by(BinTrendComb) %>%
              dplyr::summarize(N = n_distinct(unique(SS)) )
  ) %>% 
  # Split
  separate(BinTrendComb,c("A","B"),"_") %>% 
  mutate(A = factor(A,c("Largest NTC","Large NTC","Small NTC","ND","Small PTC","Large PTC","Largest PTC"),c("---","--","-","0","+","++","+++") ),
         B = factor(B,c("Largest NTC","Large NTC","Small NTC","ND","Small PTC","Large PTC","Largest PTC"),c("---","--","-","0","+","++","+++") )
  )
stopifnot(max(out_d$N) < n_distinct(ol$SS))  # Security check

# Rescale all values relative tothe comparison between stable sites
for(s1 in unique(out_d$A) ){
  for(s2 in unique(out_d$B) )
    if(s1 != "0" | s2 != "0" ){
      out_d$avg[which(out_d$A== s1 & out_d$B == s2)] <- (out_d$avg[which(out_d$A == s1 & out_d$B == s2)] - out_d$avg[which(out_d$A == "0" & out_d$B == "0")] )
    }
}
out_d$avg[which(out_d$A=="0" & out_d$B == "0")] <- 0

out_d$A <- fct_recode(out_d$A, "UC" = "0")
out_d$B <- fct_recode(out_d$B, "UC" = "0")

# Convert to matrix
m <- reshape2::acast(out_d,A~B,value.var = "avg",fun.aggregate = mean)

hc <- hclust(dist(m,"man"),method = "complete")
plot(hc)

hc <- rotate(hc, order = c("---","--","+++","-","++","+","UC") )

gtree <- ggdendrogram(hc, rotate = F, size = 4, theme_dendro = FALSE) + 
  theme_tufte(base_size = 16,base_family="Arial",ticks = F) + 
  labs(x= "",y="Manhattan distance") +
  scale_y_continuous(breaks = pretty_breaks(5)) +
  # Add Y -axis again
  theme(axis.line.y = element_line(size = .5)) 
gtree  
ggsave("F3_SorClust_Trend.png",plot=gtree,width=5,height=3)

# Construct ggplot manually with sample size inserted
m2 <- m
m2[lower.tri(m2)] <- NA
# Do the same for label
m.lab <- reshape2::acast(out_d,A~B,value.var = "N",fun.aggregate = sum)
m.lab[lower.tri(m.lab)] <- NA

m.lab.col <- ifelse( abs(reshape2::melt(m2)[,"value"]) >= 0.05,"white","black")

cols <- c("#d73027","#5A3F37","#834d9b","#4575b4")

gd <- reshape2::melt(m2) %>% subset(.,complete.cases(.)) %>% 
  # Reorder 
  ggplot(.,aes(x=Var1,y=fct_rev(Var2),fill=value)) + 
  theme_tufte(base_size = 20,base_family="Arial",ticks = F) + 
  geom_tile() + coord_equal() +# Tiles
  suppressWarnings( geom_text(data=reshape2::melt(m.lab),aes(x=Var1,y=Var2,label=value),inherit.aes = F,size=5,fontface = "bold") ) +
  scale_fill_gradient2(low = "#601200", mid = "white",high = "#001260",
                       na.value = "white",midpoint = 0,breaks = pretty_breaks(5),
                       guide = guide_colorbar(title = "",title.position="top",
                                              direction = "horizontal",nbin=100,ticks = FALSE, barwidth = 18, barheight = 1)) +
  guides(color="none") + theme(legend.position = "bottom") +
  scale_x_discrete(position = "top") +
  labs(x = "", y = "",title = "") +
  #theme(line = element_blank(),aspect.ratio = .75) +
  theme(legend.margin=margin(t = -.75, unit='cm'))#  theme(legend.margin=margin(t=0, r=0, b=0, l=0, unit="cm"))
gd
ggsave("F3_SorGrid_Trend.png",plot=gd+ theme(plot.margin=unit(c(0,0,0,0), "mm")),width=6,height=6)

# --------------------------------------------- #
#                Time Bins                      # 
# --------------------------------------------- #
# Version 2 with time bins and direction
o <- out
ss <- readRDS("resSaves/Out_MatricesSor.rds")
ss <- ss[which(names(ss) %in% unique(o$SS)) ]
o$trendchange <- trendchange
o$Trendchange_direction = ifelse(o$trendchange < 0,"N","P");o$Trendchange_direction[which(is.na(o$trendchange))] <- "S"
o$Trendchange_direction <- factor(o$Trendchange_direction,c("S","N","P"))
o$Inter <- interaction(o$Break_direction,o$Trendchange_direction)
# Combine all
ol <- data.frame()
for(study in names(ss)){
  print(study)
  sub <- ss[[study]] %>% reshape2::melt() %>% 
    # Filter to only abrupt changes present in comparison
    dplyr::filter(Var1 %in% o$SSBS,Var2  %in% o$SSBS) %>% 
    left_join(., 
              o %>% dplyr::select(SSBS,LCLU,BinTime,Trendchange_direction) %>% 
                mutate(NewComb = str_c(Trendchange_direction,"_",BinTime) ) %>% mutate(NewComb = fct_relevel(NewComb,"S_ND")) %>% dplyr::select(-BinTime,-Trendchange_direction) %>% 
                # Rename
                dplyr::rename(Var1 = SSBS),
              by = "Var1") %>% rename(NewComb_1 = NewComb) %>% dplyr::select(-Var1) %>% 
    left_join(., 
              o %>% dplyr::select(SSBS,LCLU,BinTime,Trendchange_direction) %>% 
                mutate(NewComb = str_c(Trendchange_direction,"_",BinTime) ) %>% mutate(NewComb = fct_relevel(NewComb,"S_ND")) %>% dplyr::select(-BinTime,-Trendchange_direction) %>% 
                # Rename
                dplyr::rename(Var2 = SSBS),
              by = "Var2") %>% rename(NewComb_2 = NewComb) %>% dplyr::select(-Var2) %>% 
    # Make a combinated column
    mutate(NewCombFull = paste0(NewComb_1,"__",NewComb_2)) %>% dplyr::select(-NewComb_1,-NewComb_2) %>% 
    subset(.,complete.cases(.)) %>% mutate(SS = study)
  # Save
  ol <- rbind(ol, sub)
}
# Now summarize
out_d <- ol %>% mutate(LCLUCc = paste0(LCLU.x,"_",LCLU.y)) %>% 
  dplyr::group_by(NewCombFull,LCLUCc) %>% 
  dplyr::summarize(avg = mean(value,na.rm=T)) %>% 
  # Only consider within same land use
  dplyr::filter(LCLUCc %in% c("PV_PV","SV_SV","HDV_HDV")) %>% 
  # Group and summarise
  dplyr::group_by(NewCombFull) %>% 
  dplyr::summarize(avg = mean(avg,na.rm=T)) %>% # Maximal value to get the number of studies contributing to a combination
  # Join sample sizes back in 
  left_join(.,ol %>% mutate(LCLUCc = paste0(LCLU.x,"_",LCLU.y)) %>%
              dplyr::filter(LCLUCc %in% c("PV_PV","SV_SV","HDV_HDV")) %>%
              dplyr::group_by(NewCombFull) %>%
              dplyr::summarize(N = n_distinct(unique(SS)) )
  ) %>% 
  # Split
  separate(NewCombFull,c("A","B"),"__") %>% 
  mutate(A = factor(A,levels = c("N_>10y","N_5-10y","N_>0-5y","S_ND", "P_>0-5y","P_5-10y","P_>10y"),labels = c("> 10\n (\u2B07)","5-10\n (\u2B07)","0-5\n (\u2B07)","0","0-5\n (\u2B06)","5-10\n (\u2B06)","> 10\n (\u2B06)")  ),
         B = factor(B,levels = c("N_>10y","N_5-10y","N_>0-5y","S_ND", "P_>0-5y","P_5-10y","P_>10y"),labels = c("> 10\n (\u2B07)","5-10\n (\u2B07)","0-5\n (\u2B07)","0","0-5\n (\u2B06)","5-10\n (\u2B06)","> 10\n (\u2B06)")  ) )
stopifnot(max(out_d$N) < n_distinct(ol$SS))  # Security check

# Rescale all values relative tothe comparison between stable sites
for(s1 in unique(out_d$A) ){
  for(s2 in unique(out_d$B) )
    if(s1 != "0" | s2 != "0" ){
      out_d$avg[which(out_d$A== s1 & out_d$B == s2)] <- (out_d$avg[which(out_d$A == s1 & out_d$B == s2)] - out_d$avg[which(out_d$A == "0" & out_d$B == "0")] )
    }
}
out_d$avg[which(out_d$A=="0" & out_d$B == "0")] <- 0

# Convert to matrix
m <- reshape2::acast(out_d,A~B,value.var = "avg")

gtree <- ggdendrogram(hc, rotate = F, size = 4, theme_dendro = FALSE) + 
  theme_tufte(base_size = 16,base_family="Gill Sans MT",ticks = T) + 
  labs(x= "",y="Manhattan distance") +
  scale_y_continuous(breaks = pretty_breaks(5)) +
  theme(axis.text.x = element_text(family = "Cambria"))
gtree  
ggsave("SF3_SorClust_TimeDir_Trend.png",plot=gtree,width=5,height=3)

# Construct ggplot manually with sample size inserted
m2 <- m
m2[lower.tri(m2)] <- NA
# Do the same for label
m.lab <- reshape2::acast(out_d,A~B,value.var = "N")
m.lab[lower.tri(m.lab)] <- NA

m.lab.col <- ifelse( abs(reshape2::melt(m2)[,"value"]) >= 0.05,"white","black")

gd <- reshape2::melt(m2) %>% subset(.,complete.cases(.)) %>% 
  # Reorder based on clustering
  ggplot(.,aes(x=Var1,y=fct_rev(Var2),fill=value)) + 
  theme_tufte(base_size = 20,base_family="Gill Sans MT",ticks = T) + 
  geom_tile() + # Tiles
  suppressWarnings( geom_text(data=reshape2::melt(m.lab),aes(x=Var1,y=Var2,label=value),inherit.aes = F,size=5,color = m.lab.col,fontface = "bold") ) +
  scale_fill_gradient2(low = "#601200", mid = "white",high = "#001260",
                       na.value = "white",midpoint = 0,breaks = pretty_breaks(5),
                       guide = guide_colorbar(title = "",title.position="top",
                                              direction = "horizontal",nbin=100,ticks = FALSE, barwidth = 18, barheight = 1)) +
  guides(color="none") + theme(legend.position = "bottom") +
  scale_x_discrete(position = "top") +
  labs(x = "", y = "",title = "") +
  theme(line = element_blank(),aspect.ratio = .75) +
  theme(axis.text = element_text(family = "Cambria")) +
  theme(legend.margin=margin(t = -.75, unit='cm'))#  theme(legend.margin=margin(t=0, r=0, b=0, l=0, unit="cm"))
gd
ggsave("SF3_SorGrid_TimeDir_Trend.png",plot=gd,width=5,height=5)

#### Figure 4  ####
# Idea:
# Bar plot with prediction error for each taxonomic group
assign("last.warning", NULL, envir = baseenv())
o <- out
# Filter only to extremes for both 
before = out$largest_trendbef*12
after = out$largest_trendaft*12
o$trendchange <- (after-before)
o$Trendchange_direction = ifelse(o$trendchange < 0,"N","P");o$Trendchange_direction[which(is.na(o$trendchange))] <- "S"
o$Trendchange_direction <- factor(o$Trendchange_direction,c("S","N","P"))
o$Inter <- interaction(o$Break_direction,o$Trendchange_direction)

full <- data.frame()
for(group in unique(o$TGrouping) ){
  print(group)
  sub <- subset(o,TGrouping == group)

  sfit <- glmer(Species_richness ~ Break_direction + (1|SS) + (1|SSB) + (1|LCLU) +(1|SSBS),data=o,family = "poisson",
                subset = TGrouping == group,na.action = na.exclude)
  sfit.overall.x <- format.results(sfit,"Break_direction",maxlevels = length(fixef(sfit))-1 ) %>% mutate(model = "SR") %>% dplyr::select(-group)
  # LA
  afit <- glmer(logabund ~  Break_direction + (1|SS) + (1|SSB) + (1|LCLU),data=o,family = "gaussian",
                subset = TGrouping == group,na.action = na.exclude)
  afit.overall.x <- format.results(afit,"Break_direction",maxlevels = length(fixef(afit))-1) %>% mutate(model = "LA") %>% dplyr::select(-group)
  # PIE
  pfit <- glmer(asPIE ~  Break_direction + (1|SS) + (1|SSB) + (1|LCLU),data=o,family = "gaussian",
                subset = TGrouping == group,na.action = na.exclude)
  pfit.overall.x <- format.results(pfit,"Break_direction",maxlevels = length(fixef(pfit))-1) %>% mutate(model = "PIE") %>% dplyr::select(-group)
  
  # Save | Combine
  full <- rbind(full,
                rbind(sfit.overall.x,afit.overall.x,pfit.overall.x) %>% dplyr::mutate(TGrouping = group,Type = "Magnitude")  
  )
  rm(sfit.overall.x,afit.overall.x,pfit.overall.x)
  # --- #
  # Trend change direction
  sfit <- glmer(Species_richness ~ Trendchange_direction + (1|SS) + (1|SSB) + (1|LCLU) +(1|SSBS),data=o,family = "poisson",
                subset = TGrouping == group,na.action = na.exclude)
  sfit.overall.x <- format.results(sfit,"Trendchange_direction",maxlevels = length(fixef(sfit))-1 ) %>% mutate(model = "SR") %>% dplyr::select(-group)
  # LA
  afit <- glmer(logabund ~  Trendchange_direction + (1|SS) + (1|SSB) + (1|LCLU),data=o,family = "gaussian",
                subset = TGrouping == group,na.action = na.exclude)
  afit.overall.x <- format.results(afit,"Trendchange_direction",maxlevels = length(fixef(afit))-1) %>% mutate(model = "LA") %>% dplyr::select(-group)
  # PIE
  pfit <- glmer(asPIE ~  Trendchange_direction + (1|SS) + (1|SSB) + (1|LCLU),data=o,family = "gaussian",
                subset = TGrouping == group,na.action = na.exclude)
  pfit.overall.x <- format.results(pfit,"Trendchange_direction",maxlevels = length(fixef(pfit))-1) %>% mutate(model = "PIE") %>% dplyr::select(-group)
  
  # Save | Combine
  full <- rbind(full,
                rbind(sfit.overall.x,afit.overall.x,pfit.overall.x) %>% dplyr::mutate(TGrouping = group,Type = "Trend")  
  )
  
}

# ------- #
# Reformat and convert
full$term <- str_replace(full$term,"Trendchange_direction","Break_direction")
full$term <- str_split(full$term,"Break_direction",simplify = T)[,2]

full <- full[which(!is.na(full$term)),] # Remove Intercepts
full <- full[which(full$term != ""),]
full$term <- factor(full$term,c("N","P"), labels = c("<",">") ) #c("\u&uarr;","\u&darr;") )

# Metric
full$model <- factor(full$model,levels=c("SR","LA","PIE"),
                     labels = c("Species richness", "Total abundance", "Probability\n of interspecific encounter")
                     )
# Reorder taxonomic grouping
full$TGrouping <- factor(full$TGrouping, levels = c(
  "Plants","Fungi","Invertebrates.ground","Invertebrates.flying","Amphibians","Reptiles","Birds","Mammals"
  ),
  labels = c("Plants","Fungi","Invertebrates (ground)", "Invertebrates (flying)","Amphibians","Reptiles","Birds","Mammals")
)

# colors for taxonomic groups
cols <- c("#7fbf7b","wheat1","#FF7F0E","#ffc425","#5D959A","tomato2","tan4","grey50")
full$Type_term <- paste(full$Type,full$term)
cols <- c("Magnitude >" ="#1a2a6c","Magnitude <"= "#b21f1f","Trend >" = "white","Trend <" = "white") # Blue and darkred
cols2 <- c("Magnitude >" ="#1a2a6c","Magnitude <"= "#b21f1f","Trend >" = "#1a2a6c","Trend <" = "#b21f1f") # Blue and darkred for borders
# Relative to intercept
full[,c("estimate","se.low","se.high")] <- (full[,c("estimate","se.low","se.high")]  - 1)*100

# For label adjustment
full$change <- paste0( round(full$estimate,2), "\n", full$p.stars )
full$ypos <- ifelse((full$se.low) > 0, -4.5, full$se.low -4.5  )
# For Study size
full$nSS.y <- rep(c(-42,-42,-32,-32,-17,-17), length.out =nrow(full))

# For reordering directions
full$Type_term <- factor(full$Type_term,c("Magnitude <","Trend <","Magnitude >","Trend >"))

# Save data for export
Figure4 <- full

# Plot - Barplot colored by tax. group
gB <- ggplot(full,
             aes(x=TGrouping,y=estimate, ymin=se.low, ymax=se.high,group=Type_term,fill=Type_term,shape=Type,color=term) ) + 
  theme_tufte(base_size = 20,base_family="Gill Sans MT",ticks = T) + theme(axis.ticks.x = element_blank()) +
  theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank() ) +
  # barplots with errorbars
  geom_hline(yintercept = 0,linetype="solid",color = "lightgrey") +
  # Point model
  geom_pointrange(fatten = 3,position = position_dodge(.75),size=1.5) + 
  geom_point(aes(x=0.5,y=0),size=3.5,color="black",inherit.aes = F) +
  scale_shape_manual(values = c(15,18)) + theme(legend.position = "bottom",legend.text = element_text(size=12)) +
  geom_vline(xintercept =  seq(1.5,7.5),size=.5,linetype="dotted" ) +
  scale_y_continuous(breaks=pretty_breaks(5)) +
  scale_color_manual(values = c(redBlue[1],redBlue[7]),guide = guide_legend("",ncol=4)) +  theme(legend.position = "bottom",legend.text = element_text(size=12)) +
  facet_wrap(~model,scales = "free_y",as.table = T,ncol=1, strip.position = "left") + 
  theme(strip.text = element_text(size = 18), strip.background = element_blank(), strip.placement = "outside") +#, ,panel.spacing.y = unit(-.25, "lines"),panel.spacing.x = unit(-.5, "lines")) +
#  annotate("text", -Inf, Inf, label = "Top-left", ) +
  # Add number of studies
  geom_text(aes(label=nSS, y = nSS.y ), colour="black", size=6,hjust=0.5,nudge_x = .1) + #y=min(estimate) - 0.035*min(estimate)
  # Add significance symbols
  geom_text(aes(x = TGrouping, y=ypos, label = p.stars,group=Type_term),position = position_dodge(.75), size = 5, colour = "black") +
  labs(x="",y="Difference (% \u00B1 1 SE)") + 
  theme(axis.text.x = element_text(size = 20,angle=90,hjust=1,vjust=.5),axis.ticks.x = element_blank()) + 
  theme(axis.text.x = element_blank()) +
  theme(panel.spacing.x =  unit(1, "lines")) +
  guides(color = "none", shape = "none",fill="none") +
  # Add Y -axis again
  theme(axis.line.y = element_line(size = .5)) 
gB
ggplot2::ggsave("Figure4.png",plot = gB,width = 21,height = 23,units = "cm",dpi = 400)


#### Data export for repository ####
full_data <- o
save(full_data,
     Figure2_part1,Figure2_part2,
     Figure3_part1,Figure3_part2,
     Figure4,
     file = 'C:/Users/Martin/PhD/Projects/P5_MagnitudeBreakpoints/writeup/Submission/revision_nr2/final/SupportingData.RData')



#### Figure 1  ####
# Make the map
library(tmap);library(sp)
s <- out
s <- s %>% dplyr::filter(Break_binom == 1)
s$BinMagn <- factor(s$BinMagn, levels = c("< -50%","-50% <> -25%","-25% <> 0%","ND","0% <> 25%", "25% <> 50%","> 50%"),
                    labels = c("---","--","-","0","+","++","+++") )
s$BinMagn <- droplevels(s$BinMagn)
s$BinTrend <- fct_collapse(s$BinTrend,
                           "---" = c("< -50%","Largest NTC"),
                           "--" = c("-50% <> -25%","Large NTC"),
                           "-" = c("-25% <> 0%","Small NTC"),
                           "UC" = c("ND"),
                           "+" = c("0% <> 25%","Small PTC"),
                           "++" = c("25% <> 50%","Large PTC"),
                           "+++" = c("> 50%","Largest PTC")
)
s$Trend_direction <- factor( ifelse(s$BinTrend %in% c("---","--","-"),"-","+" ) )
s$Break_direction <- factor(s$Break_direction, levels = c("N","P"),labels = c("-","+") )
# Unified legend for both
s$MagTrend <- factor( paste0(s$Break_direction, "|",s$Trend_direction) )

# ------------------------------- #
s$BinMagn <- droplevels(s$BinMagn)
s$BinTrend <- droplevels(s$BinTrend)
redBlue2 <- c("#d73027","#fc8d59","#fee090","#e0f3f8","#91bfdb","#4575b4")
s$New <- paste(s$Break_direction,s$BinTime,sep = "_")
s$New = factor(s$New,levels = c("N_>10y","N_5-10y","N_>0-5y","S_ND","P_>0-5y","P_5-10y","P_>10y"),labels = c(">10y (<)","5-10y (<)",">0-5y (<)","UC",">0-5y (>)","5-10y (>)",">10y (>)") ) # Backtransform to factor
s$BinTime <- factor(s$BinTime,levels = c("ND",">0-5y","5-10y",">10y"),labels = c("UC","0-5y","5-10y",">10y") ) # Backtransform to factor

coordinates(s) = ~Longitude+Latitude
proj4string(s) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
data(World)
laea <- "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

cols <- c("#d73027","#5A3F37","#834d9b","#4575b4")
  
w <- subset(World, continent != "Antarctica")
mp <- tm_shape(w, projection="eck4",col="transparent") + 
#  tm_grid(projection="longlat", y = c(-75,-50,-25,0,25,50,75), x = 179.99999,  
#          n.x = 0, n.y = 25, lwd = .5,col="black") +
  tm_polygons(col="grey70",alpha = 1,border.col = "grey70") + 
  tm_shape(s) + tm_bubbles(
    # Magnitude
    col="BinMagn",palette = redBlue2, 
    shape = 15,
    # Trend
#    col="BinTrend",palette = redBlue2, 
#    shape = 18,
                           title.col="",scale = .75,
                           border.col = "white", border.alpha = .5,alpha=.7,
                           size= .5,
                           legend.hist = F,legend.col.is.portrait=T,legend.hist.title = "") +
  tm_format("World",inner.margins=c(.04,.03, .02, .01),#outer.margins=c(.01,.01,.01,.01), 
                       title.size=2,title = paste0(""),title.bg.color="white",
                       title.position=c("right","top"),title.bg.alpha=0.75,frame=F,
                       bg.color="white",earth.boundary = F,earth.boundary.lwd =.5,space.color="transparent",
                        
                       legend.text.size = 2,legend.just = c("bottom"),
                       legend.show=F,legend.outside = T)
mp
tmap_save(mp, paste0("Figure1.png"),asp = 0, width=3000, height=3000,dpi = 400,bg="transparent")

# ------------------------------------------------------ #
