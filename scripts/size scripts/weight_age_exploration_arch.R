#=====================================================================================================
#Load, manipulate, explore data

#by both Krista and Mike
#=====================================================================================================

library(tidyverse)
library(mgcv)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(rgeos)
library(rgdal)
library(MuMIn)
library(visreg)
library(gratia)
library(ggspatial)


# set palette
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
theme_set(theme_bw())


# version with date!
dat <- read.csv(("./data/survey data/Litzow_pollock_02032021.csv"))

unique(dat$SURVEY) # EBS and NBS

# begin with fitting to EBS data only
# and only 1999-on to get good weights

dat$YEAR <- as.numeric(as.character((chron::years(dat$TOW_DATE))))

dat <- dat %>%
  filter(SURVEY == "EBS", YEAR >= 1999) 

# exploratory plots
space.plot <- dat %>%
  group_by(LATITUDE, LONGITUDE, YEAR) %>%
  summarise(size = log(n()+1))

ak <- ne_countries(scale = "medium", returnclass = "sf", continent="north america")

min <- min(space.plot$size)
max <- max(space.plot$size)

map.plot <- ggplot(ak) +
  geom_point(data=space.plot, aes(LONGITUDE, LATITUDE), size = 0.5) +
  geom_sf(fill="darkgoldenrod3", color=NA) + 
  coord_sf(xlim = c(-180, -156), ylim = c(52, 66), expand = FALSE) +
  facet_wrap(~YEAR)

map.plot # pretty light before 2006!

ggsave("./figs/age_weight_sampling_maps.png", width=7, height=9, units='in')

ggplot(filter(dat, AGE <= 15), aes(AGE)) +
  geom_histogram(fill="grey", color="black", bins=15) +
  facet_wrap(~YEAR, scales="free_y")

ggsave("./figs/age_sample_size_histograms.png", width=9, height=7, units='in')

ggplot(dat, aes(AGE, WEIGHT)) +
  geom_point() +
  facet_wrap(~YEAR)

ggplot(filter(dat, AGE <= 12), aes(YEAR, WEIGHT)) +
  geom_point() +
  facet_wrap(~AGE, scales="free_y") + 
  geom_smooth()

ggsave("./figs/age_weight_scatter.png", width=9, height=6, units='in')

# load climate data
clim.dat <- read.csv("data/climate data_arch.csv")

sst.scatter <- ggplot(clim.dat, aes(south.sst.ndjfm, south.sst.amj)) + 
  geom_point()
sst.scatter

amj.ts <- ggplot(clim.dat, aes(year, south.sst.amj)) +
  geom_line() +
  geom_point()

amj.ts + ylab("April-June SST (°C)") + xlab("Year") + theme_bw()

png("./figs/winter_spring_ebs_sst_plots.png", width=6, height = 3, units='in', res=300)
#Manuscript fig 2

ggpubr::ggarrange(sst.scatter, amj.ts, ncol=2, widths = c(0.45, 0.55))

dev.off()

dat$sst.ndjfm <- clim.dat$south.sst.ndjfm[match(dat$YEAR, clim.dat$year)] 
dat$era <- as.factor(if_else(dat$YEAR < 2014, 1, 2))
dat$sst.amj <- clim.dat$south.sst.amj[match(dat$YEAR, clim.dat$year)] 

ggplot(filter(dat, AGE <= 12), aes(sst.amj, WEIGHT)) +
  geom_point() +
  facet_wrap(~AGE, scales="free_y") + 
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 4))

## scale length by age to allow comparison / analysis across ages

scale.dat <- plyr::ddply(dat, "AGE", transform, sc.weight = scale(WEIGHT))

# check
ggplot(scale.dat, aes(WEIGHT, sc.weight)) +
  geom_point() +
  facet_wrap(~AGE, scales = "free")

# looks good!

# and examine day of year distributions
scale.dat$Date <- chron::dates(as.character(scale.dat$TOW_DATE))

scale.dat$julian <- lubridate::yday(scale.dat$Date)

ggplot(scale.dat, aes(julian)) +
  geom_histogram(color="black", fill="grey") +
  facet_wrap(~YEAR)

# distribution of weights by age/year
ggplot(filter(scale.dat, AGE <= 12), aes(sc.weight)) +
  geom_histogram(color="black", fill="grey") +
  facet_grid(AGE~YEAR, scales = "free")

ggplot(filter(scale.dat, AGE <= 12), aes(sc.weight)) +
  geom_histogram(color="black", fill="grey") +
  facet_wrap(~AGE, scales = "free")

# a bit of skew!
ggsave("./figs/age_sc.weight_hist.png", width=9, height=6, units='in')

#OLD BELOW HERE-----
#Not used in paper
#SKIP TO SECTION log weight models for R skew

## weight by year ----------------------------------------------------------

# create object to catch results and predict time evolution
mod.summary.out <- age.out <- data.frame()
new.dat <- data.frame(YEAR = 1999:2019, 
                      LATITUDE = mean(scale.dat$LATITUDE),
                      LONGITUDE = mean(scale.dat$LONGITUDE),
                      julian = mean(scale.dat$julian))

for(i in 1:12){
  
  mod <- gam(sc.weight ~ s(YEAR, k=6) + te(LATITUDE, LONGITUDE) + s(julian, k=4),
             data=filter(scale.dat, AGE==i))
  
  temp.out <- data.frame(age = i,
                         edf = summary(mod)$edf,
                         nominal_p = summary(mod)$s.table[1,4])
  
  mod.summary.out <- rbind(mod.summary.out, temp.out)
  
  temp.out <- data.frame(age = i,
                         year = 1999:2019,
                         weight_anomaly = predict.gam(mod, newdata = new.dat, type = "response"))
  
age.out <- rbind(age.out, temp.out)
  
}

age.out$age <- as.factor(age.out$age)

ggplot(age.out, aes(year, weight_anomaly, color=age)) +
  geom_line() +
  scale_color_viridis_d()

## based on the above, responses to the recent anomalies
## appear to fall into three groups:
## age 1-2, age 5-8, age 9-12

# confirm that!
ggplot(age.out, aes(year, weight_anomaly)) +
  geom_line() +
  geom_hline(yintercept = 0) +
  facet_wrap(~age)

ggsave("./figs/weight_time_series_smooths_by_age.png", width=9, height=6, units='in')

## analyze by these group
mod.summary.out <- age.out <- data.frame()

new.dat <- data.frame(YEAR = 1999:2019, 
                      LATITUDE = mean(scale.dat$LATITUDE),
                      LONGITUDE = mean(scale.dat$LONGITUDE),
                      julian = mean(scale.dat$julian))
# age 1-2
  mod <- gam(sc.weight ~ s(YEAR, k=6) + te(LATITUDE, LONGITUDE) + s(julian, k=4),
             data=filter(scale.dat, AGE %in% 1:2))
  
  summary(mod)
  
  plot(mod, se=F, resid=T, pages=1, rug=F)
  
  temp.out <- data.frame(age = "1-2",
                         edf = summary(mod)$edf,
                         nominal_p = summary(mod)$s.table[1,4])
  
  mod.summary.out <- rbind(mod.summary.out, temp.out)
  
  temp.out <- data.frame(age = "1-2",
                         year = 1999:2019,
                         weight_anomaly = predict.gam(mod, newdata = new.dat, type = "response"))
  
  age.out <- rbind(age.out, temp.out)
  

# 5-8
  mod <- gam(sc.weight ~ s(YEAR, k=6) + te(LATITUDE, LONGITUDE) + s(julian, k=4),
             data=filter(scale.dat, AGE %in% 5:8))

  summary(mod)
  plot(mod, se=F, resid=T, pages=1, rug=F)
  
  temp.out <- data.frame(age = "5-8",
                         edf = summary(mod)$edf,
                         nominal_p = summary(mod)$s.table[1,4])
  
  mod.summary.out <- rbind(mod.summary.out, temp.out)
  
  temp.out <- data.frame(age = "5-8",
                         year = 1999:2019,
                         weight_anomaly = predict.gam(mod, newdata = new.dat, type = "response"))
  
  age.out <- rbind(age.out, temp.out)
  
# 9-12
  mod <- gam(sc.weight ~ s(YEAR, k=6) + te(LATITUDE, LONGITUDE) + s(julian, k=4),
             data=filter(scale.dat, AGE %in% 9:12))
  
  summary(mod)
  plot(mod, se=F, resid=T, pages=1, rug=F)
  
  temp.out <- data.frame(age = "9-12",
                         edf = summary(mod)$edf,
                         nominal_p = summary(mod)$s.table[1,4])
  
  mod.summary.out <- rbind(mod.summary.out, temp.out)
  
  temp.out <- data.frame(age = "9-12",
                         year = 1999:2019,
                         weight_anomaly = predict.gam(mod, newdata = new.dat, type = "response"))
  
  age.out <- rbind(age.out, temp.out)
  

ggplot(age.out, aes(year, weight_anomaly, color=age)) +
  geom_line() +
  geom_hline(yintercept = 0) +
  scale_color_viridis_d()

## temperature effect -------------------------------------------------------------

# object to catch results and predict time evolution
mod.summary.out <- sst.out <- data.frame()
new.dat <- data.frame(sst.amj = seq(from = min(scale.dat$sst.amj), to = max(scale.dat$sst.amj), length.out = 100),
                      LATITUDE = mean(scale.dat$LATITUDE),
                      LONGITUDE = mean(scale.dat$LONGITUDE),
                      julian = mean(scale.dat$julian))

for(i in 1:12){
 
  mod <- gam(sc.weight ~ s(sst.amj, k=6) + te(LATITUDE, LONGITUDE) + s(julian, k = 4), data=filter(scale.dat, AGE==i))
  
  temp.out <- data.frame(age = i,
                         edf = summary(mod)$edf[1],
                         nominal_p = summary(mod)$s.table[1,4])
  
  mod.summary.out <- rbind(mod.summary.out, temp.out)
  
  temp.out <- data.frame(age = i,
                         sst.amj = seq(from = min(scale.dat$sst.amj), to = max(scale.dat$sst.amj), length.out = 100),
                         weight_anomaly = predict.gam(mod, newdata = new.dat, type = "response"))
  
  sst.out <- rbind(sst.out, temp.out)
  
}

sst.out$age <- as.factor(sst.out$age)

ggplot(sst.out, aes(sst.amj, weight_anomaly, color=age)) +
  geom_line() +
  scale_color_viridis_d()

# confirm that!
ggplot(sst.out, aes(sst.amj, weight_anomaly)) +
  geom_line() +
  geom_hline(yintercept = 0) +
  facet_wrap(~age)

ggsave("./figs/age_sc.weight_amj.sst_smooths.png", width=9, height=6, units='in')

## analyze by these group
mod.summary.out <- sst.out <- data.frame()


# age 1-2
mod <- gam(sc.weight ~ s(sst.amj, k=6) + te(LATITUDE, LONGITUDE) + s(julian, k = 4), data=filter(scale.dat, AGE %in% 1:2))

summary(mod)

temp.out <- data.frame(age = "1-2",
                       edf = summary(mod)$edf,
                       nominal_p = summary(mod)$s.table[1,4])

mod.summary.out <- rbind(mod.summary.out, temp.out)

temp.out <- data.frame(age = "1-2",
                       sst.amj = seq(from = min(scale.dat$sst.amj), to = max(scale.dat$sst.amj), length.out = 100),
                       weight_anomaly = predict.gam(mod, newdata = new.dat, type = "response"))

sst.out <- rbind(sst.out, temp.out)


# 5-8
mod <- gam(sc.weight ~ s(sst.amj, k=6) + te(LATITUDE, LONGITUDE) + s(julian, k = 4), data=filter(scale.dat, AGE %in% 5:8))

summary(mod)

temp.out <- data.frame(age = "5-8",
                       edf = summary(mod)$edf,
                       nominal_p = summary(mod)$s.table[1,4])

mod.summary.out <- rbind(mod.summary.out, temp.out)

temp.out <- data.frame(age = "5-8",
                       sst.amj = seq(from = min(scale.dat$sst.amj), to = max(scale.dat$sst.amj), length.out = 100),
                       weight_anomaly = predict.gam(mod, newdata = new.dat, type = "response"))

sst.out <- rbind(sst.out, temp.out)

# 9-12
mod <- gam(sc.weight ~ s(sst.amj, k=6) + te(LATITUDE, LONGITUDE) + s(julian, k = 4), data=filter(scale.dat, AGE %in% 9:12))

summary(mod)

temp.out <- data.frame(age = "9-12",
                       edf = summary(mod)$edf,
                       nominal_p = summary(mod)$s.table[1,4])

mod.summary.out <- rbind(mod.summary.out, temp.out)

temp.out <- data.frame(age = "9-12",
                       sst.amj = seq(from = min(scale.dat$sst.amj), to = max(scale.dat$sst.amj), length.out = 100),
                       weight_anomaly = predict.gam(mod, newdata = new.dat, type = "response"))

sst.out <- rbind(sst.out, temp.out)


ggplot(sst.out, aes(sst.amj, weight_anomaly, color=age)) +
  geom_line() +
  geom_hline(yintercept = 0) +
  scale_color_viridis_d()

ggsave("./figs/age_class_vs_sst_smooths.png", width = 6, height = 5, units = 'in')

## compare with sst*era model --------------------------------

# age 1-2
mod.null <- gam(sc.weight ~ s(sst.amj, k=6) + te(LATITUDE, LONGITUDE) + s(julian, k = 4), data=filter(scale.dat, AGE %in% 1:2))

mod.alt <- gam(sc.weight ~ sst.amj*era + te(LATITUDE, LONGITUDE) + s(julian, k = 4), data=filter(scale.dat, AGE %in% 1:2))

summary(mod.alt)

AICc_1.2 <- AICc(mod.null, mod.alt) 
AICc_1.2 # null model is better

# age 5-8
mod.null <- gam(sc.weight ~ s(sst.amj, k=6) + te(LATITUDE, LONGITUDE) + s(julian, k = 4), data=filter(scale.dat, AGE %in% 5:8))

mod.alt <- gam(sc.weight ~ sst.amj*era + te(LATITUDE, LONGITUDE) + s(julian, k = 4), data=filter(scale.dat, AGE %in% 5:8))

summary(mod.alt)

AICc_5.8 <- AICc(mod.null, mod.alt) 
AICc_5.8 # null model is better

# age 9-12
mod.null <- gam(sc.weight ~ s(sst.amj, k=6) + te(LATITUDE, LONGITUDE) + s(julian, k = 4), data=filter(scale.dat, AGE %in% 9:12))

mod.alt <- gam(sc.weight ~ sst.amj*era + te(LATITUDE, LONGITUDE) + s(julian, k = 4), data=filter(scale.dat, AGE %in% 9:12))

summary(mod.alt)

AICc_9.12 <- AICc(mod.null, mod.alt) 
AICc_9.12 # alt model is better; but nominal p-values are pretty high consider the # of obs in the model!

## add a size class:temp interaction

class.dat <- scale.dat %>%
  filter(AGE %in% c(1,2,9:12)) %>%
  mutate(age.class = as.factor(if_else(AGE %in% 1:2, "1-2", "9-12")))

mod1 <- gam(sc.weight ~ s(sst.amj, k=6) + te(LATITUDE, LONGITUDE) + s(julian, k = 4), data=class.dat)
summary(mod1)

mod2 <- gam(sc.weight ~ sst.amj*age.class + te(LATITUDE, LONGITUDE) + s(julian, k = 4), data=class.dat)
summary(mod2)

mod3 <- gam(sc.weight ~ s(sst.amj, k = 6, by = age.class) + te(LATITUDE, LONGITUDE) + s(julian, k = 4), data=class.dat)
summary(mod3)

AIC <- AICc(mod1, mod2, mod3)

AIC$model <- c("mod1", "mod2", "mod3")

AIC$delta.AICc <- AIC$AICc - min(AIC$AICc)

AIC <- AIC %>%
  arrange(delta.AICc) %>%
  select(model, df, AICc, delta.AICc)

AIC

plot(mod3, pages=1, se=F)

# so the interaction smooth model is clearly the best (mod3); 
# because the effects of date and lat, long are likely age-class specific,
# will plot models fit separately to each class

mod.1_2 <- gam(sc.weight ~ s(sst.amj, k=6) + te(LATITUDE, LONGITUDE) + s(julian, k = 4), data=filter(scale.dat, AGE %in% 1:2))
png("./figs/age_1-2_scaled_weight_ebs_gam.png", width=8, height=8, units='in', res=300)
plot(mod.1_2, pages=1, se=F)
dev.off()

mod.9_12 <- gam(sc.weight ~ s(sst.amj, k=6) + te(LATITUDE, LONGITUDE) + s(julian, k = 4), data=filter(scale.dat, AGE %in% 9:12))
png("./figs/age_9-12_scaled_weight_ebs_gam.png", width=8, height=8, units='in', res=300)
plot(mod.9_12, pages=1, se=F)
dev.off()

mod.9_12 <- gam(sc.weight ~ s(sst.amj, k=6) + te(LATITUDE, LONGITUDE) + s(julian, k = 4), data=filter(scale.dat, AGE %in% 9:12))
plot(mod.9_12, pages=1, se=F)


#log weight models for R skew====================================================================

#diagnostic plots for models above show they are very right skewed
#does analyzing log(weight) help diagnostics?
#YES diagnostic plots look MUCH BETTER


ggplot(dat[which(dat$AGE<11),], aes(WEIGHT, LENGTH)) + geom_point() + facet_wrap(~AGE)

ggplot(dat[which(dat$AGE<11),], aes(WEIGHT, LENGTH)) + geom_point() + facet_wrap(~AGE, scales="free")

#because we still want to be able to compare across ages, will now log then scale weight
scale.dat$logWEIGHT <- log(scale.dat$WEIGHT)
logscale.dat <- plyr::ddply(scale.dat, "AGE", transform, log_sc_weight = scale(logWEIGHT))

# check
ggplot(logscale.dat, aes(WEIGHT, log_sc_weight)) +
  geom_point() +
  facet_wrap(~AGE, scales = "free")

# looks good!




#repeating steps from above to compare w v w/o era*sst 

#March 4 2020 we are dropping era from the models because there just isn't enough SST variation in late era 
#the cold temps are entirely missing, seems to be leading to nonsensical models

#==
dat1_2 <- filter(logscale.dat, AGE %in% 1:2)

# age 1-2
log.null1 <- gam(log_sc_weight ~ as.factor(AGE) + s(sst.amj, k=6) + te(LATITUDE, LONGITUDE) + s(julian, k = 4), data=dat1_2)
gam.check(log.null1)

#log.alt1 <- gam(log_sc_weight ~ as.factor(AGE) + sst.amj*era + te(LATITUDE, LONGITUDE) + s(julian, k = 4), data=dat1_2)
#gam.check(log.alt1)

# log.by1 <- gam(log_sc_weight ~ as.factor(AGE) + s(sst.amj, by=era, k=6) + te(LATITUDE, LONGITUDE) + s(julian, k = 4), data=dat1_2)
# gam.check(log.by1) #

summary(log.alt1)

AICc_1.2 <- AICc(log.null1) 
AICc_1.2 # 


#==
dat5_8 <- filter(logscale.dat, AGE %in% 5:8)

# age 5-8
log.null5 <- gam(log_sc_weight ~ as.factor(AGE) + s(sst.amj, k=6) + te(LATITUDE, LONGITUDE) + s(julian, k = 4), data=dat5_8)
gam.check(log.null5) #not great

#log.alt5 <- gam(log_sc_weight ~ as.factor(AGE) + sst.amj*era + te(LATITUDE, LONGITUDE) + s(julian, k = 4), data=dat5_8)
#gam.check(log.alt5)

# log.by5 <- gam(log_sc_weight ~ as.factor(AGE) + s(sst.amj, by=era, k=6) + te(LATITUDE, LONGITUDE) + s(julian, k = 4), data=dat5_8)
# gam.check(log.by5)

summary(log.alt5)

AICc_5.8 <- AICc(log.null5) 
AICc_5.8 # 


#==
dat9_12 <- filter(logscale.dat, AGE %in% 9:12)

# age 9-12
log.null9 <- gam(log_sc_weight ~ as.factor(AGE) + s(sst.amj, k=6) + te(LATITUDE, LONGITUDE) + s(julian, k = 4), data=dat9_12)
gam.check(log.null9)

#log.alt9 <- gam(log_sc_weight ~ as.factor(AGE) + sst.amj*era + te(LATITUDE, LONGITUDE) + s(julian, k = 4), data=dat9_12)
#gam.check(log.alt9)

#log.by9 <- gam(log_sc_weight ~ as.factor(AGE) + s(sst.amj, by=era, k=6) + te(LATITUDE, LONGITUDE) + s(julian, k = 4), data=dat9_12)
#gam.check(log.by9) #bad hessian

summary(log.alt9)

AICc_9.12 <- AICc(log.null9) 
AICc_9.12 # 


#conclusions don't appear to change vs previous (non-log) models

#temp autocor?===================================================================================

#1 & 2==
#log.null1 fit above
draw(log.null1, select = 1)
draw(log.null1, select = 2)
draw(log.null1, select = 3)

#autocor?
E <- residuals(log.null1, type="deviance")
I1 <- !is.na(dat1_2$log_sc_weight)
Efull <- vector(length=length(dat1_2$log_sc_weight))
Efull <- NA
Efull[I1] <- E
plot(dat1_2$YEAR, Efull) #

res1 <- dat1_2
res1$Efull <- Efull

ggplot(res1, aes(YEAR, Efull)) + geom_point() + geom_smooth()

ggplot(res1, aes(as.factor(YEAR), Efull))  + geom_boxplot()

ggplot(res1, aes(as.factor(julian), Efull))  + geom_boxplot()

ggplot(res1, aes(sst.amj, Efull))  + geom_point()

ggplot(res1, aes(AGE, Efull))  + geom_point()



#5 - 8==
#log.null5 fit above
draw(log.null5, select = 1)
draw(log.null5, select = 2)
draw(log.null5, select = 3)

#autocor?
E <- residuals(log.null5, type="deviance")
I1 <- !is.na(dat5_8$log_sc_weight)
Efull <- vector(length=length(dat5_8$log_sc_weight))
Efull <- NA
Efull[I1] <- E
plot(dat5_8$YEAR, Efull) #

res5 <- dat5_8
res5$Efull <- Efull

ggplot(res5, aes(YEAR, Efull)) + geom_point() + geom_smooth()

ggplot(res5, aes(as.factor(YEAR), Efull))  + geom_boxplot()

ggplot(res5, aes(as.factor(julian), Efull))  + geom_boxplot()

ggplot(res5, aes(sst.amj, Efull))  + geom_point()

ggplot(res5, aes(AGE, Efull))  + geom_point()






#9 - 12==
#log.null5 fit above
draw(log.null9, select = 1)
draw(log.null9, select = 2)
draw(log.null9, select = 3)

#autocor?
E <- residuals(log.null9, type="deviance")
I1 <- !is.na(dat9_12$log_sc_weight)
Efull <- vector(length=length(dat9_12$log_sc_weight))
Efull <- NA
Efull[I1] <- E
plot(dat9_12$YEAR, Efull) #

res9 <- dat9_12
res9$Efull <- Efull

ggplot(res9, aes(YEAR, Efull)) + geom_point() + geom_smooth()

ggplot(res9, aes(as.factor(YEAR), Efull))  + geom_boxplot()

ggplot(res9, aes(as.factor(julian), Efull))  + geom_boxplot()

ggplot(res9, aes(sst.amj, Efull))  + geom_point()

ggplot(res9, aes(AGE, Efull))  + geom_point()


#last yrs weight anomaly===========================================================================

wmeans <- scale.dat %>% group_by(YEAR, AGE) %>% summarize(mean_annual_weight_at_age=mean(WEIGHT, na.rm=TRUE))

wtotalmeans <- scale.dat %>% group_by(AGE) %>% summarize(mean_overall_weight_at_age=mean(WEIGHT, na.rm=TRUE))

bothmeans <- left_join(wmeans, wtotalmeans)

#lag to same age in previous year
bothmeans$sameage_lastyr_weight_anom <- NA
bothmeans$sameage_lastyr_weight_anom <- bothmeans$mean_annual_weight_at_age - bothmeans$mean_overall_weight_at_age 

bothmeans$lag_year <- NA
bothmeans$lag_year <- bothmeans$YEAR + 1

#wmeans$temp_1yr_lag <- wmeans$mean_annual_sstamj

mergemeans <- bothmeans[,c(2,5:6)]

lagdat1 <- left_join(logscale.dat, mergemeans, by = c("YEAR" = "lag_year", "AGE"="AGE"))

#lag to previous age in previous year
lagdat1$prev_age <- NA
lagdat1$prev_age <- lagdat1$AGE - 1

mergemeans2 <- bothmeans[,c(2,5:6)]
mergemeans2$prevage_lastyr_weight_anom <- mergemeans2$sameage_lastyr_weight_anom
mergemeans2 <- mergemeans2[,c(1,3:4)]

lagdat <- left_join(lagdat1, mergemeans2, by = c("YEAR" = "lag_year", "prev_age"="AGE"))
wd <- getwd()
write.csv(lagdat, file=paste(wd,"/data/lagdat.csv", sep=""))

haultable <- table(lagdat$HAUL, lagdat$YEAR)
write.csv(haultable, file=paste(wd,"/data/haul_table_for_revision.csv", sep=""))

lag12 <- lagdat[which(lagdat$AGE<3),]
lag34 <- lagdat[which(lagdat$AGE<5 & lagdat$AGE>2),]
lag58 <- lagdat[which(lagdat$AGE<9 & lagdat$AGE>4),]
lag912 <- lagdat[which(lagdat$AGE<13 & lagdat$AGE>8),]
lag1315 <- lagdat[which(lagdat$AGE<16 & lagdat$AGE>12),]

lag1 <- lag12[which(lag12$AGE==1),]
lag2 <- lag12[which(lag12$AGE==2),]

lag3only <- lagdat[which(lagdat$AGE==3),]
lag4only <- lagdat[which(lagdat$AGE==4),]
lag5only <- lagdat[which(lagdat$AGE==5),]
lag6only <- lagdat[which(lagdat$AGE==6),]
lag7only <- lagdat[which(lagdat$AGE==7),]
lag8only <- lagdat[which(lagdat$AGE==8),]
lag9only <- lagdat[which(lagdat$AGE==9),]
lag10only <- lagdat[which(lagdat$AGE==10),]
lag11only <- lagdat[which(lagdat$AGE==11),]
lag12only <- lagdat[which(lagdat$AGE==12),]
lag13only <- lagdat[which(lagdat$AGE==13),]
lag14only <- lagdat[which(lagdat$AGE==14),]
lag15only <- lagdat[which(lagdat$AGE==15),]



#map=====

world <- ne_countries(scale = "medium", returnclass = "sf")

ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-178, -155), ylim = c(53, 63), expand = TRUE) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  geom_point(aes(LONGITUDE, LATITUDE), size=0.1, data=scale.dat) + theme_bw() + 
  theme( legend.position = c(0.87, 0.85), legend.key = element_blank(),
         legend.background=element_blank()) 


#

#sample size table----------------------------------------------

tab1 <- table(lag1$YEAR)
tab2 <- table(lag2$YEAR)
tab3 <- table(lag3only$YEAR)
tab4 <- table(lag4only$YEAR)
tab5 <- table(lag5only$YEAR)
tab6 <- table(lag6only$YEAR)
tab7 <- table(lag7only$YEAR)
tab8 <- table(lag8only$YEAR)
tab9 <- table(lag9only$YEAR)
tab10 <- table(lag10only$YEAR)
tab11 <- table(lag11only$YEAR)
tab12 <- table(lag12only$YEAR)
tab13 <- table(lag13only$YEAR)
tab14 <- table(lag14only$YEAR)
tab15 <- table(lag15only$YEAR)

all_age_table <- cbind(tab1, tab2, tab3, tab4, tab5,
                       tab6, tab7, tab8, tab9, tab10,
                       tab11, tab12, tab13, tab14, tab15)

table(lag1$STATIONID)

wd <- getwd()
write.csv(all_age_table, file=paste(wd,"/data/2022-02-03-updated_sample_size_table.csv", sep=""))

