#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#In response to reviews, trying to run models with all ages and age specific smoothers
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Notes:
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#load packages and data----

library(tidymv)
library(mgcv)
library(gamm4)
library(cAIC4)

wd <- getwd()
lagdat <- read.csv(file=paste(wd,"/data/lagdat.csv", sep=""), row.names=1)

lagdat$cohort <- lagdat$YEAR - lagdat$AGE

table(lagdat$AGE)
lagdat <- lagdat[which(lagdat$AGE<15),] #was asked if borrowing info across yrs would allow
#more ages to be included

lagdat$AGE <- as.factor(lagdat$AGE)
lagdat$cohort <- as.factor(lagdat$cohort)



#with the random effect-----
rall1 <- gamm4(log_sc_weight ~  s(sst.amj, by=AGE, k=4) + t2(LONGITUDE, LATITUDE) + s(julian, k = 4),
                        random=~(1|YEAR/HAUL), data=lagdat) 
gam.check(rall1$gam)
summary(rall1$gam)
summary(rall1$mer)
AIC(rall1$mer) #
#
saveRDS(rall1, file=paste(wd,"/scripts/size scripts/model_output_all-ages_no-cohort.csv", sep=""))


vb1 <- visreg(rall1$gam, "sst.amj", scale="response",ylab="Scaled log(weight-at-age)", xlab="April-June SST", rug=1)

rall1$gam$data <- lagdat
vr1 <- visreg(rall1$gam, "sst.amj", scale="response",ylab="Scaled log(weight-at-age)", xlab="April-June SST", rug=1)
plot(rall1$gam)


plot_smooths(
  model = rall1$gam,
  series = sst.amj,
  comparison = AGE
) 


#add cohort to the all ages model??----

#cohort as random

cor_all1 <- gamm4(log_sc_weight ~  s(sst.amj, by=AGE, k=4) + t2(LONGITUDE, LATITUDE) + s(julian, k = 4),
                 random=~(1|YEAR/HAUL) + (1|cohort), data=lagdat) 
gam.check(cor_all1$gam)
summary(cor_all1$gam)
summary(cor_all1$mer)
AIC(cor_all1$mer) 
plot(cor_all1$mer)


saveRDS(cor_all1, file=paste(wd,"/scripts/size scripts/model_output_all-ages_random-cohort.csv", sep=""))

vba <- visreg(cor_all1$gam, "sst.amj", "AGE", scale="response",ylab="Partial effect on scaled log(weight-at-age)", xlab="April-June SST", rug=1,
              data=lagdat)



plot_smooths(
  model = cor_all1$gam,
  series = sst.amj,
  comparison = AGE,
  facet_terms = AGE
) 

glob_all1 <- gamm4(log_sc_weight ~  s(sst.amj, k=4) + t2(LONGITUDE, LATITUDE) + s(julian, k = 4),
                  random=~(1|YEAR/HAUL) + (1|cohort), data=lagdat) 

cor_all1ML <- gamm4(log_sc_weight ~  s(sst.amj, by=AGE, k=4) + t2(LONGITUDE, LATITUDE) + s(julian, k = 4),
                  random=~(1|YEAR/HAUL) + (1|cohort), data=lagdat, REML=FALSE) 
saveRDS(cor_all1ML, file=paste(wd,"/scripts/size scripts/model_output_all-ages_random-cohort_wML.csv", sep=""))


glob_all1ML <- gamm4(log_sc_weight ~  s(sst.amj, k=4) + t2(LONGITUDE, LATITUDE) + s(julian, k = 4),
                   random=~(1|YEAR/HAUL) + (1|cohort), data=lagdat, REML=FALSE) 

#no limit on k
cor_freek <- gamm4(log_sc_weight ~  s(sst.amj, by=AGE) + t2(LONGITUDE, LATITUDE) + s(julian, k=4),
                  random=~(1|YEAR/HAUL) + (1|cohort), data=lagdat) 
gam.check(cor_freek$gam) #is edf close to k?

#double k
dblk <- gamm4(log_sc_weight ~  s(sst.amj, by=AGE, k=8) + t2(LONGITUDE, LATITUDE) + s(julian, k=8),
                   random=~(1|YEAR/HAUL) + (1|cohort), data=lagdat) 
saveRDS(dblk, file=paste(wd,"/scripts/size scripts/model_output_all-ages_no-cohort_k-doubled.csv", sep=""))


#no limit on k, no cohort
#running without cohort quickly to look at k because cohort models 
#are SUPER SLOW and memory intensive
freek <- gamm4(log_sc_weight ~  s(sst.amj, by=AGE) + t2(LONGITUDE, LATITUDE) + s(julian),
                   random=~(1|YEAR/HAUL), data=lagdat) 
gam.check(freek$gam) #is edf close to k?
saveRDS(freek, file=paste(wd,"/scripts/size scripts/model_output_all-ages_no-cohort_no-k-limit.csv", sep=""))

#alternative: fs instead of by for large # of levles?
freekFS <- gamm4(log_sc_weight ~  s(sst.amj, AGE, bs="fs") + t2(LONGITUDE, LATITUDE) + s(julian),
               random=~(1|YEAR/HAUL), data=lagdat) 
gam.check(freekFS$gam)

freekml <- gamm4(log_sc_weight ~  s(sst.amj, by=AGE) + t2(LONGITUDE, LATITUDE) + s(julian),
               random=~(1|YEAR/HAUL), data=lagdat, REML=FALSE) 
gam.check(freekml$gam) 
summary(freekml$gam)
saveRDS(freekml, file=paste(wd,"/scripts/size scripts/model_output_all-ages_no-cohort_no-k-limit_ML.csv", sep=""))
cAIC(freekml)

freekml_sel <- gamm4(log_sc_weight ~  s(sst.amj, by=AGE) + t2(LONGITUDE, LATITUDE) + s(julian),
                 random=~(1|YEAR/HAUL), data=lagdat, REML=FALSE, select=TRUE) 

#compared to no SST

nosst_all1 <- gamm4(log_sc_weight ~  t2(LONGITUDE, LATITUDE) + s(julian, k = 4) +
                      cohort,
                    random=~(1|YEAR/HAUL), data=lagdat) #needs updating
gam.check(nosst_all1$gam)
summary(nosst_all1$gam)
summary(nosst_all1$mer)
AIC(nosst_all1$mer) 


#can we get cAIC if only og up to age 11?
subdat11 <- lagdat[which(lagdat$AGE!="12" &
                           lagdat$AGE!="13" & 
                           lagdat$AGE!="14" &
                           lagdat$AGE!="15"),]
cor_11 <- gamm4(log_sc_weight ~  s(sst.amj, by=AGE, k=4) + t2(LONGITUDE, LATITUDE) + s(julian, k = 4),
                  random=~(1|YEAR/HAUL) + (1|cohort), data=subdat11) 
gam.check(cor_11$gam)
summary(cor_11$gam)
summary(cor_11$mer)
cAIC(cor_11) 

subdat9 <- subdat11[which(subdat11$AGE!="11" &
                           subdat11$AGE!="10"),]
cor_9 <- gamm4(log_sc_weight ~  s(sst.amj, by=AGE, k=4) + t2(LONGITUDE, LATITUDE) + s(julian, k = 4),
                random=~(1|YEAR/HAUL) + (1|cohort), data=subdat9) 
