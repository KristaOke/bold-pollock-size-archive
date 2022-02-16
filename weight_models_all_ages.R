#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#In response to reviews, trying to run models with all ages and age specific smoothers
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Notes:
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#load packages and data----

library(tidymv)
library(mgcv)
library(gamm4)

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

#no limit on k
cor_freek <- gamm4(log_sc_weight ~  s(sst.amj, by=AGE) + t2(LONGITUDE, LATITUDE) + s(julian),
                  random=~(1|YEAR/HAUL) + (1|cohort), data=lagdat) 
gam.check(cor_freek$gam) #is edf close to k?





#compared to no SST

nosst_all1 <- gamm4(log_sc_weight ~  t2(LONGITUDE, LATITUDE) + s(julian, k = 4) +
                      cohort,
                    random=~(1|YEAR/HAUL), data=lagdat) 
gam.check(nosst_all1$gam)
summary(nosst_all1$gam)
summary(nosst_all1$mer)
AIC(nosst_all1$mer) 