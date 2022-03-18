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
lagdat <- lagdat[which(lagdat$AGE<16),] #was asked if borrowing info across yrs would allow
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
saveRDS(rall1, file=paste(wd,"/scripts/size scripts/model_output_all-ages_no-cohort.rds", sep=""))


vb1 <- visreg(rall1$gam, "sst.amj", scale="response",ylab="Scaled log(weight-at-age)", xlab="April-June SST", rug=1)

rall1$gam$data <- lagdat
vr1 <- visreg(rall1$gam, "sst.amj", scale="response",ylab="Scaled log(weight-at-age)", xlab="April-June SST", rug=1)
plot(rall1$gam)


plot_smooths(
  model = rall1$gam,
  series = sst.amj,
  comparison = AGE
) 


#add cohort to the all ages model----

#cohort as lmer random

cor_all1 <- gamm4(log_sc_weight ~  s(sst.amj, by=AGE, k=4) + t2(LONGITUDE, LATITUDE) + s(julian, k = 4),
                 random=~(1|YEAR/HAUL) + (1|cohort), data=lagdat) 
gam.check(cor_all1$gam)
summary(cor_all1$gam)
summary(cor_all1$mer)
AIC(cor_all1$mer) 
plot(cor_all1$mer)


saveRDS(cor_all1, file=paste(wd,"/scripts/size scripts/model_output_all-ages_random-cohort.rds", sep=""))

vba <- visreg(cor_all1$gam, "sst.amj", "AGE", scale="response",ylab="Partial effect on scaled log(weight-at-age)", xlab="April-June SST", rug=1,
              data=lagdat)



plot_smooths(
  model = cor_all1$gam,
  series = sst.amj,
  comparison = AGE,
  facet_terms = AGE
) 

#cohort as mgcv random----

#using the s(,bs="re) is apparently an equivalent but much more efficient way to fit these 
#which helps because I'm running out of memory otherwise
cre_all1 <- gamm4(log_sc_weight ~  s(sst.amj, by=AGE, k=4) + t2(LONGITUDE, LATITUDE) + s(julian, k = 4) +
                    s(cohort, bs="re"),
                  random=~(1|YEAR/HAUL), data=lagdat) 
saveRDS(cre_all1, file=paste(wd,"/scripts/size scripts/model_output_all-ages_cohort-as-re_age15.rds", sep=""))

cre_all1k8 <- gamm4(log_sc_weight ~  s(sst.amj, by=AGE, k=8) + t2(LONGITUDE, LATITUDE) + s(julian, k = 8) +
                    s(cohort, bs="re"),
                  random=~(1|YEAR/HAUL), data=lagdat) 
saveRDS(cre_all1k8, file=paste(wd,"/scripts/size scripts/model_output_all-ages_cohort-as-re_k8_age15.rds", sep=""))

#old
# glob_all1 <- gamm4(log_sc_weight ~  s(sst.amj, k=4) + t2(LONGITUDE, LATITUDE) + s(julian, k = 4),
#                   random=~(1|YEAR/HAUL) + (1|cohort), data=lagdat) 
# 
# cor_all1ML <- gamm4(log_sc_weight ~  s(sst.amj, by=AGE, k=4) + t2(LONGITUDE, LATITUDE) + s(julian, k = 4),
#                   random=~(1|YEAR/HAUL) + (1|cohort), data=lagdat, REML=FALSE) 
# saveRDS(cor_all1ML, file=paste(wd,"/scripts/size scripts/model_output_all-ages_random-cohort_wML_age15.rds", sep=""))


# glob_all1ML <- gamm4(log_sc_weight ~  s(sst.amj, k=4) + t2(LONGITUDE, LATITUDE) + s(julian, k = 4),
#                    random=~(1|YEAR/HAUL) + (1|cohort), data=lagdat, REML=FALSE) 

#no limit on k
# cor_freek <- gamm4(log_sc_weight ~  s(sst.amj, by=AGE) + t2(LONGITUDE, LATITUDE) + s(julian, k=4),
#                   random=~(1|YEAR/HAUL) + (1|cohort), data=lagdat) 
# gam.check(cor_freek$gam) #is edf close to k?

#double k
# dblk <- gamm4(log_sc_weight ~  s(sst.amj, by=AGE, k=8) + t2(LONGITUDE, LATITUDE) + s(julian, k=8),
#                    random=~(1|YEAR/HAUL) + (1|cohort), data=lagdat) 
# saveRDS(dblk, file=paste(wd,"/scripts/size scripts/model_output_all-ages_no-cohort_k-doubled_age15.rds", sep=""))


#no limit on k, no cohort
#running without cohort quickly to look at k because cohort models 
#are SUPER SLOW and memory intensive
# freek <- gamm4(log_sc_weight ~  s(sst.amj, by=AGE) + t2(LONGITUDE, LATITUDE) + s(julian),
#                    random=~(1|YEAR/HAUL), data=lagdat) 
# gam.check(freek$gam) #is edf close to k?
# saveRDS(freek, file=paste(wd,"/scripts/size scripts/model_output_all-ages_no-cohort_no-k-limit_age15.rds", sep=""))

#full model all ages w global----

# G_all1 <- gamm4(log_sc_weight ~  s(sst.amj, k=4) + s(sst.amj, by=AGE, k=4) + t2(LONGITUDE, LATITUDE) + s(julian, k = 4),
#                   random=~(1|YEAR/HAUL) + (1|cohort), data=lagdat) 
# saveRDS(G_all1, file=paste(wd,"/scripts/size scripts/model_output_all-ages_random-cohort_wglobal.rds", sep=""))
G_all1_cAIC <- cAIC(G_all1)
saveRDS(G_all1_cAIC, file=paste(wd,"/scripts/size scripts/model_output_cluster/cAIC_all-ages_random-cohort_wglobal.rds", sep=""))

G_all1 <- gamm4(log_sc_weight ~  s(sst.amj, k=4) + s(sst.amj, by=AGE, k=4) + t2(LONGITUDE, LATITUDE) + s(julian, k = 4) +
                  s(cohort, bs="re"),
                random=~(1|YEAR/HAUL), data=lagdat) 
saveRDS(G_all1, file=paste(wd,"/scripts/size scripts/model_output_all-ages_cohort-as-re_wglobal_age15.rds", sep=""))

#no sst
nosst <- gamm4(log_sc_weight ~  t2(LONGITUDE, LATITUDE) + s(julian, k = 4) +
                   s(cohort, bs="re"),
                 random=~(1|YEAR/HAUL), data=lagdat) 
saveRDS(nosst, file=paste(wd,"/scripts/size scripts/model_output_all-ages_re-cohort_nosst_age15.rds", sep=""))

#refit w ML----

#full model all ages with cohort

cre_all1ML <- gamm4(log_sc_weight ~  s(sst.amj, by=AGE, k=4) + t2(LONGITUDE, LATITUDE) + s(julian, k = 4) +
                      s(cohort, bs="re"),
                    random=~(1|YEAR/HAUL) , data=lagdat, REML=FALSE) 
saveRDS(cre_all1ML, file=paste(wd,"/scripts/size scripts/model_output_all-ages_re-cohort_ML_age15.rds", sep=""))
cre_all1ML_cAIC <- cAIC(cre_all1ML)
saveRDS(cre_all1ML_cAIC, file=paste(wd,"/scripts/size scripts/model_output_cluster/cAIC_all-ages_re-cohort_ML.rds", sep=""))
MuMIn::AICc(cre_all1ML$mer)

#full model all ages w global
#DOES NOT run I think global smoother isn't estimable with the age specific ones included, not surprising
# G_all1ML <- gamm4(log_sc_weight ~  s(sst.amj, k=4) + s(sst.amj, by=AGE, k=4) + t2(LONGITUDE, LATITUDE) + s(julian, k = 4) +
#                   s(cohort, bs="re"),
#                 random=~(1|YEAR/HAUL), data=lagdat, REML=FALSE) 
# saveRDS(G_all1ML, file=paste(wd,"/scripts/size scripts/model_output_all-ages_re-cohort_wglobal_ML.rds", sep=""))
# G_all1ML_cAIC <- cAIC(G_all1ML)
# saveRDS(G_all1ML_cAIC, file=paste(wd,"/scripts/size scripts/model_output_cluster/cAIC_all-ages_re-cohort_wglobal_ML.rds", sep=""))

#global only model
GonlyML <- gamm4(log_sc_weight ~  s(sst.amj, k=4) +  t2(LONGITUDE, LATITUDE) + s(julian, k = 4) +
                   s(cohort, bs="re"),
                 random=~(1|YEAR/HAUL), data=lagdat, REML=FALSE) 
saveRDS(GonlyML, file=paste(wd,"/scripts/size scripts/model_output_all-ages_re-cohort_onlyglobal_ML_age15.rds", sep=""))
GonlyML_cAIC <- cAIC(GonlyML)
saveRDS(GonlyML_cAIC, file=paste(wd,"/scripts/size scripts/model_output_cluster/cAIC_all-ages_re-cohort_onlyglobal_ML.rds", sep=""))
MuMIn::AICc(GonlyML$mer)



#no sst ML
nosstML <- gamm4(log_sc_weight ~  t2(LONGITUDE, LATITUDE) + s(julian, k = 4) +
                   s(cohort, bs="re"),
                 random=~(1|YEAR/HAUL), data=lagdat, REML=FALSE) 
saveRDS(nosstML, file=paste(wd,"/scripts/size scripts/model_output_all-ages_re-cohort_nosst_ML_age15.rds", sep=""))
MuMIn::AICc(nosstML$mer)

#allow slightly higher k

cre_all6ML <- gamm4(log_sc_weight ~  s(sst.amj, by=AGE, k=6) + t2(LONGITUDE, LATITUDE) + s(julian, k = 6) +
                      s(cohort, bs="re"),
                    random=~(1|YEAR/HAUL) , data=lagdat, REML=FALSE) 
saveRDS(cre_all6ML, file=paste(wd,"/scripts/size scripts/model_output_all-ages_re-cohort_ML_k6_age15.rds", sep=""))
# cre_all6ML_cAIC <- cAIC(cre_all6ML)
# saveRDS(cre_all6ML_cAIC, file=paste(wd,"/scripts/size scripts/model_output_cluster/cAIC_all-ages_re-cohort_ML_k6.rds", sep=""))
MuMIn::AICc(cre_all6ML$mer)
