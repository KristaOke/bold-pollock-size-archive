#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#All ages modeled together

#March-April 2022
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Notes:
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#load packages and data----

library(tidymv)
library(mgcv)
library(gamm4)
library(cAIC4)
library(visreg)
library(tidymv)
library(mgcViz)

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
#saveRDS(cre_all1, file=paste(wd,"/scripts/size scripts/model_output_all-ages_cohort-as-re_age15.rds", sep=""))
cre_all1 <- readRDS(file=paste(wd,"/scripts/size scripts/model_output_all-ages_cohort-as-re_age15.rds", sep=""))
gam.check(cre_all1$gam)
R2.cre_all1 <- 1-var(residuals(cre_all1$gam))/(var(model.response(model.frame(cre_all1$gam))))

cre_all1k8 <- gamm4(log_sc_weight ~  s(sst.amj, by=AGE, k=8) + t2(LONGITUDE, LATITUDE) + s(julian, k = 8) +
                    s(cohort, bs="re"),
                  random=~(1|YEAR/HAUL), data=lagdat) 
#saveRDS(cre_all1k8, file=paste(wd,"/scripts/size scripts/model_output_all-ages_cohort-as-re_k8_age15.rds", sep=""))
cre_all1k8 <- readRDS(file=paste(wd,"/scripts/size scripts/model_output_all-ages_cohort-as-re_k8_age15.rds", sep=""))
gam.check(cre_all1k8$gam)
R2.cre_all8 <- 1-var(residuals(cre_all1k8$gam))/(var(model.response(model.frame(cre_all1k8$gam))))


#old models------
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
# G_all1_cAIC <- cAIC(G_all1)
# saveRDS(G_all1_cAIC, file=paste(wd,"/scripts/size scripts/model_output_cluster/cAIC_all-ages_random-cohort_wglobal.rds", sep=""))

G_all1 <- gamm4(log_sc_weight ~  s(sst.amj, k=4) + s(sst.amj, by=AGE, k=4) + t2(LONGITUDE, LATITUDE) + s(julian, k = 4) +
                  s(cohort, bs="re"),
                random=~(1|YEAR/HAUL), data=lagdat) 
#saveRDS(G_all1, file=paste(wd,"/scripts/size scripts/model_output_all-ages_cohort-as-re_wglobal_age15.rds", sep=""))
G_all1 <- readRDS(file=paste(wd,"/scripts/size scripts/model_output_all-ages_cohort-as-re_wglobal_age15.rds", sep=""))
gam.check(G_all1$gam)

#global only model
Gonly <- gamm4(log_sc_weight ~  s(sst.amj, k=4) +  t2(LONGITUDE, LATITUDE) + s(julian, k = 4) +
                   s(cohort, bs="re"),
                 random=~(1|YEAR/HAUL), data=lagdat) 
#saveRDS(Gonly, file=paste(wd,"/scripts/size scripts/model_output_all-ages_re-cohort_onlyglobal_age15.rds", sep=""))
Gonly <- readRDS(file=paste(wd,"/scripts/size scripts/model_output_all-ages_re-cohort_onlyglobal_age15.rds", sep=""))
gam.check(Gonly$gam)
R2.gonly <- 1-var(residuals(Gonly$gam))/(var(model.response(model.frame(Gonly$gam))))

#no sst
nosst <- gamm4(log_sc_weight ~  t2(LONGITUDE, LATITUDE) + s(julian, k = 4) +
                   s(cohort, bs="re"),
                 random=~(1|YEAR/HAUL), data=lagdat) 
#saveRDS(nosst, file=paste(wd,"/scripts/size scripts/model_output_all-ages_re-cohort_nosst_age15.rds", sep=""))
nosst <- readRDS(file=paste(wd,"/scripts/size scripts/model_output_all-ages_re-cohort_nosst_age15.rds", sep=""))
R2.nosst <- 1-var(residuals(nosst$gam))/(var(model.response(model.frame(nosst$gam))))


#refit w ML----

#full model all ages with cohort

cre_all1ML <- gamm4(log_sc_weight ~  s(sst.amj, by=AGE, k=4) + t2(LONGITUDE, LATITUDE) + s(julian, k = 4) +
                      s(cohort, bs="re"),
                    random=~(1|YEAR/HAUL) , data=lagdat, REML=FALSE) 
#saveRDS(cre_all1ML, file=paste(wd,"/scripts/size scripts/model_output_all-ages_re-cohort_ML_age15.rds", sep=""))
cre_all1ML <- readRDS(file=paste(wd,"/scripts/size scripts/model_output_all-ages_re-cohort_ML_age15.rds", sep=""))
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
#saveRDS(GonlyML, file=paste(wd,"/scripts/size scripts/model_output_all-ages_re-cohort_onlyglobal_ML_age15.rds", sep=""))
GonlyML <- readRDS(file=paste(wd,"/scripts/size scripts/model_output_all-ages_re-cohort_onlyglobal_ML_age15.rds", sep=""))
GonlyML_cAIC <- cAIC(GonlyML)
saveRDS(GonlyML_cAIC, file=paste(wd,"/scripts/size scripts/model_output_cluster/cAIC_all-ages_re-cohort_onlyglobal_ML.rds", sep=""))
MuMIn::AICc(GonlyML$mer)



#no sst ML
nosstML <- gamm4(log_sc_weight ~  t2(LONGITUDE, LATITUDE) + s(julian, k = 4) +
                   s(cohort, bs="re"),
                 random=~(1|YEAR/HAUL), data=lagdat, REML=FALSE) 
#saveRDS(nosstML, file=paste(wd,"/scripts/size scripts/model_output_all-ages_re-cohort_nosst_ML_age15.rds", sep=""))
nosstML <- readRDS(file=paste(wd,"/scripts/size scripts/model_output_all-ages_re-cohort_nosst_ML_age15.rds", sep=""))
MuMIn::AICc(nosstML$mer)


#best model-----

#cre_all1 is best model, was already fit using REML so can interpret not that we've compared fixed
#effects using ML and decided on it

summary(cre_all1$gam)
summary(cre_all1$mer)

plot_smooths(
      model = cre_all1$gam,
    series = sst.amj,
    # comparison = AGE,
    facet_terms = AGE) 

#fig 3 in manuscript----
cre_all1$gam$data <- lagdat
va1 <- visreg(cre_all1$gam, "sst.amj", by="AGE", #scale="response",
       partial=FALSE, layout=c(5,3), index.cond=list(c(11,12,13,14,15,
                                                       6,7,8,9,10,
                                                       1,2,3,4,5)),
       ylab="Partial effect on scaled log(weight-at-age)", xlab="April-June SST", rug=1)

#try to center around zero
cre_all1$gam$data <- lagdat
va2 <- visreg(cre_all1$gam, "sst.amj", by="AGE", cond=list(cohort=2003),
              partial=FALSE, layout=c(5,3), index.cond=list(c(11,12,13,14,15,
                                                              6,7,8,9,10,
                                                              1,2,3,4,5)),
              ylab="Partial effect on scaled log(weight-at-age)", xlab="April-June SST", rug=1, data=lagdat)

# visreg(cre_all1$gam, "sst.amj", by="AGE", #scale="response",
#        partial=FALSE,
#        ylab="Partial effect on scaled log(weight-at-age)", xlab="April-June SST", rug=2, gg=TRUE) + 
#   facet_wrap(~AGE, nrow=3)
#order is right but CI disappeared and y-axis rug makes it hard to see


#allow slightly higher k-----

cre_all6ML <- gamm4(log_sc_weight ~  s(sst.amj, by=AGE, k=6) + t2(LONGITUDE, LATITUDE) + s(julian, k = 6) +
                      s(cohort, bs="re"),
                    random=~(1|YEAR/HAUL) , data=lagdat, REML=FALSE) 
#saveRDS(cre_all6ML, file=paste(wd,"/scripts/size scripts/model_output_all-ages_re-cohort_ML_k6_age15.rds", sep=""))
cre_all6ML <- readRDS(file=paste(wd,"/scripts/size scripts/model_output_all-ages_re-cohort_ML_k6_age15.rds", sep=""))
# cre_all6ML_cAIC <- cAIC(cre_all6ML)
# saveRDS(cre_all6ML_cAIC, file=paste(wd,"/scripts/size scripts/model_output_cluster/cAIC_all-ages_re-cohort_ML_k6.rds", sep=""))
MuMIn::AICc(cre_all6ML$mer)


vb6 <- visreg(cre_all6ML$gam, "sst.amj", scale="response",ylab="Scaled log(weight-at-age)", xlab="April-June SST", rug=1)

cre_all6ML$gam$data <- lagdat
vr6 <- visreg(cre_all6ML$gam, "sst.amj", scale="response",ylab="Scaled log(weight-at-age)", xlab="April-June SST", rug=1)
plot(cre_all6ML$gam)


plot_smooths(
  model = cre_all6ML$gam,
  series = sst.amj,
  comparison = AGE
) 

#look at under v overfitting k=6 -------

b6 <- getViz(cre_all6ML$gam)
o1 <- plot( sm(b6, 1) )
o1 + l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
  l_ciLine(mul = 5, colour = "blue", linetype = 2) + 
  l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()

o2 <- plot( sm(b6, 2) )
o2 + l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
  l_ciLine(mul = 5, colour = "blue", linetype = 2) + 
  l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()

o3 <- plot( sm(b6, 3) )
o3 + l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
  l_ciLine(mul = 5, colour = "blue", linetype = 2) + 
  l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()


o4 <- plot( sm(b6, 4) )
o4 + l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
  l_ciLine(mul = 5, colour = "blue", linetype = 2) + 
  l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()


o5 <- plot( sm(b6, 5) )
o5 + l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
  l_ciLine(mul = 5, colour = "blue", linetype = 2) + 
  l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()


o6 <- plot( sm(b6, 6) )
o6 + l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
  l_ciLine(mul = 5, colour = "blue", linetype = 2) + 
  l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()


o7 <- plot( sm(b6, 7) )
o7 + l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
  l_ciLine(mul = 5, colour = "blue", linetype = 2) + 
  l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()


o8 <- plot( sm(b6, 8) )
o8 + l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
  l_ciLine(mul = 5, colour = "blue", linetype = 2) + 
  l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()


o9 <- plot( sm(b6, 9) )
o9 + l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
  l_ciLine(mul = 5, colour = "blue", linetype = 2) + 
  l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()


o10 <- plot( sm(b6, 10) )
o10 + l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
  l_ciLine(mul = 5, colour = "blue", linetype = 2) + 
  l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()


o11 <- plot( sm(b6, 11) )
o11 + l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
  l_ciLine(mul = 5, colour = "blue", linetype = 2) + 
  l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()


o12 <- plot( sm(b6, 12) )
o12 + l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
  l_ciLine(mul = 5, colour = "blue", linetype = 2) + 
  l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()


o13 <- plot( sm(b6, 13) )
o13 + l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
  l_ciLine(mul = 5, colour = "blue", linetype = 2) + 
  l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()


o14 <- plot( sm(b6, 14) )
o14 + l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
  l_ciLine(mul = 5, colour = "blue", linetype = 2) + 
  l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()


o15 <- plot( sm(b6, 15) )
o15 + l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
  l_ciLine(mul = 5, colour = "blue", linetype = 2) + 
  l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()



#compare to k = 8 -------

b8 <- getViz(cre_all1k8$gam)
e1 <- plot( sm(b8, 1) )
e1 + l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
  l_ciLine(mul = 5, colour = "blue", linetype = 2) + 
  l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()

e2 <- plot( sm(b8, 2) )
e2 + l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
  l_ciLine(mul = 5, colour = "blue", linetype = 2) + 
  l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()

e3 <- plot( sm(b8, 3) )
e3 + l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
  l_ciLine(mul = 5, colour = "blue", linetype = 2) + 
  l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()


e4 <- plot( sm(b8, 4) )
e4 + l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
  l_ciLine(mul = 5, colour = "blue", linetype = 2) + 
  l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()


e5 <- plot( sm(b8, 5) )
e5 + l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
  l_ciLine(mul = 5, colour = "blue", linetype = 2) + 
  l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()


e6 <- plot( sm(b8, 6) )
e6 + l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
  l_ciLine(mul = 5, colour = "blue", linetype = 2) + 
  l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()


e7 <- plot( sm(b8, 7) )
e7 + l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
  l_ciLine(mul = 5, colour = "blue", linetype = 2) + 
  l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()


e8 <- plot( sm(b8, 8) )
e8 + l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
  l_ciLine(mul = 5, colour = "blue", linetype = 2) + 
  l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()


e9 <- plot( sm(b8, 9) )
e9 + l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
  l_ciLine(mul = 5, colour = "blue", linetype = 2) + 
  l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()


e10 <- plot( sm(b8, 10) )
e10 + l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
  l_ciLine(mul = 5, colour = "blue", linetype = 2) + 
  l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()


e11 <- plot( sm(b8, 11) )
e11 + l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
  l_ciLine(mul = 5, colour = "blue", linetype = 2) + 
  l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()


e12 <- plot( sm(b8, 12) )
e12 + l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
  l_ciLine(mul = 5, colour = "blue", linetype = 2) + 
  l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()


e13 <- plot( sm(b8, 13) )
e13 + l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
  l_ciLine(mul = 5, colour = "blue", linetype = 2) + 
  l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()


e14 <- plot( sm(b8, 14) )
e14 + l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
  l_ciLine(mul = 5, colour = "blue", linetype = 2) + 
  l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()


e15 <- plot( sm(b8, 15) )
e15 + l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
  l_ciLine(mul = 5, colour = "blue", linetype = 2) + 
  l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()


#replot spatial surface-----------------------------------------------------

s1 <- getViz(cre_all1$gam)
psb1 <- plot(sm(s1, 16)) + l_fitRaster() + l_fitContour() + labs(title = NULL) + #l_points() +
  geom_polygon(data = map_data ("world"), 
               aes(x=long, y = lat,group=group),fill="white",color="black",
               inherit.aes = F)+coord_sf(xlim = c(-177, -158.5), ylim = c(54.5, 62), expand = FALSE) +
  scale_fill_distiller(palette = "Spectral", type = "div") +
 # theme(legend.position = "none")+ 
  theme(plot.margin = unit(c(0, 0, 0, 0.1), "cm"), 
                                       #  axis.title.x = element_blank(),
                                       #  axis.title.y = element_blank(),
                                          plot.title = element_blank()) + xlab("Longitude") +
                                        ylab("Latitude") + labs(fill = "Effect")

psb1 

#plot DOY
vdoy <- visreg(cre_all1$gam, "julian", ylab="Partial effect on 
scaled log(weight-at-age)", xlab="Day of year", rug=T, partial=F, data=lagdat)
