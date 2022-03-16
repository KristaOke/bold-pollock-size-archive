#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#models to run on cluster
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Notes:
#-these models are too big to run on Krista's computer, Mike is running them on a cluster
#-they rely on data and models saved elsewhere
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#packages
library(mgcv)
library(gamm4)
library(cAIC4)


#data
wd <- getwd()
lagdat <- read.csv(file=paste(wd,"/data/lagdat.csv", sep=""), row.names=1)

lagdat$cohort <- lagdat$YEAR - lagdat$AGE

lagdat <- lagdat[which(lagdat$AGE<15),] 
lagdat$AGE <- as.factor(lagdat$AGE)
lagdat$cohort <- as.factor(lagdat$cohort)

#models-------
#full model all ages
# rall1 <- gamm4(log_sc_weight ~  s(sst.amj, by=AGE, k=4) + t2(LONGITUDE, LATITUDE) + s(julian, k = 4),
#                random=~(1|YEAR/HAUL), data=lagdat) 
# saveRDS(rall1, file=paste(wd,"/scripts/size scripts/model_output_cluster/model_output_all-ages_no-cohort.rds", sep=""))
# rall1_cAIC <- cAIC(rall1)
# saveRDS(rall1_cAIC, file=paste(wd,"/scripts/size scripts/model_output_cluster/cAIC_all-ages_no-cohort.rds", sep=""))

#full model all ages with cohort

# cor_all1 <- gamm4(log_sc_weight ~  s(sst.amj, by=AGE, k=4) + t2(LONGITUDE, LATITUDE) + s(julian, k = 4),
#                   random=~(1|YEAR/HAUL) + (1|cohort), data=lagdat) 
# saveRDS(cor_all1, file=paste(wd,"/scripts/size scripts/model_output_all-ages_random-cohort.rds", sep=""))
# cor_all1_cAIC <- cAIC(cor_all1)
# saveRDS(cor_all1_cAIC, file=paste(wd,"/scripts/size scripts/model_output_cluster/cAIC_all-ages_random-cohort.rds", sep=""))

#full model all ages with cohort no limits on k

cor_freek <- gamm4(log_sc_weight ~  s(sst.amj, by=AGE) + t2(LONGITUDE, LATITUDE) + s(julian),
                  random=~(1|YEAR/HAUL) + (1|cohort), data=lagdat) 
saveRDS(cor_freek, file=paste(wd,"/scripts/size scripts/model_output_all-ages_random-cohort_free-k.rds", sep=""))
cor_freek_cAIC <- cAIC(cor_freek)
saveRDS(cor_freek_cAIC, file=paste(wd,"/scripts/size scripts/model_output_cluster/cAIC_all-ages_random-cohort_free-k.rds", sep=""))

#full model all ages w global

# G_all1 <- gamm4(log_sc_weight ~  s(sst.amj, k=4) + s(sst.amj, by=AGE, k=4) + t2(LONGITUDE, LATITUDE) + s(julian, k = 4),
#                   random=~(1|YEAR/HAUL) + (1|cohort), data=lagdat) 
# saveRDS(G_all1, file=paste(wd,"/scripts/size scripts/model_output_all-ages_random-cohort_wglobal.rds", sep=""))
G_all1_cAIC <- cAIC(G_all1)
saveRDS(G_all1_cAIC, file=paste(wd,"/scripts/size scripts/model_output_cluster/cAIC_all-ages_random-cohort_wglobal.rds", sep=""))

G_all1 <- gamm4(log_sc_weight ~  s(sst.amj, k=4) + s(sst.amj, by=AGE, k=4) + t2(LONGITUDE, LATITUDE) + s(julian, k = 4) +
                  s(cohort, bs="re"),
                random=~(1|YEAR/HAUL), data=lagdat) 
saveRDS(G_all1, file=paste(wd,"/scripts/size scripts/model_output_all-ages_cohort-as-re_wglobal.rds", sep=""))



#refit w ML----

#full model all ages with cohort

cre_all1ML <- gamm4(log_sc_weight ~  s(sst.amj, by=AGE, k=4) + t2(LONGITUDE, LATITUDE) + s(julian, k = 4) +
                    s(cohort, bs="re"),
                  random=~(1|YEAR/HAUL) , data=lagdat, REML=FALSE) 
saveRDS(cre_all1ML, file=paste(wd,"/scripts/size scripts/model_output_all-ages_re-cohort_ML.rds", sep=""))
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
saveRDS(GonlyML, file=paste(wd,"/scripts/size scripts/model_output_all-ages_re-cohort_onlyglobal_ML.rds", sep=""))
GonlyML_cAIC <- cAIC(GonlyML)
saveRDS(GonlyML_cAIC, file=paste(wd,"/scripts/size scripts/model_output_cluster/cAIC_all-ages_re-cohort_onlyglobal_ML.rds", sep=""))
MuMIn::AICc(GonlyML$mer)

#how does global only look without limit on k?
# GfreekML <- gamm4(log_sc_weight ~  s(sst.amj) +  t2(LONGITUDE, LATITUDE) + s(julian),
#                  random=~(1|YEAR/HAUL) + (1|cohort), data=lagdat, REML=FALSE) 
# saveRDS(GfreekML, file=paste(wd,"/scripts/size scripts/model_output_all-ages_random-cohort_onlyglobalfreek_ML.rds", sep=""))
# GfreekML_cAIC <- cAIC(GofreekML)
# saveRDS(GfreekML_cAIC, file=paste(wd,"/scripts/size scripts/model_output_cluster/cAIC_all-ages_random-cohort_onlyglobalfreek_ML.rds", sep=""))

#no sst
nosstML <- gamm4(log_sc_weight ~  t2(LONGITUDE, LATITUDE) + s(julian, k = 4) +
                   s(cohort, bs="re"),
                 random=~(1|YEAR/HAUL), data=lagdat, REML=FALSE) 
saveRDS(nosstML, file=paste(wd,"/scripts/size scripts/model_output_all-ages_re-cohort_nosst_ML.rds", sep=""))
MuMIn::AICc(nosstML$mer)

#allow slightly higher k

cre_all6ML <- gamm4(log_sc_weight ~  s(sst.amj, by=AGE, k=6) + t2(LONGITUDE, LATITUDE) + s(julian, k = 6) +
                      s(cohort, bs="re"),
                    random=~(1|YEAR/HAUL) , data=lagdat, REML=FALSE) 
saveRDS(cre_all6ML, file=paste(wd,"/scripts/size scripts/model_output_all-ages_re-cohort_ML_k6.rds", sep=""))
# cre_all6ML_cAIC <- cAIC(cre_all6ML)
# saveRDS(cre_all6ML_cAIC, file=paste(wd,"/scripts/size scripts/model_output_cluster/cAIC_all-ages_re-cohort_ML_k6.rds", sep=""))
MuMIn::AICc(cre_all6ML$mer)


#fewer ages models-----

subdat11 <- lagdat[which(lagdat$AGE!="12" &
                           lagdat$AGE!="13" & 
                           lagdat$AGE!="14" &
                           lagdat$AGE!="15"),]
cre_11 <- gamm4(log_sc_weight ~  s(sst.amj, by=AGE, k=4) + t2(LONGITUDE, LATITUDE) + s(julian, k = 4),
                random=~(1|YEAR/HAUL) + (1|cohort), data=subdat11) 
saveRDS(cre_11, file=paste(wd,"/scripts/size scripts/model_output_age11under_re-cohort.rds", sep=""))
# cor_11_cAIC <- cAIC(cor_11)
# saveRDS(cor_11_cAIC, file=paste(wd,"/scripts/size scripts/model_output_cluster/cAIC_age11under_random-cohort.rds", sep=""))


subdat9 <- subdat11[which(subdat11$AGE!="11" &
                            subdat11$AGE!="10"),]
cre_9 <- gamm4(log_sc_weight ~  s(sst.amj, by=AGE, k=4) + t2(LONGITUDE, LATITUDE) + s(julian, k = 4),
               random=~(1|YEAR/HAUL) + (1|cohort), data=subdat9) 
saveRDS(cre_9, file=paste(wd,"/scripts/size scripts/model_output_age9under_re-cohort.rds", sep=""))
# cor_9_cAIC <- cAIC(cor_9)
# saveRDS(cor_9_cAIC, file=paste(wd,"/scripts/size scripts/model_output_cluster/cAIC_age9under_random-cohort.rds", sep=""))




