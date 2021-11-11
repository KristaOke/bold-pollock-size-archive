#re-doing gamm models with lag added back in

#models==
library(fRegression)
library(gratia)
library(visreg)

#datasets loaded in weight-age-exploration.R

#can't do age 1, no age 0 data

# age 2====


lag.null2_ran <- gamm4(log_sc_weight ~  s(sst.amj, k=4) + t2(LONGITUDE, LATITUDE) + s(julian, k = 4) +
                          s(prevage_lastyr_weight_anom, k=4),
                        random=~(1|YEAR/HAUL), data=lag2complete)
summary(lag.null2_ran$gam) #lag not sig
summary(lag.null2_ran$mer) #R-sq.(adj) =  0.283   
AIC(lag.null2_ran$mer) #4057.229

laglin.null2_ran <- gamm4(log_sc_weight ~  s(sst.amj, k=4) + t2(LONGITUDE, LATITUDE) + s(julian, k = 4) +
                         prevage_lastyr_weight_anom,
                       random=~(1|YEAR/HAUL), data=lag2complete)
summary(laglin.null2_ran$gam) #lag not sig
summary(laglin.null2_ran$mer) #R-sq.(adj) =  0.283   
AIC(laglin.null2_ran$mer) 


AICc_2lag <- AICc(laglin.null2_ran$mer) #4058.66
AICc_2lag 

AICc(base.null2_ran$mer) #

AICc_2base_ran <- AICc(base.null2_ran$mer) #
AICc_2base_ran #4053.251

AICc_2base_ran - AICc_2lag
#-5.4091
draw(lag.null2_ran$gam, select=4)
draw(lag.null2_ran$gam, select=1)

lag.null2_ran$gam$data <- lag2complete
lr2 <- visreg(lag.null2_ran$gam, "sst.amj", scale="response",ylab="Scaled log(weight-at-age)", xlab="April-June SST", rug=1)

tt1 <- getViz(laglin.null2_ran$gam)
plot(sm(tt1, 2)) + l_fitRaster() + l_fitContour() + labs(title = NULL) + #l_points() +
     geom_polygon(data = map_data ("world"), 
                                   aes(x=long, y = lat,group=group),fill="white",color="black",
                                  inherit.aes = F)+coord_sf(xlim = c(-177, -158.5), ylim = c(54.5, 62), expand = FALSE) +
     scale_fill_distiller(palette = "Spectral", type = "div", limits = c(-3,3)) +
     theme(legend.position = "none")+ theme(plot.margin = unit(c(0, 0, 0, 0.1), "cm"), plot.title = element_blank(),
                                                                                       axis.title.x = element_blank(),
                                                                                      axis.title.y = element_blank())


#drop terms age 2----
AICc_2lag 

R2.l2 <- 1-var(residuals(laglin.null2_ran$gam))/(var(model.response(model.frame(laglin.null2_ran$gam))))

nosst.lag2_ran <- gamm4(log_sc_weight ~   t2(LONGITUDE, LATITUDE) + s(julian, k = 4)+
                          prevage_lastyr_weight_anom,
                        random=~(1|YEAR/HAUL), data=lag2complete) 
R2.l2.nosst <- 1-var(residuals(nosst.lag2_ran$gam))/(var(model.response(model.frame(nosst.lag2_ran$gam))))
R2.l2 - R2.l2.nosst 
AICc_2lagran_nosst <- AICc(nosst.lag2_ran$mer) #
AICc_2lagran_nosst - AICc_2lag

nolatlong.lag2_ran <- gamm4(log_sc_weight ~   s(sst.amj, k=4) + s(julian, k = 4)+
                              prevage_lastyr_weight_anom,
                            random=~(1|YEAR/HAUL), data=lag2complete) 
R2.l2.nolatlong <- 1-var(residuals(nolatlong.lag2_ran$gam))/(var(model.response(model.frame(nolatlong.lag2_ran$gam))))
R2.l2 - R2.l2.nolatlong 
AICc_2lagran_nolatlong <- AICc(nolatlong.lag2_ran$mer) #
AICc_2lagran_nolatlong - AICc_2lag

nojul.lag2_ran <- gamm4(log_sc_weight ~   s(sst.amj, k=4) + t2(LONGITUDE, LATITUDE)+
                          prevage_lastyr_weight_anom,
                        random=~(1|YEAR/HAUL), data=lag2complete) 
R2.l2.nojul <- 1-var(residuals(nojul.lag2_ran$gam))/(var(model.response(model.frame(nojul.lag2_ran$gam))))
R2.l2 - R2.l2.nojul 
AICc_2lagran_nojul <- AICc(nojul.lag2_ran$mer) #
AICc_2lagran_nojul - AICc_2lag


nolag.lag2_ran <- gamm4(log_sc_weight ~   s(sst.amj, k=4) + t2(LONGITUDE, LATITUDE),
                        random=~(1|YEAR/HAUL), data=lag2complete) 
R2.l2.nolag <- 1-var(residuals(nolag.lag2_ran$gam))/(var(model.response(model.frame(nolag.lag2_ran$gam))))
R2.l2 - R2.l2.nolag 
AICc_2lagran_nolag <- AICc(nolag.lag2_ran$mer) #
AICc_2lagran_nolag - AICc_2lag






# age 3====


lag.null3_ran <- gamm4(log_sc_weight ~  s(sst.amj, k=4) + t2(LONGITUDE, LATITUDE) + s(julian, k = 4) +
                         s(prevage_lastyr_weight_anom, k=4),
                       random=~(1|YEAR/HAUL), data=lag3complete)
summary(lag.null3_ran$gam) #lag not sig
summary(lag.null3_ran$mer)
AIC(lag.null3_ran$mer) #
#R-sq.(adj) = 0.237

laglin.null3_ran <- gamm4(log_sc_weight ~  s(sst.amj, k=4) + t2(LONGITUDE, LATITUDE) + s(julian, k = 4) +
                         prevage_lastyr_weight_anom,
                       random=~(1|YEAR/HAUL), data=lag3complete)
summary(laglin.null3_ran$gam) #lag not sig
summary(laglin.null3_ran$mer)
AIC(laglin.null3_ran$mer) #
#R-sq.(adj) = 0.237

AICc_3lag <- AICc(laglin.null3_ran$mer) #
AICc_3lag #4289.492

AICc(base.null3_ran$mer) #

AICc_3base_ran <- AICc(base.null3_ran$mer) #
AICc_3base_ran #4278.452

AICc_3base_ran - AICc_3lag
#-11.039
draw(lag.null3_ran$gam, select=4)
draw(lag.null3_ran$gam, select=1)

lag.null3_ran$gam$data <- lag3complete
lr3 <- visreg(lag.null3_ran$gam, "sst.amj", scale="response",ylab="Scaled log(weight-at-age)", xlab="April-June SST", rug=1)


#drop terms age 3----
AICc_3lag 

R2.l3 <- 1-var(residuals(laglin.null3_ran$gam))/(var(model.response(model.frame(laglin.null3_ran$gam))))

nosst.lag3_ran <- gamm4(log_sc_weight ~   t2(LONGITUDE, LATITUDE) + s(julian, k = 4)+
                          prevage_lastyr_weight_anom,
                        random=~(1|YEAR/HAUL), data=lag3complete) 
R2.l3.nosst <- 1-var(residuals(nosst.lag3_ran$gam))/(var(model.response(model.frame(nosst.lag3_ran$gam))))
R2.l3 - R2.l3.nosst 
AICc_3lagran_nosst <- AICc(nosst.lag3_ran$mer) #
AICc_3lagran_nosst - AICc_3lag

nolatlong.lag3_ran <- gamm4(log_sc_weight ~   s(sst.amj, k=4) + s(julian, k = 4)+
                              prevage_lastyr_weight_anom,
                            random=~(1|YEAR/HAUL), data=lag3complete) 
R2.l3.nolatlong <- 1-var(residuals(nolatlong.lag3_ran$gam))/(var(model.response(model.frame(nolatlong.lag3_ran$gam))))
R2.l3 - R2.l3.nolatlong 
AICc_3lagran_nolatlong <- AICc(nolatlong.lag3_ran$mer) #
AICc_3lagran_nolatlong - AICc_3lag

nojul.lag3_ran <- gamm4(log_sc_weight ~   s(sst.amj, k=4) + t2(LONGITUDE, LATITUDE)+
                          prevage_lastyr_weight_anom,
                        random=~(1|YEAR/HAUL), data=lag3complete) 
R2.l3.nojul <- 1-var(residuals(nojul.lag3_ran$gam))/(var(model.response(model.frame(nojul.lag3_ran$gam))))
R2.l3 - R2.l3.nojul 
AICc_3lagran_nojul <- AICc(nojul.lag3_ran$mer) #
AICc_3lagran_nojul - AICc_3lag


nolag.lag3_ran <- gamm4(log_sc_weight ~   s(sst.amj, k=4) + t2(LONGITUDE, LATITUDE),
                        random=~(1|YEAR/HAUL), data=lag3complete) 
R2.l3.nolag <- 1-var(residuals(nolag.lag3_ran$gam))/(var(model.response(model.frame(nolag.lag3_ran$gam))))
R2.l3 - R2.l3.nolag 
AICc_3lagran_nolag <- AICc(nolag.lag3_ran$mer) #
AICc_3lagran_nolag - AICc_3lag





# age 4====


lag.null4_ran <- gamm4(log_sc_weight ~  s(sst.amj, k=4) + t2(LONGITUDE, LATITUDE) + s(julian, k = 4) +
                         s(prevage_lastyr_weight_anom, k=4),
                       random=~(1|YEAR/HAUL), data=lag4complete)
summary(lag.null4_ran$gam) #lag not sig
summary(lag.null4_ran$mer)
AIC(lag.null4_ran$mer) #
#R-sq.(adj) = 0.257

laglin.null4_ran <- gamm4(log_sc_weight ~  s(sst.amj, k=4) + t2(LONGITUDE, LATITUDE) + s(julian, k = 4) +
                         prevage_lastyr_weight_anom,
                       random=~(1|YEAR/HAUL), data=lag4complete)
summary(laglin.null4_ran$gam) #lag not sig
summary(laglin.null4_ran$mer)
AIC(laglin.null4_ran$mer) #
#R-sq.(adj) = 0.257

AICc_4lag <- AICc(laglin.null4_ran$mer) #
AICc_4lag #5549.182

AICc(base.null4_ran$mer) #

AICc_4base_ran <- AICc(base.null4_ran$mer) #
AICc_4base_ran #5537.379

AICc_4base_ran - AICc_4lag
#-11.80
draw(lag.null4_ran$gam, select=4)
draw(lag.null4_ran$gam, select=1)

lag.null4_ran$gam$data <- lag4complete
lr4 <- visreg(lag.null4_ran$gam, "sst.amj", scale="response",ylab="Scaled log(weight-at-age)", xlab="April-June SST", rug=1)


#drop terms age 4----
AICc_4lag 

R2.l4 <- 1-var(residuals(laglin.null4_ran$gam))/(var(model.response(model.frame(laglin.null4_ran$gam))))

nosst.lag4_ran <- gamm4(log_sc_weight ~   t2(LONGITUDE, LATITUDE) + s(julian, k = 4)+
                          prevage_lastyr_weight_anom,
                        random=~(1|YEAR/HAUL), data=lag4complete) 
R2.l4.nosst <- 1-var(residuals(nosst.lag4_ran$gam))/(var(model.response(model.frame(nosst.lag4_ran$gam))))
R2.l4 - R2.l4.nosst 
AICc_4lagran_nosst <- AICc(nosst.lag4_ran$mer) #
AICc_4lagran_nosst - AICc_4lag

nolatlong.lag4_ran <- gamm4(log_sc_weight ~   s(sst.amj, k=4) + s(julian, k = 4)+
                              prevage_lastyr_weight_anom,
                            random=~(1|YEAR/HAUL), data=lag4complete) 
R2.l4.nolatlong <- 1-var(residuals(nolatlong.lag4_ran$gam))/(var(model.response(model.frame(nolatlong.lag4_ran$gam))))
R2.l4 - R2.l4.nolatlong 
AICc_4lagran_nolatlong <- AICc(nolatlong.lag4_ran$mer) #
AICc_4lagran_nolatlong - AICc_4lag

nojul.lag4_ran <- gamm4(log_sc_weight ~   s(sst.amj, k=4) + t2(LONGITUDE, LATITUDE)+
                          prevage_lastyr_weight_anom,
                        random=~(1|YEAR/HAUL), data=lag4complete) 
R2.l4.nojul <- 1-var(residuals(nojul.lag4_ran$gam))/(var(model.response(model.frame(nojul.lag4_ran$gam))))
R2.l4 - R2.l4.nojul 
AICc_4lagran_nojul <- AICc(nojul.lag4_ran$mer) #
AICc_4lagran_nojul - AICc_4lag


nolag.lag4_ran <- gamm4(log_sc_weight ~   s(sst.amj, k=4) + t2(LONGITUDE, LATITUDE),
                        random=~(1|YEAR/HAUL), data=lag4complete) 
R2.l4.nolag <- 1-var(residuals(nolag.lag4_ran$gam))/(var(model.response(model.frame(nolag.lag4_ran$gam))))
R2.l4 - R2.l4.nolag 
AICc_4lagran_nolag <- AICc(nolag.lag4_ran$mer) #
AICc_4lagran_nolag - AICc_4lag





# age 5====


lag.null5_ran <- gamm4(log_sc_weight ~  s(sst.amj, k=4) + t2(LONGITUDE, LATITUDE) + s(julian, k = 4) +
                         s(prevage_lastyr_weight_anom, k=4),
                       random=~(1|YEAR/HAUL), data=lag5complete)
summary(lag.null5_ran$gam) #lag not sig
summary(lag.null5_ran$mer)
AIC(lag.null5_ran$mer) #8095.1
#R-sq.(adj) = 0.176

laglin.null5_ran <- gamm4(log_sc_weight ~  s(sst.amj, k=4) + t2(LONGITUDE, LATITUDE) + s(julian, k = 4) +
                         prevage_lastyr_weight_anom,
                       random=~(1|YEAR/HAUL), data=lag5complete)
summary(laglin.null5_ran$gam) #lag not sig
summary(laglin.null5_ran$mer)
AIC(laglin.null5_ran$mer) #8095.1
#R-sq.(adj) = 0.176

AICc_5lag <- AICc(laglin.null5_ran$mer) #8101.141
AICc_5lag 

AICc(base.null5_ran$mer) #8091.201

AICc_5base_ran <- AICc(base.null5_ran$mer) #
AICc_5base_ran #

AICc_5base_ran - AICc_5lag
#-9.9395
draw(lag.null5_ran$gam, select=4)
draw(lag.null5_ran$gam, select=1)

lag.null5_ran$gam$data <- lag5complete
lr5 <- visreg(lag.null5_ran$gam, "sst.amj", scale="response",ylab="Scaled log(weight-at-age)", xlab="April-June SST", rug=1)


#drop terms age 5----
AICc_5lag 

R2.l5 <- 1-var(residuals(laglin.null5_ran$gam))/(var(model.response(model.frame(laglin.null5_ran$gam))))

nosst.lag5_ran <- gamm4(log_sc_weight ~   t2(LONGITUDE, LATITUDE) + s(julian, k = 4)+
                          prevage_lastyr_weight_anom,
                        random=~(1|YEAR/HAUL), data=lag5complete) 
R2.l5.nosst <- 1-var(residuals(nosst.lag5_ran$gam))/(var(model.response(model.frame(nosst.lag5_ran$gam))))
R2.l5 - R2.l5.nosst 
AICc_5lagran_nosst <- AICc(nosst.lag5_ran$mer) #
AICc_5lagran_nosst - AICc_5lag

nolatlong.lag5_ran <- gamm4(log_sc_weight ~   s(sst.amj, k=4) + s(julian, k = 4)+
                              prevage_lastyr_weight_anom,
                            random=~(1|YEAR/HAUL), data=lag5complete) 
R2.l5.nolatlong <- 1-var(residuals(nolatlong.lag5_ran$gam))/(var(model.response(model.frame(nolatlong.lag5_ran$gam))))
R2.l5 - R2.l5.nolatlong 
AICc_5lagran_nolatlong <- AICc(nolatlong.lag5_ran$mer) #
AICc_5lagran_nolatlong - AICc_5lag

nojul.lag5_ran <- gamm4(log_sc_weight ~   s(sst.amj, k=4) + t2(LONGITUDE, LATITUDE)+
                          prevage_lastyr_weight_anom,
                        random=~(1|YEAR/HAUL), data=lag5complete) 
R2.l5.nojul <- 1-var(residuals(nojul.lag5_ran$gam))/(var(model.response(model.frame(nojul.lag5_ran$gam))))
R2.l5 - R2.l5.nojul 
AICc_5lagran_nojul <- AICc(nojul.lag5_ran$mer) #
AICc_5lagran_nojul - AICc_5lag


nolag.lag5_ran <- gamm4(log_sc_weight ~   s(sst.amj, k=4) + t2(LONGITUDE, LATITUDE),
                        random=~(1|YEAR/HAUL), data=lag5complete) 
R2.l5.nolag <- 1-var(residuals(nolag.lag5_ran$gam))/(var(model.response(model.frame(nolag.lag5_ran$gam))))
R2.l5 - R2.l5.nolag 
AICc_5lagran_nolag <- AICc(nolag.lag5_ran$mer) #
AICc_5lagran_nolag - AICc_5lag





# age 6====


lag.null6_ran <- gamm4(log_sc_weight ~  s(sst.amj, k=4) + t2(LONGITUDE, LATITUDE) + s(julian, k = 4) +
                         s(prevage_lastyr_weight_anom, k=4),
                       random=~(1|YEAR/HAUL), data=lag6complete)
summary(lag.null6_ran$gam) #lag IS significant 0.0003
summary(lag.null6_ran$mer)
AIC(lag.null6_ran$mer) #9817.93
#R-sq.(adj) = 0.186

laglin.null6_ran <- gamm4(log_sc_weight ~  s(sst.amj, k=4) + t2(LONGITUDE, LATITUDE) + s(julian, k = 4) +
                         prevage_lastyr_weight_anom,
                       random=~(1|YEAR/HAUL), data=lag6complete)
summary(laglin.null6_ran$gam) #lag IS significant 0.0003
summary(laglin.null6_ran$mer)
AIC(laglin.null6_ran$mer) #9817.93
#R-sq.(adj) = 0.186

AICc_6lag <- AICc(laglin.null6_ran$mer) #9823.533
AICc_6lag 

AICc(base.null6_ran$mer) #

AICc_6base_ran <- AICc(base.null6_ran$mer) #
AICc_6base_ran #9820.437

AICc_6base_ran - AICc_6lag
#-3.095
draw(lag.null6_ran$gam, select=4)
draw(lag.null6_ran$gam, select=1)

lag.null6_ran$gam$data <- lag6complete
lr6 <- visreg(lag.null6_ran$gam, "sst.amj", scale="response",ylab="Scaled log(weight-at-age)", xlab="April-June SST", rug=1)


#drop terms age 6----
AICc_6lag 

R2.l6 <- 1-var(residuals(laglin.null6_ran$gam))/(var(model.response(model.frame(laglin.null6_ran$gam))))

nosst.lag6_ran <- gamm4(log_sc_weight ~   t2(LONGITUDE, LATITUDE) + s(julian, k = 4)+
                          prevage_lastyr_weight_anom,
                        random=~(1|YEAR/HAUL), data=lag6complete) 
R2.l6.nosst <- 1-var(residuals(nosst.lag6_ran$gam))/(var(model.response(model.frame(nosst.lag6_ran$gam))))
R2.l6 - R2.l6.nosst 
AICc_6lagran_nosst <- AICc(nosst.lag6_ran$mer) #
AICc_6lagran_nosst - AICc_6lag

nolatlong.lag6_ran <- gamm4(log_sc_weight ~   s(sst.amj, k=4) + s(julian, k = 4)+
                              prevage_lastyr_weight_anom,
                            random=~(1|YEAR/HAUL), data=lag6complete) 
R2.l6.nolatlong <- 1-var(residuals(nolatlong.lag6_ran$gam))/(var(model.response(model.frame(nolatlong.lag6_ran$gam))))
R2.l6 - R2.l6.nolatlong 
AICc_6lagran_nolatlong <- AICc(nolatlong.lag6_ran$mer) #
AICc_6lagran_nolatlong - AICc_6lag

nojul.lag6_ran <- gamm4(log_sc_weight ~   s(sst.amj, k=4) + t2(LONGITUDE, LATITUDE)+
                          prevage_lastyr_weight_anom,
                        random=~(1|YEAR/HAUL), data=lag6complete) 
R2.l6.nojul <- 1-var(residuals(nojul.lag6_ran$gam))/(var(model.response(model.frame(nojul.lag6_ran$gam))))
R2.l6 - R2.l6.nojul 
AICc_6lagran_nojul <- AICc(nojul.lag6_ran$mer) #
AICc_6lagran_nojul - AICc_6lag


nolag.lag6_ran <- gamm4(log_sc_weight ~   s(sst.amj, k=4) + t2(LONGITUDE, LATITUDE),
                        random=~(1|YEAR/HAUL), data=lag6complete) 
R2.l6.nolag <- 1-var(residuals(nolag.lag6_ran$gam))/(var(model.response(model.frame(nolag.lag6_ran$gam))))
R2.l6 - R2.l6.nolag 
AICc_6lagran_nolag <- AICc(nolag.lag6_ran$mer) #
AICc_6lagran_nolag - AICc_6lag





# age 7====


lag.null7_ran <- gamm4(log_sc_weight ~  s(sst.amj, k=4) + t2(LONGITUDE, LATITUDE) + s(julian, k = 4) +
                         s(prevage_lastyr_weight_anom, k=4),
                       random=~(1|YEAR/HAUL), data=lag7complete)
summary(lag.null7_ran$gam) #lag sig 0.031
summary(lag.null7_ran$mer)
AIC(lag.null7_ran$mer) #
#R-sq.(adj) = 0.231

laglin.null7_ran <- gamm4(log_sc_weight ~  s(sst.amj, k=4) + t2(LONGITUDE, LATITUDE) + s(julian, k = 4) +
                         prevage_lastyr_weight_anom,
                       random=~(1|YEAR/HAUL), data=lag7complete)
summary(laglin.null7_ran$gam) #lag sig 0.031
summary(laglin.null7_ran$mer)
AIC(laglin.null7_ran$mer) #
#R-sq.(adj) = 0.231

AICc_7lag <- AICc(laglin.null7_ran$mer) #
AICc_7lag #8414.382

AICc(base.null7_ran$mer) #

AICc_7base_ran <- AICc(base.null7_ran$mer) #
AICc_7base_ran #8404.605

AICc_7base_ran - AICc_7lag
#-9.7746
draw(lag.null7_ran$gam, select=4)
draw(lag.null7_ran$gam, select=1)

lag.null7_ran$gam$data <- lag7complete
lr7 <- visreg(lag.null7_ran$gam, "sst.amj", scale="response",ylab="Scaled log(weight-at-age)", xlab="April-June SST", rug=1)


#drop terms age 7----
AICc_7lag 

R2.l7 <- 1-var(residuals(laglin.null7_ran$gam))/(var(model.response(model.frame(laglin.null7_ran$gam))))

nosst.lag7_ran <- gamm4(log_sc_weight ~   t2(LONGITUDE, LATITUDE) + s(julian, k = 4)+
                          prevage_lastyr_weight_anom,
                        random=~(1|YEAR/HAUL), data=lag7complete) 
R2.l7.nosst <- 1-var(residuals(nosst.lag7_ran$gam))/(var(model.response(model.frame(nosst.lag7_ran$gam))))
R2.l7 - R2.l7.nosst 
AICc_7lagran_nosst <- AICc(nosst.lag7_ran$mer) #
AICc_7lagran_nosst - AICc_7lag

nolatlong.lag7_ran <- gamm4(log_sc_weight ~   s(sst.amj, k=4) + s(julian, k = 4)+
                              prevage_lastyr_weight_anom,
                            random=~(1|YEAR/HAUL), data=lag7complete) 
R2.l7.nolatlong <- 1-var(residuals(nolatlong.lag7_ran$gam))/(var(model.response(model.frame(nolatlong.lag7_ran$gam))))
R2.l7 - R2.l7.nolatlong 
AICc_7lagran_nolatlong <- AICc(nolatlong.lag7_ran$mer) #
AICc_7lagran_nolatlong - AICc_7lag

nojul.lag7_ran <- gamm4(log_sc_weight ~   s(sst.amj, k=4) + t2(LONGITUDE, LATITUDE)+
                          prevage_lastyr_weight_anom,
                        random=~(1|YEAR/HAUL), data=lag7complete) 
R2.l7.nojul <- 1-var(residuals(nojul.lag7_ran$gam))/(var(model.response(model.frame(nojul.lag7_ran$gam))))
R2.l7 - R2.l7.nojul 
AICc_7lagran_nojul <- AICc(nojul.lag7_ran$mer) #
AICc_7lagran_nojul - AICc_7lag


nolag.lag7_ran <- gamm4(log_sc_weight ~   s(sst.amj, k=4) + t2(LONGITUDE, LATITUDE),
                        random=~(1|YEAR/HAUL), data=lag7complete) 
R2.l7.nolag <- 1-var(residuals(nolag.lag7_ran$gam))/(var(model.response(model.frame(nolag.lag7_ran$gam))))
R2.l7 - R2.l7.nolag 
AICc_7lagran_nolag <- AICc(nolag.lag7_ran$mer) #
AICc_7lagran_nolag - AICc_7lag





# age 8====


lag.null8_ran <- gamm4(log_sc_weight ~  s(sst.amj, k=4) + t2(LONGITUDE, LATITUDE) + s(julian, k = 4) +
                         s(prevage_lastyr_weight_anom, k=4),
                       random=~(1|YEAR/HAUL), data=lag8complete)
summary(lag.null8_ran$gam) #lag is significant 1.83e-05
summary(lag.null8_ran$mer)
AIC(lag.null8_ran$mer) #
#R-sq.(adj) = 0.266

laglin.null8_ran <- gamm4(log_sc_weight ~  s(sst.amj, k=4) + t2(LONGITUDE, LATITUDE) + s(julian, k = 4) +
                         prevage_lastyr_weight_anom,
                       random=~(1|YEAR/HAUL), data=lag8complete)
summary(laglin.null8_ran$gam) #lag is significant 1.76e-05
summary(laglin.null8_ran$mer)
AIC(laglin.null8_ran$mer) #
#R-sq.(adj) = 0.266

AICc_8lag <- AICc(laglin.null8_ran$mer) #6116.906
AICc_8lag 

AICc(base.null8_ran$mer) #

AICc_8base_ran <- AICc(base.null8_ran$mer) #6114.534
AICc_8base_ran #

AICc_8base_ran - AICc_8lag
#-2.37
draw(lag.null8_ran$gam, select=4)
draw(lag.null8_ran$gam, select=1)

lag.null8_ran$gam$data <- lag8complete
lr8 <- visreg(lag.null8_ran$gam, "sst.amj", scale="response",ylab="Scaled log(weight-at-age)", xlab="April-June SST", rug=1)


tt8 <- getViz(laglin.null8_ran$gam)
plot(sm(tt8, 2)) + l_fitRaster() + l_fitContour() + labs(title = NULL) + #l_points() +
  geom_polygon(data = map_data ("world"), 
               aes(x=long, y = lat,group=group),fill="white",color="black",
               inherit.aes = F)+coord_sf(xlim = c(-177, -158.5), ylim = c(54.5, 62), expand = FALSE) +
  scale_fill_distiller(palette = "Spectral", type = "div", limits = c(-3,3)) +
  theme(legend.position = "none")+ theme(plot.margin = unit(c(0, 0, 0, 0.1), "cm"), plot.title = element_blank(),
                                         axis.title.x = element_blank(),
                                         axis.title.y = element_blank())



#drop terms age 8----
AICc_8lag 

R2.l8 <- 1-var(residuals(laglin.null8_ran$gam))/(var(model.response(model.frame(laglin.null8_ran$gam))))

nosst.lag8_ran <- gamm4(log_sc_weight ~   t2(LONGITUDE, LATITUDE) + s(julian, k = 4)+
                          prevage_lastyr_weight_anom,
                        random=~(1|YEAR/HAUL), data=lag8complete) 
R2.l8.nosst <- 1-var(residuals(nosst.lag8_ran$gam))/(var(model.response(model.frame(nosst.lag8_ran$gam))))
R2.l8 - R2.l8.nosst 
AICc_8lagran_nosst <- AICc(nosst.lag8_ran$mer) #
AICc_8lagran_nosst - AICc_8lag

nolatlong.lag8_ran <- gamm4(log_sc_weight ~   s(sst.amj, k=4) + s(julian, k = 4)+
                              prevage_lastyr_weight_anom,
                            random=~(1|YEAR/HAUL), data=lag8complete) 
R2.l8.nolatlong <- 1-var(residuals(nolatlong.lag8_ran$gam))/(var(model.response(model.frame(nolatlong.lag8_ran$gam))))
R2.l8 - R2.l8.nolatlong 
AICc_8lagran_nolatlong <- AICc(nolatlong.lag8_ran$mer) #
AICc_8lagran_nolatlong - AICc_8lag

nojul.lag8_ran <- gamm4(log_sc_weight ~   s(sst.amj, k=4) + t2(LONGITUDE, LATITUDE)+
                          prevage_lastyr_weight_anom,
                        random=~(1|YEAR/HAUL), data=lag8complete) 
R2.l8.nojul <- 1-var(residuals(nojul.lag8_ran$gam))/(var(model.response(model.frame(nojul.lag8_ran$gam))))
R2.l8 - R2.l8.nojul 
AICc_8lagran_nojul <- AICc(nojul.lag8_ran$mer) #
AICc_8lagran_nojul - AICc_8lag


nolag.lag8_ran <- gamm4(log_sc_weight ~   s(sst.amj, k=4) + t2(LONGITUDE, LATITUDE),
                        random=~(1|YEAR/HAUL), data=lag8complete) 
R2.l8.nolag <- 1-var(residuals(nolag.lag8_ran$gam))/(var(model.response(model.frame(nolag.lag8_ran$gam))))
R2.l8 - R2.l8.nolag 
AICc_8lagran_nolag <- AICc(nolag.lag8_ran$mer) #
AICc_8lagran_nolag - AICc_8lag





# age 9====


lag.null9_ran <- gamm4(log_sc_weight ~  s(sst.amj, k=4) + t2(LONGITUDE, LATITUDE) + s(julian, k = 4) +
                         s(prevage_lastyr_weight_anom, k=4),
                       random=~(1|YEAR/HAUL), data=lag9complete)
summary(lag.null9_ran$gam) #lag significant 1.11e-07
summary(lag.null9_ran$mer)
AIC(lag.null9_ran$mer) #
#R-sq.(adj) = 0.251

laglin.null9_ran <- gamm4(log_sc_weight ~  s(sst.amj, k=4) + t2(LONGITUDE, LATITUDE) + s(julian, k = 4) +
                         prevage_lastyr_weight_anom,
                       random=~(1|YEAR/HAUL), data=lag9complete)
summary(laglin.null9_ran$gam) #lag significant 1.14e-07
summary(laglin.null9_ran$mer)
AIC(laglin.null9_ran$mer) #
#R-sq.(adj) = 0.251

AICc_9lag <- AICc(laglin.null9_ran$mer) #4920.083
AICc_9lag 

AICc(base.null9_ran$mer) #

AICc_9base_ran <- AICc(base.null9_ran$mer) #4921.382
AICc_9base_ran #

AICc_9base_ran - AICc_9lag
#1.3
draw(lag.null9_ran$gam, select=4)
draw(lag.null9_ran$gam, select=1)

lag.null9_ran$gam$data <- lag9complete
lr9 <- visreg(lag.null9_ran$gam, "sst.amj", scale="response",ylab="Scaled log(weight-at-age)", xlab="April-June SST", rug=1)


tt9 <- getViz(laglin.null9_ran$gam)
plot(sm(tt9, 2)) + l_fitRaster() + l_fitContour() + labs(title = NULL) + #l_points() +
  geom_polygon(data = map_data ("world"), 
               aes(x=long, y = lat,group=group),fill="white",color="black",
               inherit.aes = F)+coord_sf(xlim = c(-177, -158.5), ylim = c(54.5, 62), expand = FALSE) +
  scale_fill_distiller(palette = "Spectral", type = "div", limits = c(-3,3)) +
  theme(legend.position = "none")+ theme(plot.margin = unit(c(0, 0, 0, 0.1), "cm"), plot.title = element_blank(),
                                         axis.title.x = element_blank(),
                                         axis.title.y = element_blank())



#drop terms age 9----
AICc_9lag 

R2.l9 <- 1-var(residuals(laglin.null9_ran$gam))/(var(model.response(model.frame(laglin.null9_ran$gam))))

nosst.lag9_ran <- gamm4(log_sc_weight ~   t2(LONGITUDE, LATITUDE) + s(julian, k = 4)+
                           prevage_lastyr_weight_anom,
                         random=~(1|YEAR/HAUL), data=lag9complete) 
R2.l9.nosst <- 1-var(residuals(nosst.lag9_ran$gam))/(var(model.response(model.frame(nosst.lag9_ran$gam))))
R2.l9 - R2.l9.nosst 
AICc_9lagran_nosst <- AICc(nosst.lag9_ran$mer) #
AICc_9lagran_nosst - AICc_9lag

nolatlong.lag9_ran <- gamm4(log_sc_weight ~   s(sst.amj, k=4) + s(julian, k = 4)+
                               prevage_lastyr_weight_anom,
                             random=~(1|YEAR/HAUL), data=lag9complete) 
R2.l9.nolatlong <- 1-var(residuals(nolatlong.lag9_ran$gam))/(var(model.response(model.frame(nolatlong.lag9_ran$gam))))
R2.l9 - R2.l9.nolatlong 
AICc_9lagran_nolatlong <- AICc(nolatlong.lag9_ran$mer) #
AICc_9lagran_nolatlong - AICc_9lag

nojul.lag9_ran <- gamm4(log_sc_weight ~   s(sst.amj, k=4) + t2(LONGITUDE, LATITUDE)+
                           prevage_lastyr_weight_anom,
                         random=~(1|YEAR/HAUL), data=lag9complete) 
R2.l9.nojul <- 1-var(residuals(nojul.lag9_ran$gam))/(var(model.response(model.frame(nojul.lag9_ran$gam))))
R2.l9 - R2.l9.nojul 
AICc_9lagran_nojul <- AICc(nojul.lag9_ran$mer) #
AICc_9lagran_nojul - AICc_9lag


nolag.lag9_ran <- gamm4(log_sc_weight ~   s(sst.amj, k=4) + t2(LONGITUDE, LATITUDE),
                         random=~(1|YEAR/HAUL), data=lag9complete) 
R2.l9.nolag <- 1-var(residuals(nolag.lag9_ran$gam))/(var(model.response(model.frame(nolag.lag9_ran$gam))))
R2.l9 - R2.l9.nolag 
AICc_9lagran_nolag <- AICc(nolag.lag9_ran$mer) #
AICc_9lagran_nolag - AICc_9lag



# age 10====


lag.null10_ran <- gamm4(log_sc_weight ~  s(sst.amj, k=4) + t2(LONGITUDE, LATITUDE) + s(julian, k = 4) +
                         s(prevage_lastyr_weight_anom, k=4),
                       random=~(1|YEAR/HAUL), data=lag10complete)
summary(lag.null10_ran$gam) #lag is sig
summary(lag.null10_ran$mer)
AIC(lag.null10_ran$mer) #
#R-sq.(adj) = 0.213

laglin.null10_ran <- gamm4(log_sc_weight ~  s(sst.amj, k=4) + t2(LONGITUDE, LATITUDE) + s(julian, k = 4) +
                          prevage_lastyr_weight_anom,
                        random=~(1|YEAR/HAUL), data=lag10complete)
summary(laglin.null10_ran$gam) #lag is sig
summary(laglin.null10_ran$mer)
AIC(laglin.null10_ran$mer) #4165.847
#R-sq.(adj) =0.213

AICc_10lag <- AICc(laglin.null10_ran$mer) #4166.148
AICc_10lag 

AICc(base.null10_ran$mer) #

AICc_10base_ran <- AICc(base.null10_ran$mer) #
AICc_10base_ran #4161.175

AICc_10base_ran - AICc_10lag
#-4.97
draw(lag.null10_ran$gam, select=4)
draw(lag.null10_ran$gam, select=1)

lag.null10_ran$gam$data <- lag10complete
lr10 <- visreg(lag.null10_ran$gam, "sst.amj", scale="response",ylab="Scaled log(weight-at-age)", xlab="April-June SST", rug=1)


tt10 <- getViz(laglin.null10_ran$gam)
plot(sm(tt10, 2)) + l_fitRaster() + l_fitContour() + labs(title = NULL) + #l_points() +
  geom_polygon(data = map_data ("world"), 
               aes(x=long, y = lat,group=group),fill="white",color="black",
               inherit.aes = F)+coord_sf(xlim = c(-177, -158.5), ylim = c(54.5, 62), expand = FALSE) +
  scale_fill_distiller(palette = "Spectral", type = "div", limits = c(-3,3)) +
  theme(legend.position = "none")+ theme(plot.margin = unit(c(0, 0, 0, 0.1), "cm"), plot.title = element_blank(),
                                         axis.title.x = element_blank(),
                                         axis.title.y = element_blank())



#drop terms age 10----
AICc_10lag 

R2.l10 <- 1-var(residuals(laglin.null10_ran$gam))/(var(model.response(model.frame(laglin.null10_ran$gam))))

nosst.lag10_ran <- gamm4(log_sc_weight ~   t2(LONGITUDE, LATITUDE) + s(julian, k = 4)+
                           prevage_lastyr_weight_anom,
                         random=~(1|YEAR/HAUL), data=lag10complete) 
R2.l10.nosst <- 1-var(residuals(nosst.lag10_ran$gam))/(var(model.response(model.frame(nosst.lag10_ran$gam))))
R2.l10 - R2.l10.nosst 
AICc_10lagran_nosst <- AICc(nosst.lag10_ran$mer) #
AICc_10lagran_nosst - AICc_10lag

nolatlong.lag10_ran <- gamm4(log_sc_weight ~   s(sst.amj, k=4) + s(julian, k = 4)+
                               prevage_lastyr_weight_anom,
                             random=~(1|YEAR/HAUL), data=lag10complete) 
R2.l10.nolatlong <- 1-var(residuals(nolatlong.lag10_ran$gam))/(var(model.response(model.frame(nolatlong.lag10_ran$gam))))
R2.l10 - R2.l10.nolatlong 
AICc_10lagran_nolatlong <- AICc(nolatlong.lag10_ran$mer) #
AICc_10lagran_nolatlong - AICc_10lag

nojul.lag10_ran <- gamm4(log_sc_weight ~   s(sst.amj, k=4) + t2(LONGITUDE, LATITUDE)+
                           prevage_lastyr_weight_anom,
                         random=~(1|YEAR/HAUL), data=lag10complete) 
R2.l10.nojul <- 1-var(residuals(nojul.lag10_ran$gam))/(var(model.response(model.frame(nojul.lag10_ran$gam))))
R2.l10 - R2.l10.nojul 
AICc_10lagran_nojul <- AICc(nojul.lag10_ran$mer) #
AICc_10lagran_nojul - AICc_10lag



nolag.lag10_ran <- gamm4(log_sc_weight ~   s(sst.amj, k=4) + t2(LONGITUDE, LATITUDE),
                         random=~(1|YEAR/HAUL), data=lag10complete) 
R2.l10.nolag <- 1-var(residuals(nolag.lag10_ran$gam))/(var(model.response(model.frame(nolag.lag10_ran$gam))))
R2.l10 - R2.l10.nolag 
AICc_10lagran_nolag <- AICc(nolag.lag10_ran$mer) #
AICc_10lagran_nolag - AICc_10lag



# age 11====


lag.null11_ran <- gamm4(log_sc_weight ~  s(sst.amj, k=4) + t2(LONGITUDE, LATITUDE) + s(julian, k = 4) +
                          s(prevage_lastyr_weight_anom, k=4),
                        random=~(1|YEAR/HAUL), data=lag11complete)
summary(lag.null11_ran$gam) #lag is sig and significantly nonlinear
summary(lag.null11_ran$mer)
AIC(lag.null11_ran$mer) #
#R-sq.(adj) = 0.233

laglin.null11_ran <- gamm4(log_sc_weight ~  s(sst.amj, k=4) + t2(LONGITUDE, LATITUDE) + s(julian, k = 4) +
                             prevage_lastyr_weight_anom,
                           random=~(1|YEAR/HAUL), data=lag11complete)
summary(laglin.null11_ran$gam) #
summary(laglin.null11_ran$mer)
AIC(laglin.null11_ran$mer) #
#R-sq.(adj) =

AICc_11lag <- AICc(laglin.null11_ran$mer) #3415.326
AICc_11lag 

AICc(base.null11_ran$mer) #

AICc_11base_ran <- AICc(base.null11_ran$mer) #3426.908
AICc_11base_ran #


AICc_11base_ran - AICc_11lag
#11.58
draw(lag.null11_ran$gam, select=4)
draw(lag.null11_ran$gam, select=1)

lag.null11_ran$gam$data <- lag11complete
lr11 <- visreg(lag.null11_ran$gam, "sst.amj", scale="response",ylab="Scaled log(weight-at-age)", xlab="April-June SST", rug=1)


tt11 <- getViz(laglin.null11_ran$gam)
plot(sm(tt11, 2)) + l_fitRaster() + l_fitContour() + labs(title = NULL) + #l_points() +
  geom_polygon(data = map_data ("world"), 
               aes(x=long, y = lat,group=group),fill="white",color="black",
               inherit.aes = F)+coord_sf(xlim = c(-177, -158.5), ylim = c(54.5, 62), expand = FALSE) +
  scale_fill_distiller(palette = "Spectral", type = "div", limits = c(-3,3)) +
  theme(legend.position = "none")+ theme(plot.margin = unit(c(0, 0, 0, 0.1), "cm"), plot.title = element_blank(),
                                         axis.title.x = element_blank(),
                                         axis.title.y = element_blank())



#drop terms age 11----
AICc_11lag 

R2.l11 <- 1-var(residuals(laglin.null11_ran$gam))/(var(model.response(model.frame(laglin.null11_ran$gam))))

nosst.lag11_ran <- gamm4(log_sc_weight ~   t2(LONGITUDE, LATITUDE) + s(julian, k = 4)+
                            prevage_lastyr_weight_anom,
                          random=~(1|YEAR/HAUL), data=lag11complete) 
R2.l11.nosst <- 1-var(residuals(nosst.lag11_ran$gam))/(var(model.response(model.frame(nosst.lag11_ran$gam))))
R2.l11 - R2.l11.nosst 
AICc_11lagran_nosst <- AICc(nosst.lag11_ran$mer) #
AICc_11lagran_nosst - AICc_11lag

nolatlong.lag11_ran <- gamm4(log_sc_weight ~   s(sst.amj, k=4) + s(julian, k = 4)+
                               prevage_lastyr_weight_anom,
                              random=~(1|YEAR/HAUL), data=lag11complete) 
R2.l11.nolatlong <- 1-var(residuals(nolatlong.lag11_ran$gam))/(var(model.response(model.frame(nolatlong.lag11_ran$gam))))
R2.l11 - R2.l11.nolatlong 
AICc_11lagran_nolatlong <- AICc(nolatlong.lag11_ran$mer) #
AICc_11lagran_nolatlong - AICc_11lag

nojul.lag11_ran <- gamm4(log_sc_weight ~   s(sst.amj, k=4) + t2(LONGITUDE, LATITUDE)+
                            prevage_lastyr_weight_anom,
                          random=~(1|YEAR/HAUL), data=lag11complete) 
R2.l11.nojul <- 1-var(residuals(nojul.lag11_ran$gam))/(var(model.response(model.frame(nojul.lag11_ran$gam))))
R2.l11 - R2.l11.nojul 
AICc_11lagran_nojul <- AICc(nojul.lag11_ran$mer) #
AICc_11lagran_nojul - AICc_11lag


nolag.lag11_ran <- gamm4(log_sc_weight ~   s(sst.amj, k=4) + t2(LONGITUDE, LATITUDE),
                         random=~(1|YEAR/HAUL), data=lag11complete) 
R2.l11.nolag <- 1-var(residuals(nolag.lag11_ran$gam))/(var(model.response(model.frame(nolag.lag11_ran$gam))))
R2.l11 - R2.l11.nolag 
AICc_11lagran_nolag <- AICc(nolag.lag11_ran$mer) #
AICc_11lagran_nolag - AICc_11lag


#plot-----


plr3 <- visregList(lr8, lr9, lr10, lr11, 
                   lr4, lr5, lr6,  lr7, 
                   lr2, lr3, 
                   collapse=TRUE,
                   labels=c("Age 8", "Age 9",  "Age 10", "Age 11", 
                            "Age 4", "Age 5", "Age 6", "Age 7",  
                            "Age 2", "Age 3"))

plot(plr3,
     ylab="Scaled log (weight-at-age)",
     xlab="April-June SST", rug=1) #for manuscript but with random






