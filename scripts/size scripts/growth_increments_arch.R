#============================================================================================================
#growth increments

#Krista
#Oct-13-2021
#============================================================================================================
#Notes: Based on suggestion from Jim I during internal review
#not currently in manuscript, just bonus!
#============================================================================================================

library(lme4)
library(mgcv)
library(gamm4)
library(cowplot)

#===

#loaded in weight_age_exploration.R
names(lagdat)

w2means <- lagdat %>% group_by(YEAR, AGE) %>% summarize(mean_annual_weight_at_age=mean(WEIGHT, na.rm=TRUE))

#loop?

w2means$prev_weight_at_age <- NA
w2means$prev_weight_at_age <- as.numeric(w2means$prev_weight_at_age)

table(w2means$AGE, w2means$YEAR)
#data in all yrs up to age 14

w2means <- w2means[which(w2means$AGE<14),]

counter <- 1
i<-1
for(i in 1:length(w2means$YEAR)){
  prev_WAA <- NA
  yr_i <- w2means$YEAR[i]
    age_i <- w2means$AGE[i]
    
    prev_yr <- yr_i - 1
    prev_age <- age_i - 1
    
   if(prev_age>0 & prev_yr>1998) { 
      prev_WAA <- w2means$mean_annual_weight_at_age[which(w2means$YEAR==prev_yr & w2means$AGE==prev_age)]
     }
  
   w2means$prev_weight_at_age[i] <- prev_WAA
     
    
  counter <- counter + 1
}

w2means$growth_increment <- w2means$mean_annual_weight_at_age - w2means$prev_weight_at_age

#all looks good, now match back to metadata

metadat <- lagdat[,c(1,4,12,17,20,21,22,23)] #select just cols w metadata, not individual data
metadat <- metadat[!duplicated(metadat),]

meanjoin <- left_join(w2means, metadat)

#take a look
ggplot(meanjoin, aes(YEAR, growth_increment, col=as.factor(AGE))) +
  geom_point() + geom_smooth() + facet_wrap(~AGE, scales="free")


ggplot(meanjoin, aes(sst.amj, growth_increment, col=as.factor(AGE))) +
  geom_point() + geom_smooth() + facet_wrap(~AGE, scales="free")

#try models==========

#cannot include the year|haul random effect

# age 2----
means2 <- meanjoin[which(meanjoin$AGE==2),]
incr.null2 <- gam(growth_increment ~ s(sst.amj, k=4), 
                  data=means2)
gam.check(incr.null2) #
summary(incr.null2) #smooth not sig
plot(incr.null2) #
draw(incr.null2, select=1)

incr.lm2 <- lm(growth_increment ~ sst.amj, 
                  data=means2)
summary(incr.lm2) #sst not significant

pg2 <- ggplot(means2, aes(sst.amj, growth_increment)) +
  geom_point() +  geom_abline(slope = coef(incr.lm2 )[[2]], intercept = coef(incr.lm2)[[1]])+
  xlab("April-June SST") + ylab("Growth increment")

# age 3----
means3 <- meanjoin[which(meanjoin$AGE==3),]
incr.null3 <- gam(growth_increment ~ s(sst.amj, k=4), 
                  data=means3)
gam.check(incr.null3) #
summary(incr.null3)#smooth not sig
plot(incr.null3) #
draw(incr.null3, select=1)

incr.lm3 <- lm(growth_increment ~ sst.amj, 
               data=means3)
summary(incr.lm3) 

pg3 <- ggplot(means3, aes(sst.amj, growth_increment)) +
  geom_point() +  geom_abline(slope = coef(incr.lm3 )[[2]], intercept = coef(incr.lm3)[[1]])+
  xlab("April-June SST") + ylab("Growth increment")


# age 4----
means4 <- meanjoin[which(meanjoin$AGE==4),]
incr.null4 <- gam(growth_increment ~ s(sst.amj, k=4), 
                  data=means4)
gam.check(incr.null4) #
summary(incr.null4)
plot(incr.null4) #smooth significant and wiggly
#pg4 <- draw(incr.null4, select=1)
pg4 <- visreg(incr.null4, "sst.amj", gg=TRUE)+
      xlab("April-June SST") + ylab("Growth increment")

incr.lm4 <- lm(growth_increment ~ sst.amj, 
               data=means4)
summary(incr.lm4) 

# pg4 <- ggplot(means4, aes(sst.amj, growth_increment)) +
#   geom_point() +  geom_abline(slope = coef(incr.lm4 )[[2]], intercept = coef(incr.lm4)[[1]])


# age 5----
means5 <- meanjoin[which(meanjoin$AGE==5),]
incr.null5 <- gam(growth_increment ~ s(sst.amj, k=4), 
                  data=means5)
gam.check(incr.null5) #
summary(incr.null5)#smooth not sig
plot(incr.null5) #
draw(incr.null5, select=1)

incr.lm5 <- lm(growth_increment ~ sst.amj, 
               data=means5)
summary(incr.lm5) #sst not sig

pg5 <- ggplot(means5, aes(sst.amj, growth_increment)) +
  geom_point() +  geom_abline(slope = coef(incr.lm5 )[[2]], intercept = coef(incr.lm5)[[1]])+
  xlab("April-June SST") + ylab("Growth increment")



# age 6----
means6 <- meanjoin[which(meanjoin$AGE==6),]
incr.null6 <- gam(growth_increment ~ s(sst.amj, k=4), 
                  data=means6)
gam.check(incr.null6) #
summary(incr.null6)#smooth not sig
plot(incr.null6) #
draw(incr.null6, select=1)

incr.lm6 <- lm(growth_increment ~ sst.amj, 
               data=means6)
summary(incr.lm6) #sst not sig

pg6 <- ggplot(means6, aes(sst.amj, growth_increment)) +
  geom_point() +  geom_abline(slope = coef(incr.lm6 )[[2]], intercept = coef(incr.lm6)[[1]])+
  xlab("April-June SST") + ylab("Growth increment")




# age 7----
means7 <- meanjoin[which(meanjoin$AGE==7),]
incr.null7 <- gam(growth_increment ~ s(sst.amj, k=4), 
                  data=means7)
gam.check(incr.null7) #
summary(incr.null7) #smooth is significant but edf=1
plot(incr.null7) #
draw(incr.null7, select=1)

incr.lm7 <- lm(growth_increment ~ sst.amj, 
               data=means7)
summary(incr.lm7) #

pg7 <- ggplot(means7, aes(sst.amj, growth_increment)) +
  geom_point() +  geom_abline(slope = coef(incr.lm7 )[[2]], intercept = coef(incr.lm7)[[1]])+
  xlab("April-June SST") + ylab("Growth increment")





# age 8----
means8 <- meanjoin[which(meanjoin$AGE==8),]
incr.null8 <- gam(growth_increment ~ s(sst.amj, k=4), 
                  data=means8)
gam.check(incr.null8) #
summary(incr.null8)#smooth not sig
plot(incr.null8) #
draw(incr.null8, select=1)


incr.lm8 <- lm(growth_increment ~ sst.amj, 
               data=means8)
summary(incr.lm8) #sst not sig

pg8 <- ggplot(means8, aes(sst.amj, growth_increment)) +
  geom_point()+  geom_abline(slope = coef(incr.lm8 )[[2]], intercept = coef(incr.lm8)[[1]])+
  xlab("April-June SST") + ylab("Growth increment")


# age 9----
means9 <- meanjoin[which(meanjoin$AGE==9),]
incr.null9 <- gam(growth_increment ~ s(sst.amj, k=4), 
                  data=means9)
gam.check(incr.null9) #
summary(incr.null9) #smooth significant but edf=1
plot(incr.null9) #
draw(incr.null9, select=1)

incr.lm9 <- lm(growth_increment ~ sst.amj, 
               data=means9)
summary(incr.lm9) #

pg9 <- ggplot(means9, aes(sst.amj, growth_increment)) +
  geom_point()+  geom_abline(slope = coef(incr.lm9 )[[2]], intercept = coef(incr.lm9)[[1]])+
  xlab("April-June SST") + ylab("Growth increment")


# age 10----
means10 <- meanjoin[which(meanjoin$AGE==10),]
incr.null10 <- gam(growth_increment ~ s(sst.amj, k=4), 
                  data=means10)
gam.check(incr.null10) #
summary(incr.null10) #smooth significant but edf=1
plot(incr.null10) #
draw(incr.null10, select=1)


incr.lm10 <- lm(growth_increment ~ sst.amj, 
               data=means10)
summary(incr.lm10) #

pg10 <- ggplot(means10, aes(sst.amj, growth_increment)) +
  geom_point() +  geom_abline(slope = coef(incr.lm10 )[[2]], intercept = coef(incr.lm10)[[1]])+
  xlab("April-June SST") + ylab("Growth increment")


# age 11----
means11 <- meanjoin[which(meanjoin$AGE==11),]
incr.null11 <- gam(growth_increment ~ s(sst.amj, k=4), 
                   data=means11)
gam.check(incr.null11) #
summary(incr.null11) #not significant
plot(incr.null11) #
draw(incr.null11, select=1)


incr.lm11 <- lm(growth_increment ~ sst.amj, 
                data=means11)
summary(incr.lm11) #

pg11 <- ggplot(means11, aes(sst.amj, growth_increment)) +
  geom_point() +  geom_abline(slope = coef(incr.lm11 )[[2]], intercept = coef(incr.lm11)[[1]]) +
  xlab("April-June SST") + ylab("Growth increment")


#plot all=====

plot_grid(pg2, pg3, pg4, pg5, pg6, pg7, pg8, pg9, pg10,pg11, ncol=3,
          labels = c('2', '3', '4', '5', '6', '7', '8', '9', '10', '11'))






