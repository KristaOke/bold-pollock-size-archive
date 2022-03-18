#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#making plot of size thru time

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Notes:
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

library(ggplot2)
library(tidyverse)

wd <- getwd()
lagdat <- read.csv(file=paste(wd,"/data/lagdat.csv", sep=""), row.names=1)

lagdat$cohort <- lagdat$YEAR - lagdat$AGE

table(lagdat$AGE)
lagdat <- lagdat[which(lagdat$AGE<16),] #was asked if borrowing info across yrs would allow
#more ages to be included

lagdat$AGE <- as.factor(lagdat$AGE)
lagdat$cohort <- as.factor(lagdat$cohort)

p1 <- ggplot(lagdat, aes(YEAR, log_sc_weight, colour=AGE)) + geom_point() + facet_wrap(~AGE)
p1

waa_avgs <- lagdat %>% group_by(AGE, YEAR) %>% summarise(mean_log_sc_weight = mean(log_sc_weight, na.rm=TRUE))

#scale data w avgs
p2 <- ggplot(lagdat, aes(YEAR, log_sc_weight, colour=AGE)) + geom_point(alpha=0.2) + 
  geom_point(aes(YEAR, mean_log_sc_weight), colour="black", data=waa_avgs) + 
  geom_line(aes(YEAR, mean_log_sc_weight), colour="black", data=waa_avgs) +
  facet_wrap(~AGE)
p2 

#just avgs
p3 <- ggplot(waa_avgs, aes(YEAR, mean_log_sc_weight, colour=AGE)) + geom_line(aes(linetype=AGE)) + 
  theme_bw()
p3

#raw data w avgs

raw_avgs <- lagdat %>% group_by(AGE, YEAR) %>% summarise(mean_raw_weight = mean(WEIGHT, na.rm=TRUE), n=n())

p4 <- ggplot(lagdat, aes(YEAR, WEIGHT)) + geom_point(colour="dark grey", alpha=0.5) + 
   geom_point(aes(YEAR, mean_raw_weight), colour="black", data=raw_avgs) + 
   geom_line(aes(YEAR, mean_raw_weight), colour="black", data=raw_avgs) +
  facet_wrap(~AGE, scales="free") + theme_bw()
p4 

#same but with samples less than 5 excluded
p5 <- ggplot(lagdat, aes(YEAR, WEIGHT)) + geom_point(colour="dark grey", alpha=0.5) + 
  geom_point(aes(YEAR, mean_raw_weight), colour="black", data=raw_avgs[which(raw_avgs$n>4),]) + 
  geom_line(aes(YEAR, mean_raw_weight), colour="black", data=raw_avgs[which(raw_avgs$n>4),]) +
  facet_wrap(~AGE, scales="free") + theme_bw()
p5 

#want broken lines so will have to create a new column with NA for any avgs with sample sizes too small
raw_avgs$mean_raw_weight_n_over_5 <- raw_avgs$mean_raw_weight
raw_avgs$mean_raw_weight_n_over_5[raw_avgs$n<5] <- NA
#double check yep looks nice
#only happens 7 times, all since 2016, all ages greater than 12

p6 <- ggplot(lagdat, aes(YEAR, WEIGHT)) + geom_point(colour="dark grey", alpha=0.5) + 
  geom_point(aes(YEAR, mean_raw_weight_n_over_5), colour="black", data=raw_avgs) + 
  geom_line(aes(YEAR, mean_raw_weight_n_over_5), colour="black", data=raw_avgs) +
  facet_wrap(~AGE, scales="free") + theme_bw()
p6 


#what if I colour points by cool/warm?
raw_avgs$period <- NA
raw_avgs$period[raw_avgs$YEAR==1999] <- "cool"
raw_avgs$period[raw_avgs$YEAR>1999 & raw_avgs$YEAR<2006] <- "warm"
raw_avgs$period[raw_avgs$YEAR>2005 & raw_avgs$YEAR<2014] <- "cool"
raw_avgs$period[raw_avgs$YEAR>2013] <- "warm"

p7 <- ggplot(lagdat, aes(YEAR, WEIGHT)) + geom_point(colour="dark grey", alpha=0.5) + 
  geom_point(aes(YEAR, mean_raw_weight_n_over_5, colour=period), data=raw_avgs) +
  scale_color_manual(values=c("blue","red"))+
  geom_line(aes(YEAR, mean_raw_weight_n_over_5), colour="black", data=raw_avgs) +
  facet_wrap(~AGE, scales="free") + theme_bw() +
  ylab("Weight-at-age (g)") + xlab("Year") + theme(legend.position=c(0.9, 0.1))
p7 


#table of hauls---------------------------------------------------------------------------

#was also a request for a table of hauls let's try to do that here too

table(lagdat$YEAR, lagdat$HAUL)

#hauldat <- lagdat %>% group_by(YEAR, HAUL) %>% summarize(n=n()) #this gets n in haul not quite

hauldat <- lagdat %>% group_by(YEAR) %>% summarize(n_hauls=n_distinct(HAUL)) 
#double checked looks right









