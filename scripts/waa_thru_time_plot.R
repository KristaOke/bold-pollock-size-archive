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
lagdat <- lagdat[which(lagdat$AGE<15),] #was asked if borrowing info across yrs would allow
#more ages to be included

lagdat$AGE <- as.factor(lagdat$AGE)
lagdat$cohort <- as.factor(lagdat$cohort)

p1 <- ggplot(lagdat, aes(YEAR, log_sc_weight, colour=AGE)) + geom_point() + facet_wrap(~AGE)
p1

waa_avgs <- lagdat %>% group_by(AGE, YEAR) %>% summarise(mean_log_sc_weight = mean(log_sc_weight))

p2 <- ggplot(lagdat, aes(YEAR, log_sc_weight, colour=AGE)) + geom_point(alpha=0.2) + 
  geom_point(aes(YEAR, mean_log_sc_weight), colour=black, data=waa_avgs) + facet_wrap(~AGE)
p2 #fix colour issue


