

median_julian <- median(lagdat$julian)

#get catch COG to use in predict--------

#will use only 2019 since most recent
lag2019 <- lagdat[which(lagdat$YEAR==2019),]

#going to need the other survey data I think? With abundance?
#cog_lat <- sum()/sum()

#for now let's try with median values which look VERY similar
#to the COG estimated by Thorson et al 2017 for individuals 50cm and up
median_lat <- median(lagdat$LATITUDE)
median_long <- median(lagdat$LONGITUDE)

wd <- getwd()
cre_all1 <- readRDS(file=paste(wd,"/scripts/size scripts/model_output_all-ages_cohort-as-re_age15.rds", sep=""))

#start by creating a df with the observed temps to predict
tempseries <- unique(lagdat$sst.amj)
outputdf <- NULL
outputdf <- data.frame(matrix(ncol = 6, nrow = 0))
coln <- c("sst.amj", "AGE", "LONGITUDE", "LATITUDE", "julian", "cohort")
colnames(outputdf) <- coln

i <- 1
for(i in 1:length(tempseries)){
newdf <- data.frame(matrix(ncol = 6, nrow = 15))
coln <- c("sst.amj", "AGE", "LONGITUDE", "LATITUDE", "julian", "cohort")
colnames(newdf) <- coln
newdf$AGE <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)
newdf$AGE <- as.factor(newdf$AGE)
newdf$LATITUDE <- median_lat
newdf$LONGITUDE <- median_long
newdf$julian <- median_julian
newdf$cohort <- as.factor("2018") #setting to an arbitrary value
#newdf$cohort <- as.factor(newdf$cohort)
#above stays the same
newdf$sst.amj <- tempseries[i]

outputdf <- rbind(outputdf, newdf)
} #looks right


outputdf$predicted_values <- predict(cre_all1$gam, newdata=outputdf)
predicted_waa_obs_temps <- outputdf

#get the mean for each age
age_mean_predicted <- predicted_waa_obs_temps %>% group_by(AGE) %>%
  summarize(mean_predicted = mean(predicted_waa_obs_temps$predicted_values, na.rm=TRUE))

ggplot(predicted_waa_obs_temps, aes(sst.amj, predicted_values)) + geom_point() +
  facet_wrap(~AGE)

write_csv(predicted_waa_df, file=paste(wd,"/predicted_waa_df.csv", sep=""))
