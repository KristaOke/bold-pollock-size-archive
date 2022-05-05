

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
  facet_wrap(~AGE) #looks right!

#the original scaling was done using scale, which subtracts mean and divides by sd
#I am going to make sure that I get the same answer doing that manually 
#to make sure that I'm getting the same sd as scale used


test_scaling <- lagdat %>% group_by(AGE) %>%
  summarize(sd_logW = sd(logWEIGHT, na.rm=TRUE), mean_logW = mean(logWEIGHT, na.rm=TRUE),
            scaling = (logWEIGHT - mean_logW)/sd_logW)
test_scaling$scaling == lagdat$log_sc_weight #look same

#so now we will use the age specific sd to multiply

age_sds <- lagdat %>% group_by(AGE) %>%
  summarize(sd_logW = sd(logWEIGHT, na.rm=TRUE))

#first get predicted values for new temps
#very similar to the above but now predict for series of temps increasing at 0.1 degree

#start by creating a df with the observed temps to predict
range(lagdat$sst.amj)
newtempseries <- seq(from=1.8, to=5.4, by=0.1)
newoutputdf <- NULL
newoutputdf <- data.frame(matrix(ncol = 6, nrow = 0))
coln <- c("sst.amj", "AGE", "LONGITUDE", "LATITUDE", "julian", "cohort")
colnames(newoutputdf) <- coln

i <- 1
for(i in 1:length(newtempseries)){
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
  newdf$sst.amj <- newtempseries[i]
  
  newoutputdf <- rbind(newoutputdf, newdf)
} #looks right


newoutputdf$predicted_values <- predict(cre_all1$gam, newdata=newoutputdf)
new_predicted_waa_temps <- newoutputdf

#now subtract the mean predicted valuye from observed temps
age_mean_predicted
#the mean is the same for all ages, we can grab age 1 assuming that remains true (check!)
mean_predicted <- as.numeric(age_mean_predicted[1,2])

new_predicted_waa_temps$predicted_minus_mean <- 
  new_predicted_waa_temps$predicted_values - mean_predicted

#and multiply by sd by age
new_predicted_waa_temps$corrected_predicted_values <- NA
len_values <- length(new_predicted_waa_temps$predicted_minus_mean)
h <- 1
for(h in 1:len_values){
  new_predicted_waa_temps$corrected_predicted_values[h] <- 
    new_predicted_waa_temps$predicted_minus_mean[h] * 
    age_sds$sd_logW[which(age_sds$AGE==new_predicted_waa_temps$AGE[h])]
} #looks good

View(new_predicted_waa_temps)

#WEIRD columns are arrays, fix that
new_predicted_waa_temps$predicted_values <- as.numeric(new_predicted_waa_temps$predicted_values)
new_predicted_waa_temps$corrected_predicted_values <- as.numeric(new_predicted_waa_temps$corrected_predicted_values)
new_predicted_waa_temps$predicted_minus_mean <- as.numeric(new_predicted_waa_temps$predicted_minus_mean)

#need to overwrite
write_csv(new_predicted_waa_temps, file=paste(wd,"/predicted_waa_df.csv", sep="")) 

#and another output for ONLY temps of interest
pred_interest <- new_predicted_waa_temps[which(new_predicted_waa_temps$sst.amj=="2"|
                                                 new_predicted_waa_temps$sst.amj=="2.5"|
                                               new_predicted_waa_temps$sst.amj=="3"|
                                               new_predicted_waa_temps$sst.amj=="3.5"|
                                               new_predicted_waa_temps$sst.amj=="4"|
                                               new_predicted_waa_temps$sst.amj=="4.5"|
                                               new_predicted_waa_temps$sst.amj=="5"),]
age_temp_value <- pred_interest[,c("sst.amj", "AGE", "corrected_predicted_values")]
#transpose to make it easier to enter into excel
Tage_temp_value <- as.data.frame(t(age_temp_value))


write_csv(Tage_temp_value, file=paste(wd,"/transposed_predictions_interest_temps.csv", sep="")) 
