

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


newdf <- data.frame(matrix(ncol = 6, nrow = 15))
coln <- c("sst.amj", "AGE", "LONGITUDE", "LATITUDE", "julian", "cohort")
colnames(newdf) <- coln
newdf$AGE <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)
newdf$AGE <- as.factor(newdf$AGE)
newdf$LATITUDE <- median_lat
newdf$LONGITUDE <- median_long
newdf$julian <- median_julian
# newdf$cohort <- c(2019-1, 2019-2, 2019-3, 2019-4, 2019-5,
#                   2019-6, 2019-7, 2019-8, 2019-9, 2019-10,
#                   2019-11, 2019-12, 2019-13, 2019-14, 2019-15)
newdf$cohort <- as.factor("2019") #setting to an arbitary value
newdf$cohort <- as.factor(newdf$cohort)
#above stays the same
newdf$sst.amj <- 2 #EDIT HERE BASED ON TEMP OF INTEREST

predicted_waa_2deg <- predict(cre_all1$gam, newdata=newdf)
predicted_waa_2.5deg <- predict(cre_all1$gam, newdata=newdf)
predicted_waa_3deg <- predict(cre_all1$gam, newdata=newdf)
predicted_waa_3.5deg <- predict(cre_all1$gam, newdata=newdf)
predicted_waa_4deg <- predict(cre_all1$gam, newdata=newdf)
predicted_waa_4.5deg <- predict(cre_all1$gam, newdata=newdf)
predicted_waa_5deg <- predict(cre_all1$gam, newdata=newdf)

predicted_waa_df <- as.data.frame(rbind(predicted_waa_2deg,
                          predicted_waa_2.5deg,
                          predicted_waa_3deg,
                          predicted_waa_3.5deg,
                          predicted_waa_4deg,
                          predicted_waa_4.5deg,
                          predicted_waa_5deg))

write_csv(predicted_waa_df, file=paste(wd,"/predicted_waa_df.csv", sep=""))
