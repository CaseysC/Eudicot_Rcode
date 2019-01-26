#####################
## Celine Caseys
## Plant Science dept. UC davis
## January 2019
###################

source('~/Desktop/Github_Eudicot/1.Creating_Dataset.R')
source('~/Desktop/Github_Eudicot/2.Summary_RawData.R')

Data_cleaning <-  merge(Full.Data.7sp,Summary_IsoGen, by="Key_GenoIsol")

splt_Data_cleaning <- split(Data_cleaning, Data_cleaning$Taxon)
splt_Data_cleaning <- as.array(splt_Data_cleaning)

dimnames(splt_Data_cleaning) <- list(c("Brapa_Isogen", "Cendivia_IsoGen", "Cintybus_IsoGen", "Glycine_IsoGen", "Helianthus_IsoGen", "Lactuca_IsoGen", "Solanum_IsoGen"))

####################################### Data cleaning 

## To partition biological from technical failures, the lesion area distribution was analyzed 
## for each species and empirical thresholds were fixed (see below, median and size thresholds).
## A lesion below the size threshold threshold was considered a technical error only if the median 
## of lesion area for a plant genotype - strain pair was larger than the threshold. 
## The rationale is the following: when most lesions are of small size, the likelihood of biological 
## reasons for such small lesion areas is high, while when the majority of lesion areas are large, the likelihood of technical 
## error is high.

splt_Data_cleaning_filterOutliers <- list()
Failed_lesions <- list()

########### Species 1: Brassica rapa

splt_Data_cleaning[[1]] <- mutate(splt_Data_cleaning[[1]],
                                  median_logic=median.isogen > 0.15,
                                  threshold_logic=Scale.LS<0.1,
                                  TBC_param=as.factor(paste(median_logic, threshold_logic, sep="")),
                                  count_TBC=length(Scale.LS[TBC_param=="TRUETRUE"]))

# Select all data points that are not 'TRUETRUE'
splt_Data_cleaning_filterOutliers[[1]] <- splt_Data_cleaning[[1]] %>% filter(TBC_param!="TRUETRUE")

# Select all data points that are 'TRUETRUE', e.g. the data removed from the dataset
Failed_lesions[[1]] <- splt_Data_cleaning[[1]] %>% filter(TBC_param=="TRUETRUE")

################### Species 2: C.endivia

splt_Data_cleaning[[2]] <- mutate(splt_Data_cleaning[[2]], 
                                  median_logic=median.isogen > 0.2,
                                  threshold_logic=Scale.LS<0.15, 
                                  TBC_param=as.factor(paste(median_logic,threshold_logic, sep="")),
                                  count_TBC=length(Scale.LS[TBC_param=="TRUETRUE"]))

splt_Data_cleaning_filterOutliers[[2]] <- splt_Data_cleaning[[2]] %>% filter(TBC_param!="TRUETRUE")
Failed_lesions[[2]] <- splt_Data_cleaning[[2]] %>% filter(TBC_param=="TRUETRUE")

################### Species 3: C.intybus

splt_Data_cleaning[[3]] <- mutate(splt_Data_cleaning[[3]], 
                                  median_logic=median.isogen > 0.2,
                                  threshold_logic=Scale.LS<0.15, 
                                  TBC_param=as.factor(paste(median_logic,threshold_logic, sep="")),
                                  count_TBC=length(Scale.LS[TBC_param=="TRUETRUE"]))

splt_Data_cleaning_filterOutliers[[3]] <- splt_Data_cleaning[[3]] %>% filter(TBC_param!="TRUETRUE")
Failed_lesions[[3]] <- splt_Data_cleaning[[3]] %>% filter(TBC_param=="TRUETRUE")


################### Species 4: Glycine max

splt_Data_cleaning[[4]] <- mutate(splt_Data_cleaning[[4]], 
                                  median_logic=median.isogen > 0.15,
                                  threshold_logic=Scale.LS<0.1, 
                                  TBC_param=as.factor(paste(median_logic,threshold_logic, sep="")),
                                  count_TBC=length(Scale.LS[TBC_param=="TRUETRUE"]))

splt_Data_cleaning_filterOutliers[[4]] <- splt_Data_cleaning[[4]] %>% filter(TBC_param!="TRUETRUE")
Failed_lesions[[4]] <- splt_Data_cleaning[[4]] %>% filter(TBC_param=="TRUETRUE")

################### Species 5: Helianthus annuus

splt_Data_cleaning[[5]] <- mutate(splt_Data_cleaning[[5]], 
                                  median_logic=median.isogen > 0.3,
                                  threshold_logic=Scale.LS<0.25, 
                                  TBC_param=as.factor(paste(median_logic,threshold_logic, sep="")),
                                  count_TBC=length(Scale.LS[TBC_param=="TRUETRUE"]))

splt_Data_cleaning_filterOutliers[[5]] <- splt_Data_cleaning[[5]] %>% filter(TBC_param!="TRUETRUE")
Failed_lesions[[5]] <- splt_Data_cleaning[[5]] %>% filter(TBC_param=="TRUETRUE")

################### Species 6: Lactuca 

splt_Data_cleaning[[6]] <- mutate(splt_Data_cleaning[[6]], 
                                  median_logic=median.isogen > 0.15,
                                  threshold_logic=Scale.LS<0.1, 
                                  TBC_param=as.factor(paste(median_logic,threshold_logic, sep="")),
                                  count_TBC=length(Scale.LS[TBC_param=="TRUETRUE"]))

splt_Data_cleaning_filterOutliers[[6]] <- splt_Data_cleaning[[6]] %>% filter(TBC_param!="TRUETRUE")
Failed_lesions[[6]] <- splt_Data_cleaning[[6]] %>% filter(TBC_param=="TRUETRUE")


################### Species 7: Solanum

splt_Data_cleaning[[7]] <- mutate(splt_Data_cleaning[[7]], 
                                  median_logic=median.isogen > 0.15,
                                  threshold_logic=Scale.LS<0.1, 
                                  TBC_param=as.factor(paste(median_logic,threshold_logic, sep="")),
                                  count_TBC=length(Scale.LS[TBC_param=="TRUETRUE"]))

splt_Data_cleaning_filterOutliers[[7]] <- splt_Data_cleaning[[7]] %>% filter(TBC_param!="TRUETRUE")
Failed_lesions[[7]] <- splt_Data_cleaning[[7]] %>% filter(TBC_param=="TRUETRUE")



########################## Reassembling Data

Eudicot7sp_clean_median <- rbind(splt_Data_cleaning_filterOutliers[[1]], splt_Data_cleaning_filterOutliers[[2]], splt_Data_cleaning_filterOutliers[[3]], splt_Data_cleaning_filterOutliers[[4]], splt_Data_cleaning_filterOutliers[[5]], splt_Data_cleaning_filterOutliers[[6]], splt_Data_cleaning_filterOutliers[[7]])
#write.table(Eudicot7sp_clean_median, file="Eudicot7sp_clean_median.txt")

Names <- c("Brapa", "Cendivia", "Cintybus", "Glycine", "Helianthus", "Lactuca", "Solanum")

FailedLesions_7sp_median <- rbind(Failed_lesions[[1]], Failed_lesions[[2]], Failed_lesions[[3]], Failed_lesions[[4]],Failed_lesions[[5]],Failed_lesions[[6]],Failed_lesions[[7]])
#write.table(FailedLesions_7sp_median, file="FailedLesions_7sp.txt")




