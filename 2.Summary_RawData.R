#####################
## Author: Celine Caseys
## Plant Science dept. UC davis
## January 2019
###################

source('~/Desktop/Github_Eudicot/1.Creating_Dataset.R')
library(tidyverse)

TaxOrder <- Accessions.7sp[, c(2,7)]

myplots <- list()

## Calculate the mean, median, standard deviation, minimum, maximum of the lesions at different levels of the dataset. 
## The length is the number of datapoints considered for the summary statistics. 
## This is done at the taxon level, for each Botrytis isolate, for each plant genotype. 
## It is also calculated for pairs of plant genotypes - isolates and taxon-isolates. 

############### Data at the taxon level (Taxon = plant species)

Summary_species <- summarise(group_by(Full.Data.7sp, Taxon), 
                             mean.sp=mean(Scale.LS), 
                             median.sp=median(Scale.LS), 
                             sd.sp=sd(Scale.LS), 
                             min.sp=min(Scale.LS), 
                             max.sp=max(Scale.LS), 
                             length.sp=length(Scale.LS), 
                             cv.sp=sd.sp/mean.sp)

Summary_species <-  merge(Summary_species, TaxOrder, by="Taxon")

# write.table(Summary_species, file="Raw_Summary_species.txt")

############### Data at the Isolate level

Summary_Isolate <- summarise(group_by(Full.Data.7sp, IsolateID), 
                             mean.iso=mean(Scale.LS), 
                             median.iso=median(Scale.LS), 
                             sd.iso=sd(Scale.LS), 
                             min.iso=min(Scale.LS), 
                             max.iso=max(Scale.LS), 
                             cv.iso=sd.iso/mean.iso)

Summary_Isolate <-  merge(Summary_Isolate, B_host, by="IsolateID")

 # write.table(Summary_Isolate, file="Raw_Summary_Isolate.txt")
  
############### Data at the Genotype level

Summary_PGeno <- summarise(group_by(Full.Data.7sp, PlantGeno), 
                             mean.geno=mean(Scale.LS), 
                             median.geno=median(Scale.LS), 
                             sd.geno=sd(Scale.LS), 
                             min.geno=min(Scale.LS), 
                             max.geno=max(Scale.LS), 
                             length.geno=length(Scale.LS), 
                             cv.geno=sd.geno/mean.geno)

Summary_PGeno <-  merge(Summary_PGeno, Accessions.7sp[,c(2,3,6,7)], by="PlantGeno")
# write.table(Summary_PGeno, file="Raw_Summary_PGeno.txt")

############### Data at the Genotype x isolate level

Summary_IsoGen <- summarise(group_by(Full.Data.7sp, Key_GenoIsol), 
                             mean.isogen=mean(Scale.LS), 
                             median.isogen=median(Scale.LS), 
                             sd.isogen=sd(Scale.LS), 
                             min.isogen=min(Scale.LS), 
                             max.isogen=max(Scale.LS), 
                             length.isogen=length(Scale.LS), 
                             cv.isogen=sd.isogen/mean.isogen)

Summary_IsoGen$Key_GenoIsol1 <- Summary_IsoGen$Key_GenoIsol
Summary_IsoGen = separate(data=Summary_IsoGen, col=Key_GenoIsol1, into=c("PGeno", "Isolate"), sep="\\-")

# # write.table(Summary_IsoGen, file="Raw_Summary_IsoGen.txt")

############### Data at the Isolate x taxon  level

Summary_IsoTax <- summarise(group_by(Full.Data.7sp, Isolate_taxon_key), 
                            mean.isotax=mean(Scale.LS), 
                            median.isotax=median(Scale.LS), 
                            sd.isotax=sd(Scale.LS), 
                            min.isotax=min(Scale.LS), 
                            max.isotax=max(Scale.LS), 
                            length.isotax=length(Scale.LS), 
                            cv.isotax=sd.isotax/mean.isotax)

Summary_IsoTax$Isolate_taxon_key1 <- Summary_IsoTax$Isolate_taxon_key
Summary_IsoTax = separate(data=Summary_IsoTax, col=Isolate_taxon_key1, into=c("Isolate", "Taxon"), sep="\\-")

# write.table(Summary_IsoTax, file="Raw_Summary_IsoTax.txt")
