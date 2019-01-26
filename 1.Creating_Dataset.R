#####################
## Author: Celine Caseys
## Plant Science dept. UC davis
## January 2019
###################

library(tidyverse)
library(stringr)

Full.Data.7sp <- read.table(file="FullMetaDat_7sp.txt", header=T)
Accessions.7sp <- read.table(file="accession.list7sp.txt", header=T)

B_host <-  read.table(file="Isolate_Host_origin_CORRECTED.txt", header=T)
colnames(B_host)[1] <-  c("IsolateID")

#Replace all . with _ to avoid confusion of numbers with dates
Full.Data.7sp$IsolateID <- as.factor(gsub('\\.', '_', Full.Data.7sp$IsolateID, fixed=F))
# Typos as Gallo3 instead of Gallo2
levels(Full.Data.7sp$IsolateID)[76] <- c("Gallo2")

# gsub perform replacement on all matches
# Correct isolate names (no dot to avoid date format issues)
Full.Data.7sp$PlantGeno <- as.factor(gsub('\\-', '.', Full.Data.7sp$PlantGeno, fixed=F))

#Merging datasets

Full.Data.7sp <- merge(Full.Data.7sp, Accessions.7sp[, -c(1,2)], by="PlantGeno")
Full.Data.7sp <- merge(Full.Data.7sp, B_host, by="IsolateID")


#create column Taxon_genotype (to use as key)
Full.Data.7sp$Key_GenoTax <- as.factor(paste(Full.Data.7sp$Taxon, Full.Data.7sp$PlantGeno, sep="-"))

#create column Isolate_genotype (to use as key)
Full.Data.7sp$Key_GenoIsol <- as.factor(paste(Full.Data.7sp$PlantGeno,Full.Data.7sp$IsolateID, sep="-"))

#create column Isolate_Taxon (to use as key)
Full.Data.7sp$Isolate_taxon_key <- as.factor(paste(Full.Data.7sp$IsolateID, Full.Data.7sp$Taxon, sep="-"))

#Create column with Taxon and Domestication

levels(Full.Data.7sp$Domest) <- c("D", "D", "D", "L", "W", "W", "W")
Full.Data.7sp$Improvement <- Full.Data.7sp$Domest
levels(Full.Data.7sp$Improvement) <- c("High", "Low", "Low")
Full.Data.7sp$ImprovS <- Full.Data.7sp$Improvement
levels(Full.Data.7sp$ImprovS) <- c("H", "L")

Full.Data.7sp$TaxonS <- Full.Data.7sp$Taxon
levels(Full.Data.7sp$TaxonS) <- c("Bra", "Cen", "Cin", "Gly", "Hel", "Lac", "Sol")
### Short names for plotting

Full.Data.7sp$Wi_Do <- as.factor(paste(Full.Data.7sp$TaxonS, Full.Data.7sp$Domest, sep="-"))
Full.Data.7sp$ImpTax <- as.factor(paste(Full.Data.7sp$TaxonS, Full.Data.7sp$ImprovS, sep="-"))

#create column Isolate_Taxon_WiDo (to use as key)
Full.Data.7sp$Isolate_taxon_Wido_key <- as.factor(paste(Full.Data.7sp$IsolateID, Full.Data.7sp$ImpTax, sep="-"))

#create column Isolate x order

Full.Data.7sp$OrderS <- Full.Data.7sp$Order
levels(Full.Data.7sp$OrderS) <- c("Aster", "Brass", "Fab", "Solan")
Full.Data.7sp$IsORder <- as.factor(paste(Full.Data.7sp$IsolateID, Full.Data.7sp$OrderS, sep="-"))

#write.table(Full.Data.7sp, file="Full.Data.7sp.txt")

#------------------------------------ footnotes
# accession.list <- as.data.frame(unique(Full.Data.7sp$Key_GenoTax))
# creating the original plant accession list
# write.table(accession.list, file="accession.list7sp.txt")
