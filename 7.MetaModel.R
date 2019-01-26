library(MASS)
library(lme4)
library(lmerTest)
library(car)
library(lsmeans)

source('~/Desktop/Github_Eudicot/1.Creating_Dataset.R')

Eudicot7sp_lsmeans <- read.table(file="LSmeanIndP_7spData_row_Geno.txt", header=T)
colnames(Eudicot7sp_lsmeans)[7] <- c("PlantGeno")

Eudicot7sp_lsmeans <- merge(Eudicot7sp_lsmeans, Accessions.7sp, by="PlantGeno")

Eudicot7sp_lsmeans$Domest <- Eudicot7sp_lsmeans$Wi_Do
levels(Eudicot7sp_lsmeans$Domest) <- c("D", "L", "W")
Eudicot7sp_lsmeans$Improvement <- Eudicot7sp_lsmeans$Domest
levels(Eudicot7sp_lsmeans$Improvement) <- c("High", "Low", "Low")
Eudicot7sp_lsmeans$ImprovS <- Eudicot7sp_lsmeans$Improvement
levels(Eudicot7sp_lsmeans$ImprovS) <- c("H", "L")

Eudicot7sp_lsmeans <- merge(Eudicot7sp_lsmeans, B_host, by="IsolateID")

#Nested meta model
# To run the model, uncomment the following lines.
# ! This model takes a few days of computation to completely run !

#Meta3<- lm(emmean ~ IsolateID + Taxon + Taxon/Improvement + Taxon/PlantGeno + IsolateID*Taxon + IsolateID*Taxon/Improvement + IsolateID*Taxon/PlantGeno , data=Eudicot7sp_lsmeans)
#Anova_meta <- anova(Meta3)
