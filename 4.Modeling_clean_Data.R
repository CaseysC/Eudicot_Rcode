#####################
## Celine Caseys
## Plant Science dept. UC davis
## January 2019
###################


library(MASS)
library(lme4)
library(lmerTest)
library(car)
library(lsmeans)


Eudicot_clean <- read.table(file="Eudicot7sp_clean_median.txt", header=T)

Eudicot_clean$IndPlant1 <- as.factor(paste(Eudicot_clean$PlantGeno,Eudicot_clean$Exp,Eudicot_clean$IndPlant, sep="-"))

EudicotC_data_Sp <- split(Eudicot_clean, Eudicot_clean$Taxon)
EudicotC_data_Sp  <- as.array(EudicotC_data_Sp)
list2env(EudicotC_data_Sp , envir=.GlobalEnv)
rm(EudicotC_data_Sp)


##################_____________________________________ Lactuca

############ Correction for the way experimental flats and replicates are coded.
Lactuca <- droplevels.data.frame(Lactuca)
Lactuca$Rep <- as.factor(paste(Lactuca$Exp, Lactuca$Rep, sep="-"))
Lactuca$Flat <- as.factor(paste(Lactuca$Rep, Lactuca$Flat, sep="-"))

print(paste("starting Lactuca", Sys.time()))
Lactuca_Cmodel <- lm(Scale.LS ~ IsolateID + Improvement/PlantGeno+ IsolateID*Improvement + IsolateID*Improvement/PlantGeno + Exp + Flat + IndPlant1, data=Lactuca)
Lactuca_Cmodel_anova<- Anova(Lactuca_Cmodel)

Lactuca_Cmodel <- lmer(Scale.LS ~ IsolateID + Improvement/PlantGeno+ IsolateID*Improvement + IsolateID*Improvement/PlantGeno + (1|Exp) + (1|Flat) + (1|IndPlant1), data=Lactuca)
Lactuca_Cmodel_anova<- Anova(Lactuca_Cmodel)
write.table(Lactuca_Cmodel_anova, file="Lactuca_Cmodel_anova.txt")

print(paste("starting Lactuca random", Sys.time()))
Lactuca_Cmodel_random <- rand(Lactuca_Cmodel)
Lactuca_Cmodel_random  <- as.data.frame(print(Lactuca_Cmodel_random))
write.table(Lactuca_Cmodel_random, file="Lactuca_Cmodel_random.txt")

print(paste("starting Lactuca random", Sys.time()))
Lactuca_Cmodel_residual <- residuals(Lactuca_Cmodel)
write.table(Lactuca_Cmodel_residual, file="Lactuca_Cmodel_residual.txt")

print(paste("starting Lactuca lsmeans", Sys.time()))
Lactuca_Cmodel_iso <- lmer(Scale.LS ~ IsolateID + (1|Exp) + (1|Flat) + (1|IndPlant1), data=Lactuca)
Lactuca_Cmodel_iso_anova <- anova(Lactuca_Cmodel_iso)  
write.table(Lactuca_Cmodel_iso_anova, file="Lactuca_Cmodel_iso_anova.txt")
LesionIso_lsm_Lactuca <- lsmeans(Lactuca_Cmodel_iso, "IsolateID")
write.table(LesionIso_lsm_Lactuca, file="LesionIso_lsm_Lactuca.txt")

Lactuca_Cmodel_Geno <- lmer(Scale.LS ~ PlantGeno + (1|Exp) + (1|Flat) + (1|IndPlant1), data=Lactuca) 
Lactuca_Cmodel_Geno_anova <- anova(  Lactuca_Cmodel_Geno)  
write.table(Lactuca_Cmodel_Geno_anova, file="Lactuca_Cmodel_Geno_anova.txt")
LesionGeno_lsm_Lactuca <- lsmeans(Lactuca_Cmodel_Geno, "PlantGeno")
LesionGeno_lsm_Lactuca <- as.data.frame(print(LesionGeno_lsm_Lactuca))
write.table(LesionGeno_lsm_Lactuca, file="LesionGeno_lsm_Lactuca.txt")

Lactuca_Cmodel_improv <- lmer(Scale.LS ~ Improvement + (1|Exp) + (1|Flat) + (1|IndPlant1), data=Lactuca) 
Lactuca_Cmodel_improv_anova <- anova(Lactuca_Cmodel_improv)  
write.table(Lactuca_Cmodel_improv_anova, file="Lactuca_Cmodel_improv_anova.txt")
LesionDom_lsm_Lactuca <- lsmeans(Lactuca_Cmodel_improv, "Improvement")
LesionDom_lsm_Lactuca  <- as.data.frame(print(LesionDom_lsm_Lactuca))
write.table(LesionDom_lsm_Lactuca, file="LesionDom_lsm_Lactuca.txt")

################################################## Helianthus
rm(list=ls())

Eudicot_clean <- read.table(file="~/Documents/Science_projects/Davis_project/Eudicot project/DataForCeline/7species_data/Eudicot7sp_clean_median.txt", header=T)
Eudicot_clean$IndPlant1 <- as.factor(paste(Eudicot_clean$PlantGeno,Eudicot_clean$Exp,Eudicot_clean$IndPlant, sep="-"))

EudicotC_data_Sp <- split(Eudicot_clean, Eudicot_clean$Taxon)
EudicotC_data_Sp  <- as.array(EudicotC_data_Sp)
list2env(EudicotC_data_Sp , envir=.GlobalEnv)
rm(EudicotC_data_Sp)

Helianthus <- droplevels.data.frame(Helianthus)
print(paste("starting Helianthus", Sys.time()))

Helianthus_Cmodel <- lmer(Scale.LS ~ IsolateID + Improvement/PlantGeno+ IsolateID*Improvement + IsolateID*Improvement/PlantGeno + (1|Exp) + (1|Flat) + (1|IndPlant1), data=Helianthus)
Helianthus_Cmodel_anova<- Anova(Helianthus_Cmodel)
write.table(Helianthus_Cmodel_anova, file="Helianthus_Cmodel_anova.txt")

print(paste("starting Helianthus random", Sys.time()))
Helianthus_Cmodel_random <- rand(Helianthus_Cmodel)
Helianthus_Cmodel_random  <- as.data.frame(print(Helianthus_Cmodel_random))
write.table(Helianthus_Cmodel_random, file="Helianthus_Cmodel_random.txt")

print(paste("starting Helianthus residual", Sys.time()))
Helianthus_Cmodel_residual <- residuals(Helianthus_Cmodel)
Helianthus_Cmodel_residual  <- as.data.frame(print(Helianthus_Cmodel_residual))
write.table(Helianthus_Cmodel_residual, file="Helianthus_Cmodel_residual.txt")

print(paste("starting Helianthus lsmeans", Sys.time()))
Helianthus_Cmodel_iso <- lmer(Scale.LS ~ IsolateID + (1|Exp) + (1|Flat) + (1|IndPlant1), data=Helianthus)
Helianthus_Cmodel_iso_anova <- anova(Helianthus_Cmodel_iso)  
write.table(Helianthus_Cmodel_iso_anova, file="Helianthus_Cmodel_iso_anova.txt")
LesionIso_lsm_Helianthus <- lsmeans(Helianthus_Cmodel_iso, "IsolateID")
write.table(LesionIso_lsm_Helianthus, file="LesionIso_lsm_Helianthus.txt")

Helianthus_Cmodel_Geno <- lmer(Scale.LS ~ PlantGeno + (1|Exp) + (1|Flat) + (1|IndPlant1), data=Helianthus) 
Helianthus_Cmodel_Geno_anova <- anova(  Helianthus_Cmodel_Geno)  
write.table(Helianthus_Cmodel_Geno_anova, file="Helianthus_Cmodel_Geno_anova.txt")
LesionGeno_lsm_Helianthus <- lsmeans(Helianthus_Cmodel_Geno, "PlantGeno")
LesionGeno_lsm_Helianthus <- as.data.frame(print(LesionGeno_lsm_Helianthus))
write.table(LesionGeno_lsm_Helianthus, file="LesionGeno_lsm_Helianthus.txt")

Helianthus_Cmodel_improv <- lmer(Scale.LS ~ Improvement + (1|Exp) + (1|Flat) + (1|IndPlant1), data=Helianthus) 
Helianthus_Cmodel_improv_anova <- anova(Helianthus_Cmodel_improv)  
write.table(Helianthus_Cmodel_improv_anova, file="Helianthus_Cmodel_improv_anova.txt")
LesionDom_lsm_Helianthus <- lsmeans(Helianthus_Cmodel_improv, "Improvement")
LesionDom_lsm_Helianthus  <- as.data.frame(print(LesionDom_lsm_Helianthus))
write.table(LesionDom_lsm_Helianthus, file="LesionDom_lsm_Helianthus.txt")

################------------------------------ Brapa

rm(list=ls())

Eudicot_clean <- read.table(file="~/Documents/Science_projects/Davis_project/Eudicot project/DataForCeline/7species_data/Eudicot7sp_clean_median.txt", header=T)
Eudicot_clean$IndPlant1 <- as.factor(paste(Eudicot_clean$PlantGeno,Eudicot_clean$Exp,Eudicot_clean$IndPlant, sep="-"))

EudicotC_data_Sp <- split(Eudicot_clean, Eudicot_clean$Taxon)
EudicotC_data_Sp  <- as.array(EudicotC_data_Sp)
list2env(EudicotC_data_Sp , envir=.GlobalEnv)
rm(EudicotC_data_Sp)

Brapa <- droplevels.data.frame(Brapa)

print(paste("starting Brapa", Sys.time()))

Brapa_Cmodel <- lmer(Scale.LS ~ IsolateID + Improvement/PlantGeno+ IsolateID*Improvement + IsolateID*Improvement/PlantGeno + (1|Exp) + (1|Flat) + (1|IndPlant1), data=Brapa)

Brapa_Cmodel_anova<- Anova(Brapa_Cmodel)
write.table(Brapa_Cmodel_anova, file="Brapa_Cmodel_anova.txt")

print(paste("starting Brapa random", Sys.time()))
Brapa_Cmodel_random <- rand(Brapa_Cmodel)
Brapa_Cmodel_random  <- as.data.frame(print(Brapa_Cmodel_random))
write.table(Brapa_Cmodel_random, file="Brapa_Cmodel_random.txt")

print(paste("starting Brapa residual", Sys.time()))
Brapa_Cmodel_residual <- residuals(Brapa_Cmodel)
Brapa_Cmodel_residual  <- as.data.frame(print(Brapa_Cmodel_residual))
write.table(Brapa_Cmodel_residual, file="Brapa_Cmodel_residual.txt")

print(paste("starting Brapa lsmeans", Sys.time()))

Brapa_Cmodel_iso <- lmer(Scale.LS ~ IsolateID + (1|Exp) + (1|Flat) + (1|IndPlant1), data=Brapa)
Brapa_Cmodel_iso_anova <- anova(Brapa_Cmodel_iso) 
Brapa_Cmodel_iso_rand <- rand(Brapa_Cmodel_iso) 

write.table(Brapa_Cmodel_iso_anova, file="Brapa_Cmodel_iso_anova.txt")
LesionIso_lsm_Brapa <- lsmeans(Brapa_Cmodel_iso, "IsolateID")
write.table(LesionIso_lsm_Brapa, file="LesionIso_lsm_Brapa.txt")

Brapa_Cmodel_Geno <- lmer(Scale.LS ~ PlantGeno + (1|Exp) + (1|Flat) + (1|IndPlant1), data=Brapa) 
Brapa_Cmodel_Geno_anova <- anova(  Brapa_Cmodel_Geno)  
write.table(Brapa_Cmodel_Geno_anova, file="Brapa_Cmodel_Geno_anova.txt")
LesionGeno_lsm_Brapa <- lsmeans(Brapa_Cmodel_Geno, "PlantGeno")
LesionGeno_lsm_Brapa <- as.data.frame(print(LesionGeno_lsm_Brapa))
write.table(LesionGeno_lsm_Brapa, file="LesionGeno_lsm_Brapa.txt")

Brapa_Cmodel_improv <- lmer(Scale.LS ~ Improvement + (1|Exp) + (1|Flat) + (1|IndPlant1), data=Brapa) 
Brapa_Cmodel_improv_anova <- anova(Brapa_Cmodel_improv)  
write.table(Brapa_Cmodel_improv_anova, file="Brapa_Cmodel_improv_anova.txt")
LesionDom_lsm_Brapa <- lsmeans(Brapa_Cmodel_improv, "Improvement")
LesionDom_lsm_Brapa  <- as.data.frame(print(LesionDom_lsm_Brapa))
write.table(LesionDom_lsm_Brapa, file="LesionDom_lsm_Brapa.txt")

#############____________________________________ Cendivia

rm(list=ls())

Eudicot_clean <- read.table(file="~/Documents/Science_projects/Davis_project/Eudicot project/DataForCeline/7species_data/Eudicot7sp_clean_median.txt", header=T)
Eudicot_clean$IndPlant1 <- as.factor(paste(Eudicot_clean$PlantGeno,Eudicot_clean$Exp,Eudicot_clean$IndPlant, sep="-"))

EudicotC_data_Sp <- split(Eudicot_clean, Eudicot_clean$Taxon)
EudicotC_data_Sp  <- as.array(EudicotC_data_Sp)
list2env(EudicotC_data_Sp , envir=.GlobalEnv)
rm(EudicotC_data_Sp)

Cendivia <- droplevels.data.frame(Cendivia)

print(paste("starting Cendivia", Sys.time()))

Cendivia_Cmodel <- lmer(Scale.LS ~ IsolateID + Improvement/PlantGeno+ IsolateID*Improvement + IsolateID*Improvement/PlantGeno + (1|Exp) + (1|Flat) + (1|IndPlant1), data=Cendivia)


Cendivia_Cmodel_anova<- Anova(Cendivia_Cmodel)
write.table(Cendivia_Cmodel_anova, file="Cendivia_Cmodel_anova.txt")

print(paste("starting Cendivia random", Sys.time()))
Cendivia_Cmodel_random <- rand(Cendivia_Cmodel)
Cendivia_Cmodel_random  <- as.data.frame(print(Cendivia_Cmodel_random))
write.table(Cendivia_Cmodel_random, file="Cendivia_Cmodel_random.txt")

print(paste("starting Cendivia residual", Sys.time()))
Cendivia_Cmodel_residual <- residuals(Cendivia_Cmodel)
Cendivia_Cmodel_residual  <- as.data.frame(print(Cendivia_Cmodel_residual))
write.table(Cendivia_Cmodel_residual, file="Cendivia_Cmodel_residual.txt")

print(paste("starting Cendivia lsmean", Sys.time()))

Cendivia_Cmodel_iso <- lmer(Scale.LS ~ IsolateID + (1|Exp) + (1|Flat) + (1|IndPlant1), data=Cendivia)
Cendivia_Cmodel_iso_anova <- anova(Cendivia_Cmodel_iso)  
write.table(Cendivia_Cmodel_iso_anova, file="Cendivia_Cmodel_iso_anova.txt")
LesionIso_lsm_Cendivia <- lsmeans(Cendivia_Cmodel_iso, "IsolateID")
write.table(LesionIso_lsm_Cendivia, file="LesionIso_lsm_Cendivia.txt")

Cendivia_Cmodel_Geno <- lmer(Scale.LS ~ PlantGeno + (1|Exp) + (1|Flat) + (1|IndPlant1), data=Cendivia) 
Cendivia_Cmodel_Geno_anova <- anova(  Cendivia_Cmodel_Geno)  
write.table(Cendivia_Cmodel_Geno_anova, file="Cendivia_Cmodel_Geno_anova.txt")
LesionGeno_lsm_Cendivia <- lsmeans(Cendivia_Cmodel_Geno, "PlantGeno")
LesionGeno_lsm_Cendivia <- as.data.frame(print(LesionGeno_lsm_Cendivia))
write.table(LesionGeno_lsm_Cendivia, file="LesionGeno_lsm_Cendivia.txt")

Cendivia_Cmodel_improv <- lmer(Scale.LS ~ Improvement + (1|Exp) + (1|Flat) + (1|IndPlant1), data=Cendivia) 
Cendivia_Cmodel_improv_anova <- anova(Cendivia_Cmodel_improv)  
write.table(Cendivia_Cmodel_improv_anova, file="Cendivia_Cmodel_improv_anova.txt")
LesionDom_lsm_Cendivia <- lsmeans(Cendivia_Cmodel_improv, "Improvement")
LesionDom_lsm_Cendivia  <- as.data.frame(print(LesionDom_lsm_Cendivia))
write.table(LesionDom_lsm_Cendivia, file="LesionDom_lsm_Cendivia.txt")

#############____________________________________________ Cintybus

rm(list=ls())

Eudicot_clean <- read.table(file="~/Documents/Science_projects/Davis_project/Eudicot project/DataForCeline/7species_data/Eudicot7sp_clean_median.txt", header=T)
Eudicot_clean$IndPlant1 <- as.factor(paste(Eudicot_clean$PlantGeno,Eudicot_clean$Exp,Eudicot_clean$IndPlant, sep="-"))

EudicotC_data_Sp <- split(Eudicot_clean, Eudicot_clean$Taxon)
EudicotC_data_Sp  <- as.array(EudicotC_data_Sp)
list2env(EudicotC_data_Sp , envir=.GlobalEnv)
rm(EudicotC_data_Sp)

Cintybus <- droplevels.data.frame(Cintybus)

print(paste("starting Cintybus", Sys.time()))

Cintybus_Cmodel <- lmer(Scale.LS ~ IsolateID + Improvement/PlantGeno+ IsolateID*Improvement + IsolateID*Improvement/PlantGeno + (1|Exp) + (1|Flat) + (1|IndPlant1), data=Cintybus)

Cintybus_Cmodel_anova<- Anova(Cintybus_Cmodel)
write.table(Cintybus_Cmodel_anova, file="Cintybus_Cmodel_anova.txt")

print(paste("starting Cintybus random", Sys.time()))

Cintybus_Cmodel_random <- rand(Cintybus_Cmodel)
Cintybus_Cmodel_random  <- as.data.frame(print(Cintybus_Cmodel_random))
write.table(Cintybus_Cmodel_random, file="Cintybus_Cmodel_random.txt")

print(paste("starting Cintybus residual", Sys.time()))
Cintybus_Cmodel_residual <- residuals(Cintybus_Cmodel)
Cintybus_Cmodel_residual  <- as.data.frame(print(Cintybus_Cmodel_residual))
write.table(Cintybus_Cmodel_residual, file="Cintybus_Cmodel_residual.txt")

print(paste("starting cintybus ls means", Sys.time()))

Cintybus_Cmodel_iso <- lmer(Scale.LS ~ IsolateID + (1|Exp) + (1|Flat) + (1|IndPlant1), data=Cintybus)
Cintybus_Cmodel_iso_anova <- anova(Cintybus_Cmodel_iso)  
write.table(Cintybus_Cmodel_iso_anova, file="Cintybus_Cmodel_iso_anova.txt")
LesionIso_lsm_Cintybus <- lsmeans(Cintybus_Cmodel_iso, "IsolateID")
write.table(LesionIso_lsm_Cintybus, file="LesionIso_lsm_Cintybus.txt")

Cintybus_Cmodel_Geno <- lmer(Scale.LS ~ PlantGeno + (1|Exp) + (1|Flat) + (1|IndPlant1), data=Cintybus) 
Cintybus_Cmodel_Geno_anova <- anova(  Cintybus_Cmodel_Geno)  
write.table(Cintybus_Cmodel_Geno_anova, file="Cintybus_Cmodel_Geno_anova.txt")
LesionGeno_lsm_Cintybus <- lsmeans(Cintybus_Cmodel_Geno, "PlantGeno")
LesionGeno_lsm_Cintybus <- as.data.frame(print(LesionGeno_lsm_Cintybus))
write.table(LesionGeno_lsm_Cintybus, file="LesionGeno_lsm_Cintybus.txt")

Cintybus_Cmodel_improv <- lmer(Scale.LS ~ Improvement + (1|Exp) + (1|Flat) + (1|IndPlant1), data=Cintybus) 
Cintybus_Cmodel_improv_anova <- anova(Cintybus_Cmodel_improv)  
write.table(Cintybus_Cmodel_improv_anova, file="Cintybus_Cmodel_improv_anova.txt")
LesionDom_lsm_Cintybus <- lsmeans(Cintybus_Cmodel_improv, "Improvement")
LesionDom_lsm_Cintybus  <- as.data.frame(print(LesionDom_lsm_Cintybus))
write.table(LesionDom_lsm_Cintybus, file="LesionDom_lsm_Cintybus.txt")

#############_________________________________________ Glycine

rm(list=ls())

Eudicot_clean <- read.table(file="~/Documents/Science_projects/Davis_project/Eudicot project/DataForCeline/7species_data/Eudicot7sp_clean_median.txt", header=T)
Eudicot_clean$IndPlant1 <- as.factor(paste(Eudicot_clean$PlantGeno,Eudicot_clean$Exp,Eudicot_clean$IndPlant, sep="-"))

EudicotC_data_Sp <- split(Eudicot_clean, Eudicot_clean$Taxon)
EudicotC_data_Sp  <- as.array(EudicotC_data_Sp)
list2env(EudicotC_data_Sp , envir=.GlobalEnv)
rm(EudicotC_data_Sp)

Glycine <- droplevels.data.frame(Glycine)

print(paste("starting glycine", Sys.time()))

Glycine_Cmodel <- lmer(Scale.LS ~ IsolateID + Improvement/PlantGeno+ IsolateID*Improvement + IsolateID*Improvement/PlantGeno + (1|Exp) + (1|Flat) + (1|IndPlant1), data=Glycine)

Glycine_Cmodel_anova<- Anova(Glycine_Cmodel)
write.table(Glycine_Cmodel_anova, file="Glycine_Cmodel_anova.txt")

print(paste("starting glycine random", Sys.time()))

Glycine_Cmodel_random <- rand(Glycine_Cmodel)
Glycine_Cmodel_random  <- as.data.frame(print(Glycine_Cmodel_random))
write.table(Glycine_Cmodel_random, file="Glycine_Cmodel_random.txt")

print(paste("starting Glycine residual", Sys.time()))
Glycine_Cmodel_residual <- residuals(Glycine_Cmodel)
Glycine_Cmodel_residual  <- as.data.frame(print(Glycine_Cmodel_residual))
write.table(Glycine_Cmodel_residual, file="Glycine_Cmodel_residual.txt")

print(paste("starting glycine lsmeans", Sys.time()))

Glycine_Cmodel_iso <- lmer(Scale.LS ~ IsolateID + (1|Exp) + (1|Flat) + (1|IndPlant1), data=Glycine)
Glycine_Cmodel_iso_anova <- anova(Glycine_Cmodel_iso)  
write.table(Glycine_Cmodel_iso_anova, file="Glycine_Cmodel_iso_anova.txt")
LesionIso_lsm_Glycine <- lsmeans(Glycine_Cmodel_iso, "IsolateID")
write.table(LesionIso_lsm_Glycine, file="LesionIso_lsm_Glycine.txt")

Glycine_Cmodel_Geno <- lmer(Scale.LS ~ PlantGeno + (1|Exp) + (1|Flat) + (1|IndPlant1), data=Glycine) 
Glycine_Cmodel_Geno_anova <- anova(  Glycine_Cmodel_Geno)  
write.table(Glycine_Cmodel_Geno_anova, file="Glycine_Cmodel_Geno_anova.txt")
LesionGeno_lsm_Glycine <- lsmeans(Glycine_Cmodel_Geno, "PlantGeno")
LesionGeno_lsm_Glycine <- as.data.frame(print(LesionGeno_lsm_Glycine))
write.table(LesionGeno_lsm_Glycine, file="LesionGeno_lsm_Glycine.txt")

Glycine_Cmodel_improv <- lmer(Scale.LS ~ Improvement + (1|Exp) + (1|Flat) + (1|IndPlant1), data=Glycine) 
Glycine_Cmodel_improv_anova <- anova(Glycine_Cmodel_improv)  
write.table(Glycine_Cmodel_improv_anova, file="Glycine_Cmodel_improv_anova.txt")
LesionDom_lsm_Glycine <- lsmeans(Glycine_Cmodel_improv, "Improvement")
LesionDom_lsm_Glycine  <- as.data.frame(print(LesionDom_lsm_Glycine))
write.table(LesionDom_lsm_Glycine, file="LesionDom_lsm_Glycine.txt")


#############________________________________________________________ Solanum

rm(list=ls())

Eudicot_clean <- read.table(file="~/Documents/Science_projects/Davis_project/Eudicot project/DataForCeline/7species_data/Eudicot7sp_clean_median.txt", header=T)
Eudicot_clean$IndPlant1 <- as.factor(paste(Eudicot_clean$PlantGeno,Eudicot_clean$Exp,Eudicot_clean$IndPlant, sep="-"))

EudicotC_data_Sp <- split(Eudicot_clean, Eudicot_clean$Taxon)
EudicotC_data_Sp  <- as.array(EudicotC_data_Sp)
list2env(EudicotC_data_Sp , envir=.GlobalEnv)
rm(EudicotC_data_Sp)

Solanum <- droplevels.data.frame(Solanum)

print(paste("starting Solanum", Sys.time()))

Solanum_Cmodel <- lmer(Scale.LS ~ IsolateID + Improvement/PlantGeno+ IsolateID*Improvement + IsolateID*Improvement/PlantGeno + (1|Exp) + (1|Flat)+ (1|IndPlant1), data=Solanum)

Solanum_Cmodel_anova<- Anova(Solanum_Cmodel)
write.table(Solanum_Cmodel_anova, file="Solanum_Cmodel_anova.txt")

print(paste("starting solanum random", Sys.time()))
Solanum_Cmodel_random <- rand(Solanum_Cmodel)
Solanum_Cmodel_random  <- as.data.frame(print(Solanum_Cmodel_random))
write.table(Solanum_Cmodel_random, file="Solanum_Cmodel_random.txt")

print(paste("starting Solanum residual", Sys.time()))
Solanum_Cmodel_residual <- residuals(Solanum_Cmodel)
Solanum_Cmodel_residual  <- as.data.frame(print(Solanum_Cmodel_residual))
write.table(Solanum_Cmodel_residual, file="Solanum_Cmodel_residual.txt")


print(paste("starting solanum lsmean", Sys.time()))

Solanum_Cmodel_iso <- lmer(Scale.LS ~ IsolateID + (1|Exp) + (1|Flat) + (1|IndPlant1), data=Solanum)
Solanum_Cmodel_iso_anova <- anova(Solanum_Cmodel_iso)  
write.table(Solanum_Cmodel_iso_anova, file="Solanum_Cmodel_iso_anova.txt")
LesionIso_lsm_Solanum <- lsmeans(Solanum_Cmodel_iso, "IsolateID")
write.table(LesionIso_lsm_Solanum, file="LesionIso_lsm_Solanum.txt")

Solanum_Cmodel_Geno <- lmer(Scale.LS ~ PlantGeno + (1|Exp) + (1|Flat) + (1|IndPlant1), data=Solanum) 
Solanum_Cmodel_Geno_anova <- anova(  Solanum_Cmodel_Geno)  
write.table(Solanum_Cmodel_Geno_anova, file="Solanum_Cmodel_Geno_anova.txt")
LesionGeno_lsm_Solanum <- lsmeans(Solanum_Cmodel_Geno, "PlantGeno")
LesionGeno_lsm_Solanum <- as.data.frame(print(LesionGeno_lsm_Solanum))
write.table(LesionGeno_lsm_Solanum, file="LesionGeno_lsm_Solanum.txt")

Solanum_Cmodel_improv <- lmer(Scale.LS ~ Improvement + (1|Exp) + (1|Flat) + (1|IndPlant1), data=Solanum) 
Solanum_Cmodel_improv_anova <- anova(Solanum_Cmodel_improv)  
write.table(Solanum_Cmodel_improv_anova, file="Solanum_Cmodel_improv_anova.txt")
LesionDom_lsm_Solanum <- lsmeans(Solanum_Cmodel_improv, "Improvement")
LesionDom_lsm_Solanum  <- as.data.frame(print(LesionDom_lsm_Solanum))
write.table(LesionDom_lsm_Solanum, file="LesionDom_lsm_Solanum.txt")
