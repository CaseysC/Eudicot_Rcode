#####################
## Celine Caseys
## Plant Science dept. UC davis
## January 2019
###################


library(MASS)
library(lme4)
library(lmerTest)
library(car)
library(emmeans)

Eudicot_clean <- read.table(file="Eudicot7sp_clean_median.txt", header=T)
Eudicot_clean$IndPlant1 <- as.factor(paste(Eudicot_clean$PlantGeno,Eudicot_clean$Exp,Eudicot_clean$IndPlant, sep="-"))
EudicotC_data_PGeno <- split(Eudicot_clean, Eudicot_clean$PlantGeno)
EudicotC_data_PGeno <-  as.array(EudicotC_data_PGeno)

Geno_list <- Eudicot_clean$PlantGeno

Colnames1 <-  colnames(Eudicot_clean)
Geno <- names(EudicotC_data_PGeno)

Genolist <- as.data.frame(dimnames(EudicotC_data_PGeno))
colnames(Genolist)[1] <- c("PlantGeno")
Genolist$NumID <- rownames(Genolist)

Genolist <- merge(Genolist, Accessions.7sp, by='PlantGeno')

## Updated to include only Individual plants

for (f in c(1:length(EudicotC_data_PGeno))) {
  temp <- droplevels.data.frame(data.frame(EudicotC_data_PGeno[f]))
  colnames(temp) <- Colnames1
  D1 <- lmer(Scale.LS ~ IsolateID + (1|IndPlant) + (1|Exp) + (1|Flat), data=temp)
  Anova_temp <- Anova(D1)
  #write.table(Anova_temp, file=paste(Geno[f], "_anova.txt", sep=""))
  #write.table(Rand_temp, file=paste(Geno[f], "_rand.txt", sep=""))
  D1A_Isolate_lsm <- emmeans(D1, ~IsolateID, lmer.df='satterthwaite')
  lsm_temp <- as.data.frame(print(D1A_Isolate_lsm))
  lsm_temp$Geno <- rep(Geno[f], dim(lsm_temp)[1])
  write.table(lsm_temp, file=paste(Geno[f], "_lsm.txt", sep=""))
  rm(D1,temp,Anova_temp, Rand_temp, D1A_Isolate_lsm, lsm_temp)
  print(f)
}

################################ Compiling the Least square Means

#### In rows format
f <-1
LSmean_Data_Geno <- as.data.frame(matrix(ncol=7))
temp1 <- read.table(file=paste(Geno[f], "_lsm.txt", sep=""), header=T)
temp2 <- read.table(file=paste(Geno[f+1], "_lsm.txt", sep=""), header=T)
LSmean_Data_Geno <-rbind(temp1, temp2)

for (f in c(3:length(EudicotC_data_PGeno))) {
  temp1 <- read.table(file=paste(Geno[f], "_lsm.txt", sep=""), header=T)
  LSmean_Data_Geno <-rbind(LSmean_Data_Geno, temp1)
}

#write.table(LSmean_Data_Geno, file="LSmeanIndP_7spData_raw_Geno.txt")

#### In columns format
f <-1
temp1 <- read.table(file=paste(Geno[f], "_lsm.txt", sep=""), header=T)
temp1.1 <- temp1[,1:2]
colnames(temp1.1) <- c("IsolateID", levels(temp1$Geno))

temp2 <- read.table(file=paste(Geno[f+1], "_lsm.txt", sep=""), header=T)
temp2.1 <- temp2[,1:2]
colnames(temp2.1) <- c("IsolateID", levels(temp2$Geno))

LSmean_Data_col_Geno <-merge(temp2.1, temp1.1, by="IsolateID", all=T)

for (f in c(3:length(EudicotC_data_PGeno))) {
  temp1 <- read.table(file=paste(Geno[f], "_lsm.txt", sep=""), header=T)
  temp1.1 <- temp1[,1:2]
  colnames(temp1.1) <- c("IsolateID", levels(temp1$Geno))
  
  LSmean_Data_col_Geno <-merge(LSmean_Data_col_Geno, temp1.1, by="IsolateID", all=T)
}

#write.table(LSmean_Data_col_Geno, file="LSmeanIndP_7spData_col_Geno.txt")

rownames(LSmean_Data_col_Geno) <- LSmean_Data_col_Geno[,1]
LSmean_Data_col_GenoT <- as.data.frame(t(LSmean_Data_col_Geno))
LSmean_Data_col_GenoT$PlantGeno <- rownames(LSmean_Data_col_GenoT)
LSmean_Data_col_GenoT <-  LSmean_Data_col_GenoT[-1,]
LSmean_Data_col_GenoT <-  merge(LSmean_Data_col_GenoT, Accessions.7sp, by='PlantGeno')

#write.table(LSmean_Data_col_GenoT, file="LSmeanIndP_7spData_col_Geno_All.txt")

