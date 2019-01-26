#####################
## Celine Caseys
## Plant Science dept. UC davis
## January 2019
###################

library(tidyverse)

Eudicot8sp_lsmeans <- read.table(file="Eudi_Ara_Lsmeans1_8sp_Col_Geno_All.txt", header=T)
Eudicot_info <- Eudicot8sp_lsmeans[, c(99:105)]

Eudicot_info$Domest <- Eudicot_info$Wi_Do
levels(Eudicot_info$Domest) <-  c("High", "Low", "Other", "Low")

Eudicot_info$clade <- Eudicot_info$Order
levels(Eudicot_info$clade) <- c("Asterids", "Rosids", "Rosids", "Asterids")

Eudicot8sp_lsmeans <- arrange(Eudicot8sp_lsmeans, PGeno)
rownames(Eudicot8sp_lsmeans) <-  Eudicot8sp_lsmeans$PGeno
Eudicot7sp_heatmap <- Eudicot8sp_lsmeans[, -c(99:105)]

Eudicot_info <- arrange(Eudicot_info, PGeno)

Isolate_HO <- read.table(file="Isolate_Host_origin_CORRECTED.txt", header=T)
Isolate_HO <- arrange(Isolate_HO, Isolate)
Host <- as.factor(Isolate_HO$Host)
Origin <- as.factor(Isolate_HO$Origin)

# to have the species as column
Eudicot7sp_heatmapt1 <- as.data.frame(t(Eudicot7sp_heatmap))
Eudicot7sp_heatmapt <- Eudicot7sp_heatmapt1[complete.cases(Eudicot7sp_heatmapt1), ]

## Z-scoring the data
Eudicot7sp_heatmapStand <-  as.data.frame(matrix(ncol=90, nrow=91))
Eudicot7sp_heatmapStand1 <-  as.data.frame(matrix(ncol=90, nrow=98))

for ( i in c(1:90)){
Eudicot7sp_heatmapStand[,i] <- scale(Eudicot7sp_heatmapt[, i], center=T, scale=T)  
}
colnames(Eudicot7sp_heatmapStand) <-  colnames(Eudicot7sp_heatmapt)
rownames(Eudicot7sp_heatmapStand) <-  rownames(Eudicot7sp_heatmapt)

for ( i in c(1:90)){
  Eudicot7sp_heatmapStand1[,i] <- scale(Eudicot7sp_heatmapt1[, i], center=T, scale=T)  
}
colnames(Eudicot7sp_heatmapStand1) <-  colnames(Eudicot7sp_heatmapt1)
rownames(Eudicot7sp_heatmapStand1) <-  rownames(Eudicot7sp_heatmapt1)

