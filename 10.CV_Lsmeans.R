#####################
## Celine Caseys
## Plant Science dept. UC davis
## January 2019
###################


source('~/Desktop/Github_Eudicot/6.ZcoreLSmeans.R')

Virulence_general <- rowMeans(Eudicot7sp_heatmapStand1, na.rm=T)
Virulence_general <- as.data.frame(Virulence_general)
Virulence_general$IsolateID <- rownames(Virulence_general)
Virulence_general$IsolateID <-  as.factor(gsub('X0', '0', Virulence_general$IsolateID, fixed=F))
Virulence_general$IsolateID <-  as.factor(gsub('X2004', '2004', Virulence_general$IsolateID, fixed=F))
Virulence_general$IsolateID <-  as.factor(gsub('X94', '94', Virulence_general$IsolateID, fixed=F))
Virulence_general$Eudicot <-  as.factor(rep(c("Eudicot"), 98))

source('~/Desktop/Github_Eudicot/1.Creating_Dataset.R')
library(car)
library(tidyverse)
Eudicot7sp_lsmeans <- read.table(file="Eudi_Ara_Lsmeans1_8sp_row_Geno_All_plot.txt", header=T)
Eudicot7sp_lsmeans$IsoTax <- as.factor(paste(Eudicot7sp_lsmeans$Isolate, Eudicot7sp_lsmeans$Taxon, sep="-"))


# Summary statistics at the isolate x taxon level
LSm_Summary_IsoTax <- summarise(group_by(Eudicot7sp_lsmeans, IsoTax), 
                            mean.isotax=mean(LsMeans, na.rm=T),
                            median.isotax=median(LsMeans), 
                            sd.isotax=sd(LsMeans), 
                            min.isotax=min(LsMeans), 
                            max.isotax=max(LsMeans), 
                            length.isotax=length(LsMeans), 
                            cv.isotax=sd.isotax/mean.isotax,
                            SE.mean.isotax=sd.isotax/sqrt(length.isotax))

LSm_Summary_IsoTax$IsoTax1 <- LSm_Summary_IsoTax$IsoTax
LSm_Summary_IsoTax = separate(data=LSm_Summary_IsoTax, col=IsoTax1, into=c("IsolateID", "Taxon"), sep="\\-")
LSm_Summary_IsoTax$IsolateID <-  as.factor(LSm_Summary_IsoTax$IsolateID)
LSm_Summary_IsoTax$Taxon <-  as.factor(LSm_Summary_IsoTax$Taxon)

LSm_Summary_IsoTax <-  LSm_Summary_IsoTax[complete.cases(LSm_Summary_IsoTax), ]
  
### Calculating the host specificity 
## as the coefficient of variation of means at the isotax level

CV_btwSp <- summarise(group_by(LSm_Summary_IsoTax, IsolateID), 
                      mean.Iso=mean(mean.isotax, na.rm=T), 
                      median.Iso=median(mean.isotax), 
                      sd.Iso=sd(mean.isotax), 
                      min.Iso=min(mean.isotax), 
                      max.Iso=max(mean.isotax), 
                      Present_in=length(mean.isotax),                                      
                      cv_betweenSp=sd.Iso/mean.Iso)


CV_btwSp1 <- arrange(CV_btwSp, cv_betweenSp)
CV_btwSp1 <- CV_btwSp1[, -c(3,5,6)]


SpecVir <- merge(Virulence_general, CV_btwSp, by="IsolateID", all=T)
SpecVir <- merge(SpecVir, B_host, by="IsolateID")


############## Modeling

# linear model
model1 <-  lm(cv_betweenSp ~ Virulence_general, data=SpecVir)
Anova(model1)
summary(model1)
SpecVir$linearResiduals <- model1$residuals
err1 <- predict(model1, se.fit = TRUE)

SpecVir$Lin_lci <- err1$fit - 1.96 * err1$se.fit
SpecVir$Lin_fit <- err1$fit
SpecVir$Lin_uci <- err1$fit + 1.96 * err1$se.fit

ggplot(SpecVir, aes(x = Virulence_general, y = Lin_fit)) +
  theme_bw() +
  geom_line() +
  geom_smooth(aes(ymin = Lin_lci, ymax = Lin_uci), stat = "identity") +
  geom_point(data = SpecVir, aes(x = Virulence_general, y = cv_betweenSp))

# Quadratic model
SpecVir$Vir2 <- (SpecVir$Virulence_general)^2
model_quad <-  lm(cv_betweenSp ~ Virulence_general + Vir2, data=SpecVir)
summary(model_quad)
Anova(model_quad)
SpecVir$QuadResiduals <- model_quad$residuals
err_quad <- predict(model_quad, se.fit = TRUE)

SpecVir$Quad_lci <- err_quad$fit - 1.96 * err_quad$se.fit
SpecVir$Quad_fit <- err_quad$fit
SpecVir$Quad_uci <- err_quad$fit + 1.96 * err_quad$se.fit

ggplot(SpecVir, aes(x = Virulence_general, y = Quad_fit)) +
  theme_bw() +
  geom_line() +
  geom_smooth(aes(ymin = Quad_lci, ymax = Quad_uci), stat = "identity") +
  geom_point(data = SpecVir, aes(x = Virulence_general, y = cv_betweenSp))



### Fig. 4 based on the quadratic model

Host_color1 <- c("#8dd3c7", "#b2df8a", "#984ea3","#ccebc5","#ffed6f","#e7298a", "#c994c7","#fe9929", "#4eb3d3", "#f768a1", "#fcbba1", "#e31a1c")

Spe_vir <-  ggplot(SpecVir, aes(Virulence_general, cv_betweenSp)) +
  geom_point(aes(color=Host2), size=2.5)+
  scale_color_manual(values=Host_color1)+
  theme_bw() +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), se=TRUE, colour=I("black"))+
  xlab("Virulence")+
  ylab("Specificity")
# theme(legend.position="none")


# Isolates highlight: B05.10, UKrazz, DavisNavel

B05 <- LSm_Summary_IsoTax[which(LSm_Summary_IsoTax$IsolateID=='B05_10'),]
Dav <- LSm_Summary_IsoTax[which(LSm_Summary_IsoTax$IsolateID=='DavisNavel'),]
KT <- LSm_Summary_IsoTax[which(LSm_Summary_IsoTax$IsolateID=='KatieTomato'),]
UKR <- LSm_Summary_IsoTax[which(LSm_Summary_IsoTax$IsolateID=='UKRazz'),]

Iso3Plot <- rbind(B05, Dav, KT, UKR)

CV_btwSp$Eudicot <-  as.factor(rep(c("Eudicot"), 98))

B05_plot <- ggplot(B05, aes(Taxon, mean.isotax)) +
  geom_bar(stat="identity", fill="skyblue")+ 
  geom_errorbar(aes(ymin=mean.isotax-SE.mean.isotax, ymax=mean.isotax+SE.mean.isotax), width=.2,
                position=position_dodge(.9)) +
  ggtitle("B05.10")+
  ylab("Lesion Area [cm2]")+
  theme_bw()

Dav_plot <- ggplot(Dav, aes(Taxon, mean.isotax)) +
  geom_bar(stat="identity", fill="Navy") +
  geom_errorbar(aes(ymin=mean.isotax-SE.mean.isotax, ymax=mean.isotax+SE.mean.isotax), width=.2,
                position=position_dodge(.9)) +
  ggtitle("Davis Navel")+
  ylab("Lesion Area [cm2]")+
  theme_bw()

UKR_plot <- ggplot(UKR, aes(Taxon, mean.isotax)) +
  geom_bar(stat="identity", fill="steelblue")+
  geom_errorbar(aes(ymin=mean.isotax-SE.mean.isotax, ymax=mean.isotax+SE.mean.isotax), width=.2,
                position=position_dodge(.9)) +
  ggtitle("UKRazz") +
  ylab("Lesion Area [cm2]")+
  theme_bw()


KT_Plot <- ggplot(KT, aes(Taxon, mean.isotax)) +
  geom_bar(stat="identity", fill="blue")+
  geom_errorbar(aes(ymin=mean.isotax-SE.mean.isotax, ymax=mean.isotax+SE.mean.isotax), width=.2,
                position=position_dodge(.9)) +
  ggtitle("Katie Tomato") +
  ylab("Lesion Area [cm2]")+
  theme_bw()
  
