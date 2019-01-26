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
Eudicot7sp_heatmap <- as.data.frame(t(Eudicot7sp_heatmap))

Eudicot7sp_heatmap <- Eudicot7sp_heatmap[complete.cases(Eudicot7sp_heatmap), ]


Eudicot7sp_heatmapStand <-  as.data.frame(matrix(ncol=90, nrow=91))

for ( i in c(1:90)){
Eudicot7sp_heatmapStand[,i] <- scale(Eudicot7sp_heatmap[, i], center=T, scale=T)  
}
colnames(Eudicot7sp_heatmapStand) <-  colnames(Eudicot7sp_heatmap)
rownames(Eudicot7sp_heatmapStand) <-  rownames(Eudicot7sp_heatmap)

Eudicot7sp_heatmapSMean <-  as.data.frame(matrix(ncol=90, nrow=91))

for ( i in c(1:90)){
  Eudicot7sp_heatmapSMean[,i] <- scale(Eudicot7sp_heatmap[, i], center=T, scale=F)  
}
colnames(Eudicot7sp_heatmapSMean) <-  colnames(Eudicot7sp_heatmap)
rownames(Eudicot7sp_heatmapSMean) <-  rownames(Eudicot7sp_heatmap)

Eudicot7sp_heatmapS <- as.matrix(t(Eudicot7sp_heatmapStand))

library(iheatmapr)

library(RColorBrewer)
jBrewColors <- brewer.pal(n = 8, name = "Dark2")
jBrewColors1 <- brewer.pal(n = 6, name = "Pastel1")

Host_color <- c("#8dd3c7", "#b2df8a", "#984ea3","#ccebc5","#ffed6f","#e7298a", "#f1eef6","#fe9929", "#4eb3d3", "#f768a1", "#fcbba1", "#e31a1c")
Origin_color <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f")
jBrewColor <- c(jBrewColors,jBrewColors1)

Sp_color <-  c("#c7eae5","#4575b4", "#feb24c","#fd8d3c", "#74add1","#fee090", "#b8e186" ,"#d73027")
Domest_color <- c("Black","Grey", "white")

Order_color <- c("#ffffbf", "#d8daeb", "#b2abd2","#fdb863")
#Order_color <- c("#E0FFFF", "#CD5555", "#F08080","#00BFFF")
Clade_color <- c("#b35806", "#542788")
#Clade_color <- c("#3A5FCD", "#FF6347")




Eudicot_heatmap_stand <- main_heatmap(Eudicot7sp_heatmapS, name = "Standardized lesion size") %>%
  add_col_clustering() %>%
  add_row_clustering() %>%
  #add_col_clustering() %>%
  #add_col_labels() %>%
  #add_row_labels() %>%
  add_row_title("Plant Genotypes") %>%
  add_col_title("Botrytis Isolates") %>%
  add_col_annotation(data.frame("Host" = Isolate_HO$Host2,
                                "Origin"= Isolate_HO$Origin), side="top",
                     color=list("Host"=Host_color,
                                "Origin"=Origin_color))  %>%
  add_row_annotation(data.frame("Improv."=Eudicot_info$Domest,
                                "Species" = Eudicot_info$Taxon,
                                "Order"= Eudicot_info$Order,
                                "Clade"= Eudicot_info$clade), side="left",
                     color=list("Improv."=Domest_color,
                                "Species"=Sp_color,
                                "Order"=Order_color,
                                "Clade"=Clade_color))

Eudicot_heatmap_stand
Eudicot_heatmap_stand %>% save_iheatmap("Eudicot8sp_heatmap_stand.pdf", vwidth=2300, vheight=1700) 

library(MASS)
library(pvclust)
rownames(Eudicot7sp_heatmapS) <- Eudicot_info$SNames
PGenost_iso_pv <- pvclust(Eudicot7sp_heatmapS, nboot=200, method.hclust = "centroid")
plot(PGenost_iso_pv, main="Iso HClust 'complete'")
pvrect(PGenost_iso_pv, alpha = .95)

tNoAdmixt_heatmapS <-  t(Eudicot7sp_heatmapS)
Noadmixt_genoo_pv <- pvclust(tNoAdmixt_heatmapS, nboot=200, method.hclust = "centroid")
plot(Noadmixt_genoo_pv, main="Plant Genotypes Hierarchical Clustering method='complete'")
pvrect(Noadmixt_genoo_pv, alpha = .95)



Eudicot7sp_heatmap <- as.matrix(t(Eudicot7sp_heatmap))

Eudicot_heatmap <- main_heatmap(Eudicot7sp_heatmap, name = "Standardized lesion size") %>%
  add_col_clustering() %>%
  add_row_clustering() %>%
  add_col_labels()%>%
  add_row_labels()%>%
  add_row_title("Genotypes") %>%
  add_col_title("Isolates") %>%
  add_row_annotation(data.frame("Order"= Eudicot_info$Order,
                                "Taxon" = Eudicot_info$Taxon), side="left")

htmlwidgets::saveWidget(as_widget(Eudicot_heatmap), "Eudicotheatmapimp.html")

############################################## Plotting Violin plot from LSmeans


Eudi8sp <- read.table(file="Eudi_Ara_Lsmeans1_8sp_row_Geno_All_plot.txt", header=T)

Eudi8sp$Taxon1 <- Eudi8sp$Taxon
levels(Eudi8sp$Taxon1) <- c("RosB_Ath", "RosB_Bra", "AstA_Cen", "AstA_Cin", "RosF_Glyc", "AstA_Han", "AstA_Lact", "AstS_Sol" )
Eudi8sp$clade <- Eudi8sp$Order
levels(Eudi8sp$clade) <- c("Asterids", "Rosids", "Rosids", "Asterids")

Eudi8sp$IsoTax <- paste(Eudi8sp$Isolate, Eudi8sp$Taxon, sep="-")

Eudi8sp$PlotOrder <- factor(Eudi8sp$Taxon)
levels(Eudi8sp$PlotOrder) <- c(6, 7, 4, 5,8,2, 3,1)
Eudi8sp$PlotOrder <- as.numeric(Eudi8sp$PlotOrder)

Eudi8sp$Taxon <- factor(Eudi8sp$Taxon, levels = Eudi8sp$Taxon[order(Eudi8sp$PlotOrder)])


ggplot(Eudi8sp, aes(Taxon, LsMeans, fill=Wi_Do)) +
  geom_violin()

#### Geom Split violin from https://stackoverflow.com/questions/35717353/split-violin-plot-with-ggplot2
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, draw_group = function(self, data, ..., draw_quantiles = NULL){
  data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
  grp <- data[1,'group']
  newdata <- plyr::arrange(transform(data, x = if(grp%%2==1) xminv else xmaxv), if(grp%%2==1) y else -y)
  newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
  newdata[c(1,nrow(newdata)-1,nrow(newdata)), 'x'] <- round(newdata[1, 'x']) 
  if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
    stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 
                                              1))
    quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
    aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
    aesthetics$alpha <- rep(1, nrow(quantiles))
    both <- cbind(quantiles, aesthetics)
    quantile_grob <- GeomPath$draw_panel(both, ...)
    ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
  }
  else {
    ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
  }
})


geom_split_violin <- function (mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, position = position, show.legend = show.legend, inherit.aes = inherit.aes, params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

ggplot(Eudi8sp, aes(Taxon, LsMeans, fill=Wi_Do)) +
  geom_split_violin() 

levels(Eudi8sp$Wi_Do) <- c("coi1", "Col-0", "High",  "Low" ,"Low")


ggplot(Eudi8sp, aes(Taxon1, LsMeans, fill=Wi_Do)) +
  geom_split_violin(scale="width") +
  geom_boxplot(width=0.15)+ 
  scale_fill_manual(values=c("#80cdc1", "#c7eae5", "#000000", "#878787"))+
  theme_minimal()




### Getting the summaries for writing the results. 

Summary_species_Eudi8 <- summarise(group_by(Eudi8sp, Taxon), 
                             mean.sp=mean(LsMeans, na.rm = TRUE), 
                             median.sp=median(LsMeans, na.rm = TRUE), 
                             sd.sp=sd(LsMeans, na.rm = TRUE), 
                             min.sp=min(LsMeans, na.rm = TRUE), 
                             max.sp=max(LsMeans, na.rm = TRUE), 
                             length.sp=length(LsMeans), 
                             cv.sp=sd.sp/mean.sp)

Summary_IsoTax_Eudi8 <- summarise(group_by(Eudi8sp, IsoTax), 
                                   mean.sp=mean(LsMeans, na.rm = TRUE), 
                                   median.sp=median(LsMeans, na.rm = TRUE), 
                                   sd.sp=sd(LsMeans, na.rm = TRUE), 
                                   min.sp=min(LsMeans, na.rm = TRUE), 
                                   max.sp=max(LsMeans, na.rm = TRUE), 
                                   length.sp=length(LsMeans), 
                                   cv.sp=sd.sp/mean.sp)

Summary_IsoTax_Eudi8$IsoTax1 <- Summary_IsoTax_Eudi8$IsoTax
Summary_IsoTax_Eudi8 <- separate(Summary_IsoTax_Eudi8,IsoTax1, into=c("Isolate", "Taxon"), sep="-" )
Summary_IsoTax_Eudi8$Isolate <-  as.factor(Summary_IsoTax_Eudi8$Isolate)
Summary_IsoTax_Eudi8$Taxon <-  as.factor(Summary_IsoTax_Eudi8$Taxon)

#write.table(Summary_IsoTax_Eudi8, file="SummaryLSmean_IsoTax_Eudi8.txt")

Summary_IsoTax_Eudi8_splt <- split(Summary_IsoTax_Eudi8, Summary_IsoTax_Eudi8$Taxon) 
TaxNames <- names(Summary_IsoTax_Eudi8_splt)

temp1 <- Summary_IsoTax_Eudi8_splt[[1]][, c(8,9)]
colnames(temp1) <- c(paste("CV", TaxNames[1], sep=""), "IsolateID")

temp2 <- Summary_IsoTax_Eudi8_splt[[2]][, c(8,9)]
colnames(temp2) <- c(paste("CV", TaxNames[2], sep=""), "IsolateID")

TaxCVcol_8sp <- merge(temp1, temp2, by="IsolateID", all=T)

for (f in c(3:length(Summary_IsoTax_Eudi8_splt))) {
  temp1 <- Summary_IsoTax_Eudi8_splt[[f]][, c(8,9)]
  colnames(temp1) <- c(paste("CV", TaxNames[f], sep=""), "IsolateID")
  
  TaxCVcol_8sp <-merge(TaxCVcol_8sp, temp1, by="IsolateID", all=T)
}

rm(Summary_IsoTax_Eudi8_splt, temp1, temp2)

#SpeciesLevel_Mapping <- merge(TaxSum_col_8sp, TaxSD_col_8sp, by="IsolateID")
#SpeciesLevel_Mapping <- merge(SpeciesLevel_Mapping, TaxCVcol_8sp, by="IsolateID")

#write.table( SpeciesLevel_Mapping, file="SpeciesLevel_Mapping_stat.txt")
SpeciesLevel_Mapping <- read.table(file="SpeciesLevel_Mapping_stat.txt", header=T)


Summary_TaxWido <-  summarise(group_by(Eudi8sp, Tax_Wido), 
                             mean=mean(LsMeans, na.rm = TRUE), 
                             median=median(LsMeans, na.rm = TRUE), 
                             sd=sd(LsMeans, na.rm = TRUE), 
                             min=min(LsMeans, na.rm = TRUE), 
                             max=max(LsMeans, na.rm = TRUE), 
                             length=length(LsMeans), 
                             cv=sd/mean)

#write.table(Summary_TaxWido, file="LSmean_Summary_SpeciesImprov.txt")

#######################################################################

Virulence_general <- rowMeans(Eudicot7sp_heatmapStand, na.rm=T)
Virulence_general_1 <- rowMeans(Eudicot7sp_heatmap, na.rm=T)

Virulence_general <- as.data.frame(Virulence_general_1)
Virulence_general$IsolateID <- rownames(Virulence_general)
Virulence_general$FromZero <- Virulence_general$Virulence_general + 1.7010967927
Virulence_general$FromZero <-  as.numeric(Virulence_general$FromZero)
Virulence_general$IsolateID <-  as.factor(gsub('X0', '0', Virulence_general$IsolateID, fixed=F))

Virulence_mean <- rowMeans(Eudicot7sp_heatmapSMean, na.rm=T)
Virulence_mean <- as.data.frame(Virulence_mean)
Virulence_mean$IsolateID <- rownames(Virulence_general)
Virulence_mean$IsolateID <-  as.factor(gsub('X0', '0', Virulence_mean$IsolateID, fixed=F))

Virulence_general$Eudicot <-  as.factor(rep(c("Eudicot"), 98))

ggplot(Virulence_general, aes(FromZero)) +
  geom_density()+
  theme_bw()

ggplot(Virulence_general, aes(Eudicot, FromZero)) +
  geom_boxplot()+
  theme_bw() +
  coord_flip()

source('~/Documents/Science_projects/Davis_project/Eudicot project/DataForCeline/FinalRcodes/9.CV_Lsmeans.R')

Virulence_general$IsolateID <-  as.factor(gsub('X0', '0', Virulence_general$IsolateID, fixed=F))

SpecVir <- merge(Virulence_general, Virulence_mean)
SpecVir <- merge(SpecVir, CV_btwSp)
SpecVir <- merge(SpecVir, B_host)

library(car)


model1 <-  lm(cv_betweenSp ~ FromZero, data=SpecVir)
Anova(model1)
summary(model1)
SpecVir$linearResiduals <- model1$residuals
err1 <- predict(model1, se.fit = TRUE)

SpecVir$Lin_lci <- err1$fit - 1.96 * err1$se.fit
SpecVir$Lin_fit <- err1$fit
SpecVir$Lin_uci <- err1$fit + 1.96 * err1$se.fit

ggplot(SpecVir, aes(x = FromZero, y = Lin_fit)) +
  theme_bw() +
  geom_line() +
  geom_smooth(aes(ymin = Lin_lci, ymax = Lin_uci), stat = "identity") +
  geom_point(data = SpecVir, aes(x = FromZero, y = cv_betweenSp))


SpecVir$Vir2 <- (SpecVir$FromZero)^2
model_quad <-  lm(cv_betweenSp ~ FromZero + Vir2, data=SpecVir)
summary(model_quad)
Anova(model_quad)
SpecVir$QuadResiduals <- model_quad$residuals
err_quad <- predict(model_quad, se.fit = TRUE)

SpecVir$Quad_lci <- err_quad$fit - 1.96 * err_quad$se.fit
SpecVir$Quad_fit <- err_quad$fit
SpecVir$Quad_uci <- err_quad$fit + 1.96 * err_quad$se.fit

ggplot(SpecVir, aes(x = FromZero, y = Quad_fit)) +
  theme_bw() +
  geom_line() +
  geom_smooth(aes(ymin = Quad_lci, ymax = Quad_uci), stat = "identity") +
  geom_point(data = SpecVir, aes(x = FromZero, y = cv_betweenSp))


write.table(SpecVir, file="Specializ_virulence_data.txt")


ggplot(SpecVir, aes(FromZero, cv_betweenSp, color=Host_color)) +
  geom_point()+
  theme_bw() +
  #geom_smooth(method=lm, se=TRUE)+
  xlab("Virulence")+
  ylab("Specificity")

ggplot(SpecVir, aes(FromZero, cv_betweenSp, color=Admixed5)) +
  geom_point()+
  theme_bw() +
  #geom_smooth(method=lm, se=TRUE)+
  xlab("Virulence")+
  ylab("Specificity")


########## Publication plot Here
jBrewColor

library(gridExtra)
library(ggpubr)

Host_color1 <- c("#8dd3c7", "#b2df8a", "#984ea3","#ccebc5","#ffed6f","#e7298a", "#c994c7","#fe9929", "#4eb3d3", "#f768a1", "#fcbba1", "#e31a1c")

Spe_vir <-  ggplot(SpecVir, aes(Virulence_general, cv_betweenSp)) +
  geom_point(aes(color=Host2), size=2.5)+
  scale_color_manual(values=Host_color1)+
  theme_bw() +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), se=TRUE, colour=I("black"))+
  xlab("Virulence")+
  ylab("Specificity")
  # theme(legend.position="none")
 

Vio_vir <-  ggplot(SpecVir, aes(Eudicot, Virulence_general)) +
  geom_boxplot(width=0.15)+ 
  theme_bw() +
  coord_flip()+
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

Vio_spe <- ggplot(SpecVir, aes(Eudicot, cv_betweenSp)) +
  geom_boxplot(aes(Eudicot, cv_betweenSp), width=0.15)+ 
  theme_bw() +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())


pdf(file="SpeVir_boxPlot1.pdf")
ggarrange(Vio_spe, Spe_vir, NULL, Vio_vir, align = "v" )
dev.off()

pdf(file="SpeVir_legend.pdf")
print(ggplot(SpecVir, aes(Virulence_general, cv_betweenSp)) +
  geom_point(aes(color=Host2), size=2)+
  scale_color_manual(values=Host_color1)+
  theme_bw() +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), se=TRUE, colour=I("black"))+
  xlab("Virulence")+
  ylab("Specificity"))
dev.off()



##################








ggplot(SpecVir, aes(Virulence_mean, cv_betweenSp)) +
  geom_point()+
  theme_bw() +
  geom_smooth(method=lm, se=TRUE)+
  xlab("Virulence")+
  ylab("Specificity")+
  ggtitle("mean centered")


cor.test(SpecVir$Virulence_general, SpecVir$cv_betweenSp)

####################

require(vioplot)
require(devtools)
require(digest)

vioplot2 <- function (x, ..., range = 1.5, h = NULL, ylim = NULL, names = NULL, 
                      horizontal = FALSE, col = "magenta", border = "black", lty = 1, 
                      lwd = 1, rectCol = "black", colMed = "white", pchMed = 19, 
                      at, add = FALSE, wex = 1, drawRect = TRUE, side="both") 
{
  datas <- list(x, ...)
  n <- length(datas)
  if (missing(at)) 
    at <- 1:n
  upper <- vector(mode = "numeric", length = n)
  lower <- vector(mode = "numeric", length = n)
  q1 <- vector(mode = "numeric", length = n)
  q2 <- vector(mode = "numeric", length = n)
  q3 <- vector(mode = "numeric", length = n)
  med <- vector(mode = "numeric", length = n)
  base <- vector(mode = "list", length = n)
  height <- vector(mode = "list", length = n)
  baserange <- c(Inf, -Inf)
  args <- list(display = "none")
  radj <- ifelse(side == "right", 0, 1)
  ladj <- ifelse(side == "left", 0, 1)
  if (!(is.null(h))) 
    args <- c(args, h = h)
  med.dens <- rep(NA, n)
  for (i in 1:n) {
    data <- datas[[i]]
    data.min <- min(data)
    data.max <- max(data)
    q1[i] <- quantile(data, 0.25)
    q2[i] <- quantile(data, 0.5)
    q3[i] <- quantile(data, 0.75)
    med[i] <- median(data)
    iqd <- q3[i] - q1[i]
    upper[i] <- min(q3[i] + range * iqd, data.max)
    lower[i] <- max(q1[i] - range * iqd, data.min)
    est.xlim <- c(min(lower[i], data.min), max(upper[i], 
                                               data.max))
    smout <- do.call("sm.density", c(list(data, xlim = est.xlim), 
                                     args))
    med.dat <- do.call("sm.density", 
                       c(list(data, xlim=est.xlim,
                              eval.points=med[i], display = "none")))
    med.dens[i] <- med.dat$estimate
    hscale <- 0.4/max(smout$estimate) * wex
    base[[i]] <- smout$eval.points
    height[[i]] <- smout$estimate * hscale
    med.dens[i] <- med.dens[i] * hscale
    t <- range(base[[i]])
    baserange[1] <- min(baserange[1], t[1])
    baserange[2] <- max(baserange[2], t[2])
  }
  if (!add) {
    xlim <- if (n == 1) 
      at + c(-0.5, 0.5)
    else range(at) + min(diff(at))/2 * c(-1, 1)
    if (is.null(ylim)) {
      ylim <- baserange
    }
  }
  if (is.null(names)) {
    label <- 1:n
  }
  else {
    label <- names
  }
  boxwidth <- 0.05 * wex
  if (!add) 
    plot.new()
  if (!horizontal) {
    if (!add) {
      plot.window(xlim = xlim, ylim = ylim)
      axis(2)
      axis(1, at = at, label = label)
    }
    box()
    for (i in 1:n) {
      polygon(x = c(at[i] - radj*height[[i]], rev(at[i] + ladj*height[[i]])), 
              y = c(base[[i]], rev(base[[i]])), 
              col = col, border = border, 
              lty = lty, lwd = lwd)
      if (drawRect) {
        lines(at[c(i, i)], c(lower[i], upper[i]), lwd = lwd, 
              lty = lty)
        rect(at[i] - radj*boxwidth/2, 
             q1[i], 
             at[i] + ladj*boxwidth/2, 
             q3[i], col = rectCol)
        # median line segment
        lines(x = c(at[i] - radj*med.dens[i], 
                    at[i], 
                    at[i] + ladj*med.dens[i]),
              y = rep(med[i],3))
      }
    }
  }
  else {
    if (!add) {
      plot.window(xlim = ylim, ylim = xlim)
      axis(1)
      axis(2, at = at, label = label)
    }
    box()
    for (i in 1:n) {
      polygon(c(base[[i]], rev(base[[i]])), 
              c(at[i] - radj*height[[i]], rev(at[i] + ladj*height[[i]])), 
              col = col, border = border, 
              lty = lty, lwd = lwd)
      if (drawRect) {
        lines(c(lower[i], upper[i]), at[c(i, i)], lwd = lwd, 
              lty = lty)
        rect(q1[i], at[i] - radj*boxwidth/2, q3[i], at[i] + 
               ladj*boxwidth/2, col = rectCol)
        lines(y = c(at[i] - radj*med.dens[i], 
                    at[i], 
                    at[i] + ladj*med.dens[i]),
              x = rep(med[i],3))
      }
    }
  }
  invisible(list(upper = upper, lower = lower, median = med, 
                 q1 = q1, q3 = q3))
}


plot(x=NULL, y=NULL,
     xlim = c(0, 8), ylim=c(0, 5),
     type="n", ann=FALSE, axes=F)

axis(1, at=c(1:8),  labels=c("Solanum", "H.annuus", "Lactuca", "C.endivia", "C.intybus", "A.thaliana", "B.rapa", "G.max"))
axis(2)

for (i in unique(Eudi8sp$Taxon1)) {
  for (j in unique(Eudi8sp$Wi_Do)){
    vioplot2(Eudi8sp$LsMeans[which(Eudi8sp$Taxon1 == i & Eudi8sp$Wi_Do == j)],
             at = ifelse(i == c(1:8) "A", 1, 2),
             side = ifelse(j == 1, "left", "right"),
             col = ifelse(j == 1, "purple", "blue", "lightblue"),
             add = T)
  }
}
title("Violin plot", xlab="Treatment")
legend("bottomright", fill = c("purple", "lightblue"),
       legend = c("Group 1", "Group 2"), box.lty=0)



