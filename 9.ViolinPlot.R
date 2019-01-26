#####################
## Celine Caseys
## Plant Science dept. UC davis
## January 2019
###################


########### Plotting Violin plot from LSmeans

Eudi8sp <- read.table(file="Eudi_Ara_Lsmeans1_8sp_row_Geno_All_plot.txt", header=T)

Eudi8sp$Taxon1 <- Eudi8sp$Taxon
levels(Eudi8sp$Taxon1) <- c("RosB_Ath", "RosB_Bra", "AstA_Cen", "AstA_Cin", "RosF_Glyc", "AstA_Han", "AstA_Lact", "AstS_Sol" )
Eudi8sp$clade <- Eudi8sp$Order
levels(Eudi8sp$clade) <- c("Asterids", "Rosids", "Rosids", "Asterids")

Eudi8sp$IsoTax <- paste(Eudi8sp$Isolate, Eudi8sp$Taxon, sep="-")

Eudi8sp$PlotOrder <- factor(Eudi8sp$Taxon)
levels(Eudi8sp$PlotOrder) <- c(6, 7, 4, 5,8,2, 3,1)
Eudi8sp$PlotOrder <- as.numeric(Eudi8sp$PlotOrder)

Eudi8sp$Taxon <- factor(Eudi8sp$Taxon, ordered=is.ordered(Eudi8sp$PlotOrder))


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


levels(Eudi8sp$Wi_Do) <- c("Col-0","High","Low" ,"Low", "coi1" )

# ggplot(Eudi8sp, aes(Taxon1, LsMeans, fill=Wi_Do)) +
#   geom_split_violin(scale="width") +
#   geom_boxplot(width=0.15)+ 
#   scale_fill_manual(values=c("#c7eae5", "#000000", "#878787", "#80cdc1"))+
#   theme_minimal()




