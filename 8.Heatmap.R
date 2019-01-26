#####################
## Celine Caseys
## Plant Science dept. UC davis
## January 2019
###################

source('~/Desktop/Github_Eudicot/6.ZcoreLSmeans.R')

### Heatmap

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
Clade_color <- c("#b35806", "#542788")

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

#Eudicot_heatmap_stand
#Eudicot_heatmap_stand %>% save_iheatmap("Eudicot8sp_heatmap_stand.pdf", vwidth=2300, vheight=1700) 

library(MASS)
library(pvclust)
rownames(Eudicot7sp_heatmapS) <- Eudicot_info$SNames
# Put back nboot to 20'000 for final results. Put here at 200 for a quick demo run
PGenost_iso_pv <- pvclust(Eudicot7sp_heatmapS, nboot=200, method.hclust = "complete")
#plot(PGenost_iso_pv, main="Iso HClust 'complete'")
#pvrect(PGenost_iso_pv, alpha = .95)

tNoAdmixt_heatmapS <-  t(Eudicot7sp_heatmapS)
# Put back nboot to 20'000 for final results. Put here at 200 for a quick demo run
Noadmixt_genoo_pv <- pvclust(tNoAdmixt_heatmapS, nboot=200, method.hclust = "complete")
#plot(Noadmixt_genoo_pv, main="Plant Genotypes Hierarchical Clustering method='complete'")
#pvrect(Noadmixt_genoo_pv, alpha = .95)

