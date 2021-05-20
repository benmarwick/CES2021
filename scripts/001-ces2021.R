# from https://github.com/yesdavid/designspace_culttax_article_2021/blob/main/2_script/3_code_late_neolithic_early_bronze_age.R
# with minor modifications

library(cluster)
library(ggplot2)
library(readr)
library(Momocs)
library(dplyr)
library(maptree)
library(ggimage)
library(phangorn)
library(rworldmap)
library(raster)
library(NbClust)

if (!requireNamespace("ggtree", quietly = TRUE)){
  if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
    BiocManager::install("ggtree")
    library(ggtree)
  } else {
    BiocManager::install("ggtree")
    library(ggtree)
  }
} else {library(ggtree)}


outlines_combined_nicolas_2016 <- readRDS(file = file.path("data/outlines_combined_nicholas_2016.RDS"))


outlines_combined_nicolas_2016_centered <- 
  Momocs::coo_centre(outlines_combined_nicolas_2016) # center
outlines_combined_nicolas_2016_centered_scaled <- 
  Momocs::coo_scale(outlines_combined_nicolas_2016_centered) # scale
outlines_combined_nicolas_2016_centered_scaled <- 
  Momocs::coo_slidedirection(outlines_combined_nicolas_2016_centered_scaled, 
                             direction = "up")
# Momocs::panel(outlines_combined_nicolas_2016_centered_scaled)


# unification of outlines with catalogue-dataframe
nicolas_fleches_2016_catalog_ids_coordinates <- 
  readr::read_csv(file = "data/nicolas_fleches_2016_catalog_ids_with_coordinates.csv")

outlines_combined_nicolas_2016_names <- 
  names(outlines_combined_nicolas_2016_centered_scaled)
outlines_combined_nicolas_2016_names_splitted <- 
  strsplit(outlines_combined_nicolas_2016_names, split = "_")

ID_and_artefact_ID_list <- list()
for (name_index in 1:length(outlines_combined_nicolas_2016_names)){
  
  ID_and_artefact_ID_interrim_df <- data.frame(ID = paste0(outlines_combined_nicolas_2016_names_splitted[[name_index]][1], "-", outlines_combined_nicolas_2016_names_splitted[[name_index]][2]),
                                               ID_artefact <- outlines_combined_nicolas_2016_names[[name_index]])
  names(ID_and_artefact_ID_interrim_df) <- c("ID", "ID_artefact")
  ID_and_artefact_ID_list[[name_index]] <- ID_and_artefact_ID_interrim_df
  
}

ID_and_artefact_ID_df <- 
  do.call("rbind", ID_and_artefact_ID_list)

nicolas_fleches_2016_catalog_ids_coordinates_artefact_ID <- 
  dplyr::inner_join(ID_and_artefact_ID_df, 
                    nicolas_fleches_2016_catalog_ids_coordinates, 
                    by = "ID")


outlines_combined_nicolas_2016_centered_scaled <- 
  Momocs::Out(outlines_combined_nicolas_2016_centered_scaled$coo,
              fac = nicolas_fleches_2016_catalog_ids_coordinates_artefact_ID)

outlines_combined_nicolas_2016_centered_scaled <- 
  Momocs::filter(outlines_combined_nicolas_2016_centered_scaled, 
                 !ID_artefact %in% c("UK_60_XX_pseudo_no_10", "UK_15_XX_pseudo_no_4")) 
# remove fragmented outliers

length(which(nicolas_fleches_2016_catalog_ids_coordinates_artefact_ID$country == "Denmark")) # number of artefacts from denmark
length(which(nicolas_fleches_2016_catalog_ids_coordinates_artefact_ID$country == "France")) # number of artefacts from france
length(which(nicolas_fleches_2016_catalog_ids_coordinates_artefact_ID$country == "United Kingdom")) # number of artefacts from the uk

# harmonic calibration
outlines_combined_nicolas_2016_centered_scaled_harmonics <- 
  Momocs::calibrate_harmonicpower_efourier(outlines_combined_nicolas_2016_centered_scaled, 
                                           plot = F)  # Estimates the number of harmonics required for the Fourier methods implemented in Momocs.
# outlines_combined_nicolas_2016_centered_scaled_harmonics

# efourier
outlines_combined_nicolas_2016_centered_scaled_efourier <- 
  Momocs::efourier(outlines_combined_nicolas_2016_centered_scaled,
                   nb.h = as.matrix(outlines_combined_nicolas_2016_centered_scaled_harmonics[["minh"]])[[4,1]], # harmonics for 99.9%
                   norm = F) # you selected `norm=TRUE`, which is not recommended. See ?efourier --> probably no problem in our case


# PCA
outlines_combined_nicolas_2016_centered_scaled_PCA <- 
  Momocs::PCA(outlines_combined_nicolas_2016_centered_scaled_efourier) # PCA on Coe objects, using prcomp.

minimum_no_of_pcs_nicolas <- 
  Momocs::scree_min(outlines_combined_nicolas_2016_centered_scaled_PCA,
                    prop = 0.95) # minimum number of axis to use to retain a given proportion (i.e. prop = 0.99 to describe 99% of the variation) -- reduces computing time in the phylogeny estimation step:

Momocs::scree_plot(outlines_combined_nicolas_2016_centered_scaled_PCA)

# nice PCA plot 
# Create groups
pch.group <- c(rep(21, times=length(which(outlines_combined_nicolas_2016_centered_scaled_PCA$fac$country == "Denmark"))), 
               rep(22, times=length(which(outlines_combined_nicolas_2016_centered_scaled_PCA$fac$country == "France"))),
               rep(23, times=length(which(outlines_combined_nicolas_2016_centered_scaled_PCA$fac$country == "United Kingdom"))))
col.group <- c(rep("skyblue2", times=length(which(outlines_combined_nicolas_2016_centered_scaled_PCA$fac$country == "Denmark"))), 
               rep("gold", times=length(which(outlines_combined_nicolas_2016_centered_scaled_PCA$fac$country == "France"))),
               rep("green2", times=length(which(outlines_combined_nicolas_2016_centered_scaled_PCA$fac$country == "United Kingdom"))))

# plot
plot(outlines_combined_nicolas_2016_centered_scaled_PCA$x[,1],
     outlines_combined_nicolas_2016_centered_scaled_PCA$x[,2],
     xlab=paste("PCA 1 (", round(summary(outlines_combined_nicolas_2016_centered_scaled_PCA)$importance[2]*100, 1), "%)", sep = ""),
     ylab=paste("PCA 2 (", round(summary(outlines_combined_nicolas_2016_centered_scaled_PCA)$importance[5]*100, 1), "%)", sep = ""),
     pch=pch.group,
     col="black",
     bg=col.group,
     cex=1,
     las=1,
     asp=1,
     main = "Nicolas 2016; PCA of Bell Beaker arrowhead outline shapes")

# Add grid lines
abline(v=0, lty=2, col="grey50")
abline(h=0, lty=2, col="grey50")

# Add legend
legend("bottomleft", 
       legend=c("Denmark", "France", "United Kingdom"), 
       col="black", 
       pt.bg=c("skyblue2", "gold", "green2"), 
       pch=c(21, 22, 24), 
       pt.cex=1.5)



#### remove outliers
nicolas_outliers_db <- 
  fpc::dbscan(outlines_combined_nicolas_2016_centered_scaled_PCA$x, 
              eps = 0.3, 
              MinPts = 3)


plot(nicolas_outliers_db, 
     outlines_combined_nicolas_2016_centered_scaled_PCA$x, 
     main = "DBSCAN", 
     frame = FALSE)

nicolas_outliers_cluster <- 
  data.frame(name = row.names(outlines_combined_nicolas_2016_centered_scaled_PCA$x), 
             value = nicolas_outliers_db$cluster, 
             row.names = NULL)

nicolas_outliers_cluster_outlier_names <- 
  subset(nicolas_outliers_cluster, value != 1)

outlines_combined_nicolas_2016_centered_scaled_PCA$fac$outlier_names <- NA
outlines_combined_nicolas_2016_centered_scaled$fac$cluster <- as.factor(nicolas_outliers_cluster$value)

for (vector_index in 1:length(match(nicolas_outliers_cluster_outlier_names$name, 
                                    outlines_combined_nicolas_2016_centered_scaled_PCA$fac$ID_artefact))){
  
  current_index <- match(nicolas_outliers_cluster_outlier_names$name, 
                         outlines_combined_nicolas_2016_centered_scaled_PCA$fac$ID_artefact)[vector_index]
  
  outlines_combined_nicolas_2016_centered_scaled_PCA$fac$outlier_names[current_index] <- nicolas_outliers_cluster_outlier_names$name[vector_index]
}


# look at the shape of the outliers
nicolas_2016_with_outliers <- 
  Momocs::slice(outlines_combined_nicolas_2016_centered_scaled, 
                match(nicolas_outliers_cluster_outlier_names$name, 
                      outlines_combined_nicolas_2016_centered_scaled_PCA$fac$ID_artefact))

Momocs::panel(nicolas_2016_with_outliers,
              fac = "country_code",
              names = T,
              col = "grey")

# analysis without outliers
nicolas_2016_without_outliers <- 
  Momocs::slice(outlines_combined_nicolas_2016_centered_scaled, 
                -match(nicolas_outliers_cluster_outlier_names$name, 
                       outlines_combined_nicolas_2016_centered_scaled_PCA$fac$ID_artefact))


# harmonic calibration
nicolas_2016_without_outliers_harmonics <- 
  Momocs::calibrate_harmonicpower_efourier(nicolas_2016_without_outliers, 
                                           plot = F)  # Estimates the number of harmonics required for the Fourier methods implemented in Momocs. This is the only step in this section that produces data we need in the subsequent step.

# efourier
nicolas_2016_without_outliers_efourier <- 
  Momocs::efourier(nicolas_2016_without_outliers,
                   nb.h = as.matrix(nicolas_2016_without_outliers_harmonics[["minh"]])[[4,1]], # harmonics for 99.9%
                   norm = F) 
# PCA
nicolas_2016_without_outliers_PCA <- 
  Momocs::PCA(nicolas_2016_without_outliers_efourier) # PCA on Coe objects, using prcomp.

nicolas_2016_screeplot_wo_outliers <- 
  Momocs::scree_plot(nicolas_2016_without_outliers_PCA, 
                     nax = 1:50) 

# see how many PCs will give 99.9%  39 PCs gets there
summary(nicolas_2016_without_outliers_PCA)

minimum_no_of_pcs_nicolas_without_outliers <- 39
  # ncol(nicolas_2016_without_outliers_PCA$x) #use all PC axes


gg_nicolas <- Momocs::PCcontrib(nicolas_2016_without_outliers_PCA,
                                nax = 1:5,
                                sd.r = c(-2,-1,0,1,2)) 
gg_nicolas$gg + 
  geom_polygon(fill="slategrey", col="black") + 
  theme_minimal(base_size = 20) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

ggsave("figures/nicolas_arrowheads_pca_contrib.png",
       h = 10, w = 8)


library(ggplot2)
gg_nicolas <- gg_nicolas$gg +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

# nice PCA plot 
# Create groups
pch.group <- c(rep(21, times=length(which(nicolas_2016_without_outliers_PCA$fac$country == "Denmark"))), 
               rep(22, times=length(which(nicolas_2016_without_outliers_PCA$fac$country == "France"))),
               rep(23, times=length(which(nicolas_2016_without_outliers_PCA$fac$country == "United Kingdom"))))
col.group <- c(rep("skyblue2", times=length(which(nicolas_2016_without_outliers_PCA$fac$country == "Denmark"))), 
               rep("gold", times=length(which(nicolas_2016_without_outliers_PCA$fac$country == "France"))),
               rep("green2", times=length(which(nicolas_2016_without_outliers_PCA$fac$country == "United Kingdom"))))

# plot
plot(nicolas_2016_without_outliers_PCA$x[,1],
     nicolas_2016_without_outliers_PCA$x[,2],
     xlab=paste("PCA 1 (", round(summary(nicolas_2016_without_outliers_PCA)$importance[2]*100, 1), "%)", sep = ""),
     ylab=paste("PCA 2 (", round(summary(nicolas_2016_without_outliers_PCA)$importance[5]*100, 1), "%)", sep = ""),
     pch=pch.group,
     col="black",
     bg=col.group,
     cex=1,
     las=1,
     asp=1,
     main = "Nicolas 2016; PCA of Bell Beaker arrowhead outline shapes without outliers")
# Add grid lines
abline(v=0, lty=2, col="grey50")
abline(h=0, lty=2, col="grey50")
# Add legend
legend("bottomleft", 
       legend=c("Denmark", "France", "United Kingdom"), 
       col="black", 
       pt.bg=c("skyblue2", "gold", "green2"), 
       pch=c(21, 22, 24), 
       pt.cex=1.5)

############################# NJ typochronology ############################# 


typochronologie_csv <- readr::read_csv("data/nicolas_2017_typochronologie.csv")

typochronologie_csv <- dplyr::distinct(typochronologie_csv, ID_country, .keep_all = T)

typochronologie_csv_UK <- subset(typochronologie_csv, Country == "UK")
typochronologie_csv_FR <- subset(typochronologie_csv, Country == "FR")


nicolas_2016_without_outliers_PCA_as_df <- 
  as.data.frame(nicolas_2016_without_outliers_PCA$x[,1:minimum_no_of_pcs_nicolas_without_outliers])


nicolas_2016_without_outliers_PCA_as_df$ID_country <- sapply(rownames(nicolas_2016_without_outliers_PCA_as_df),
                                                             function(x){
                                                               strsplit(x, split = "(?<=.{5})", perl = TRUE)[[1]][1] # RegEx "(?<=.{5})" means: split into chunks of 5 characters long
                                                             })  

########## FR
nicolas_2016_without_outliers_PCA_as_df_subset_typochron_FR <- subset(nicolas_2016_without_outliers_PCA_as_df, ID_country %in% typochronologie_csv_FR$ID_country)

names_artefacts_ID <- data.frame(artefact_ID = rownames(nicolas_2016_without_outliers_PCA_as_df_subset_typochron_FR),
                                 ID_country = nicolas_2016_without_outliers_PCA_as_df_subset_typochron_FR$ID_country)
names_artefacts_ID_and_period <- dplyr::left_join(names_artefacts_ID, typochronologie_csv_FR[,c("ID_country", "Period")], by = "ID_country")
names_artefacts_ID_and_period$Period <- as.factor(names_artefacts_ID_and_period$Period)

# plot
plot(nicolas_2016_without_outliers_PCA_as_df_subset_typochron_FR$PC1,
     nicolas_2016_without_outliers_PCA_as_df_subset_typochron_FR$PC2,
     xlab="PC1", 
     ylab="PC2",
     cex=1,
     las=1,
     asp=1,
     main = "Nicolas 2016; PCA of French Bell Beaker arrowhead outline shapes without outliers")
# Add grid lines
abline(v=0, lty=2, col="grey50")
abline(h=0, lty=2, col="grey50")

nicolas_2016_without_outliers_PCA_as_df_subset_typochron_FR$ID_country <- NULL

nicolas_typochron_FR_dist <- dist(nicolas_2016_without_outliers_PCA_as_df_subset_typochron_FR, method = "euclidean")

nicolas_typochron_FR_NJ <- phangorn::NJ(nicolas_typochron_FR_dist)

additional_information_period <- data.frame(Period = names_artefacts_ID_and_period[,c("Period")])
rownames(additional_information_period) <- names_artefacts_ID_and_period$artefact_ID


# need to run to at least here before moving to the next script file

NJ_nicolas_ggtree <- ggtree(nicolas_typochron_FR_NJ) %<+% names_artefacts_ID_and_period +
  geom_tiplab(size=3, 
              aes(label = ID_country,
                  color = ID_country)) +
  geom_treescale() + 
  scale_colour_discrete(na.translate = F) + 
  guides(color=guide_legend(title="Site ID"))

NJ_nicolas_ggtree



