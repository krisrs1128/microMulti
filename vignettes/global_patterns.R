
################################################################################
# R script accompanying global_patterns.Rmd
################################################################################

## ---- libraries ----
library("phyloseq")
library("plyr")
library("dplyr")
library("microMulti")
library("mvarVis")
library("RCurl")

## ---- get-data ----
data_path <- file.path(tempdir(), "microbiomeData.RData")
download.file("https://dl.dropboxusercontent.com/s/5raqtzi9qmhl9ot/microbiomeData.RData",
              data_path, "libcurl")
global_patterns <- get(load(data_path))$global_patterns %>%
  subset_taxa(Phylum == "Bacteroidetes")

## ---- DCA ----
gp_dca <- ordinate(global_patterns)
gp_dca_mvar <- convert_to_mvar(gp_dca)
gp_dca_mvar <- supp_annotation(gp_dca_mvar, global_patterns)
plot_mvar_d3(gp_dca_mvar, c("arrow", "point"), height = 800)

## ---- PCA ----
gp_pca <- ordi(global_patterns@otu_table, method = "ade4_pca", scannf = F, nf = 5)
gp_pca_mvar <- supp_annotation(gp_pca, global_patterns, "li", "co")
plot_mvar_d3(gp_pca_mvar, c("point", "arrow"), height = 500)

## ---- preprocess-pca ----
X <- otu_table(global_patterns)@.Data
hist(log(1 + X), 100)
hist(log(X), 100)

gp_pca <- ordi(scale(log(1 + X)), method = "ade4_pca", scannf = F, nf = 10)
gp_pca_mvar <- supp_annotation(gp_pca, global_patterns, "li", "co")
plot_mvar_d3(gp_pca_mvar, c("point", "arrow"), height = 500)

################################################################################
# multitable methods
################################################################################

## ---- correspondence ----
gp_cca <- ordinate(global_patterns, formula = GP1 ~ ., , method = "CCA")
gp_cca_mvar <- convert_to_mvar(gp_cca)
plot_mvar_d3(gp_cca_mvar, c("arrow", "point"), height = 500)
