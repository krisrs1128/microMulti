
################################################################################
# R script accompanying global_patterns.Rmd
################################################################################

## ---- libraries ----
library("phyloseq")
library("plyr")
library("dplyr")
library("mvarVis")
library("RCurl")

## ---- utils ----
supp_annotation <- function(mvar_object, phyloseq_object, n1 = "species",
                            n2 = "site") {
  mvar_object@table[[n1]]@annotation <- data.frame(mvar_object@table[[n1]]@annotation,
                                                   phyloseq_object@tax_table)
  mvar_object@table[[n2]]@annotation <- data.frame(mvar_object@table[[n2]]@annotation,
                                                   phyloseq_object@sam_data)
  mvar_object
}

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

