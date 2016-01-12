
################################################################################
# R script to accompany pregnancy.Rmd
################################################################################

## ---- libraries ----
library("caret")
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
preg <- get(load(data_path))$pregnancy
preg

## ---- prelim ----
otu_table(preg) <- transform_sample_counts(otu_table(preg), function(OTU) OTU/sum(OTU))


## ---- MDS ----
# this is in the paper (they color everything by a k-medoids clustering, though)
braydist <- phyloseq::distance(otu_table(preg), method="bray")
ord <- ordinate(otu_table(preg), method = "MDS", distance = braydist)
plot_ordination(preg, ord, col = "Preterm")

## ---- LDA ----
y <- as.factor(sample_data(preg)$Preterm)
sample_pred_X <- sample_data(preg)[, c(10:57)] # don't use obviously correlated variables
X <- cbind(otu_table(preg)@.Data, sample_pred_X)
train(x = X, y = y, method = "glmnet", trControl = trainControl(verbose = T, number = 4))

## ---- KCCA ----
