
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
preg <- get(load(data_path))$pregnancy$PSPost$Vaginal_Swab
preg

## ---- prelim ----
otu_table(preg) <- transform_sample_counts(otu_table(preg), function(OTU) OTU/sum(OTU))

## ---- MDS ----
# this is in the paper (they color everything by a k-medoids clustering, though)
braydist <- phyloseq::distance(otu_table(preg), method="bray")
ord <- ordinate(otu_table(preg), method = "MDS", distance = braydist)
plot_ordination(preg, ord, col = "Preterm")
rm(list(ord, braydist))

## ---- glmnet ----
keep <- genefilter_sample(otu_table(preg), filterfun_sample(function(x) x > 0),
                          A = 0.25 * nsamples(preg))
otu_table(preg) <- otu_table(preg)[, keep]

y <- as.factor(sample_data(preg)$Preterm)
#sample_pred_X <- sample_data(preg)[, c()] # don't use obviously correlated variables
#X <- cbind(otu_table(preg)@.Data, sample_pred_X)
X <- otu_table(preg)@.Data
glmnet_model <- train(x = X, y = y, method = "glmnet",
                      trControl = trainControl(verbose = T))
glmnet_model
impt <- varImp(glmnet_model)$importance
impt$id <- rownames(impt)
sorted_vars <- impt %>%
  arrange(desc(Overall)) %>%
  select(id) %>% unlist()
preg@tax_table[sorted_vars[1:10], 4:7] # top ten most predictive species

## ---- KCCA ----
