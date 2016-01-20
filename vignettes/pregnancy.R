
################################################################################
# R script to accompany pregnancy.Rmd
################################################################################

## ---- libraries ----
library("caret")
library("gflasso")
library("phyloseq")
library("plyr")
library("dplyr")
library("microMulti")
library("mvarVis")
library("RCurl")
library("data.table")
library("ggplot2")
theme_set(theme_bw())

## ---- get-data ----
data_path <- file.path(tempdir(), "microbiomeData.RData")
download.file("https://dl.dropboxusercontent.com/s/5raqtzi9qmhl9ot/microbiomeData.RData",
              data_path, "libcurl")
preg <- get(load(data_path))$pregnancy$PSPost$Vaginal_Swab
preg

## ---- prelim ----
raw_otu <- otu_table(preg)
p_otu <- transform_sample_counts(otu_table(preg), function(OTU) OTU/sum(OTU))

## ---- MDS ----
# this is in the paper (they color everything by a k-medoids clustering, though)
braydist <- phyloseq::distance(p_otu, method="bray")
ord <- ordinate(p_otu, method = "MDS", distance = braydist)
plot_ordination(preg, ord, col = "Preterm")
rm(ord, braydist)

## ---- filter-otus ----
keep <- genefilter_sample(raw_otu, filterfun_sample(function(x) x > 0),
                          A = 0.1 * nsamples(preg))
raw_otu <- raw_otu[, keep]

## ---- get-dates ----
master <- sample_data(preg)@.Data %>%
  as.data.frame() %>%
  as.data.table()
setnames(master, colnames(sample_data(preg)))
master <- master %>%
  group_by(SubjectID) %>%
  mutate(relative_day = diff_dates(DateColl))
master %>%
  select(SubjectID, relative_day, DateColl) %>%
  arrange(SubjectID, relative_day)

diff_dates <- function(x) {
  dates <- strptime(x, "%m/%d/%y %H:%M")
  difftime(dates, min(dates), units = "days") %>%
    as.numeric()
}

## ---- ts-plots ----
for(j in seq_len(ncol(raw_otu))) {
  X <- data.table(master, raw_otu[, j])
  mX <- melt(X, id.vars = colnames(master))
  print(ggplot(mX) +
    geom_point(aes(x = relative_day, y = value, col = SubjectNice), alpha = 0.3, size = 1) +
    geom_line(aes(x = relative_day, y = value, col = SubjectNice), alpha = 0.2))
  Sys.sleep(.4)
}

## ---- glmnet ----
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

## ---- multitask-model ----
Y <- data.frame(sample_data(preg))[, c("Marginal", "Preterm")]
table(Y)
Y <- colwise(as.numeric)(Y)
# Not really applicable, since response is binary
gflasso(Y, X, matrix(c(1, .1, .1, 1), 2, 2),
        list(lambda = .05, delta_conv = 1e-4))
