#Install Packages
install.packages("pacman")
pacman::p_load("BiocManager","glmnet", "caret","ellipsis", "vctrs", "survival", "survminer", "readr", "grpreg","testthat", "pkgload", "devtools", "ROCit", "car", "ggpubr", "gplots", "testthat", "pkgload")
install.packages("tidyr")
suppressPackageStartupMessages(library("survival","survminer", "grpreg"))
suppressPackageStartupMessages(library(glmnet))
remotes::install_github("cran/plotmo", force = TRUE)
suppressPackageStartupMessages(library(plotmo))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(caret, compareGroups))
suppressPackageStartupMessages(library(devtools))
suppressPackageStartupMessages(library(ROCit))

install.packages("varhandle")
suppressPackageStartupMessages(library("varhandle"))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library("car"))
suppressPackageStartupMessages(library(gridExtra, grid, "ggpubr"))
suppressPackageStartupMessages(library(lattice, forcats, ellipsis))
suppressPackageStartupMessages(library(vctrs))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(survival))
library(dplyr)
library(survminer)

setwd("C:\\Users\\lizec\\Desktop\\CI_Project")

#Read radiomics data
data <- read.csv("final_output.csv")
data <- as.data.frame(data)

data <- data %>%
  filter(ID != "LRAD89_CT")

mask_tumour <- data[grepl("CT\\.nii\\.gz$", data$Mask), ]
mask_CTs <- data[grepl("CTs\\.nii\\.gz$", data$Mask), ]
mask_CTp <- data[grepl("CTp\\.nii\\.gz$", data$Mask), ]

#adjust ID
mask_tumour$ID <- gsub("_CT", "", mask_tumour$ID)
mask_CTs$ID <- gsub("_CT", "", mask_CTs$ID)
mask_CTp$ID <- gsub("_CT", "", mask_CTp$ID)

difference_dataset1 <- anti_join(mask_tumour, mask_CTp, by = "ID")
#596 does not have ctp

#remove unwanted data
tumour_filtered <- mask_tumour[, -c(2:33)]
mask_CTs_filtered <- mask_CTs[, -c(2:33)]
mask_CTp_filtered <- mask_CTp[, -c(2:33)]

tumour_filtered <- mask_tumour[, -c(4:5)]
mask_CTs_filtered <- mask_CTs[, -c(4:5)]
mask_CTp_filtered <- mask_CTp[, -c(4:5)]


#Read clinical data
clinical_data <- read_excel("Cardiac_clinical.xlsx")
clinical_data <- as.data.frame(clinical_data)
clinical_data <- clinical_data[-101, ]
clinical_data <- clinical_data[, -1]
clinical_data$ID <- paste("LRAD", clinical_data$ID, sep = "")
clinical_data <- subset(clinical_data, select = c(ID, OS.event, Overall.survival..days., Age, Stage, Performance, CAC, all_calc))
clinical_data <- clinical_data[!is.na(clinical_data$OS.event), ]
#merge data and remove missing clinical data

merged_data <- merge(clinical_data,tumour_filtered,by = "ID")

#adjust stage
for (i in 1:nrow(merged_data)) {
  if (is.na(merged_data$Stage[i])) {
    merged_data$Stage[i] <- 0
  } else if (merged_data$Stage[i] == "1a") {
    merged_data$Stage[i] <- 1
  } else if (merged_data$Stage[i] == "1b") {
    merged_data$Stage[i] <- 2
  } else if (merged_data$Stage[i] == "2a") {
    merged_data$Stage[i] <- 3
  } else if (merged_data$Stage[i] == "2b") {
    merged_data$Stage[i] <- 4
  } else if (merged_data$Stage[i] == "3") {
    merged_data$Stage[i] <- 5
  } else if (merged_data$Stage[i] == "4") {
    merged_data$Stage[i] <- 6
  } else {
    merged_data$Stage[i] <- 7
  } 
}

merged_df <- merged_data %>%
  inner_join(mask_CTp_filtered, by = "ID") %>%
  inner_join(mask_CTs_filtered, by = "ID")


merged_df$Stage <- as.numeric(merged_df$Stage)
class(merged_df$Performance)
class(merged_df$CAC)
class(merged_df$all_calc)

#create test and training data frames
set.seed(32)

train.index <- createDataPartition(merged_df$Age | merged_df$Overall.survival..days. | merged_df$Stage | merged_df$Performance | merged_df$CAC |
                                   merged_df$all_calc, p = 2/3, list = FALSE) #age, stage, survival 
train_df <- merged_df[ train.index,]
test_df <- merged_df[-train.index,]


#create covariates
#covariates <- c(colnames(merged_data))
#covariates <- covariates[6:ncol(merged_data)]

#standardisation
train_df <- train_df[, -c(9:38), drop = FALSE]
test_df <- test_df[, -c(9:38), drop = FALSE]

train_df <- train_df[, -c(11:12), drop = FALSE]
test_df <- test_df[, -c(11:12), drop = FALSE]


covariates <- c(colnames(train_df))
covariates <- covariates[9:ncol(train_df)]

train_sd <- train_df
for (i in 6:ncol(train_df)){
  if (sd((train_df[[covariates[i-6+1]]]), na.rm = TRUE) != 0){
    train_sd[[covariates[i-6+1]]]<- scale(train_df[[covariates[i-6+1]]])
  } else{
    train_sd[[covariates[i-6+1]]] <- 0
  }
}

train_sd[is.na(train_sd)] <- 0

test_sd <- test_df
for (i in 6:ncol(test_df)){
  if (sd((test_df[[covariates[i-6+1]]]), na.rm = TRUE) != 0){
    test_sd[[covariates[i-6+1]]]<- scale(test_df[[covariates[i-6+1]]])
  } else{
    test_sd[[covariates[i-6+1]]] <- 0
  }
}

test_sd[is.na(test_sd)] <- 0
#univariate cox regression

#one_level_cols <- which(sapply(train_sd, function(x) length(unique(x))) == 1)
#one_level_cols
#removed_one_levels <- train_sd[,-one_level_cols]

#univ_train_sd <- train_sd[, -c(1,4:5)]

#univ_covariates <- colnames(univ_train_sd)
#univ_formulas <- sapply(univ_covariates,
                        #function(x) as.formula(paste('Surv(Overall.survival..days., OS.event)~', x)))
#univ_models <- lapply( univ_formulas, function(x){coxph(x, data = univ_train_sd )})


# Identify columns with only one unique value
one_level_cols <- which(sapply(train_sd, function(x) length(unique(x))) == 1)

# Remove columns with only one unique value
removed_one_levels <- train_sd[,-one_level_cols]
removed_one_levels

# Remove the dependent variable and other unwanted columns from the data
univ_train_sd <- removed_one_levels[, -c(1, 4:8)]

# Get the names of the remaining covariates
univ_covariates <- colnames(univ_train_sd)

# Create formulas for Cox regression models for each covariate
univ_formulas <- sapply(univ_covariates,
                        function(x) as.formula(paste('Surv(Overall.survival..days., OS.event)~', x)))

# Fit Cox regression models for each covariate
univ_models <- lapply(univ_formulas, function(x) {coxph(x, data = train_sd)})



univ_results <- lapply(univ_models, function(x) {
  x <- summary(x)
  p.value <- signif(x$wald["pvalue"], digits = 2)
  wald.test <- signif(x$wald["test"], digits = 2)
  beta <- signif(x$coef[1], digits = 2)
  HR <- signif(x$coef[2], digits = 2)
  HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
  HR.confint.upper <- signif(x$conf.int[,"upper .95"], 2)
  HR <- paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")")
  res <- c(beta, HR, wald.test, p.value)
  names(res) <- c("beta", "HR (95% CI for HR)", "wald.test", "p.value")
  return(res)
})

univ_results <- univ_results[-c(1, 2)]

# Get the row lengths for each element in univ_results
row_lengths <- sapply(univ_results, NROW)

# Find the row length that differs from the majority
unique_length <- unique(row_lengths)
majority_length <- unique_length[which.max(tabulate(match(row_lengths, unique_length)))]

# Identify the elements with the differing row length
differing_elements <- names(row_lengths[row_lengths != majority_length])

# Remove the elements with differing row lengths from univ_results
univ_results_cleaned <- univ_results[setdiff(names(univ_results), differing_elements)]



res <- t(as.data.frame(univ_results_cleaned, check.names = FALSE))
coxresult <- as.data.frame(res)
coxresult$FDR<-p.adjust(coxresult$p.value,method="fdr")
coxresult<-coxresult[order(coxresult$FDR, decreasing=F),]

#selection according to p.value 
filtered_coxresult <- coxresult["p.value"]

features <- rownames(filtered_coxresult)
filtered_coxresult <- cbind(features, filtered_coxresult)
rownames(filtered_coxresult) <- NULL
filtered_coxresult$p.value <- as.numeric(filtered_coxresult$p.value)
class(filtered_coxresult$p.value)
regular_pval <- format(filtered_coxresult$p.value, scientific = FALSE)
filtered_coxresult$p.value <- regular_pval

p_thres <- 0.05

coxresult_significant <- filtered_coxresult %>% 
  filter(p.value<p_thres)



#Lasso

train_features <- coxresult_significant$features

train_matching <- names(train_sd)[names(train_sd) %in% train_features]
train_subset <- subset(train_sd, select = train_matching)
selected_columns <- c("ID", "Overall.survival..days.", "OS.event", "diagnostics_Mask.interpolated_Minimum.y")
train_raw <- train_sd[selected_columns]

#test_subset <- subset(test_sd, select = train_matching)

train_subset <- na.omit(train_subset)
removed_rows <- attr(train_subset, "na.action")
removed_rows

nFolds <- 50
foldid <- sample(rep(seq(nFolds), length.out = nrow(train_subset)))

forlasso <- train_subset
fit <- glmnet(x=as.matrix(forlasso), y= Surv(train_raw$"Overall.survival..days.", train_raw$"OS.event"), family="cox", alpha=1, nfolds = nFolds,foldid = foldid)
plot_glmnet(fit,xvar='lambda',label=7)

cvfit <- cv.glmnet(x=as.matrix(forlasso), y= Surv(train_raw$"Overall.survival..days.", train_raw$"OS.event"), family="cox",alpha=1,  nfolds = nFolds, foldid = foldid)
plot(cvfit)

fit <- glmnet(x=as.matrix(forlasso), y= Surv(train_raw$"Overall.survival..days.", train_raw$"OS.event"), family="cox",alpha=1,nfolds = nFolds, lambda=cvfit$lambda.min, foldid = foldid)
fit$beta[,1]

# Extract the coefficient values
coefficients <- fit$beta[,1]

# Get the features with non-zero coefficients
non_zero_features <- coefficients[coefficients != 0]
non_zero_features

#Predicting RPV
#training
my_prediction_model_train <- predict(fit, newx = as.matrix(train_subset[,train_features]), cvfit$lambda.min)


# K-mean clustering

predict.kmeans <- function(object, newdata){
  centers <- object$centers
  n_centers <- nrow(centers)
  dist_mat <- as.matrix(dist(rbind(centers, newdata)))
  dist_mat <- dist_mat[-seq(n_centers), seq(n_centers)]
  max.col(-dist_mat)
}

set.seed(5)
no_partition <- 5
k_mean_lasso <- cbind(my_prediction_model_train,forlasso)
colnames(k_mean_lasso)[1] <- "RPV"
km.res <- kmeans(k_mean_lasso, no_partition, nstart = 1)

k_mean_lasso$grouping <- km.res$cluster
k_mean_lasso$grouping

class(k_mean_lasso$grouping)

k_mean_lasso$grouping <- as.numeric(k_mean_lasso$grouping)

k_mean_lasso <- cbind(k_mean_lasso,train_raw$OS.event)
k_mean_lasso <- cbind(k_mean_lasso,train_raw$Overall.survival..days.)


colnames(k_mean_lasso)[604] <- "OS.event"
colnames(k_mean_lasso)[605] <- "Overall.survival..days."

for (i in 1:nrow(k_mean_lasso)) {
  if (k_mean_lasso$grouping[i] %in% c(3,4)) {
    k_mean_lasso$group_final[i] <- 1
  } else {
    k_mean_lasso$group_final[i] <- 2
  }
}

km_fit <- survfit(Surv(Overall.survival..days.,OS.event) ~ group_final, data= k_mean_lasso)

theme <- theme(axis.line = element_line(colour = "black"),
               panel.grid.major = element_line(colour = "white"),
               panel.grid.minor = element_line(colour = "white"),
               panel.border = element_blank(),
               panel.background = element_blank()
) 

a <- ggsurvplot(
  km_fit,     # survfit object with calculated statistics.
  data = k_mean_lasso,               # data used to fit survival curves. 
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,1100),        # present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 100,     # break X axis in time intervals by 100.
  ggtheme = theme,         # customize plot and risk table with a theme.
  risk.table.y.text.col = F, # colour risk table text annotations.
  risk.table.y.text = T, # show bars instead of names in text annotations
  # in legend of risk table
  legend.labs=c("high risk","low risk"),
  risk.table = T,
  palette="jco",
  tables.theme = theme_survminer(font.main = 12),
  title = "Kaplan-Meier Plot Stratified Based on RPV \n (Training set, Tumour radiomics)",
  xlab="Time in Days",
  ylab="Probability of Overall Survival",
  #surv.median.line = "v",
  ylim=c(0,1),
  cumevents=F,
  surv.scale="percent",
  font.main = c(14, "bold"),
  font.x = c(12), 
  font.y = c(12), font.tickslab = c(12),
  font.legend = c(12), risk.table.fontsize = 4
)

a

library(survminer)


cate_dis <- train_raw$OS.event
cate_dis

for (i in 1:length(cate_dis)) {
  if (cate_dis[i] == "1"){
    cate_dis[i] <- '+'
  } else {
    cate_dis[i] <- '-'
  }
}



score_dis <- unname(my_prediction_model_train)



ROCit_obj <- rocit(score=score_dis, class=cate_dis, negref = "+", method ="bin")


AUC_obj <- ciAUC(ROCit_obj, level = 0.95)
p <- plot(ROCit_obj)
text(0.95, 0.2, paste0("AUC=", round(AUC_obj$AUC, 2), ", 95% CI [", round(AUC_obj$lower, 2), ",", round(AUC_obj$upper, 2), "]"), adj = 1, font = 4, cex=1.0)
title("Prediction of 3-year survival (Training Set - Tumour)")




#testing

test_subset <- subset(test_sd, select = train_matching)

train_subset <- na.omit(train_subset)
removed_rows <- attr(train_subset, "na.action")
removed_rows

my_prediction_model_test <- predict(fit, newx = as.matrix(test_subset[,train_features]), cvfit$lambda.min)

set.seed(14)
k_mean_test <- cbind(my_prediction_model_test,test_subset)
no_partition <-  10
colnames(k_mean_test)[1] <- "RPV"
km.test <- kmeans(k_mean_test, no_partition, nstart = 1)
k_mean_test$grouping <- km.test$cluster
k_mean_test$grouping

k_mean_test$grouping <- as.numeric(k_mean_test$grouping)
selected_columns_test <- c("ID", "Overall.survival..days.", "OS.event", "diagnostics_Mask.interpolated_VolumeNum.y")
test_raw <- test_sd[selected_columns_test]

k_mean_test <- cbind(k_mean_test,test_raw$OS.event)
k_mean_test <- cbind(k_mean_test,test_raw$Overall.survival..days.)
colnames(k_mean_test)[604] <- "OS.event"
colnames(k_mean_test)[605] <- "Overall.survival..days."

for (i in 1:nrow(k_mean_test)){
  if (k_mean_test$grouping[i] %in% c(1,2,10)) {
    k_mean_test$group_final[i] <- 1
  } else {
    k_mean_test$group_final[i] <- 2
  }
}

k_mean_test$group_final
km_fit_test <- survfit(Surv(Overall.survival..days.,OS.event) ~ group_final, data= k_mean_test)


b <- ggsurvplot(
  km_fit_test,     # survfit object with calculated statistics.
  data = k_mean_test,               # data used to fit survival curves. 
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,1100),        # present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 100,     # break X axis in time intervals by 100.
  ggtheme = theme,         # customize plot and risk table with a theme.
  risk.table.y.text.col = F, # colour risk table text annotations.
  risk.table.y.text = T, # show bars instead of names in text annotations
  # in legend of risk table
  legend.labs=c("high risk","low risk"),
  risk.table = T,
  palette="jco",
  tables.theme = theme_survminer(font.main = 12),
  title = "Kaplan-Meier Plot Stratified Based on RPV \n (Testing set, Tumour radiomics)",
  xlab="Time in Days",
  ylab="Probability of Overall Survival",
  #surv.median.line = "v",
  ylim=c(0,1),
  cumevents=F,
  surv.scale="percent",
  font.main = c(14, "bold"),
  font.x = c(12), 
  font.y = c(12), font.tickslab = c(12),
  font.legend = c(12), risk.table.fontsize = 4
)

b


cate_test <- test_raw$OS.event
cate_test


for (i in 1:length(cate_test)) {
  if (cate_test[i] == "1"){
    cate_test[i] <- '+'
  } else {
    cate_test[i] <- '-'
  }
}


score_test <- unname(my_prediction_model_test)



ROCit_obj_test <- rocit(score=score_test, class=cate_test, negref = "+", method ="bin")


AUC_obj_test <- ciAUC(ROCit_obj_test, level = 0.95)
q <- plot(ROCit_obj_test)
text(0.95, 0.2, paste0("AUC=", round(AUC_obj_test$AUC, 2), ", 95% CI [", round(AUC_obj_test$lower, 2), ",", round(AUC_obj_test$upper, 2), "]"), adj = 1, font = 4, cex=1.0)
title("Prediction of 3-year survival (Testing Set - Tumour)")


# Calculate True Positives (TP), False Positives (FP), and False Negatives (FN)
train_pred <- k_mean_lasso$group_final
train_true <- train_raw$OS.event

TP_train <- sum(train_pred == 1 & train_true == 1)
FP_train <- sum(train_pred == 1 & train_true == 0)
FN_train <- sum(train_pred == 2 & train_true == 1)

# Calculate Precision, Recall, and F1 score
precision_train <- TP_train / (TP_train + FP_train)
recall_train <- TP_train / (TP_train + FN_train)
f1_score_train <- 2 * (precision_train * recall_train) / (precision_train + recall_train)

print(paste("Training Set - Precision:", precision_train))
print(paste("Training Set - Recall:", recall_train))
print(paste("Training Set - F1 Score:", f1_score_train))


# Calculate True Positives (TP), False Positives (FP), and False Negatives (FN)
test_pred <- k_mean_test$group_final
test_true <- test_raw$OS.event

TP_test <- sum(test_pred == 1 & test_true == 1)
FP_test <- sum(test_pred == 1 & test_true == 0)
FN_test <- sum(test_pred == 2 & test_true == 1)

# Calculate Precision, Recall, and F1 score
precision_test <- TP_test / (TP_test + FP_test)
recall_test <- TP_test / (TP_test + FN_test)
f1_score_test <- 2 * (precision_test * recall_test) / (precision_test + recall_test)

print(paste("Testing Set - Precision:", precision_test))
print(paste("Testing Set - Recall:", recall_test))
print(paste("Testing Set - F1 Score:", f1_score_test))

export_train <- data.frame(ID = train_sd[, 1])
write.csv(export_train, "train_sd.csv", row.names = FALSE)
write.csv(test_sd, "test_sd.csv", row.names = FALSE)
export_train
