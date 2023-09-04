#set working directory  
setwd("C:\\Users\\lizec\\Desktop\\CI_Project")
pacman::p_load("BiocManager","glmnet", "caret","ellipsis", "vctrs", "survival", "survminer", "readr", "grpreg","testthat", "pkgload", "devtools", "ROCit", "car", "ggpubr", "gplots", "testthat", "pkgload")

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


#Read radiomics data
data <- read.csv("final_output.csv")
data <- as.data.frame(data)

#filter labels
calc_masks <- data[grepl("label", data$Mask), ]

#adjust ID
calc_masks$ID <- gsub("_CT", "", calc_masks$ID)

#remove empty inputs
calc_masks <- calc_masks %>%
  filter(diagnostics_Versions_PyRadiomics != "")

#Load clinical data
clinical_data <- read_excel("Cardiac_clinical.xlsx")
clinical_data <- as.data.frame(clinical_data)
clinical_data <- clinical_data[-100, ]
clinical_data <- clinical_data[-1, ]
clinical_data <- clinical_data[, -1]
clinical_data$ID <- paste("LRAD", clinical_data$ID, sep = "")
clinical_data <- subset(clinical_data, select = c(ID, OS.event, Overall.survival..days., Age, Stage, Performance, CAC, all_calc))

merged_data <- merge(clinical_data,calc_masks,by = "ID")

#Coronary artery calcificaiton
coronary_calc <- merged_data[grepl("label_2|label_3|label_4", merged_data$Mask), ]

calc_RCA <- merged_data[grepl("label_2", merged_data$Mask), ]
calc_RCA <- calc_RCA[,-(2:40) ]
calc_RCA <- calc_RCA[,-(4:5) ]

calc_LMS_LAD <- merged_data[grepl("label_3", merged_data$Mask), ]
calc_LMS_LAD <- calc_LMS_LAD[,-(2:40) ]
calc_LMS_LAD <- calc_LMS_LAD[,-(4:5) ]

calc_LCX <- merged_data[grepl("label_4", merged_data$Mask), ]
calc_LCX <- calc_LCX[,-(2:40) ]
calc_LCX <- calc_LCX[,-(4:5) ]

aorta <- merged_data[grepl("label_1", merged_data$Mask), ]
aorta <- aorta[,-(2:40) ]
aorta <- aorta[,-(4:5) ]

av <- merged_data[grepl("label_5", merged_data$Mask), ]
av <- av[,-(2:40) ]
av <- av[,-(4:5) ]

mv <- merged_data[grepl("label_6", merged_data$Mask), ]
mv <- mv[,-(2:40) ]
mv <- mv[,-(4:5) ]

#tumour mask
mask_tumour <- data[grepl("CT\\.nii\\.gz$", data$Mask), ]

#adjust ID
mask_tumour$ID <- gsub("_CT", "", mask_tumour$ID)


#merge tumour radiomics with coronary artery radiomics
merged_df <- mask_tumour %>%
  left_join(aorta, by = "ID") %>% 
  left_join(calc_RCA, by = "ID") %>%
  left_join(calc_LMS_LAD, by = "ID") %>%
  left_join(calc_LCX, by = "ID") %>% 
  left_join(av, by = "ID") %>% 
  left_join(mv, by = "ID")

merged_df <- merged_df[,-(2:33) ]
merged_df <- merged_df[,-(4:5) ]

merged_df <- clinical_data %>% 
  inner_join(merged_df, by = "ID")

#adjust stage
for (i in 1:nrow(merged_df)) {
  if (is.na(merged_df$Stage[i])) {
    merged_df$Stage[i] <- 0
  } else if (merged_df$Stage[i] == "1a") {
    merged_df$Stage[i] <- 1
  } else if (merged_df$Stage[i] == "1b") {
    merged_df$Stage[i] <- 2
  } else if (merged_df$Stage[i] == "2a") {
    merged_df$Stage[i] <- 3
  } else if (merged_df$Stage[i] == "2b") {
    merged_df$Stage[i] <- 4
  } else if (merged_df$Stage[i] == "3") {
    merged_df$Stage[i] <- 5
  } else if (merged_df$Stage[i] == "4") {
    merged_df$Stage[i] <- 6
  } else {
    merged_df$Stage[i] <- 7
  } 
}
merged_df$Stage <- as.numeric(merged_df$Stage)

#create test and training data frames
set.seed(32)

#balance for agatston score
train.index <- createDataPartition(merged_df$Age | merged_df$Overall.survival..days. | merged_df$Stage | merged_df$Performance | merged_df$CAC |
                                     merged_df$all_calc, p = 2/3, list = FALSE) #age, stage, survival 
#train.index_old <- createDataPartition(merged_df$Age | merged_df$Overall.survival..days. | merged_df$Stage, p = 2/3, list = FALSE)
#train_df_old <- merged_df[ train.index_old,]
train_df <- merged_df[ train.index,]
test_df <- merged_df[-train.index,]


covariates <- c(colnames(merged_df))
covariates
covariates <- covariates[9:ncol(merged_df)]


#standardisation
train_sd <- train_df
for (i in 9:ncol(train_df)){
  if (sd((train_df[[covariates[i-9+1]]]), na.rm = TRUE) != 0){
    train_sd[[covariates[i-9+1]]]<- scale(train_df[[covariates[i-9+1]]])
  } else{
    train_sd[[covariates[i-9+1]]] <- 0
  }
}

train_sd[is.na(train_sd)] <- 0

test_sd <- test_df
for (i in 9:ncol(test_df)){
  if (sd((test_df[[covariates[i-9+1]]]), na.rm = TRUE) != 0){
    test_sd[[covariates[i-9+1]]]<- scale(test_df[[covariates[i-9+1]]])
  } else{
    test_sd[[covariates[i-9+1]]] <- 0
  }
}

test_sd[is.na(test_sd)] <- 0


one_level_cols <- which(sapply(train_sd, function(x) length(unique(x))) == 1)
one_level_cols
removed_one_levels <- train_sd[,-one_level_cols]

univ_train_sd <- train_sd[, -c(1,4:8)]

univ_covariates <- colnames(univ_train_sd)
univ_formulas <- sapply(univ_covariates,
                        function(x) as.formula(paste('Surv(Overall.survival..days., OS.event)~', x)))
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = univ_train_sd )})

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

res <- t(as.data.frame(univ_results, check.names = FALSE))
coxresult <- as.data.frame(res)
coxresult$FDR<-p.adjust(coxresult$p.value,method="fdr")
coxresult<-coxresult[order(coxresult$FDR, decreasing=F),]


#selection according to p.value 
filtered_coxresult <- coxresult[order(coxresult$p.value), ]

features <- rownames(filtered_coxresult)
filtered_coxresult <- cbind(features, filtered_coxresult)
rownames(filtered_coxresult) <- NULL

p_thres <- 0.05

coxresult_significant <- filtered_coxresult %>% 
  filter(p.value<p_thres)


#Lasso
train_features <- coxresult_significant$features

train_matching <- names(train_sd)[names(train_sd) %in% train_features]
train_subset <- subset(train_sd, select = train_matching)
selected_columns <- c("ID", "Overall.survival..days.", "OS.event", "diagnostics_Mask.interpolated_Minimum.y")
train_raw <- train_sd[selected_columns]

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
non_zero_features <- names(coefficients[coefficients != 0])
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

set.seed(16)
no_partition <- 10
k_mean_lasso <- cbind(my_prediction_model_train,forlasso)
colnames(k_mean_lasso)[1] <- "RPV"
km.res <- kmeans(k_mean_lasso, no_partition, nstart = 1)

k_mean_lasso$grouping <- km.res$cluster
k_mean_lasso$grouping

class(k_mean_lasso$grouping)

k_mean_lasso$grouping <- as.numeric(k_mean_lasso$grouping)

k_mean_lasso <- cbind(k_mean_lasso,train_raw$OS.event)
k_mean_lasso <- cbind(k_mean_lasso,train_raw$Overall.survival..days.)


colnames(k_mean_lasso)[524] <- "OS.event"
colnames(k_mean_lasso)[525] <- "Overall.survival..days."

for (i in 1:nrow(k_mean_lasso)) {
  if (k_mean_lasso$grouping[i] %in% c(1,3,4,5,6)) {
    k_mean_lasso$group_final[i] <- 1
  } else {
    k_mean_lasso$group_final[i] <- 2
  }
}
k_mean_lasso$group_final
k_mean_lasso$grouping

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
  title = "Kaplan-Meier Plot Stratified Based on RPV \n (Training set, Tumour + All calcification radiomics)",
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




cate_dis <- train_raw$OS.event
cate_dis

for (i in 1:length(cate_dis)) {
  if (cate_dis[i] == "1"){
    cate_dis[i] <- '+'
  } else {
    cate_dis[i] <- '-'
  }
}

library(ROCit)

score_dis <- unname(my_prediction_model_train)

score_dis

ROCit_obj <- rocit(score=score_dis, class=cate_dis, negref = "+", method ="bin")


AUC_obj <- ciAUC(ROCit_obj, level = 0.95)
p <- plot(ROCit_obj)
text(0.95, 0.2, paste0("AUC=", round(AUC_obj$AUC, 2), ", 95% CI [", round(AUC_obj$lower, 2), ",", round(AUC_obj$upper, 2), "]"), adj = 1, font = 4, cex=1.0)
title("Prediction of 3-year survival \n (Training Set - Tumour + All calcification)")




#testing

test_subset <- subset(test_sd, select = train_matching)
train_subset <- na.omit(train_subset)
removed_rows <- attr(train_subset, "na.action")
removed_rows

my_prediction_model_test <- predict(fit, newx = as.matrix(test_subset[,train_features]), cvfit$lambda.min)

set.seed(35)
k_mean_test <- cbind(my_prediction_model_test,test_subset)
no_partition <-  10
colnames(k_mean_test)[1] <- "RPV"
k_mean_test[is.na(k_mean_test)] <- 0
km.test <- kmeans(k_mean_test, no_partition, nstart = 1)


k_mean_test$grouping <- km.test$cluster
k_mean_test$grouping

k_mean_test$grouping <- as.numeric(k_mean_test$grouping)
selected_columns <- c("ID", "Overall.survival..days.", "OS.event", "diagnostics_Mask.interpolated_VolumeNum.y")
test_raw <- test_sd[selected_columns]

k_mean_test <- cbind(k_mean_test,test_raw$OS.event)
k_mean_test <- cbind(k_mean_test,test_raw$Overall.survival..days.)
colnames(k_mean_test)[524] <- "OS.event"
colnames(k_mean_test)[525] <- "Overall.survival..days."

for (i in 1:nrow(k_mean_test)){
  if (k_mean_test$grouping[i] %in% c(1,4,7,8,9)){
    k_mean_test$group_final[i] <- 1
  }
  else{
    k_mean_test$group_final[i] <- 2
  }
}

k_mean_test$group_final
km_fit_test <- survfit(Surv(Overall.survival..days.,OS.event) ~ group_final, data= k_mean_test)


b <- ggsurvplot(
  km_fit_test,     # survfit object with calculated statistics.
  data = k_mean_test,               # data used to fit survival curves. http://127.0.0.1:31261/graphics/plot_zoom_png?width=860&height=900
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
  title = "Kaplan-Meier Plot Stratified Based on RPV \n (Testing set, Tumour + All calcification radiomics)",
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

my_prediction_model_test

ROCit_obj_test <- rocit(score=score_test, class=cate_test, negref = "+", method ="bin")


AUC_obj_test <- ciAUC(ROCit_obj_test, level = 0.95)
q <- plot(ROCit_obj_test)
text(0.95, 0.2, paste0("AUC=", round(AUC_obj_test$AUC, 2), ", 95% CI [", round(AUC_obj_test$lower, 2), ",", round(AUC_obj_test$upper, 2), "]"), adj = 1, font = 4, cex=1.0)
title("Prediction of 3-year survival \n (Testing set, Tumour + all calc)")


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

write.csv(score_dis, "all_calc_radiomic_vector.csv", row.names = FALSE)
write.csv(score_test, "all_calc_radiomic_vector_2.csv", row.names = FALSE)

