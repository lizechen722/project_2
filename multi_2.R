library(gtsummary)
library(survival)
update.packages()
install.packages("xfun")
library(xfun)
xfun:::xfun_deps("yaml")
library(dplyr)

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
clinical_data <- subset(clinical_data, select = c(ID, OS.event, Overall.survival..days., Age, Gender, Stage, Performance, CAC, all_calc))
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

class(merged_df$Stage)

merged_df <- merged_data %>%
  inner_join(mask_CTp_filtered, by = "ID") %>%
  inner_join(mask_CTs_filtered, by = "ID")

merged_df$Stage <- as.numeric(merged_df$Stage)

# build survival model table
t1 <- coxph(Surv(Overall.survival..days.,OS.event) ~ Age + Performance + Gender + Stage + CAC + all_calc, data = merged_df) %>% 
  gtsummary::tbl_regression(exp = TRUE)



#ALL Calc + tumour

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
clinical_data <- subset(clinical_data, select = c(ID, OS.event, Overall.survival..days., Age, Gender,Stage, Performance, CAC, all_calc))

merged_data <- merge(clinical_data,calc_masks,by = "ID")
#Coronary artery calcificaiton
coronary_calc <- merged_data[grepl("label_2|label_3|label_4", merged_data$Mask), ]

calc_RCA <- merged_data[grepl("label_2", merged_data$Mask), ]
calc_RCA <- calc_RCA[,-(2:41) ]
calc_RCA <- calc_RCA[,-(4:5) ]

calc_LMS_LAD <- merged_data[grepl("label_3", merged_data$Mask), ]
calc_LMS_LAD <- calc_LMS_LAD[,-(2:41) ]
calc_LMS_LAD <- calc_LMS_LAD[,-(4:5) ]

calc_LCX <- merged_data[grepl("label_4", merged_data$Mask), ]
calc_LCX <- calc_LCX[,-(2:41) ]
calc_LCX <- calc_LCX[,-(4:5) ]

aorta <- merged_data[grepl("label_1", merged_data$Mask), ]
aorta <- aorta[,-(2:41) ]
aorta <- aorta[,-(4:5) ]

av <- merged_data[grepl("label_5", merged_data$Mask), ]
av <- av[,-(2:41) ]
av <- av[,-(4:5) ]

mv <- merged_data[grepl("label_6", merged_data$Mask), ]
mv <- mv[,-(2:41) ]
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

# build survival model table
t2 <- coxph(Surv(Overall.survival..days.,OS.event) ~ Age + Performance + Gender + Stage + CAC + all_calc, data = merged_df) %>% 
  gtsummary::tbl_regression(exp = TRUE)



###All Calc 
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
clinical_data <- subset(clinical_data, select = c(ID, OS.event, Overall.survival..days., Age, Gender,Stage, Performance, CAC, all_calc))

merged_data <- merge(clinical_data,calc_masks,by = "ID")
#Coronary artery calcificaiton
coronary_calc <- merged_data[grepl("label_2|label_3|label_4", merged_data$Mask), ]

calc_RCA <- merged_data[grepl("label_2", merged_data$Mask), ]
calc_RCA <- calc_RCA[,-(2:41) ]
calc_RCA <- calc_RCA[,-(4:5) ]

calc_LMS_LAD <- merged_data[grepl("label_3", merged_data$Mask), ]
calc_LMS_LAD <- calc_LMS_LAD[,-(2:41) ]
calc_LMS_LAD <- calc_LMS_LAD[,-(4:5) ]

calc_LCX <- merged_data[grepl("label_4", merged_data$Mask), ]
calc_LCX <- calc_LCX[,-(2:41) ]
calc_LCX <- calc_LCX[,-(4:5) ]

aorta <- merged_data[grepl("label_1", merged_data$Mask), ]
aorta <- aorta[,-(2:41) ]
aorta <- aorta[,-(4:5) ]

av <- merged_data[grepl("label_5", merged_data$Mask), ]
av <- av[,-(2:41) ]
av <- av[,-(4:5) ]

mv <- merged_data[grepl("label_6", merged_data$Mask), ]
mv <- mv[,-(2:41) ]
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


merged_df <- merged_df[,-(2:1077) ]


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

# build survival model table
coxph(Surv(Overall.survival..days.,OS.event) ~ Age + Performance + Gender + Stage + CAC + all_calc, data = merged_df) %>% 
  gtsummary::tbl_regression(exp = TRUE)




#coronary calc
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
clinical_data <- subset(clinical_data, select = c(ID, OS.event, Overall.survival..days., Age, Stage, Gender, Performance, CAC, all_calc, Agatston_Score))

merged_data <- merge(clinical_data,calc_masks,by = "ID")

#Coronary artery calcificaiton
coronary_calc <- merged_data[grepl("label_2|label_3|label_4", merged_data$Mask), ]

calc_RCA <- merged_data[grepl("label_2", merged_data$Mask), ]
calc_RCA <- calc_RCA[,-(2:42) ]
calc_RCA <- calc_RCA[,-(4:5) ]

calc_LMS_LAD <- merged_data[grepl("label_3", merged_data$Mask), ]
calc_LMS_LAD <- calc_LMS_LAD[,-(2:42) ]
calc_LMS_LAD <- calc_LMS_LAD[,-(4:5) ]

calc_LCX <- merged_data[grepl("label_4", merged_data$Mask), ]
calc_LCX <- calc_LCX[,-(2:42) ]
calc_LCX <- calc_LCX[,-(4:5) ]

#tumour mask
mask_tumour <- data[grepl("CT\\.nii\\.gz$", data$Mask), ]

#adjust ID
mask_tumour$ID <- gsub("_CT", "", mask_tumour$ID)


#merge tumour radiomics with coronary artery radiomics
merged_df <- mask_tumour %>%
  left_join(calc_RCA, by = "ID") %>%
  left_join(calc_LMS_LAD, by = "ID") %>%
  left_join(calc_LCX, by = "ID")

merged_df <- merged_df[,-(2:1077) ]


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

# build survival model table
coxph(Surv(Overall.survival..days.,OS.event) ~ Age + Performance + Gender + Stage + CAC + all_calc + Agatston_Score, data = merged_df) %>% 
  gtsummary::tbl_regression(exp = TRUE)



####CAC with tumour

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
clinical_data <- subset(clinical_data, select = c(ID, OS.event, Overall.survival..days., Age, Stage, Gender, Performance, CAC, all_calc))

merged_data <- merge(clinical_data,calc_masks,by = "ID")

#Coronary artery calcificaiton
coronary_calc <- merged_data[grepl("label_2|label_3|label_4", merged_data$Mask), ]

calc_RCA <- merged_data[grepl("label_2", merged_data$Mask), ]
calc_RCA <- calc_RCA[,-(2:41) ]
calc_RCA <- calc_RCA[,-(4:5) ]

calc_LMS_LAD <- merged_data[grepl("label_3", merged_data$Mask), ]
calc_LMS_LAD <- calc_LMS_LAD[,-(2:41) ]
calc_LMS_LAD <- calc_LMS_LAD[,-(4:5) ]

calc_LCX <- merged_data[grepl("label_4", merged_data$Mask), ]
calc_LCX <- calc_LCX[,-(2:41) ]
calc_LCX <- calc_LCX[,-(4:5) ]

#tumour mask
mask_tumour <- data[grepl("CT\\.nii\\.gz$", data$Mask), ]

#adjust ID
mask_tumour$ID <- gsub("_CT", "", mask_tumour$ID)


#merge tumour radiomics with coronary artery radiomics
merged_df <- mask_tumour %>%
  left_join(calc_RCA, by = "ID") %>%
  left_join(calc_LMS_LAD, by = "ID") %>%
  left_join(calc_LCX, by = "ID")



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

coxph(Surv(Overall.survival..days.,OS.event) ~ Age + Performance + Gender + Stage + CAC + all_calc, data = merged_df) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# merge tables
tbl_merge(
  tbls = list(t1, t2, t3, t4),
  tab_spanner = c("**Tumour only**", "**Tumour + All Calc**", "**All Calc only**", "**CAC**")
)
