library(survival)
library(tidyverse)
library(readxl)
library(tidytidbits)
library(survivalAnalysis)
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

#remove unwanted data
tumour_filtered <- mask_tumour[, -c(2:33)]
mask_CTs_filtered <- mask_CTs[, -c(2:33)]
mask_CTp_filtered <- mask_CTp[, -c(2:33)]

tumour_filtered <- mask_tumour[, -c(4:5)]
mask_CTs_filtered <- mask_CTs[, -c(4:5)]
mask_CTp_filtered <- mask_CTp[, -c(4:5)]

#merge data and remove missing clinical data

merged_data <- merge(clinical_data,tumour_filtered,by = "ID")

merged_df <- merged_data %>%
  inner_join(mask_CTp_filtered, by = "ID") %>%
  inner_join(mask_CTs_filtered, by = "ID")


################
clinical_data <- read_excel("Cardiac_clinical.xlsx")
clinical_data <- as.data.frame(clinical_data)
clinical_data <- clinical_data[-100, ] 
clinical_data <- clinical_data[-1, ]
clinical_data <- clinical_data[, -1]
clinical_data$ID <- paste("LRAD", clinical_data$ID, sep = "")


covariate_names <- c(
  age = "Age at Dx",
  stage = "Stage",
  performance = "Performance"
)

coxph(Surv(Overall.survival..days.,OS.event) ~ Age + Gender + Performance, data = merged_df)
summary()




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
clinical_data <- subset(clinical_data, select = c(ID, OS.event, Overall.survival..days., Age, Stage, Gender, Performance))

merged_data <- merge(clinical_data,calc_masks,by = "ID")

#Coronary artery calcificaiton
coronary_calc <- merged_data[grepl("label_2|label_3|label_4", merged_data$Mask), ]

calc_RCA <- merged_data[grepl("label_2", merged_data$Mask), ]
calc_RCA <- calc_RCA[,-(2:38) ]
calc_RCA <- calc_RCA[,-(3:4) ]

calc_LMS_LAD <- merged_data[grepl("label_3", merged_data$Mask), ]
calc_LMS_LAD <- calc_LMS_LAD[,-(2:38) ]
calc_LMS_LAD <- calc_LMS_LAD[,-(3:4) ]

calc_LCX <- merged_data[grepl("label_4", merged_data$Mask), ]
calc_LCX <- calc_LCX[,-(2:38) ]
calc_LCX <- calc_LCX[,-(3:4) ]

aorta <- merged_data[grepl("label_1", merged_data$Mask), ]
aorta <- aorta[,-(2:38) ]
aorta <- aorta[,-(3:4) ]

av <- merged_data[grepl("label_5", merged_data$Mask), ]
av <- av[,-(2:38) ]
av <- av[,-(3:4) ]

mv <- merged_data[grepl("label_6", merged_data$Mask), ]
mv <- mv[,-(2:38) ]
mv <- mv[,-(3:4) ]

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

coxph(Surv(Overall.survival..days.,OS.event) ~ Age + Gender + Performance, data = merged_df)





################## CAC 
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
clinical_data <- subset(clinical_data, select = c(ID, OS.event, Overall.survival..days., Age, Stage, Gender, Performance))

merged_data <- merge(clinical_data,calc_masks,by = "ID")

#Coronary artery calcificaiton
coronary_calc <- merged_data[grepl("label_2|label_3|label_4", merged_data$Mask), ]

calc_RCA <- merged_data[grepl("label_2", merged_data$Mask), ]
calc_RCA <- calc_RCA[,-(2:38) ]
calc_RCA <- calc_RCA[,-(3:4) ]

calc_LMS_LAD <- merged_data[grepl("label_3", merged_data$Mask), ]
calc_LMS_LAD <- calc_LMS_LAD[,-(2:38) ]
calc_LMS_LAD <- calc_LMS_LAD[,-(3:4) ]

calc_LCX <- merged_data[grepl("label_4", merged_data$Mask), ]
calc_LCX <- calc_LCX[,-(2:38) ]
calc_LCX <- calc_LCX[,-(3:4) ]

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

coxph(Surv(Overall.survival..days.,OS.event) ~ Age + Gender + Performance, data = merged_df)



