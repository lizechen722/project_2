setwd("C:\\Users\\lizec\\Desktop\\CI_Project")
clinical_data <- read_excel("Cardiac_clinical.xlsx")
clinical_data <- as.data.frame(clinical_data)
clinical_data <- clinical_data[, -1]
clinical_data$ID <- paste("LRAD", clinical_data$ID, sep = "")


data <- read.csv("radiomic_vector.csv")
data <- as.data.frame(data)

filtered_ids <- data$X[data$CAC != 0]
clinical_data$CAC <- ifelse(clinical_data$ID %in% filtered_ids, 0, 1)

filtered_ids_all_calc <- data$X[data$all_calc != 0]
clinical_data$all_calc <- ifelse(clinical_data$ID %in% filtered_ids_all_calc, 0, 1)

write.csv(clinical_data, "clinical_data_updated.csv", row.names = FALSE)


install.packages("psych")
library(psych)

age_summary <- summarySE(clinical_data, measurevar = "Age", groupvars = NULL)

age_summary <- clinical_data %>%
  summarise(mean_age = mean(Age, na.rm = TRUE),
            se_age = sd(Age, na.rm = TRUE) / sqrt(n()))

print(age_summary)

clinical_data_cleaned <- clinical_data[!is.na(clinical_data$OS.event), ]



gender_counts <- table(clinical_data_cleaned$Gender)
print(gender_counts)

mean_age <- mean(clinical_data_cleaned$Age, na.rm = TRUE)
sd_age <- sd(clinical_data_cleaned$Age, na.rm = TRUE)

cat("Mean Age:", mean_age, "\n")
cat("Standard Deviation:", sd_age, "\n")


count_younger_than_70 <- sum(clinical_data_cleaned$Age < 70, na.rm = TRUE)
count_younger_than_70

sum(clinical_data_cleaned$Histology == 3, na.rm = TRUE)


class(clinical_data_cleaned$Stage)

sum(grepl("3", clinical_data$Performance, fixed = TRUE), na.rm = TRUE)



sum(rowSums(clinical_data[c("CAC", "all_calc")] == 0) > 0)





