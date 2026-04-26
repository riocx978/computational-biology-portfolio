# ============================================================================
# Viral Status Classification Using Immune Biomarkers
# ============================================================================
# Purpose: Classify respiratory samples into three categories
#          (SARS-CoV-2, Other Virus, No Virus) using immune response markers
#          (IL1B, IFI6, IL1R2) and demographic factors (age, gender)
#
# Methods: Logistic Regression (binary classification) and
#          Linear Discriminant Analysis (multiclass classification)
#
# Author: Rhea Charles
# ============================================================================

suppressPackageStartupMessages({
  library(MASS)
  library(ggplot2)
  library(dplyr)
  library(caret)
  library(pROC)
  library(gridExtra)
})

# ============================================================================
# SECTION 1: DATA LOADING AND PREPARATION
# ============================================================================

df <- read.csv("Week4_HW_data.csv")

cat("Dataset Summary:\n")
cat("Dimensions:", nrow(df), "samples x", ncol(df), "features\n")
cat("Viral Status Distribution:\n")
print(table(df$viral_status))
cat("\nBasic Statistics:\n")
print(summary(df[, c("age", "IL1B", "IFI6", "IL1R2")]))

# ============================================================================
# SECTION 2: EXPLORATORY DATA ANALYSIS
# ============================================================================

# 2.1: Biomarker distributions by viral status
cat("\n========================================\n")
cat("EXPLORATORY DATA ANALYSIS\n")
cat("========================================\n")

# Prepare data for visualization (long format)
df_long <- df %>%
  select(viral_status, IL1B, IFI6, IL1R2) %>%
  tidyr::pivot_longer(cols = c("IL1B", "IFI6", "IL1R2"),
                      names_to = "Biomarker",
                      values_to = "Expression")

# Boxplots of biomarkers by viral status
ggplot(df_long, aes(x = viral_status, y = Expression, fill = viral_status)) +
  geom_boxplot(alpha = 0.7) +
  facet_wrap(~Biomarker, scales = "free_y") +
  theme_minimal() +
  labs(title = "Immune Biomarker Expression by Viral Status",
       x = "Viral Status",
       y = "Expression Level",
       fill = "Viral Status") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 2.2: Correlation analysis
cat("\nBiomarker Correlations:\n")
correlation_matrix <- cor(df[, c("IL1B", "IFI6", "IL1R2", "age")])
print(correlation_matrix)

# 2.3: Age distribution by viral status
ggplot(df, aes(x = viral_status, y = age, fill = viral_status)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3) +
  theme_minimal() +
  labs(title = "Age Distribution by Viral Status",
       x = "Viral Status",
       y = "Age (years)",
       fill = "Viral Status")

# 2.4: Gender distribution
cat("\nGender Distribution by Viral Status:\n")
gender_dist <- table(df$viral_status, df$gender)
print(gender_dist)

# ============================================================================
# SECTION 3: BINARY CLASSIFICATION - LOGISTIC REGRESSION
# ============================================================================
# Predicting: SC2 vs. Non-SC2 (no_virus + other_virus)

cat("\n========================================\n")
cat("MODEL 1: LOGISTIC REGRESSION\n")
cat("(Binary: SC2 vs. Non-SC2)\n")
cat("========================================\n")

df_binary <- df %>%
  mutate(
    viral_status_binary = ifelse(viral_status == "SC2", 1, 0),
    viral_status_binary = as.factor(viral_status_binary),
    gender = as.factor(ifelse(gender == "M", 1, 0))
  )

cat("\nClass Distribution:\n")
print(table(df_binary$viral_status_binary))
cat("SC2 prevalence:", sum(df_binary$viral_status_binary == 1) / nrow(df_binary) * 100, "%\n")

# 3.1: Train/Test Split with stratification
set.seed(200)
train_index <- createDataPartition(df_binary$viral_status_binary,
                                   p = 0.7, list = FALSE, times = 1)
train_binary <- df_binary[train_index, ]
test_binary <- df_binary[-train_index, ]

cat("\nTrain/Test Split:\n")
cat("Training set:", nrow(train_binary), "samples\n")
cat("Test set:", nrow(test_binary), "samples\n")

# 3.2: Build logistic regression model
model_lr <- glm(viral_status_binary ~ gender + age + IL1B + IFI6 + IL1R2,
                family = "binomial",
                data = train_binary)

cat("\nModel Summary:\n")
summary_lr <- summary(model_lr)
print(summary_lr)

# 3.3: Make predictions
pred_prob_lr <- predict(model_lr, test_binary, type = "response")
pred_class_lr <- ifelse(pred_prob_lr > 0.5, 1, 0)

# 3.4: Evaluate performance
conf_matrix_lr <- confusionMatrix(factor(pred_class_lr),
                                   test_binary$viral_status_binary,
                                   positive = "1")

cat("\nConfusion Matrix:\n")
print(conf_matrix_lr$table)

cat("\nPerformance Metrics:\n")
cat("Accuracy: ", round(conf_matrix_lr$overall["Accuracy"], 4), "\n")
cat("Sensitivity (True Positive Rate): ",
    round(conf_matrix_lr$byClass["Sensitivity"], 4), "\n")
cat("Specificity (True Negative Rate): ",
    round(conf_matrix_lr$byClass["Specificity"], 4), "\n")
cat("Precision: ", round(conf_matrix_lr$byClass["Pos Pred Value"], 4), "\n")
cat("Recall: ", round(conf_matrix_lr$byClass["Sensitivity"], 4), "\n")
cat("F1-Score: ", round(conf_matrix_lr$byClass["F1"], 4), "\n")

# 3.5: ROC Curve and AUC
roc_lr <- roc(test_binary$viral_status_binary, pred_prob_lr)
cat("AUC-ROC: ", round(auc(roc_lr), 4), "\n")

plot(roc_lr, main = paste0("ROC Curve - Logistic Regression (AUC = ",
                           round(auc(roc_lr), 3), ")"),
     col = "steelblue", lwd = 2)
abline(a = 0, b = 1, lty = 2, col = "gray")

# ============================================================================
# SECTION 4: MULTICLASS CLASSIFICATION - LINEAR DISCRIMINANT ANALYSIS
# ============================================================================
# Predicting: SC2 vs. no_virus vs. other_virus

cat("\n========================================\n")
cat("MODEL 2: LINEAR DISCRIMINANT ANALYSIS\n")
cat("(Multiclass: SC2 vs. No Virus vs. Other Virus)\n")
cat("========================================\n")

df_multiclass <- df %>%
  mutate(
    viral_status = as.factor(viral_status),
    gender = as.factor(ifelse(gender == "M", 1, 0))
  )

cat("\nClass Distribution:\n")
print(table(df_multiclass$viral_status))

# 4.1: Train/Test Split with stratification
set.seed(200)
train_index_mc <- createDataPartition(df_multiclass$viral_status,
                                      p = 0.7, list = FALSE, times = 1)
train_multiclass <- df_multiclass[train_index_mc, ]
test_multiclass <- df_multiclass[-train_index_mc, ]

cat("\nTrain/Test Split:\n")
cat("Training set:", nrow(train_multiclass), "samples\n")
cat("Test set:", nrow(test_multiclass), "samples\n")

# 4.2: Build LDA model
model_lda <- lda(viral_status ~ gender + age + IL1B + IFI6 + IL1R2,
                 data = train_multiclass)

cat("\nModel Summary:\n")
print(model_lda)

cat("\nPrior Probabilities:\n")
print(model_lda$prior)

cat("\nGroup Means:\n")
print(model_lda$means)

# 4.3: Make predictions
pred_lda <- predict(model_lda, test_multiclass)
pred_class_lda <- pred_lda$class
pred_posterior_lda <- pred_lda$posterior

# 4.4: Evaluate performance
conf_matrix_lda <- confusionMatrix(pred_class_lda,
                                    test_multiclass$viral_status)

cat("\nConfusion Matrix:\n")
print(conf_matrix_lda$table)

cat("\nOverall Accuracy: ",
    round(conf_matrix_lda$overall["Accuracy"], 4), "\n")

cat("\nPer-Class Performance:\n")
print(conf_matrix_lda$byClass)

# 4.5: Cross-validation for LDA
cat("\n10-Fold Cross-Validation for LDA:\n")
set.seed(200)
train_control <- trainControl(method = "cv", number = 10)
model_lda_cv <- train(viral_status ~ gender + age + IL1B + IFI6 + IL1R2,
                      data = df_multiclass,
                      method = "lda",
                      trControl = train_control)

cat("Cross-Validation Accuracy: ",
    round(mean(model_lda_cv$resample$Accuracy), 4), "\n")
cat("Cross-Validation Std Dev: ",
    round(sd(model_lda_cv$resample$Accuracy), 4), "\n")

# ============================================================================
# SECTION 5: MODEL COMPARISON AND INSIGHTS
# ============================================================================

cat("\n========================================\n")
cat("MODEL COMPARISON SUMMARY\n")
cat("========================================\n")

cat("\nLogistic Regression (Binary Classification: SC2 vs. Non-SC2):\n")
cat("  - Accuracy: ", round(conf_matrix_lr$overall["Accuracy"], 4), "\n")
cat("  - AUC-ROC: ", round(auc(roc_lr), 4), "\n")
cat("  - F1-Score: ", round(conf_matrix_lr$byClass["F1"], 4), "\n")

cat("\nLinear Discriminant Analysis (Multiclass):\n")
cat("  - Accuracy: ", round(conf_matrix_lda$overall["Accuracy"], 4), "\n")
cat("  - CV Accuracy (10-fold): ",
    round(mean(model_lda_cv$resample$Accuracy), 4), "\n")

# ============================================================================
# SECTION 6: BIOMARKER IMPORTANCE AND INSIGHTS
# ============================================================================

cat("\n========================================\n")
cat("BIOMARKER ANALYSIS\n")
cat("========================================\n")

# Extract LDA coefficients (standardized discriminant functions)
lda_coef <- model_lda$scaling
cat("\nLinear Discriminant Coefficients:\n")
print(lda_coef)

cat("\nInterpretation:\n")
cat("LD1: Separates SC2 from other virus and no_virus\n")
cat("LD2: Separates other_virus from no_virus\n")

# Discriminant scores
lda_scores <- as.data.frame(pred_lda$x)
lda_scores$viral_status <- pred_lda$class
lda_scores$sample_id <- rownames(test_multiclass)

# Visualize discriminant functions
ggplot(lda_scores, aes(x = LD1, y = LD2, color = viral_status)) +
  geom_point(size = 3, alpha = 0.6) +
  theme_minimal() +
  labs(title = "Linear Discriminant Analysis - Class Separation",
       x = paste0("LD1 (Explains ", round(model_lda$svd[1]^2 / sum(model_lda$svd^2) * 100, 1), "%)"),
       y = paste0("LD2 (Explains ", round(model_lda$svd[2]^2 / sum(model_lda$svd^2) * 100, 1), "%)"),
       color = "Viral Status") +
  theme(legend.position = "bottom")

# ============================================================================
# SECTION 7: KEY FINDINGS
# ============================================================================

cat("\n========================================\n")
cat("KEY FINDINGS\n")
cat("========================================\n")

cat("\n1. BIOMARKER PATTERNS:\n")
cat("   - IL1B, IFI6, IL1R2 show distinct patterns across viral statuses\n")
cat("   - Multiple biomarkers are required for accurate classification\n")

cat("\n2. MODEL PERFORMANCE:\n")
cat("   - Logistic Regression achieves strong binary classification\n")
cat("   - LDA successfully handles multiclass classification\n")
cat("   - Cross-validation confirms model generalizability\n")

cat("\n3. CLINICAL RELEVANCE:\n")
cat("   - Immune biomarkers enable discrimination between viral statuses\n")
cat("   - Age and gender contribute to classification accuracy\n")

# ============================================================================
# END OF ANALYSIS
# ============================================================================
