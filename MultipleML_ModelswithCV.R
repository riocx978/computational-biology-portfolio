# ============================================================================
# Heart Disease Prediction: Comparative ML Models
# ============================================================================
# Purpose: Compare multiple classification models (Logistic Regression,
#          Naive Bayes, LASSO, Random Forest, Gradient Boosting) on the
#          Heart Disease dataset
# Author: Rhea Charles
# ============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(MLDataR)
  library(dplyr)
  library(tidyr)
  library(tidymodels)
  library(data.table)
  library(ConfusionTableR)
  library(OddsPlotty)
  library(ROCR)
  library(ggplot2)
  library(glmnet)
  library(caret)
  library(randomForest)
  library(caTools)
  library(gbm)
  library(pROC)
})

# ============================================================================
# SECTION 1: DATA LOADING AND PREPROCESSING
# ============================================================================

data(package = "MLDataR")
df <- heartdisease

# Encode categorical variables as binary
df$Sex <- ifelse(df$Sex == "F", 1, 0)
df$RestingECG <- ifelse(df$RestingECG == "ST", 1, 0)
df$Angina <- ifelse(df$Angina == "Y", 1, 0)

# ============================================================================
# SECTION 2: EXPLORATORY DATA ANALYSIS
# ============================================================================

# 2.1: Distribution of key continuous variables by Heart Disease status
p1 <- ggplot(df, aes(x = factor(HeartDisease), y = MaxHR,
                     fill = factor(HeartDisease))) +
  geom_boxplot() +
  labs(x = "Heart Disease Status", y = "Maximum Heart Rate (MaxHR)",
       fill = "Heart Disease") +
  theme_minimal() +
  theme(legend.position = "bottom")

p2 <- ggplot(df, aes(x = factor(HeartDisease), y = HeartPeakReading,
                     fill = factor(HeartDisease))) +
  geom_boxplot() +
  labs(x = "Heart Disease Status", y = "Heart Peak Reading",
       fill = "Heart Disease") +
  theme_minimal() +
  theme(legend.position = "bottom")

# 2.2: Feature correlation analysis
df_selected <- select(df, MaxHR, Angina, HeartPeakReading, HeartDisease)
correlation_matrix <- cor(df_selected)

# Visualize correlation matrix as heatmap
correlation_df <- as.data.frame(as.table(correlation_matrix))
colnames(correlation_df) <- c("Var1", "Var2", "Correlation")

ggplot(correlation_df, aes(x = Var1, y = Var2, fill = Correlation)) +
  geom_tile() +
  geom_text(aes(label = round(Correlation, 2)), color = "black", size = 3) +
  scale_fill_gradient(low = "orange", high = "red") +
  theme_minimal() +
  labs(title = "Feature Correlation Matrix",
       x = "Variables", y = "Variables") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

cat("\nFeature correlations with HeartDisease:\n")
print(cor(df)[, "HeartDisease"])

# ============================================================================
# SECTION 3: MODEL 1 - LOGISTIC REGRESSION
# ============================================================================

cat("\n========================================\n")
cat("MODEL 1: LOGISTIC REGRESSION\n")
cat("========================================\n")

df_lr <- heartdisease
n <- nrow(df_lr)
set.seed(200)
ntest <- trunc(n / 3)
testid <- sample(1:n, ntest)

# Fit logistic regression model
model_lr <- glm(HeartDisease ~ ., data = df_lr[-testid, ], family = binomial)

# Make predictions
lpred <- predict(model_lr, newdata = df_lr[testid, ], type = "response")
pred <- ifelse(lpred > 0.5, 1, 0)

# Calculate metrics
accuracy_lr <- mean(pred == df_lr$HeartDisease[testid])
sensitivity_lr <- sum(pred == 1 & df_lr$HeartDisease[testid] == 1) /
                  sum(df_lr$HeartDisease[testid] == 1)
specificity_lr <- sum(pred == 0 & df_lr$HeartDisease[testid] == 0) /
                  sum(df_lr$HeartDisease[testid] == 0)
precision_lr <- sum(pred == 1 & df_lr$HeartDisease[testid] == 1) /
                sum(pred == 1)
recall_lr <- sensitivity_lr
f1_lr <- 2 * (precision_lr * recall_lr) / (precision_lr + recall_lr)

# Print results
cat("Logistic Regression - Model Summary:\n")
cat("Accuracy:  ", round(accuracy_lr, 4), "\n")
cat("Sensitivity: ", round(sensitivity_lr, 4), "\n")
cat("Specificity: ", round(specificity_lr, 4), "\n")
cat("Precision:   ", round(precision_lr, 4), "\n")
cat("Recall:      ", round(recall_lr, 4), "\n")
cat("F1-Score:    ", round(f1_lr, 4), "\n")

# Confusion matrix
conf_matrix_lr <- table(Predicted = pred, Actual = df_lr$HeartDisease[testid])
cat("\nConfusion Matrix:\n")
print(conf_matrix_lr)

# ============================================================================
# SECTION 4: MODEL 2 - NAIVE BAYES
# ============================================================================

cat("\n========================================\n")
cat("MODEL 2: NAIVE BAYES\n")
cat("========================================\n")

df_nb <- na.omit(heartdisease)
df_nb$HeartDisease <- as.factor(df_nb$HeartDisease)

set.seed(100)
dt_nb <- sort(sample(nrow(df_nb), nrow(df_nb) * 0.7))
train_nb <- df_nb[dt_nb, ]
test_nb <- df_nb[-dt_nb, ]

# Train Naive Bayes model with 10-fold cross-validation
trctrl <- trainControl(method = "cv", number = 10, savePredictions = TRUE)
model_nb <- train(HeartDisease ~ .,
                  data = df_nb,
                  method = "naive_bayes",
                  trControl = trctrl,
                  tuneLength = 0)

# Make predictions
pred_nb <- predict(model_nb, test_nb)
conf_matrix_nb <- confusionMatrix(test_nb$HeartDisease, pred_nb)
accuracy_nb <- mean(pred_nb == test_nb$HeartDisease)

cat("Naive Bayes - Model Performance:\n")
cat("Accuracy: ", round(accuracy_nb, 4), "\n")
cat("\nConfusion Matrix:\n")
print(conf_matrix_nb$table)
cat("\nDetailed Metrics:\n")
print(conf_matrix_nb$byClass)

# ============================================================================
# SECTION 5: MODEL 3 - LASSO REGULARIZATION
# ============================================================================

cat("\n========================================\n")
cat("MODEL 3: LASSO REGULARIZATION\n")
cat("========================================\n")

df_lasso <- heartdisease
x <- model.matrix(HeartDisease ~ ., data = df_lasso)[, -1]
y <- df_lasso$HeartDisease

# Fit LASSO model with cross-validation
fit_lasso <- cv.glmnet(x, y, alpha = 1, family = "binomial")

# Extract best model
best_lambda <- fit_lasso$lambda.min
best_model_lasso <- glmnet(x, y, alpha = 1, lambda = best_lambda,
                           family = "binomial")

cat("LASSO - Selected Features (non-zero coefficients):\n")
lasso_coefs <- coef(best_model_lasso)
print(lasso_coefs[lasso_coefs[, 1] != 0, , drop = FALSE])

# Perform 10-fold cross-validation
set.seed(123)
folds <- sample(rep(1:10, length = nrow(df_lasso)))
test_errors <- numeric(10)

for (k in 1:10) {
  x_train <- x[folds != k, ]
  y_train <- y[folds != k]
  x_test <- x[folds == k, ]
  y_test <- y[folds == k]

  # Fit LASSO on training fold
  cv_lasso <- cv.glmnet(x_train, y_train, alpha = 1, family = "binomial")

  # Predict on test fold
  pred_lasso <- predict(cv_lasso, newx = x_test, s = "lambda.min",
                        type = "response")
  pred_binary <- ifelse(pred_lasso > 0.5, 1, 0)

  # Store test error
  test_errors[k] <- mean(pred_binary != y_test)
}

mean_test_error <- mean(test_errors)
cat("Mean Cross-Validation Error: ", round(mean_test_error, 4), "\n")
cat("Mean Accuracy: ", round(1 - mean_test_error, 4), "\n")

# ============================================================================
# SECTION 6: MODEL 4 - RANDOM FOREST
# ============================================================================

cat("\n========================================\n")
cat("MODEL 4: RANDOM FOREST\n")
cat("========================================\n")

df_rf <- heartdisease
df_rf$HeartDisease <- as.factor(df_rf$HeartDisease)
df_rf$Sex <- as.factor(df_rf$Sex)
df_rf$RestingECG <- as.factor(df_rf$RestingECG)
df_rf$Angina <- as.factor(df_rf$Angina)
df_rf <- na.omit(df_rf)

set.seed(222)
dt_rf <- sort(sample(nrow(df_rf), nrow(df_rf) * 0.7))
train_rf <- df_rf[dt_rf, ]
test_rf <- df_rf[-dt_rf, ]

# Train Random Forest
model_rf <- randomForest(HeartDisease ~ ., data = train_rf, proximity = TRUE)

# Predictions
pred_train_rf <- predict(model_rf, train_rf)
pred_test_rf <- predict(model_rf, test_rf)

cat("Random Forest - Training Set Performance:\n")
conf_train_rf <- confusionMatrix(pred_train_rf, train_rf$HeartDisease)
print(conf_train_rf$table)

cat("\nRandom Forest - Test Set Performance:\n")
conf_test_rf <- confusionMatrix(pred_test_rf, test_rf$HeartDisease)
print(conf_test_rf$table)

cat("\nTest Set Metrics:\n")
print(conf_test_rf$byClass)

# Feature importance
cat("\nTop 10 Important Features:\n")
importance_rf <- importance(model_rf)
print(head(importance_rf[order(importance_rf[, 1], decreasing = TRUE), ], 10))

# ============================================================================
# SECTION 7: MODEL 5 - GRADIENT BOOSTING
# ============================================================================

cat("\n========================================\n")
cat("MODEL 5: GRADIENT BOOSTING\n")
cat("========================================\n")

df_gb <- heartdisease
df_gb$HeartDisease <- as.numeric(df_gb$HeartDisease)
df_gb$Sex <- as.factor(df_gb$Sex)
df_gb$RestingECG <- as.factor(df_gb$RestingECG)
df_gb$Angina <- as.factor(df_gb$Angina)

set.seed(123)
dt_gb <- sample(nrow(df_gb), 0.7 * nrow(df_gb))
train_gb <- df_gb[dt_gb, ]
test_gb <- df_gb[-dt_gb, ]

# Train Gradient Boosting model
model_gb <- gbm(HeartDisease ~ .,
                data = train_gb,
                distribution = "bernoulli",
                n.trees = 10000,
                shrinkage = 0.01,
                interaction.depth = 4)

# Predictions
pred_gb <- predict(model_gb, newdata = test_gb, n.trees = 10000,
                   type = "response")
pred_binary_gb <- ifelse(pred_gb > 0.5, 1, 0)

# Performance metrics
accuracy_gb <- mean(pred_binary_gb == test_gb$HeartDisease)
rocobj_gb <- roc(test_gb$HeartDisease, pred_gb)
auc_gb <- round(auc(rocobj_gb), 4)

conf_matrix_gb <- confusionMatrix(factor(pred_binary_gb),
                                   factor(test_gb$HeartDisease))

cat("Gradient Boosting - Performance Metrics:\n")
cat("Accuracy: ", round(accuracy_gb, 4), "\n")
cat("AUC-ROC: ", auc_gb, "\n")
cat("Precision: ", round(conf_matrix_gb$byClass["Pos Pred Value"], 4), "\n")
cat("Recall: ", round(conf_matrix_gb$byClass["Sensitivity"], 4), "\n")
cat("F1-Score: ", round(conf_matrix_gb$byClass["F1"], 4), "\n")

# Feature importance
cat("\nTop 10 Important Features:\n")
summary_gb <- summary(model_gb, plotit = FALSE)
print(head(summary_gb, 10))

# ============================================================================
# SECTION 8: MODEL COMPARISON
# ============================================================================

cat("\n========================================\n")
cat("MODEL COMPARISON SUMMARY\n")
cat("========================================\n")

comparison_df <- data.frame(
  Model = c("Logistic Regression", "Naive Bayes", "LASSO",
            "Random Forest", "Gradient Boosting"),
  Accuracy = c(
    round(accuracy_lr, 4),
    round(accuracy_nb, 4),
    round(1 - mean_test_error, 4),
    round(sum(diag(conf_test_rf$table)) / sum(conf_test_rf$table), 4),
    round(accuracy_gb, 4)
  ),
  Sensitivity = c(
    round(sensitivity_lr, 4),
    round(conf_matrix_nb$byClass["Sensitivity"], 4),
    NA,
    round(conf_test_rf$byClass["Sensitivity"], 4),
    round(conf_matrix_gb$byClass["Sensitivity"], 4)
  ),
  Specificity = c(
    round(specificity_lr, 4),
    round(conf_matrix_nb$byClass["Specificity"], 4),
    NA,
    round(conf_test_rf$byClass["Specificity"], 4),
    NA
  ),
  F1_Score = c(
    round(f1_lr, 4),
    round(conf_matrix_nb$byClass["F1"], 4),
    NA,
    round(conf_test_rf$byClass["F1"], 4),
    round(conf_matrix_gb$byClass["F1"], 4)
  )
)

print(comparison_df)

# ============================================================================
# END OF ANALYSIS
# ============================================================================
