# =============================================================================
# COVID-19 vs. Other Viral Infection Classification
# Cross-Validation: Logistic Regression vs. Generalized Additive Model (GAM)
# Author: Rhea Charles | University of South Florida
#
# Description:
#   Compares logistic regression and GAM for classifying SARS-CoV-2 (SC2)
#   infection vs. other viral infections using immune gene expression markers
#   (IL1B, IFI6, IL1R2) and demographic covariates.
#
#   Both models are evaluated with 10-fold cross-validation repeated 10 times.
#
# Input:
#   Week4_HW_data.csv — patient-level data with viral status, demographics,
#                        and gene expression values
#
# Output:
#   Confusion matrices and accuracy comparison for each model
#   Bar chart comparing model accuracies
# =============================================================================


# -----------------------------------------------------------------------------
# 0. Setup
# -----------------------------------------------------------------------------

library(caret)
library(mgcv)
library(ggplot2)
library(gam)

set.seed(42)

data_path <- "data/Week4_HW_data.csv"   # <-- update path as needed


# -----------------------------------------------------------------------------
# 1. Load and prepare data
# -----------------------------------------------------------------------------

df <- read.csv(data_path) |> na.omit()

# Binary outcome: SARS-CoV-2 vs. all other viruses
df$viral_status <- factor(ifelse(df$viral_status == "SC2", "SC2", "other virus"))
df$gender       <- factor(df$gender)


# -----------------------------------------------------------------------------
# 2. Model definitions
# -----------------------------------------------------------------------------
# Predictors: gender (demographic), age, and three immune response genes
#   IL1B  — interleukin-1 beta, pro-inflammatory cytokine
#   IFI6  — interferon alpha-inducible protein 6, antiviral response marker
#   IL1R2 — interleukin-1 receptor type 2, immune regulation

PREDICTORS <- "gender + age + IL1B + IFI6 + IL1R2"
OUTCOME    <- "viral_status"
THRESHOLD  <- 0.5


fit_logistic <- function(train, test) {
  # Logistic regression — assumes linear relationship between predictors and log-odds
  formula <- as.formula(paste(OUTCOME, "~", PREDICTORS))
  fit     <- glm(formula, data = train, family = "binomial")
  probs   <- predict(fit, newdata = test, type = "response")
  preds   <- ifelse(probs > THRESHOLD, "SC2", "other virus")
  data.frame(Prediction = preds, Reference = test[[OUTCOME]])
}


fit_gam <- function(train, test) {
  # GAM — allows non-linear smoothed terms for continuous predictors (df=4 spline)
  formula <- as.formula(paste(
    OUTCOME, "~ gender + s(age, df=4) + s(IL1B, df=4) + s(IFI6, df=4) + s(IL1R2, df=4)"
  ))
  fit   <- gam(formula, data = train, family = "binomial")
  probs <- predict(fit, newdata = test, type = "response")
  preds <- ifelse(probs > THRESHOLD, "SC2", "other virus")
  data.frame(Prediction = preds, Reference = test[[OUTCOME]])
}


# -----------------------------------------------------------------------------
# 3. Cross-validation — 10-fold repeated 10 times
# -----------------------------------------------------------------------------
# Repeated CV reduces variance in performance estimates compared to a single run.
# Each iteration uses a different random seed for fold assignment.

run_cv <- function(model_fn, data, k = 10, repeats = 10) {
  all_preds <- data.frame()

  for (i in seq_len(repeats)) {
    set.seed(i)
    folds <- createFolds(data[[OUTCOME]], k = k)
    train <- data[-folds[[i]], ]
    test  <- data[ folds[[i]], ]

    preds <- model_fn(train, test)

    # Align factor levels before binding
    lvls              <- levels(data[[OUTCOME]])
    preds$Prediction  <- factor(preds$Prediction, levels = lvls)
    preds$Reference   <- factor(preds$Reference,  levels = lvls)

    all_preds <- rbind(all_preds, preds)
  }

  all_preds
}

cat("Running logistic regression CV...\n")
logistic_preds <- run_cv(fit_logistic, df)

cat("Running GAM CV...\n")
gam_preds <- run_cv(fit_gam, df)


# -----------------------------------------------------------------------------
# 4. Evaluate performance
# -----------------------------------------------------------------------------

logistic_cm <- confusionMatrix(logistic_preds$Prediction, logistic_preds$Reference)
gam_cm      <- confusionMatrix(gam_preds$Prediction,      gam_preds$Reference)

cat("\n=== Logistic Regression ===\n")
print(logistic_cm)

cat("\n=== GAM ===\n")
print(gam_cm)

logistic_accuracy <- mean(logistic_preds$Prediction == logistic_preds$Reference)
gam_accuracy      <- mean(gam_preds$Prediction      == gam_preds$Reference)

cat("\n=== Model Comparison ===\n")
if (logistic_accuracy > gam_accuracy) {
  cat("Logistic regression performs better — accuracy:", round(logistic_accuracy, 4), "\n")
} else if (gam_accuracy > logistic_accuracy) {
  cat("GAM performs better — accuracy:", round(gam_accuracy, 4), "\n")
} else {
  cat("Both models tied — accuracy:", round(logistic_accuracy, 4), "\n")
}


# -----------------------------------------------------------------------------
# 5. Visualize accuracy comparison
# -----------------------------------------------------------------------------

accuracy_df <- data.frame(
  Model    = c("Logistic Regression", "GAM"),
  Accuracy = c(logistic_accuracy, gam_accuracy)
)

ggplot(accuracy_df, aes(x = Model, y = Accuracy, fill = Model)) +
  geom_bar(stat = "identity", width = 0.5, show.legend = FALSE) +
  geom_text(aes(label = round(Accuracy, 3)), vjust = -0.5, size = 4) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  labs(
    title    = "SC2 vs. Other Viral Infection: Model Accuracy Comparison",
    subtitle = "10-fold cross-validation, 10 repeats",
    x        = "Model",
    y        = "Accuracy"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, color = "gray50"))

ggsave("output/model_accuracy_comparison.png", width = 6, height = 5, dpi = 300)
