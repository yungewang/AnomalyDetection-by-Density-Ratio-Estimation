library(pROC)
library(R.matlab)
source("brulsif.R")
data <- readMat("benchmarks.mat")
banana <- data$banana
real_anomalies <- which(banana[[2]] == -1)  # Anomaly indices
positions_banana <- which(banana[[2]] == 1)  # Normal data indices
x_train <- banana[[1]][positions_banana, ]  # Using only normal data for training

# Full dataset for testing
x_test <- banana[[1]]
predicted_anomaly_indices <- list()
auc_values <- numeric(20)
set.seed(123)  
for (i in 1:20) {
  selected_dimensions <- sample(1:2, 2)
  x_train_subset <- x_train[, selected_dimensions]  
  x_test_subset <- x_test[, selected_dimensions]    
  detection_result <- detect_anomalies(
    test_data = x_test_subset,
    train_data = x_train_subset,
    alpha = 0,
    lambdai = 0,
    sigmai = NULL, 
    b = 50, 
    fold = 5, 
    threshold_quantile = 0.2, 
    method = "pearson", 
    penalty = "l2"
  )
  outlier_indices <- detection_result$outlier_indices
  predicted_anomaly_indices[[i]] <- outlier_indices
  predicted_labels <- ifelse(1:nrow(x_test) %in% outlier_indices, 1, 0)
  true_labels <- ifelse(1:nrow(x_test) %in% real_anomalies, 1, 0)
  roc_curve <- roc(true_labels, predicted_labels)
  auc_values[i] <- auc(roc_curve)  
  cat("Iteration", i, "completed with selected dimensions:", selected_dimensions, "\n")
}
avg_auc <- mean(auc_values)
cat("AUC values for each iteration:\n")
print(auc_values)
cat("Average AUC across all iterations:", avg_auc, "\n")


