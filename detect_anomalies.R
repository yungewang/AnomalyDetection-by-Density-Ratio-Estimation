detect_anomalies<- function(test_data, train_data, alpha = 0.1, sigmai = NULL, 
                            lambdai = 10^(seq(-3, 1, length.out = 9)), 
                            b = 50, fold = 5, 
                            threshold_quantile = 0.05, 
                            method = "pearson", penalty = "l2"
) {
  # alpha
  if (alpha < 0 || alpha >= 1 || length(alpha) != 1) {
    stop("Parameter alpha must be in [0, 1).")
  }
  
  # threshold_quantile
  if (threshold_quantile <= 0 || threshold_quantile >= 1 || length(threshold_quantile) > 1) {
    stop("Parameter threshold_quantile should be a scalar in (0, 1).")
  }
  
  # b
  if (b <= 0 || b %% 1 != 0) {
    stop("Parameter b must be a positive integer.")
  }
  
  # fold
  if (fold <= 0 || fold %% 1 != 0) {
    stop("Parameter fold must be a positive integer.")
  } 
  
  result <- bRuLSIF(test_data,train_data, alpha = alpha, sigmai = sigmai, lambdai = lambdai, b = b, fold = fold, method, penalty)
  
  # Get density ratios for test data
  density_ratios <- result$r(test_data)
  
  # Set a threshold based on quantile
  threshold <- quantile(density_ratios, probs = threshold_quantile)
  
  # Identify anomalies
  is_outlier <- density_ratios < threshold
  outlier_indices <- which(density_ratios < threshold)
  
  list(density_ratios = density_ratios, threshold = threshold, is_outlier=is_outlier, outlier_indices=outlier_indices)
}


