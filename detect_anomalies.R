detect_anomalies<- function(test_set, training_set, 
                                           window_size, 
                                           step_size = NULL,
                                           alpha = 0, 
                                           sigmai = NULL, 
                                           lambdai = 10^(seq(-3, 1, l = 9)), 
                                           threshold = 0.9, 
                                           method = "pearson", 
                                           penalty = "l2") {
  
  # Ensure the inputs are matrices
  if (is.vector(training_set)) training_set <- as.matrix(training_set, ncol = 1)
  if (is.vector(test_set)) test_set <- as.matrix(test_set, ncol = 1)
  
  n_test <- nrow(test_set)
  n_train <- nrow(training_set)
  
  if (n_test < window_size) stop("Window size must be less than or equal to the number of rows in the test set.")
  
  anomaly_flags <- numeric(n_test - window_size + 1) # Store anomaly flags (0: normal, 1: anomaly)
  divergence_scores <- numeric(n_test - window_size + 1) # Store divergence scores for visualization
  
  for (start_idx in seq(1, n_test - window_size + 1, by = step_size)) {  # Use step_size to control the window shift
    # Define the test set window
    test_window <- test_set[start_idx:(start_idx + window_size - 1), , drop = FALSE]
    
    # Calculate the divergence score using the provided bRuLSIF function
    result <- bRuLSIF(test_window, training_set, alpha = alpha, sigmai = sigmai, lambdai = lambdai, method = method, penalty = penalty)
    
    # Compute divergence score (absolute value in case of negative scores)
    divergence_score <- abs(result$score)
    divergence_scores[start_idx] <- divergence_score
  }
  
  # Normalize the scores by dividing by the maximum divergence score
  max_score <- max(divergence_scores)
  normalized_scores <- divergence_scores / max_score
  
  # Identify anomalies based on the threshold
  anomaly_flags[normalized_scores > threshold] <- 1
  anomaly_points <- which(normalized_scores > threshold) + step_size + window_size
  # Return results as a list
  return(list(
    max_score = max_score,
    divergence_scores = divergence_scores,
    normalized_scores = normalized_scores,
    anomaly_flags = anomaly_flags,
    anomaly_points=anomaly_points
  ))
}


