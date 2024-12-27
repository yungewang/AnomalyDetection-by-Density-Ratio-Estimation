sliding_window_divergence_brulsif<- function(
    train_data,          # Training dataset
    test_data,           # Test dataset
    window_size = 10,    # Size of the sliding window
    step_size = 1,       # Step size for the sliding window
    threshold = 0.9,     # Threshold for anomaly detection, default=0.9
    runs = 1,            # Number of runs 
    alpha = 0,           # Hyperparameter for bRuLSIF, a number between 0 and 1
    sigmai = NULL,       # Hyperparameter for kernel bandwidth in bRuLSIF
    lambdai = 10^(seq(-3, 1, length.out = 9)), # Regularization parameter
    b = 50,              # Parameter for kernel centers in bRuLSIF
    fold = 5,            # Number of folds for cross-validation in bRuLSIF
    method = "pearson",  # Divergence estimation method, method="pearson", "bregman"
    penalty = "l2"       # Penalty type for bRuLSIF, penalty="l1", "l2"
{
   # Ensure input data is in matrix form
  if (!is.matrix(train_data)) {
    train_data <- matrix(train_data, ncol = 1)
  }
  if (!is.matrix(test_data)) {
    test_data <- matrix(test_data, ncol = 1)
  }
  # Number of test points
  n_test <- nrow(test_data)
  
  if (n_test < window_size) {
    stop("Test set is too small for the specified window size.")
  }
  # Check if the step size is valid
  if (step_size <= 0) {
    stop("step_size must be a positive integer.")
  }
  # Define the start and end indices for all sliding windows
  windows_start <- seq(1, n_test - window_size + 1, by = step_size)
  windows_end   <- windows_start + window_size - 1
  
  # Initialize tracking vectors for anomaly points
  point_frequency <- numeric(n_test)  # How often each point is marked as anomaly
  point_scores <- numeric(n_test)     # Sum of scores for each point
  point_counts <- numeric(n_test)     # How many times each point gets a score
  
  for (run_idx in seq_len(runs)) {
    # Compute window scores
    window_scores <- numeric(length(windows_start))
    
    for (w_idx in seq_along(windows_start)) {
      w_start <- windows_start[w_idx]
      w_end   <- windows_end[w_idx]
      
      test_window <- test_data[w_start:w_end, , drop = FALSE]
      # Calculate the bRuLSIF score for the current window
      res_brulsif <- bRuLSIF(
        x_de    = train_data,
        x_nu    = test_window,
        alpha   = alpha,
        sigmai  = sigmai,
        lambdai = lambdai,
        b       = b,
        fold    = fold,
        method  = method,
        penalty = penalty
      )
      window_scores[w_idx] <- res_brulsif$score # Store the score for the current window
    }
    
    # Normalize scores
    max_score <- max(window_scores)
    if (max_score == 0) {
      normalized_scores <- rep(0, length(window_scores))
    } else {
      normalized_scores <- window_scores / max_score # Scale scores to [0, 1]
    }
    
   # Identify windows with scores exceeding the threshold
    anomalous_windows <- which(normalized_scores > threshold)
    
    if(length(anomalous_windows) > 0) {
     # Find the anomaly point for each anomalous window
      anomaly_points <- windows_start[anomalous_windows] + step_size + window_size
     # Ensure anomaly points are within valid test indices
      anomaly_points <- anomaly_points[anomaly_points <= n_test]
       # Update the frequency, scores, and counts for detected points
      point_frequency[anomaly_points] <- point_frequency[anomaly_points] + 1
      point_scores[anomaly_points] <- point_scores[anomaly_points] + 
        normalized_scores[anomalous_windows]
      point_counts[anomaly_points] <- point_counts[anomaly_points] + 1
    }
  }
  
  # Calculate mean scores for points that were detected as anomalies
  mean_point_scores <- numeric(n_test)
  mean_point_scores[point_counts > 0] <- point_scores[point_counts > 0] / 
    point_counts[point_counts > 0]

  list(
    anomaly_indices = anomaly_points,        # Indices of detected anomalies
    point_frequency = point_frequency,       # Frequency of anomaly detection per point
    mean_point_scores = mean_point_scores,   # Mean scores for detected anomalies
    windows_start = windows_start,           # Start indices of sliding windows
    windows_end = windows_end                # End indices of sliding windows
  )
}


