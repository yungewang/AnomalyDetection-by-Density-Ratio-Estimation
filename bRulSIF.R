library(MASS)
library(kernlab)
 
ginv2 <- function (X, tol = sqrt(.Machine$double.eps)) 
{
  if (length(dim(X)) > 2L || !(is.numeric(X) || is.complex(X))) 
    stop("'X' must be a numeric or complex matrix")
  if (!is.matrix(X)) 
    X <- as.matrix(X)
  Xsvd <- svd(X, LINPACK=TRUE)
  if (is.complex(X)) 
    Xsvd$u <- Conj(Xsvd$u)
  Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
  if (all(Positive)) 
    Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
  else if (!any(Positive)) 
    array(0, dim(X)[2L:1L])
  else Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) * 
                                               t(Xsvd$u[, Positive, drop = FALSE]))
}
bRuLSIF <- function(x_de, x_nu, alpha = 0, sigmai = NULL, lambdai = 10^(seq(-3, 1, l = 9)), 
                    b = 50, fold = 5, method = "pearson", penalty = "l2") {
  # Ensure inputs are matrices
  if (is.vector(x_de)) x_de <- as.matrix(x_de, ncols = 1)
  if (is.vector(x_nu)) x_nu <- as.matrix(x_nu, ncols = 1)
  
  n_de <- nrow(x_de) 
  n_nu <- nrow(x_nu)
  n <- n_de + n_nu
  
  # Estimate sigmai if necessary
  if (is.null(sigmai)) {
    x <- c(x_nu, x_de)
    if (n < 500) med <- median(x)
    else med <- median(x[sample(1:length(x), 500)])
    sigmai <- med * (3:7 / 5)
  }
  
  n_min <- min(n_de, n_nu)
  n_l <- length(lambdai)
  n_s <- length(sigmai)
  
  # Choose random centers
  rand_index <- sample(n_nu)  
  b <- min(b, n_nu)  
  x_ce <- x_nu[rand_index[1:b], ]  
  score_cv <- matrix(0, n_s, n_l)
 
  if (n_s == 1 && n_l == 1) {
    sigma_chosen = sigmai
    lambda_chosen = lambdai
  } else {
    # k-fold cross-validation
    cv_index_de <- sample(n_de)
    cv_split_de <- floor((1:n_de - 1) * fold / n_de) + 1
    cv_index_nu <- sample(n_nu)
    cv_split_nu <- floor((1:n_nu - 1) * fold / n_nu) + 1
    
    for (sigma_index in 1:n_s) {
      sigma <- sigmai[sigma_index]
      K_de <- kernelMatrix(rbfdot(sigma = 1 / (2 * sigma^2)), as.matrix(x_de), as.matrix(x_ce))
      K_nu <- kernelMatrix(rbfdot(sigma = 1 / (2 * sigma^2)), as.matrix(x_nu), as.matrix(x_ce))
      score_tmp <- matrix(0, fold, n_l)
      
      for (lambda_index in 1:n_l) {
        lambda <- lambdai[lambda_index]
        for (k in 1:fold) {
          Ktmp1 <- K_de[cv_index_de[cv_split_de != k], ]
          Ktmp2 <- K_nu[cv_index_nu[cv_split_nu != k], ]
          Ktmp <- alpha / nrow(Ktmp2) * t(Ktmp2) %*% Ktmp2 + 
            (1 - alpha) / nrow(Ktmp1) * t(Ktmp1) %*% Ktmp1
          mKtmp <- colMeans(K_de[cv_index_de[cv_split_de != k], ])
          mK_de <- (1-alpha) * colMeans(K_de[cv_index_de[cv_split_de != k], ])
          mK_nu <- alpha * colMeans(K_nu[cv_index_nu[cv_split_nu != k], ])
          if (penalty == "l2") {
            thetat_cv <- try(ginv(Ktmp + lambda * diag(b)) %*% mKtmp)
          } else if (penalty == "l1") {
            thetat_cv <- 1/(mK_de + mK_nu)
          }
          if (inherits(thetat_cv, "try-error")) {
            score_tmp[k, lambda_index] <- NA
          } else {
            thetah_cv <- thetat_cv
            score_tmp[k, lambda_index] <- if (method == "pearson") {
              # Pearson divergence format
              calculate_pearson_score(K_de, K_nu, cv_index_de, cv_index_nu, cv_split_de, cv_split_nu, k, thetah_cv, alpha)
            } else if (method == "bregman") {
              # Bregman divergence format
              calculate_bregman_score(K_de, K_nu, cv_index_de, cv_index_nu, cv_split_de, cv_split_nu, k, thetah_cv, alpha)
            }
          }
        }
      }
      score_cv[sigma_index, ] <- colMeans(score_tmp, na.rm = TRUE)
    }
    best_score <- min(score_cv, na.rm = TRUE)
    best_score_pos <- which(score_cv == best_score, arr.ind = TRUE)[1, ]
    sigma_chosen_index <- best_score_pos[1]
    lambda_chosen_index <- best_score_pos[2]
    lambda_chosen <- lambdai[lambda_chosen_index]
    sigma_chosen <- sigmai[sigma_chosen_index]
  }
  
  K_de <- kernelMatrix(rbfdot(sigma = 1 / (2 * sigma_chosen^2)), as.matrix(x_de), as.matrix(x_ce))
  K_nu <- kernelMatrix(rbfdot(sigma = 1 / (2 * sigma_chosen^2)), as.matrix(x_nu), as.matrix(x_ce))
  thetat <- if (penalty == "l2") {
    ginv(alpha / n_nu * t(K_nu) %*% K_nu + (1 - alpha) / n_de * t(K_de) %*% K_de + lambda_chosen * diag(b)) %*% colMeans(K_nu)
  } else if (penalty == "l1") {
    1/((1-alpha)*colMeans(K_de) + alpha *colMeans(K_nu))
  }
  
  thetat <- pmax(thetat, 0)
  wh_x_nu <- t(K_de %*% thetat)  
  
  r <- function(x) { 
    K <- kernelMatrix(rbfdot(sigma = 1 / (2 * sigma_chosen^2)), as.matrix(x_ce), as.matrix(x))
    array(t(K) %*% thetat)
  }
  g_nu <- t(K_nu %*% thetat)
  g_de <- t(K_de %*% thetat)
  if (method == "pearson") {
    # Pearson divergence
    rPE <- -alpha / 2 * mean(g_nu^2) - 
      (1 - alpha) / 2 * mean(g_de^2) + mean(g_nu) -1/2
  } else if (method == "bregman") {
    # Pearson-like scaled Bregman divergence
    rPE <- 1/2 * mean(g_nu) - (2-alpha)/(2*(1-alpha)) * mean(g_de) +
      1/(2*(1-alpha))
  }
  
  list(wtrain = wh_x_nu, r = r, thetat = thetat, score=rPE)
  
}


calculate_pearson_score <- function(K_de, K_nu, cv_index_de, cv_index_nu, cv_split_de, cv_split_nu, k, thetah_cv, alpha) {
  -alpha / 2 * mean((K_nu[cv_index_nu[cv_split_nu == k], ] %*% thetah_cv)^2) - 
    (1 - alpha) / 2 * mean((K_de[cv_index_de[cv_split_de == k], ] %*% thetah_cv)^2) +
    mean(K_nu[cv_index_nu[cv_split_nu == k], ] %*% thetah_cv) -1/2
}


calculate_bregman_score <- function(K_de, K_nu, cv_index_de, cv_index_nu, cv_split_de, cv_split_nu, k, thetah_cv, alpha) {
  1/2 * mean((K_nu[cv_index_nu[cv_split_nu == k], ] %*% thetah_cv))- (2-alpha)/(2*(1-alpha)) * mean((K_de[cv_index_de[cv_split_de == k], ] %*% thetah_cv)) +
    1/(2*(1-alpha))
}
