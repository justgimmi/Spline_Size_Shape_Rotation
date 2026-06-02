# Sample the angle parameters -----

log_density_lambda <- function(init){
  val <- sum(init$Q_R[1,1,] * init$residual[1,1,] +
               init$Q_R[1,2,] * init$residual[2,1,] +
               init$Q_R[2,1,] * init$residual[1,2,] +
               init$Q_R[2,2,] * init$residual[2,2,])
  return(-0.5 * val)
}


sample_lambda <- function(init, hyper, X){
  # This function build the slice sampler for the angle parameters 
  n <- init$n # number of samples
  k <- dim(X)[1] # number of landmarks
  # In this way, for each unit, we are decomposing the R^TSigma^-1R in three different elements
  
  for (i in 1:n) {
    
    Q_R_11 <- init$Q_R[1,1,i] 
    Q_R_12 <- init$Q_R[1,2,i] 
    Q_R_22 <- init$Q_R[2,2,i]
    
    e <- rexp(1)
    log_thresh <- hyper$log_lambda - e
    v <- runif(1)
    left <- init$lambdas[i] - v * hyper$width_lambda
    right <- left + hyper$width_lambda
    bool <- TRUE

    v_old_x <- X[, 1, i] - init$mean_i[, 1, i] # compute the old value for the residuals xof unit i
    v_old_y <- X[, 2, i] - init$mean_i[, 2, i] # compute the old value for the residuals y of unit i
    quad_old <- sum(v_old_x^2) * Q_R_11 + 
      2 * sum(v_old_x * v_old_y) * Q_R_12 + 
      sum(v_old_y^2) * Q_R_22 # this is the contribution of the i-th landmark on the log-likelihood
    # it is very easy to evaluate
    
    while (bool == TRUE) {
      lambda_raw <- runif(n = 1, min = left, max = right)
      lambda_cand <- lambda_raw %% (2*pi)
      
      R <- matrix(c(cos(lambda_cand), -sin(lambda_cand), 
                                    sin(lambda_cand), cos(lambda_cand)), 
                                  byrow = T, nrow = 2, ncol = 2)
      eta_matrix <- matrix(init$eta[i,], nrow = length(init$thetas), ncol = 2, byrow = T)
      mean_i <- init$alphas[i]*init$mean%*%R+ eta_matrix 
      Q_R <- (1/(init$alphas[i]^2)) *t(R)%*%init$Sigma_inv%*%R
      residual <- X[,,i] - mean_i

      quad_new <- sum(residual[,1]^2) * Q_R[1,1] + 
        2 * sum(residual[,1] * residual[,2]) * Q_R[1,2] + 
        sum(residual[,2]^2) * Q_R[2,2] 
      
      log_lambda_prop <- hyper$log_lambda - 0.5 * sum(quad_new - quad_old)
      
      if (log_lambda_prop >= log_thresh) {
        init$lambdas[i] <- lambda_cand
        init$mean_i[, , i] <- mean_i
        init$R[,,i] <- R
        init$residual[, , i] <- t(residual)%*%residual
        init$Q_R[,,i] <- (1/(init$alphas[i]^2)) *t(R)%*%init$Sigma_inv%*%R
        hyper$log_lambda <- log_lambda_prop
        bool <- FALSE
      } else {
        if (lambda_raw < init$lambdas[i]) {
          left <- lambda_raw
        } else {
          right <- lambda_raw
        }
      }
    }
  }
  
  return(list(init = init, hyper = hyper))
}
