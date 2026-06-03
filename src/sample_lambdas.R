# Sample the angle parameters -----

log_density_lambda <- function(init){
  val <- sum(init$Q_R[1,1,] * init$residual[1,1,] +
               init$Q_R[1,2,] * init$residual[2,1,] +
               init$Q_R[2,1,] * init$residual[1,2,] +
               init$Q_R[2,2,] * init$residual[2,2,])
  return(-0.5 * val)
}

eval_log_lik_i <- function(lambda_candidate, hyper, init, i, X) {
  # this function is just a wrap up to compute the contribution of the new element inside the log-lik
  lambda_mod <- lambda_candidate %% (2*pi)
  R <- matrix(c(cos(lambda_mod), -sin(lambda_mod), 
                sin(lambda_mod), cos(lambda_mod)), 
              byrow = T, nrow = 2, ncol = 2)
  eta_matrix <- matrix(init$eta[i,], nrow = length(init$thetas), ncol = 2, byrow = T)
  mean_i <- init$alphas[i]*init$mean%*%R+ eta_matrix 
  Q_R <- (1/(init$alphas[i]^2)) *t(R)%*%init$Sigma_inv%*%R
  residual <- X[,,i] - mean_i
  
  quad_new <- sum(residual[,1]^2) * Q_R[1,1] + 
    2 * sum(residual[,1] * residual[,2]) * Q_R[1,2] + 
    sum(residual[,2]^2) * Q_R[2,2] 
  
  
  return(list(log_lik = -0.5 * sum(quad_new), residual  = residual, Q_R = Q_R, 
              mean_i = mean_i, R = R))
}



sample_lambda <- function(init, hyper, X){
  # This function build the slice sampler for the angle parameters 
  n <- init$n # number of samples
  k <- dim(X)[1] # number of landmarks
  # In this way, for each unit, we are decomposing the R^TSigma^-1R in three different elements
  m <- hyper$m
  w <- hyper$width_lambda
  for (i in 1:n) {
    
    Q_R_11 <- init$Q_R[1,1,i] 
    Q_R_12 <- init$Q_R[1,2,i] 
    Q_R_22 <- init$Q_R[2,2,i]
    
    e <- rexp(1)
    v <- runif(1)
    left <- init$lambdas[i] - v * w
    right <- left + w
    bool <- TRUE

    v_old_x <- X[, 1, i] - init$mean_i[, 1, i] # compute the old value for the residuals x of unit i
    v_old_y <- X[, 2, i] - init$mean_i[, 2, i] # compute the old value for the residuals y of unit i
    quad_old <- sum(v_old_x^2) * Q_R_11 + 
      2 * sum(v_old_x * v_old_y) * Q_R_12 + 
      sum(v_old_y^2) * Q_R_22 # this is the contribution of the i-th landmark on the log-likelihood
    # it is very easy to evaluate
    
    log_lik_i_old <- -0.5 * sum(quad_old)
    log_thresh <- log_lik_i_old - e
    u <- runif(1)
    J <- floor(m*u)
    up <- (m-1) - J
    
    while ((right - left) < 2*pi && eval_log_lik_i(left, hyper, init, i, X)$log_lik >= log_thresh && J > 0) {
      left <- left - w
      J <- J - 1
    }
    while ((right - left) < 2*pi && eval_log_lik_i(right, hyper, init,i, X)$log_lik >= log_thresh && up > 0) {
      
      right <- right + w
      up <- up - 1
    }
    
    if ((right - left) >= 2*pi) {
      right <- left + 2*pi
    }
    
    
    while (bool == TRUE) {
      lambda_raw <- left + runif(n = 1)*(right - left)
      
      res_cand <- eval_log_lik_i(lambda_raw, hyper, init, i, X)
      
      log_lambda_prop <- res_cand$log_lik
      
      if (log_lambda_prop >= log_thresh) {
        init$lambdas[i] <- lambda_raw  %% (2*pi)
        init$mean_i[, , i] <- res_cand$mean_i
        init$R[,,i] <- res_cand$R
        init$residual[, , i] <- t(res_cand$residual)%*%res_cand$residual
        init$Q_R[,,i] <- res_cand$Q_R
        hyper$log_lambda <- hyper$log_lambda +  res_cand$log_lik - log_lik_i_old
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
