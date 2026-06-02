# Sample the angle parameters -----

log_density_theta <- function(init){
  # log_dens <- matrix(0, nrow = 2, ncol = 2)
  # for (r in 1:init$n) {
  #   log_dens <- log_dens + init$Q_R[,,r]%*%init$residual[,,r]
  #   
  # }
  # return((-1/2) * sum(diag(log_dens )))
  
  val <- sum(init$Q_R[1,1,] * init$residual[1,1,] +
               init$Q_R[1,2,] * init$residual[2,1,] +
               init$Q_R[2,1,] * init$residual[1,2,] +
               init$Q_R[2,2,] * init$residual[2,2,])
  return(-0.5 * val)
}


sample_theta <- function(init, hyper, X){
  # This function build the slice sampler for the angle parameters 
  n <- init$n # number of samples
  k <- dim(X)[1] # number of landmarks
  # In this way, for each unit, we are decomposing the R^TSigma^-1R in three different elements
  Q_R_11 <- init$Q_R[1,1,] 
  Q_R_12 <- init$Q_R[1,2,] 
  Q_R_22 <- init$Q_R[2,2,]
  
  for (l in 1:k) {
    e <- rexp(1)
    log_thresh <- hyper$log_theta - e
    v <- runif(1)
    left <- init$thetas[l] - v * hyper$width_theta
    right <- left + hyper$width_theta
    bool <- TRUE
    theta_l <- init$thetas[l]
    v_old_x <- X[l, 1, ] - init$mean_i[l, 1, ] # compute the old value for the residuals x 
    v_old_y <- X[l, 2, ] - init$mean_i[l, 2, ] # compute the old value for the residuals x 
    quad_old <- v_old_x^2 * Q_R_11 + 
      2 * v_old_x * v_old_y * Q_R_12 + 
      v_old_y^2 * Q_R_22 # this is the contribution of the l-th landmark on the log-likelihood
    # it is very easy to evaluate
    
    while (bool == TRUE) {
      theta_raw <- runif(n = 1, min = left, max = right)
      theta_cand <- theta_raw %% (2*pi)
      B_l_raw <- Basis_Construction(theta_cand, hyper$n_basis, hyper$degree)
      r_l_raw <- as.numeric(exp(B_l_raw %*% as.vector(init$betas)))
      mean_l_raw <- c(r_l_raw * cos(theta_cand), r_l_raw * sin(theta_cand))
      out_x_raw <- mean_l_raw[1] * init$R[1, 1, ] + mean_l_raw[2] * init$R[2, 1, ]
      out_y_raw <- mean_l_raw[1] * init$R[1, 2, ] + mean_l_raw[2] * init$R[2, 2, ]
      
      mean_i_l_1_raw <- init$alphas * out_x_raw + init$eta[, 1]
      mean_i_l_2_raw <- init$alphas * out_y_raw + init$eta[, 2]
      v_new_x <- X[l, 1, ] - mean_i_l_1_raw
      v_new_y <- X[l, 2, ] - mean_i_l_2_raw
      
      quad_new <- v_new_x^2 * Q_R_11 + 
        2 * v_new_x * v_new_y * Q_R_12 + 
        v_new_y^2 * Q_R_22# this is the contribution of the new l-th theta
      
      log_theta_prop <- hyper$log_theta - 0.5 * sum(quad_new - quad_old)
      
      if (log_theta_prop >= log_thresh) {
        init$thetas[l] <- theta_cand
        init$r[l] <- r_l_raw
        init$mean[l, ] <- mean_l_raw
        init$mean_i[l, 1, ] <- mean_i_l_1_raw
        init$mean_i[l, 2, ] <- mean_i_l_2_raw
        
        init$residual[1, 1, ] <- init$residual[1, 1, ] - v_old_x^2 + v_new_x^2
        init$residual[1, 2, ] <- init$residual[1, 2, ] - v_old_x * v_old_y + v_new_x * v_new_y
        init$residual[2, 1, ] <- init$residual[2, 1, ] - v_old_x * v_old_y + v_new_x * v_new_y
        init$residual[2, 2, ] <- init$residual[2, 2, ] - v_old_y^2 + v_new_y^2
        
        hyper$log_theta <- log_theta_prop
        bool <- FALSE
      } else {
        if (theta_raw < theta_l) {
          left <- theta_raw
        } else {
          right <- theta_raw
        }
      }
    }
  }
  
  return(list(init = init, hyper = hyper))
}
