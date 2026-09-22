# Function to sample betas -----
# we have decided to implement elliptical slice sampling for the betas

eval_log_lik_beta <- function(gamma_cand, init, hyper, X, k_l){
  n <- init$n
  B_sim <- hyper$B_sim %*% hyper$C_bar
  r <- exp(B_sim %*% as.vector(gamma_cand))# compute the current value of the radius
  mu_mean_x <- r*cos(init$thetas)
  mu_mean_y <- r*sin(init$thetas)
  mean <- cbind(mu_mean_x, mu_mean_y) # average configuration 
  S_mat_cand <- S_mat_constructor(hyper$S_mat, mean)
  D_p_cand <- D_p_constructor(mean, D_p_prev = hyper$D_p)
  #D_p_cand <- hyper$D_p
  W_cand <- D_p_cand %*%init$psi
  #W_cand <- init$W
  mean_i <- array(NA, dim = c(k_l, 2, n))
  residual <- array(NA, dim = c(2, 2, n) )
  for (i in 1:n) { # compute rotation matrix 

    eta_matrix <- matrix(init$eta[i,], nrow = k_l, ncol = 2, byrow = T)
    mean_i[,,i] <- init$alphas[i]*(mean + S_mat_cand %*% W_cand)%*%init$R[,,i] + eta_matrix # compute the mean configuration for every unit
    #mean_i[,,i] <- init$alphas[i]*(mean + project_warp(S_mat_cand, W_cand, init$thetas))%*%init$R[,,i] + eta_matrix
    
    res <- backsolve(init$chol_c,X[,,i] - mean_i[,,i],  transpose = TRUE)
    residual[,,i] <- t(res)%*%res
  }
  
  val <- sum(init$Q_R[1,1,] * residual[1,1,] +
               init$Q_R[1,2,] * residual[2,1,] +
               init$Q_R[2,1,] * residual[1,2,] +
               init$Q_R[2,2,] * residual[2,2,])
  
  val_2 <- -2*init$k_l*sum(log(init$alphas)) - ((init$n*init$k_l)/2)*log(det(init$Sigma)) -
    init$n*log(det(init$C))
  
  val_3 <- -0.5 * 1/((init$sigma_2)) * tr(t(W_cand) %*% S_mat_cand %*% W_cand) + 
    determinant(t(D_p_cand) %*% S_mat_cand %*% D_p_cand, logarithm = TRUE)$modulus
  if (init$z == 1) {
    log_lik <- -0.5 * val + val_2 + val_3
  }
  else{
    log_lik <- -0.5 * val + val_2
  }
  return(list(log_lik =  log_lik, residual = residual, mean = mean, mean_i = mean_i,
              S_mat = S_mat_cand, W = W_cand, D_p = D_p_cand, val_3 = val_3))
}

sample_beta <- function(init, hyper, X, k_l){
  #gamma_samp <- rnorm(hyper$L)
  #gamma_samp <- sqrt(init$tau) * solve(t(hyper$U_lambda), rnorm(hyper$L))
  gamma_samp <- solve(t(hyper$U_lambda), rnorm(hyper$L))
  #gamma_samp <- sqrt(init$tau)*rnorm(hyper$L)
  u <- runif(n = 1)
  if (init$z == 1) {
    thresh <- hyper$log_lik + hyper$traces + log(u)
  }
  else{
  #thresh <- hyper$log_lik + hyper$traces + log(u)
    thresh <- hyper$log_lik + log(u)
  }
  theta <- runif(n = 1)*2*pi
  theta_min <- theta - 2*pi
  theta_max <- theta
  gamma_cand <- init$gammas * cos(theta) + gamma_samp*sin(theta)
  #beta_prop <- beta_cand - (hyper$P_Null %*% beta_cand) + (hyper$P_Null %*% init$betas)
  res_cand <-eval_log_lik_beta(gamma_cand, init, hyper, X, k_l) 
  # while(res_cand$log_lik + res_cand$val_3 <=   thresh){
  while(res_cand$log_lik<=   thresh){
    #print(theta)
    #print(res_cand$log_lik)
    #print(thresh)
    if (theta >=  0) {
      theta_max <- theta
    }
    else{
      theta_min <- theta
    }
  
    theta <- runif(n = 1, min = theta_min, max = theta_max)
  
    #beta_cand <- init$betas * cos(theta) + beta_samp*sin(theta)
    gamma_cand <- init$gammas * cos(theta) + gamma_samp*sin(theta)
    #beta_prop <- beta_cand - (hyper$P_Null %*% beta_cand) + (hyper$P_Null %*% init$betas)
    res_cand <-eval_log_lik_beta(gamma_cand, init, hyper, X, k_l) 
  }
  if (init$z == 1) {
    hyper$log_lik <- res_cand$log_lik - res_cand$val_3
  }
  
  else{
    hyper$log_lik <- res_cand$log_lik
    }
  init$gammas <- gamma_cand 
  init$betas <- hyper$C_bar %*% gamma_cand
  init$residual <- res_cand$residual
  init$mean_i <- res_cand$mean_i
  init$mean <- res_cand$mean
  init$W <- res_cand$W
  hyper$S_mat <- res_cand$S_mat
  hyper$D_p <- res_cand$D_p
  hyper$traces <- res_cand$val_3

}

