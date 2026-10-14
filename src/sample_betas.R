# Function to sample betas -----
# we have decided to implement elliptical slice sampling for the betas

eval_log_lik_beta <- function(gamma_cand, init, hyper, X, k_l){
  n <- init$n # number of units 
  k_free <- hyper$k_free
  B_sim <- hyper$B_sim %*% hyper$C_bar # basis function
  r <- exp(B_sim %*% as.vector(gamma_cand))# compute the current value of the radius
  mean <- cbind(r*cos(init$thetas),  r*sin(init$thetas)) # mean cand
  S_mat_cand <- S_mat_constructor(mean) # S cand
  D_p_cand <- D_p_constructor(mean, k_l) # D_p_cand
  D_P_S_cand <- S_mat_cand %*% D_p_cand # proposed S x D_p (k x k-3)
  W_cand <- tensor(D_p_cand, init$psi, 2, 1) # create the candidate value W
  SW_cand <- tensor(S_mat_cand, W_cand, 2, 1)# candidate value of SW
  app_cand <- sweep(SW_cand, MARGIN = c(1, 2), STATS = mean, FUN = "+")
  mean_i <- init$mean_i # mean_i app
  residual <- init$residual # residual app 
  
  for (i in 1:n) { # compute rotation matrix 
    mean_i[,,i] <- init$alphas[i]*app_cand[,,i]%*%init$R[,,i] + init$eta_array[,,i] # compute the mean configuration for every unit
    res <- X[,,i] - mean_i[,,i]
    residual[,,i] <- crossprod(res) # residual computation
  }
  
  val <- sum(init$Q_R[1,1,] * residual[1,1,] +
               init$Q_R[1,2,] * residual[2,1,] +
               init$Q_R[2,1,] * residual[1,2,] +
               init$Q_R[2,2,] * residual[2,2,])
  
  val_2 <- -2*init$k_l*sum(log(init$alphas)) - ((init$n*init$k_l)/2)*log(det(init$Sigma))
  
  K <- crossprod(D_p_cand, D_P_S_cand) # variance matrix
  Lk <- chol(K) # cholesky
  Kinv    <- chol2inv(Lk) # inverse of the matrix 
  V <-  median(diag(D_P_S_cand %*% Kinv %*% t(D_P_S_cand))) # V
  # lp_sig <- -log(V) - log(1 + (init$sigma_2/(V*hyper$a_sigma))^2)
  lp_sig  <- hyper$a_s2 * log(V) - hyper$b_s2 * V / init$sigma_2
  logdetK <- n*(chol_det(Lk) - k_free*log(init$sigma_2))
  val_3 <- -0.5 * (1/init$sigma_2) * trace_function(W_cand, S_mat_cand) +
      logdetK + lp_sig
    # val_3 <- -0.5 * (1/init$sigma_2) * trace_function_psi(W_cand, S_mat_cand)+
    #    logdetK + lp_sig
  log_lik <- -0.5 * val + val_2 + val_3

  return(list(log_lik =  log_lik, residual = residual, mean = mean, mean_i = mean_i,
              S_mat = S_mat_cand, W = W_cand, D_p = D_p_cand, val_3 = val_3, r = r,V = V, SW = SW_cand))
}

sample_beta <- function(init, hyper, X, k_l){

  gamma_samp <- solve(t(hyper$U_lambda), rnorm(hyper$L))
  #gamma_samp <- sqrt(init$tau)*rnorm(hyper$L)
  u <- runif(n = 1)
  thresh <- hyper$log_lik + hyper$traces + log(u)
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
  hyper$log_lik <- res_cand$log_lik - res_cand$val_3
  init$gammas <- gamma_cand 
  init$betas <- hyper$C_bar %*% gamma_cand
  init$residual <- res_cand$residual
  init$mean_i <- res_cand$mean_i
  init$mean <- res_cand$mean
  init$W <- res_cand$W
  hyper$S_mat <- res_cand$S_mat
  hyper$D_p <- res_cand$D_p
  hyper$SW <- res_cand$SW
  hyper$traces <- res_cand$val_3
  init$r <- res_cand$r
  hyper$V <- res_cand$V

}

