###### new idea ###### 



eval_log_lik_theta_joint <- function(theta_cand, init, hyper, X) {
  k_l <- init$k_l
  n <- init$n
  
  B_cand <- Basis_Construction(theta_cand, hyper$n_basis, hyper$degree)
  B_sim <- B_cand %*% hyper$C_bar
  r_cand <- exp(B_sim %*% as.vector(init$gammas))
  
  # 3. Construct the global latent mean shape
  mu_mean_x_cand <- r_cand*cos(theta_cand)
  mu_mean_y_cand <- r_cand*sin(theta_cand)
  mean_cand <- cbind(mu_mean_x_cand, mu_mean_y_cand)
  
  # 4. Evaluate the aggregate fit across all subjects
  val <- 0
  mean_i_array <- array(0, dim = c(k_l, 2, n))
  residual_array <- array(0, dim = c(2, 2, n))
  mat_cand <- covariance_mat(init$mat_dist, theta_cand)
  C_cand <- exp(-mat_cand/init$phi)
  chol_c_cand <- chol(C_cand)
  
  
  S_mat_cand <- S_mat_constructor(hyper$S_mat, mean_cand)
  D_p_cand <- D_p_constructor(mean_cand, theta_cand, init$betas, r_cand, hyper$n_basis, hyper$degree, D_p_prev = hyper$D_ref)
  #D_p_cand <- hyper$D_p
  #W_cand <- init$W
  W_cand <- D_p_cand %*%init$psi
  for (i in 1:n) {
    eta_matrix <- matrix(init$eta[i, ], nrow = k_l, ncol = 2, byrow = TRUE)
    #mean_i_cand <- init$alphas[i] * (mean_cand + S_mat_cand%*% W_cand)%*% init$R[,,i] + eta_matrix
    mean_i_cand <- init$alphas[i] * (mean_cand + project_warp(S_mat_cand, W_cand, theta_cand))%*% init$R[,,i] + eta_matrix
    res_cand <- X[,,i] - mean_i_cand
    res_cand <- backsolve(chol_c_cand,res_cand,  transpose = TRUE)
    res_mat <- t(res_cand) %*% res_cand
    
    mean_i_array[,,i] <- mean_i_cand
    residual_array[,,i] <- res_mat
    val <- val + sum(init$Q_R[,,i] * res_mat)
  }
  
  val_2 <- -2 * k_l * sum(log(init$alphas)) - ((n * k_l) / 2) * log(det(init$Sigma)) - n*log(det(C_cand))
  K <- t(D_p_cand) %*% S_mat_cand%*% D_p_cand
  Kinv <- chol2inv(chol(K))
  w <- 1 / diag(Kinv)
  w_mat <- diag(w, ncol = length(w), nrow = length(w))
  V <-  mean(diag(S_mat_cand %*% D_p_cand %*% solve(K) %*% t(D_p_cand) %*% t(S_mat_cand)))
  lp_sig <- -log(V) - log(1 + (init$sigma_2/(V*hyper$a_sigma))^2)
  logdetK <- as.numeric(determinant(K/init$sigma_2, logarithm = TRUE)$modulus)
  
  if (init$z  == 1) {
    k_free <- ncol(hyper$D_p)
    val_3 <- -0.5 * (1/init$sigma_2) * tr(t(W_cand) %*% S_mat_cand %*% W_cand) +
      logdetK  + lp_sig
    log_lik_total <- -0.5 * val + val_2 + val_3
  }
  else{
    k_free <- ncol(hyper$D_p)
    val_3 <- -0.5 * (1/init$sigma_2) * tr(t(init$psi) %*% w_mat %*%  init$psi) + sum(log(w)) - k_free*log(init$sigma_2) + lp_sig
    log_lik_total <- -0.5 * val + val_2 + val_3
  }
  
  return(list(log_lik = log_lik_total, 
              mean = mean_cand, 
              mean_i = mean_i_array, 
              residual = residual_array,
              B_sim = B_cand,
              C = C_cand, 
              chol_c = chol_c_cand,
              mat_cand = mat_cand,
              S_mat = S_mat_cand,
              W = W_cand,
              D_p = D_p_cand,
              val_3 = val_3,
              r = r_cand))
}


sample_theta_logit_normal <- function(init, hyper, X){
  # this function tries to implement the logit-normal reparametrization of the gaps  
  k_l <- init$k_l
  #beta_samp <- sqrt(init$tau) * hyper$U_lambda%*%beta_samp
  omega_samp <- c(rnorm(k_l), 0)
  u <- runif(n = 1)
  if (init$z == 0) {
    thresh <- hyper$log_lik + hyper$traces + log(u)
  }
  else {
    thresh <- hyper$log_lik + hyper$traces + log(u)
  }
  #thresh <- hyper$log_lik + hyper$traces + log(u)
  angle <- runif(n = 1)*2*pi
  angle_min <- angle - 2*pi
  angle_max <- angle
  omega_cand <- init$omega_theta * cos(angle) + omega_samp*sin(angle)
  gaps_cand <- exp(omega_cand)/(sum(exp(omega_cand))) 
  theta_cand <- cumsum(gaps_cand)[-c(k_l + 1)]*(2*pi)
  res_cand <-eval_log_lik_theta_joint(theta_cand, init, hyper, X) 
  while(res_cand$log_lik <=   thresh){
    #print(theta)
    #print(res_cand$log_lik)
    #print(thresh)
    if (angle >=  0) {
      angle_max <- angle
    }
    else{
      angle_min <- angle
    }
    
    angle <- runif(n = 1, min = angle_min, max = angle_max)
    omega_cand <- init$omega_theta * cos(angle) + omega_samp*sin(angle)
    gaps_cand <- exp(omega_cand)/(sum(exp(omega_cand))) 
    theta_cand <- cumsum(gaps_cand)[-c(k_l + 1)]*(2*pi)
    res_cand <- eval_log_lik_theta_joint(theta_cand, init, hyper, X) 
  }
  if (init$z == 1) {
    hyper$log_lik <- res_cand$log_lik - res_cand$val_3 
  }
  else{
    hyper$log_lik <- res_cand$log_lik - res_cand$val_3
  }
  #hyper$log_lik <- res_cand$log_lik - res_cand$val_3
  init$omega_theta <- omega_cand 
  init$thetas <- theta_cand
  init$residual <- res_cand$residual
  init$mean_i <- res_cand$mean_i
  init$mean <- res_cand$mean
  init$mat_dist <- res_cand$mat_cand
  init$C <- res_cand$C
  init$chol_c <- res_cand$chol_c
  hyper$B_sim <- res_cand$B_sim
  # for (i in 1:k_l) {
  #   for (j in 1:k_l) {
  #     init$mat_dist[i, j] <- min(2*pi - abs(init$thetas[i] - init$thetas[j]),
  #                                 abs(init$thetas[i] - init$thetas[j]))
  #     init$mat_dist[j, i] <- init$mat_dist[i, j]
  #   }
  # }
  # init$C <- exp(-init$mat_dist/init$phi)
  hyper$S_mat <- res_cand$S_mat
  hyper$D_p <- res_cand$D_p
  init$W <- res_cand$W
  init$r <- res_cand$r
  hyper$traces <- res_cand$val_3

}
