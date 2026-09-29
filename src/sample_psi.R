sample_psi <- function(init_param, hyper, X) {
  n <- init_param$n # number of units
  k_l <- init_param$k_l # number of landmarks
  r_curr <- exp(Basis_Construction(init_param$thetas, hyper$n_basis, hyper$degree) %*% init_param$betas) #current radius
  D_p <- D_p_constructor(init_param$mean, init_param$thetas, init_param$betas, r_curr, hyper$n_basis, hyper$degree, hyper$D_p) #current D_p
  k_free <- ncol(D_p) # k_l - 5
  S_mat <- hyper$S_mat # S_mat
  K <- t(D_p) %*% S_mat%*% D_p # precision of psi if z = 1 
  Kinv <- chol2inv(chol(K)) #variance of psi if z = 1
  w    <- 1 / diag(Kinv) # precision of psi if z = 0
  w_mat <- diag(w, nrow = k_free, ncol = k_free)
  Sigma_inv <- init_param$Sigma_inv # inverse of Sigma_e
  sigma_2 <- init_param$sigma_2 # variance parameter 
  C_inv <- chol2inv(init_param$chol_c)
  
  # 1. Compute G = D_p' * S * C^{-1} * S * D_p  ((k-3) x (k-3))
  Dp_S <- t(D_p) %*% S_mat
  G <- Dp_S %*% C_inv %*% t(Dp_S)
  
  # 2. Compute M = sum_{i=1}^n (1 / alpha_i) * D_p' * S * C^{-1} * Z_i * R_i'  ((k-3) x 2)
  M <- matrix(0, nrow = k_free, ncol = 2)
  for (i in 1:n) {
    R_i <- init_param$R[, , i]
    alpha_i <- init_param$alphas[i]
    eta_mat <- matrix(init_param$eta[i, ], nrow = k_l, ncol = 2, byrow = TRUE)
    
    # Z_i = X_i - alpha_i * mu^* * R_i - 1_k * eta_i'
    Z_i <- X[, , i] - alpha_i * (init_param$mean %*% R_i) - eta_mat
    
    M <- M + (1 / alpha_i) * (Dp_S %*% C_inv %*% Z_i %*% t(R_i))
  }
  
  # 3. Construct Precision matrix Omega_psi  (2(k-3) x 2(k-3))
  term1 <- n * kronecker(Sigma_inv, G)
  if (init$z == 1) {
    term2 <- (1 / sigma_2) * kronecker(diag(2), K)
  }
  
  else{
    term2 <- (1 / sigma_2) * kronecker(diag(2),  w_mat)
  }
  Omega_psi <- term1 + term2
  
  # 4. Construct RHS vector = (I_2 \otimes M) %*% vec(Sigma_inv)  (2(k-3) x 1)
  vec_Sigma_inv <- as.vector(Sigma_inv) # Stacked columns of Sigma^{-1}
  rhs <- kronecker(diag(2), M) %*% vec_Sigma_inv
  
  # 5. Draw vec(psi) ~ Normal_{2(k-3)}(mu_psi, Omega_psi^{-1}) via Cholesky factor
  U_omega <- chol(Omega_psi) # U_omega' %*% U_omega = Omega_psi
  
  # Solve Omega_psi * mu_psi = rhs
  mu_psi <- backsolve(U_omega, forwardsolve(t(U_omega), rhs))
  
  # Draw standardized Gaussian noise and solve U_omega * noise_sample = z
  z <- rnorm(2 * k_free)
  vec_psi_sample <- mu_psi + backsolve(U_omega, z)
  
  # 6. Reconstruct matrix psi ((k-3) x 2) and update W (k x 2)
  init_param$psi <- matrix(vec_psi_sample, nrow = k_free, ncol = 2, byrow = FALSE)
  init_param$W <- D_p %*% init_param$psi
  
  # 7. Update dependent subject-level means and residuals
  for (i in 1:n) {
    R_i <- init_param$R[, , i]
    alpha_i <- init_param$alphas[i]
    eta_mat <- matrix(init_param$eta[i, ], nrow = k_l, ncol = 2, byrow = TRUE)
    
    init_param$mean_i[, , i] <- alpha_i * (init_param$mean + S_mat %*% init_param$W) %*% R_i + eta_mat
    # init_param$mean_i[, , i] <- alpha_i * (init_param$mean + project_warp(S_mat, init_param$W, init_param$thetas)) %*% R_i + eta_mat
    #project_warp(S_mat_cand, W_cand, theta_cand)
    res <- backsolve(init_param$chol_c, X[, , i] - init_param$mean_i[, , i], transpose = TRUE)
    init_param$residual[, , i] <- t(res) %*% res
  }
  #hyper$log_lik <- log_density(init_param)
  V <-  max(diag(S_mat %*% D_p %*% Kinv%*% t(D_p) %*% t(S_mat)))
  lp_sig <- -log(V) - log(1 + (init_param$sigma_2/(V*hyper$a_sigma))^2)
  logdetK <- as.numeric(determinant(K/init_param$sigma_2, logarithm = TRUE)$modulus)

  
  
  if (init_param$z == 1) {
    hyper$traces <- -0.5 * (1/init_param$sigma_2) * tr(t(init_param$W) %*% hyper$S_mat %*% init_param$W) +
      logdetK  + lp_sig
  }
  else{
    
    hyper$traces <- -0.5 * (1/init_param$sigma_2) * tr(t(init_param$psi) %*% w_mat %*% init_param$psi) + sum(log(w)) - k_free*log(init_param$sigma_2) +
      lp_sig
  }
  
  hyper$D_p <- D_p
  init_param$r <- r_curr
  hyper$log_lik <- log_density(init)
}
