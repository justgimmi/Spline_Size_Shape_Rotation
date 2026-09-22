sample_psi <- function(init_param, hyper, X) {
  n <- init_param$n
  k_l <- init_param$k_l
  D_p <- D_p_constructor(init_param$mean)
  k_free <- ncol(D_p) # k_l - 3
  
  S_mat <- hyper$S_mat
  Sigma_inv <- init_param$Sigma_inv
  sigma_2 <- init_param$sigma_2
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
  term2 <- (1 / sigma_2) * kronecker(diag(2), t(D_p) %*% S_mat %*% D_p)
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
  hyper$log_lik <- log_density(init_param)
  hyper$traces <- -0.5 * 1/((init_param$sigma_2)) * tr(t(init_param$W) %*% hyper$S_mat %*% init_param$W) + 
    log(det(t(hyper$D_p) %*% hyper$S_mat %*%hyper$D_p  ))
  hyper$D_p <- D_p
}
