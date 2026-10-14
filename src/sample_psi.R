# sample_psi <- function(init_param, hyper, X) {
#   n <- init_param$n # number of units
#   k_l <- init_param$k_l # number of landmarks
#   r_curr <- exp(Basis_Construction(init_param$thetas, hyper$n_basis, hyper$degree) %*% init_param$betas) #current radius
#   D_p <- D_p_constructor(init_param$mean, init_param$thetas, init_param$betas, r_curr, hyper$n_basis, hyper$degree, hyper$D_ref) #current D_p
#   k_free <- ncol(D_p) # k_l - 5
#   S_mat <- hyper$S_mat # S_mat
#   K <- t(D_p) %*% S_mat%*% D_p# precision of psi if z = 1
#   Kinv <- chol2inv(chol(K)) #variance of psi if z = 1
#   w    <- 1 / diag(Kinv) # precision of psi if z = 0
#   w_mat <- diag(w, nrow = k_free, ncol = k_free)
#   Sigma_inv <- init_param$Sigma_inv # inverse of Sigma_e
#   sigma_2 <- init_param$sigma_2 # variance parameter
#   C_inv <- chol2inv(init_param$chol_c)
# 
#   # 1. Compute G = D_p' * S * C^{-1} * S * D_p  ((k-3) x (k-3))
#   Dp_S <- t(D_p) %*% S_mat
#   G <- Dp_S %*% C_inv %*% t(Dp_S)
# 
#   # 2. Compute M = sum_{i=1}^n (1 / alpha_i) * D_p' * S * C^{-1} * Z_i * R_i'  ((k-3) x 2)
#   M <- matrix(0, nrow = k_free, ncol = 2)
#   for (i in 1:n) {
#     R_i <- init_param$R[, , i]
#     alpha_i <- init_param$alphas[i]
#     eta_mat <- matrix(init_param$eta[i, ], nrow = k_l, ncol = 2, byrow = TRUE)
# 
#     # Z_i = X_i - alpha_i * mu^* * R_i - 1_k * eta_i'
#     Z_i <- X[, , i] - alpha_i * (init_param$mean %*% R_i) - eta_mat
# 
#     M <- M + (1 / alpha_i) * (Dp_S %*% C_inv %*% Z_i %*% t(R_i))
#   }
# 
#   # 3. Construct Precision matrix Omega_psi  (2(k-3) x 2(k-3))
#   term1 <- n * kronecker(Sigma_inv, G)
#   if (init_param$z == 1) {
#     term2 <- (1 / sigma_2) * kronecker(diag(2), K)
#   }
# 
#   else{
#     term2 <- (1 / sigma_2) * kronecker(diag(2),  w_mat)
#   }
#   Omega_psi <- term1 + term2
# 
#   # 4. Construct RHS vector = (I_2 \otimes M) %*% vec(Sigma_inv)  (2(k-3) x 1)
#   vec_Sigma_inv <- as.vector(Sigma_inv) # Stacked columns of Sigma^{-1}
#   rhs <- kronecker(diag(2), M) %*% vec_Sigma_inv
# 
#   # 5. Draw vec(psi) ~ Normal_{2(k-3)}(mu_psi, Omega_psi^{-1}) via Cholesky factor
#   U_omega <- chol(Omega_psi) # U_omega' %*% U_omega = Omega_psi
# 
#   # Solve Omega_psi * mu_psi = rhs
#   mu_psi <- backsolve(U_omega, forwardsolve(t(U_omega), rhs))
# 
#   # Draw standardized Gaussian noise and solve U_omega * noise_sample = z
#   z <- rnorm(2 * k_free)
#   vec_psi_sample <- mu_psi + backsolve(U_omega, z)
# 
#   # 6. Reconstruct matrix psi ((k-3) x 2) and update W (k x 2)
#   for (i in 1:n) {
#     init_param$psi[,,i] <- matrix(vec_psi_sample, nrow = k_free, ncol = 2, byrow = FALSE)
#     init_param$W[,,i] <- D_p %*% init_param$psi[,,i]
#   }
# 
#   # 7. Update dependent subject-level means and residuals
#   for (i in 1:n) {
#     R_i <- init_param$R[, , i]
#     alpha_i <- init_param$alphas[i]
#     eta_mat <- matrix(init_param$eta[i, ], nrow = k_l, ncol = 2, byrow = TRUE)
# 
#     init_param$mean_i[, , i] <- alpha_i * (init_param$mean + S_mat %*% init_param$W[,,i]) %*% R_i + eta_mat
#     # init_param$mean_i[, , i] <- alpha_i * (init_param$mean + project_warp(S_mat, init_param$W, init_param$thetas)) %*% R_i + eta_mat
#     #project_warp(S_mat_cand, W_cand, theta_cand)
#     res <- backsolve(init_param$chol_c, X[, , i] - init_param$mean_i[, , i], transpose = TRUE)
#     init_param$residual[, , i] <- t(res) %*% res
#   }
#   #hyper$log_lik <- log_density(init_param)
#   V <-  median(diag(S_mat %*% D_p %*% Kinv%*% t(D_p) %*% t(S_mat)))
#   lp_sig <- -log(V) - log(1 + (init_param$sigma_2/(V*hyper$a_sigma))^2)
#   logdetK <- as.numeric(determinant(K/init_param$sigma_2, logarithm = TRUE)$modulus)
# 
# 
# 
#   hyper$traces <- -0.5 * (1/init_param$sigma_2) * trace_function_psi(init_param$W, hyper$S_mat) +
#       logdetK  + lp_sig
#   # else{
#   #
#   #   hyper$traces <- -0.5 * (1/init_param$sigma_2) * tr(t(init_param$psi) %*% w_mat %*% init_param$psi) + sum(log(w)) - k_free*log(init_param$sigma_2) +
#   #     lp_sig
#   # }
# 
#   hyper$D_p <- D_p
#   init_param$r <- r_curr
#   hyper$log_lik <- log_density(init_param)
# }


sample_psi <- function(init, hyper, X) {
  n   <- init$n
  k_l <- init$k_l
  r_curr <- exp(Basis_Construction(init$thetas, hyper$n_basis, hyper$degree) %*% init$betas)
  D_p <- D_p_constructor(init$mean, init$thetas, init$betas, r_curr,
                         hyper$n_basis, hyper$degree, D_p_prev = hyper$D_p)
  k_free  <- ncol(D_p)
  S_mat   <- hyper$S_mat
  K       <- crossprod(D_p, S_mat %*% D_p)          # prior precision of each column of psi_i (x sigma2)
  Kinv    <- chol2inv(chol(K))
  Sigma_inv <- init$Sigma_inv
  sigma_2   <- init$sigma_2
  C_inv     <- chol2inv(init$chol_c)

  # 1. G = D' S C^-1 S D  (common to all subjects)
  Dp_S <- t(D_p) %*% S_mat
  G    <- Dp_S %*% C_inv %*% t(Dp_S)

  # 2. Precision of vec(psi_i): the same for every subject (alpha_i and R_i cancel)
  Omega_psi <- kronecker(Sigma_inv, G) + kronecker(diag(2), K / sigma_2)
  U_omega   <- chol(Omega_psi)

  # 3. One draw per subject
  for (i in 1:n) {
    R_i     <- init$R[, , i]
    alpha_i <- init$alphas[i]
    eta_mat <- matrix(init$eta[i, ], nrow = k_l, ncol = 2, byrow = TRUE)

    Z_i <- X[, , i] - alpha_i * (init$mean %*% R_i) - eta_mat
    M_i <- (1 / alpha_i) * (Dp_S %*% C_inv %*% Z_i %*% t(R_i))
    rhs <- as.vector(M_i %*% Sigma_inv)                          # = vec(M_i Sigma^-1)

    mu_psi <- backsolve(U_omega, forwardsolve(t(U_omega), rhs))
    vec_psi_sample <- mu_psi + backsolve(U_omega, rnorm(2 * k_free))

    init$psi[, , i] <- matrix(vec_psi_sample, nrow = k_free, ncol = 2)
    init$W[, , i]   <- D_p %*% init$psi[, , i]

    init$mean_i[, , i] <- alpha_i * (init$mean + S_mat %*% init$W[, , i]) %*% R_i + eta_mat
    res <- backsolve(init$chol_c, X[, , i] - init$mean_i[, , i], transpose = TRUE)
    init$residual[, , i] <- crossprod(res)
  }

  # 4. Prior term 
  V <- median(diag(S_mat %*% D_p %*% Kinv %*% t(D_p) %*% t(S_mat)))
  # V <- 1
  #lp_sig  <- -log(V) - log(1 + (sigma_2 / (V * hyper$a_sigma))^2)
  lp_sig  <- hyper$a_s2 * log(V) - hyper$b_s2 * V / init$sigma_2
  logdetK <- as.numeric(determinant(K / sigma_2, logarithm = TRUE)$modulus)
  hyper$traces <- -0.5 / sigma_2 * trace_function_psi(init$W, S_mat) + n * logdetK + lp_sig

  hyper$D_p <- D_p
  init$r   <- r_curr
  hyper$log_lik  <- log_density(init)
}
