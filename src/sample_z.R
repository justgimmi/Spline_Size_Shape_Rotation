sample_z <- function(init, hyper, X){
  n <- init$n
  k_l <- init$k_l
  D_p <- D_p_constructor(init$mean, D_p_prev = hyper$D_p)
  k_free <- ncol(D_p)
  
  S_mat <- hyper$S_mat
  Sigma_inv <- init$Sigma_inv
  sigma_2 <- init$sigma_2
  C_inv <- chol2inv(init$chol_c)
  chi <- init$chi
  
  # --- Quantità comuni ad entrambi gli stati (dipendono SOLO dai dati) ---
  Dp_S <- t(D_p) %*% S_mat
  G <- Dp_S %*% C_inv %*% t(Dp_S)
  
  M <- matrix(0, nrow = k_free, ncol = 2)
  for (i in 1:n) {
    R_i <- init$R[, , i]
    alpha_i <- init$alphas[i]
    eta_mat <- matrix(init$eta[i, ], nrow = k_l, ncol = 2, byrow = TRUE)
    Z_i <- X[, , i] - alpha_i * (init$mean %*% R_i) - eta_mat
    M <- M + (1 / alpha_i) * (Dp_S %*% C_inv %*% Z_i %*% t(R_i))
  }
  
  term1 <- n * kronecker(Sigma_inv, G)                 # contributo della likelihood, uguale in entrambi gli stati
  vec_Sigma_inv <- as.vector(Sigma_inv)
  rhs <- kronecker(diag(2), M) %*% vec_Sigma_inv        # idem
  B_gamma <- hyper$B_sim %*% hyper$C_bar   # k x Lgamma
  tangent <- cbind(init$mean[,1] * B_gamma, init$mean[,2] * B_gamma)  # k x 2Lgamma, colonne = d(mu)/d(gamma_j)
  
  # proietta le colonne di D_p sullo spazio tangente e guarda quanto "entrano"
  P_tangent <- tangent %*% solve(t(tangent) %*% tangent, t(tangent))  # proiettore su span(tangent)
  overlap <- diag(t(hyper$D_p) %*% P_tangent %*% hyper$D_p)  # per ogni colonna di D_p, quanta norma cade nello spazio tangente
  Omega_bend_base <- t(D_p) %*% S_mat %*% D_p / sigma_2 # base condivisa spike/slab
  
  # --- Slab: prior larga ---
  term2_slab <- kronecker(diag(2), Omega_bend_base)
  Omega_psi_slab <- term1 + term2_slab
  
  # --- Spike: STESSA geometria, solo più concentrata (fattore relativo, non assoluto) ---
  var_spike <- 1e-10                                   # frazione della slab: da tarare, ma è RELATIVA
  prec_spike <- 1/var_spike
  term2_spike <- kronecker(diag(2), diag(k_free))/var_spike
  Omega_psi_spike <- term1 + term2_spike
  
  U_slab  <- chol(Omega_psi_slab)
  U_spike <- chol(Omega_psi_spike)
  
  mu_psi_slab  <- backsolve(U_slab,  forwardsolve(t(U_slab),  rhs))
  mu_psi_spike <- backsolve(U_spike, forwardsolve(t(U_spike), rhs))
  
  # --- Log-determinanti (nello spazio 2*k_free) ---
  chol_bend_base <- chol(Omega_bend_base)
  logdet_prior_slab  <- 2 * (2 * sum(log(diag(chol_bend_base))))                 # |I2 (x) Omega_bend_base|
  logdet_prior_spike <- (2 * k_free) * log(prec_spike) # |I2 (x) Omega_bend_base/c_spike|
  
  logdet_post_slab  <- 2 * sum(log(diag(U_slab)))
  logdet_post_spike <- 2 * sum(log(diag(U_spike)))
  
  quad_slab  <- sum((U_slab  %*% mu_psi_slab)^2)
  quad_spike <- sum((U_spike %*% mu_psi_spike)^2)
  
  # --- log p(X|z=0) - log p(X|z=1), verso corretto per inclusion_prob = P(z=1|X) ---
  # log_ratio <- -0.5*logdet_prior_slab + 0.5*logdet_prior_spike +
  #   0.5*logdet_post_slab  - 0.5*logdet_post_spike  -
  #   0.5*quad_slab         + 0.5*quad_spike
  log_ratio <- -0.5*logdet_prior_slab +
    0.5*logdet_post_slab -
    0.5*quad_slab
  #print(inclusion_prob)
  inclusion_prob <- 1 / (1 + ((1 - chi) / chi) * exp(log_ratio))
  z <- sample(c(0, 1), size = 1, prob = c(1 - inclusion_prob, inclusion_prob))
  print(inclusion_prob)
  #print(log_ratio)
  # --- Ricampiona psi dalla conditional CORRISPONDENTE (mai dalla prior pura) ---
  z_raw <- rnorm(2 * k_free)
  if (z == 0) {
    #vec_psi <- mu_psi_spike + backsolve(U_spike, z_raw)
    vec_psi <- matrix(0, nrow = k_free, ncol = 2)
    init$sigma_2 <- rhalfcauchy(n = 1, scale = hyper$a_sigma)
  } else {
    vec_psi <- mu_psi_slab + backsolve(U_slab, z_raw)
  }
  init$psi <- matrix(vec_psi, nrow = k_free, ncol = 2)
  init$z <- z
  init$W <- D_p %*% init$psi
  
  # --- Un solo aggiornamento di mean_i/residual, con il psi appena estratto ---
  for (i in 1:n) {
    R_i <- init$R[, , i]
    alpha_i <- init$alphas[i]
    eta_mat <- matrix(init$eta[i, ], nrow = k_l, ncol = 2, byrow = TRUE)
    init$mean_i[, , i] <- alpha_i * (init$mean + S_mat %*% init$W) %*% R_i + eta_mat
    res <- backsolve(init$chol_c, X[, , i] - init$mean_i[, , i], transpose = TRUE)
    init$residual[, , i] <- t(res) %*% res
  }
  
  hyper$log_lik <- log_density(init)
  hyper$D_p <- D_p
  hyper$traces <- -0.5 * (1/init$sigma_2) * tr(t(init$W) %*% S_mat %*% init$W) + 
    as.numeric(determinant(t(D_p) %*% S_mat %*% D_p, logarithm = TRUE)$modulus)
  
}



# sample_z <- function(init, hyper, X){
#   n <- init$n
#   k_l <- init$k_l
#   D_p <- D_p_constructor(init$mean, D_p_prev = hyper$D_p)
#   k_free <- ncol(D_p) # k_l - 3
#   
#   S_mat <- hyper$S_mat
#   Sigma_inv <- init$Sigma_inv
#   sigma_2 <- init$sigma_2
#   C_inv <- chol2inv(init$chol_c)
#   
#   # 1. Compute G = D_p' * S * C^{-1} * S * D_p
#   Dp_S <- t(D_p) %*% S_mat
#   G <- Dp_S %*% C_inv %*% t(Dp_S)
#   
#   # 2. Compute M
#   M <- matrix(0, nrow = k_free, ncol = 2)
#   for (i in 1:n) {
#     R_i <- init$R[, , i]
#     alpha_i <- init$alphas[i]
#     eta_mat <- matrix(init$eta[i, ], nrow = k_l, ncol = 2, byrow = TRUE)
#     Z_i <- X[, , i] - alpha_i * (init$mean %*% R_i) - eta_mat
#     M <- M + (1 / alpha_i) * (Dp_S %*% C_inv %*% Z_i %*% t(R_i))
#   }
#   
#   # 3. Definisci Precisioni a Priori (Slab e Continuous Spike)
#   # Termine Likelihood
#   term1 <- n * kronecker(Sigma_inv, G)
#   
#   # Termine Slab Prior
#   Omega_bend_base <- t(D_p) %*% S_mat %*% D_p / sigma_2
#   term2_slab <- kronecker(diag(2), Omega_bend_base)
#   
#   # Termine Spike Prior (Varianza = 1e-5 -> Precisione = 1e5)
#   var_spike <- 1e-3
#   prec_spike <- 1 / var_spike
#   term2_spike <- prec_spike * kronecker(diag(2), diag(k_free))
#   
#   # Matrici di Precisione Posterior
#   Omega_psi_slab <- term1 + term2_slab
#   Omega_psi_spike <- term1 + term2_spike
#   
#   # 4. RHS Vector
#   vec_Sigma_inv <- as.vector(Sigma_inv)
#   rhs <- kronecker(diag(2), M) %*% vec_Sigma_inv
#   
#   # 5. Factorizzazioni di Cholesky e Medie Posterior
#   U_slab <- chol(Omega_psi_slab)
#   U_spike <- chol(Omega_psi_spike)
#   
#   mu_psi_slab <- backsolve(U_slab, forwardsolve(t(U_slab), rhs))
#   mu_psi_spike <- backsolve(U_spike, forwardsolve(t(U_spike), rhs))
#   
#   # 6. Log-Determinanti e Forme Quadratiche
#   # Prior Log-Dets (spazio 2*k_free)
#   chol_bend_base <- chol(Omega_bend_base)
#   logdet_prior_slab <- 2 * (2 * sum(log(diag(chol_bend_base))))
#   logdet_prior_spike <- (2 * k_free) * log(prec_spike)
#   
#   # Posterior Log-Dets
#   logdet_post_slab <- 2 * sum(log(diag(U_slab)))
#   logdet_post_spike <- 2 * sum(log(diag(U_spike)))
#   
#   # Forme quadratiche mu' * Omega * mu
#   quad_slab <- sum((U_slab %*% mu_psi_slab)^2)
#   quad_spike <- sum((U_spike %*% mu_psi_spike)^2)
#   
#   # 7. Log-Ratio Corretto: log( p(X | z=1) / p(X | z=0) )
#   log_ratio <- -0.5 * logdet_prior_slab + 0.5 * logdet_prior_spike - 
#     0.5 * logdet_post_slab  - 0.5 * logdet_post_spike  + 
#     0.5 * quad_slab         + 0.5 * quad_spike
#   
#   # Probabilità di inclusione ed estrazione di z
#   inclusion_prob <- 1 / (1 + ((1 - init$chi) / init$chi) * exp(log_ratio))
#   init$z <- sample(c(0, 1), size = 1, prob = c(1 - inclusion_prob, inclusion_prob))
#   
#   # 8. Campionamento di psi condizionato a z
#   z_raw <- rnorm(2 * k_free)
#   print(log_ratio)
#   if (init$z == 0) {
#     vec_psi_sample <- mu_psi_spike + backsolve(U_spike, z_raw)
#     init$psi <- matrix(vec_psi_sample, nrow = k_free, ncol = 2, byrow = FALSE)
#     init$sigma_2 <- rhalfcauchy(n = 1, scale = hyper$a_sigma)
#   } else {
#     vec_psi_sample <- mu_psi_slab + backsolve(U_slab, z_raw)
#     init$psi <- matrix(vec_psi_sample, nrow = k_free, ncol = 2, byrow = FALSE)
#   }
#   
#   init$W <- D_p %*% init$psi
#   
#   # 9. Aggiornamento medie e residui soggetti
#   for (i in 1:n) {
#     R_i <- init$R[, , i]
#     alpha_i <- init$alphas[i]
#     eta_mat <- matrix(init$eta[i, ], nrow = k_l, ncol = 2, byrow = TRUE)
#     init$mean_i[, , i] <- alpha_i * (init$mean + S_mat %*% init$W) %*% R_i + eta_mat
#     res <- backsolve(init$chol_c, X[, , i] - init$mean_i[, , i], transpose = TRUE)
#     init$residual[, , i] <- t(res) %*% res
#   }
#   
#   hyper$log_lik <- log_density(init)
#   hyper$D_p <- D_p
#   
# }

# 
# sample_z <- function(init, hyper, X){
#   n <- init$n
#   k_l <- init$k_l
#   D_p <- D_p_constructor(init$mean, D_p_prev = hyper$D_p)
#   k_free <- ncol(D_p) # k_l - 3
#   
#   S_mat <- hyper$S_mat
#   mean_i <- array(0, dim = c(k_l,2, n ))
#   residual <- array(0, dim = c(2, 2, n))
#   if (condition) {
#     
#   }
#   for (i in 1:n) {
#     R_i <- init$R[, , i]
#     alpha_i <- init$alphas[i]
#     eta_mat <- matrix(init$eta[i, ], nrow = k_l, ncol = 2, byrow = TRUE)
#     
#     mean_i[, , i] <- alpha_i * (init$mean + S_mat %*% init$W) %*% R_i + eta_mat
#     # init$mean_i[, , i] <- alpha_i * (init$mean + project_warp(S_mat, init$W, init$thetas)) %*% R_i + eta_mat
#     #project_warp(S_mat_cand, W_cand, theta_cand)
#     res <- backsolve(init$chol_c, X[, , i] - mean_i[, , i], transpose = TRUE)
#     residual[, , i] <- t(res) %*% res
#   }
#   inclusion_prob <- 1/(1 + ((1 - init$chi)/(init$chi))*exp(log_ratio))
#   init$z <- sample(c(0, 1), size = 1, prob = c(1 - inclusion_prob, inclusion_prob))
#   
#   if (init$z == 0) {
#     init$psi <- matrix(0, nrow = k_free, ncol = 2)
#     #init$sigma_2 <- rhalfcauchy(n = 1, scale = hyper$a_sigma)
#   }
#   
#   else{
#     z <- rnorm(2 * k_free)
#     vec_psi_sample <- mu_psi + backsolve(U_omega, z)
#     
#     # 6. Reconstruct matrix psi ((k-3) x 2) and update W (k x 2)
#     init$psi <- matrix(vec_psi_sample, nrow = k_free, ncol = 2, byrow = FALSE) 
#   }
#   init$W <- D_p %*% init$psi
#   
#   # 7. Update dependent subject-level means and residuals
#   for (i in 1:n) {
#     R_i <- init$R[, , i]
#     alpha_i <- init$alphas[i]
#     eta_mat <- matrix(init$eta[i, ], nrow = k_l, ncol = 2, byrow = TRUE)
#     
#     init$mean_i[, , i] <- alpha_i * (init$mean + S_mat %*% init$W) %*% R_i + eta_mat
#     # init$mean_i[, , i] <- alpha_i * (init$mean + project_warp(S_mat, init$W, init$thetas)) %*% R_i + eta_mat
#     #project_warp(S_mat_cand, W_cand, theta_cand)
#     res <- backsolve(init$chol_c, X[, , i] - init$mean_i[, , i], transpose = TRUE)
#     init$residual[, , i] <- t(res) %*% res
#   }
#   hyper$log_lik <- log_density(init)
#   hyper$traces <- -0.5 * 1/((init$sigma_2)) * tr(t(init$W) %*% hyper$S_mat %*% init$W) + 
#     log(det(t(hyper$D_p) %*% hyper$S_mat %*%hyper$D_p  ))
#   hyper$D_p <- D_p
#   cat(
#     "log_ratio =", log_ratio,
#     "chi =", init$chi,
#     "p1 =", inclusion_prob,
#     "\n"
#   )
# }

