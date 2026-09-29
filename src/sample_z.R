# sample_z <- function(init, hyper, X){
#   n <- init$n
#   k_l <- init$k_l
#   r_curr <- exp(Basis_Construction(init$thetas, hyper$n_basis, hyper$degree) %*% init$betas)
#   D_p <- D_p_constructor(init$mean, init$thetas, init$betas, r_curr, hyper$n_basis, hyper$degree, D_p_prev = hyper$D_p)
#   k_free <- ncol(D_p)
#   
#   S_mat <- hyper$S_mat
#   Sigma_inv <- init$Sigma_inv
#   sigma_2 <- init$sigma_2
#   C_inv <- chol2inv(init$chol_c)
#   chi <- init$chi
#   
#   # --- Quantità comuni ad entrambi gli stati (dipendono SOLO dai dati) ---
#   Dp_S <- t(D_p) %*% S_mat
#   G <- Dp_S %*% C_inv %*% t(Dp_S)
#   
#   M <- matrix(0, nrow = k_free, ncol = 2)
#   for (i in 1:n) {
#     R_i <- init$R[, , i]
#     alpha_i <- init$alphas[i]
#     eta_mat <- matrix(init$eta[i, ], nrow = k_l, ncol = 2, byrow = TRUE)
#     Z_i <- X[, , i] - alpha_i * (init$mean %*% R_i) - eta_mat
#     M <- M + (1 / alpha_i) * (Dp_S %*% C_inv %*% Z_i %*% t(R_i))
#   }
#   
#   term1 <- n * kronecker(Sigma_inv, G)                
#   vec_Sigma_inv <- as.vector(Sigma_inv)
#   rhs <- kronecker(diag(2), M) %*% vec_Sigma_inv      
# 
#   Omega_bend_base <- t(D_p) %*% S_mat %*% D_p / sigma_2 # base condivisa spike/slab
#   
#   # --- Slab: prior larga ---
#   term2_slab <- kronecker(diag(2), Omega_bend_base)
#   Omega_psi_slab <- term1 + term2_slab
#   
#   var_spike <- sigma_2                                  
#   prec_spike <- 1/var_spike
#   term2_spike <- kronecker(diag(2), diag(k_free))/var_spike
#   Omega_psi_spike <- term1 + term2_spike
#   
#   U_slab  <- chol(Omega_psi_slab)
#   U_spike <- chol(Omega_psi_spike)
#   
#   mu_psi_slab  <- backsolve(U_slab,  forwardsolve(t(U_slab),  rhs))
#   mu_psi_spike <- backsolve(U_spike, forwardsolve(t(U_spike), rhs))
#   
#   # --- Log-determinanti (nello spazio 2*k_free) ---
#   chol_bend_base <- chol(Omega_bend_base)
#   logdet_prior_slab  <- 2 * (2 * sum(log(diag(chol_bend_base))))                 # |I2 (x) Omega_bend_base|
#   logdet_prior_spike <- (2 * k_free) * log(prec_spike) # |I2 (x) Omega_bend_base/c_spike|
#   
#   logdet_post_slab  <- 2 * sum(log(diag(U_slab)))
#   logdet_post_spike <- 2 * sum(log(diag(U_spike)))
#   
#   quad_slab  <- sum((U_slab  %*% mu_psi_slab)^2)
#   quad_spike <- sum((U_spike %*% mu_psi_spike)^2)
#   
#   # --- log p(X|z=0) - log p(X|z=1), verso corretto per inclusion_prob = P(z=1|X) ---
#   log_ratio <- -0.5*logdet_prior_slab + 0.5*logdet_prior_spike +
#     0.5*logdet_post_slab  - 0.5*logdet_post_spike  -
#     0.5*quad_slab         + 0.5*quad_spike
#   # log_ratio <- -0.5*logdet_prior_slab +
#   #   0.5*logdet_post_slab -
#   #   0.5*quad_slab
#   #print(inclusion_prob)
#   inclusion_prob <- 1 / (1 + ((1 - chi) / chi) * exp(log_ratio))
#   z <- sample(c(0, 1), size = 1, prob = c(1 - inclusion_prob, inclusion_prob))
#   print(inclusion_prob)
#   #print(log_ratio)
#   # --- Ricampiona psi dalla conditional CORRISPONDENTE (mai dalla prior pura) ---
#   z_raw <- rnorm(2 * k_free)
#   if (z == 0) {
#     vec_psi <- mu_psi_spike + backsolve(U_spike, z_raw)
#     #vec_psi <- matrix(0, nrow = k_free, ncol = 2)
#     #init$sigma_2 <- rhalfcauchy(n = 1, scale = hyper$a_sigma)
#   } else {
#     vec_psi <- mu_psi_slab + backsolve(U_slab, z_raw)
#   }
#   init$psi <- matrix(vec_psi, nrow = k_free, ncol = 2)
#   init$z <- z
#   init$W <- D_p %*% init$psi
#   
#   # --- Un solo aggiornamento di mean_i/residual, con il psi appena estratto ---
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
#   if (init$z == 1) {
#     hyper$traces <- -0.5 * (1/init$sigma_2) * tr(t(init$W) %*% S_mat %*% init$W) + 
#       as.numeric(determinant(t(D_p) %*% S_mat %*% D_p, logarithm = TRUE)$modulus)
#   }
#   
#   else{
#     hyper$traces <- -0.5 * (1/init$sigma_2) * tr(t(init$psi) %*% init$psi) -k_free*log(init$sigma_2)
#   }
#   
# }
# 

compute_val <- function(W, init, hyper, X){ # compute the likelihood given the new value of W
  val <- 0
  n <- init$n
  k_l <- init$k_l
  S_mat <- hyper$S_mat
  for (i in 1:n) {
    R_i <- init$R[,,i]
    alpha_i <- init$alphas[i]
    eta_mat <- matrix(init$eta[i,], nrow = k_l, ncol = 2, byrow = TRUE)
    mean_i <- alpha_i * (init$mean + S_mat %*% W) %*% R_i + eta_mat
    res <- backsolve(init$chol_c, X[,,i] - mean_i, transpose = TRUE)
    val <- val + sum(init$Q_R[,,i] * (t(res) %*% res))
  }
  val_2 <- -2*init$k_l*sum(log(init$alphas)) - ((init$n*init$k_l)/2)*log(det(init$Sigma)) - ((init$n))*log(det(init$C))
  return(-0.5 *val + val_2)
}

# 
sample_z <- function(init, hyper, X){
  n <- init$n # number of units
  k_l <- init$k_l # number of landmarks
  r_curr <- exp(Basis_Construction(init$thetas, hyper$n_basis, hyper$degree) %*% init$betas) # current value of the radius
  D_p <- D_p_constructor(init$mean, init$thetas, init$betas, r_curr, hyper$n_basis, hyper$degree, D_p_prev = hyper$D_ref) # current D_p
  k_free <- ncol(D_p) # current free dimensions of psi
  q_i <- 0.5 # sampling probability  of z = 0
  z_cand <- sample(c(0, 1), size = 1, prob = c(q_i, 1 - q_i)) # propose z
  #z_cand <- 0
  # sample from the full conditional of psi|z --> only one bit is different if you recall
  S_mat <- hyper$S_mat
  Sigma_inv <- init$Sigma_inv
  sigma_2 <- init$sigma_2
  C_inv <- chol2inv(init$chol_c)
  Dp_S <- t(D_p) %*% S_mat
  G <- Dp_S %*% C_inv %*% t(Dp_S)
  M <- matrix(0, nrow = k_free, ncol = 2)

  K <- t(D_p) %*% S_mat%*% D_p
  Kinv <- chol2inv(chol(K))
  w    <- 1 / diag(Kinv)
  w_mat <- diag(w, nrow = k_free, ncol = k_free)
  V <-  mean(diag(S_mat %*% D_p %*% Kinv %*% t(D_p) %*% t(S_mat)))
  lp_sig <- -log(V) - log(1 + (init$sigma_2/(V*hyper$a_sigma))^2)
  logdetK <- as.numeric(determinant(K/init$sigma_2, logarithm = TRUE)$modulus)
  for (i in 1:n) {
    R_i <- init$R[, , i]
    alpha_i <- init$alphas[i]
    eta_mat <- matrix(init$eta[i, ], nrow = k_l, ncol = 2, byrow = TRUE)
    Z_i <- X[, , i] - alpha_i * (init$mean %*% R_i) - eta_mat
    M <- M + (1 / alpha_i) * (Dp_S %*% C_inv %*% Z_i %*% t(R_i))
  }
  term1 <- n * kronecker(Sigma_inv, G)
  rhs <- kronecker(diag(2), M) %*% as.vector(Sigma_inv)

  # --- Full conditional under z=1--
  Omega_bend_base <- t(D_p) %*% S_mat %*% D_p / sigma_2
  Omega_psi_slab <- term1 + kronecker(diag(2), Omega_bend_base)
  U_slab <- chol(Omega_psi_slab)
  mu_psi_slab <- backsolve(U_slab, forwardsolve(t(U_slab), rhs))

  # --- Full conditional under z=0 ---
  Omega_psi_iso <- term1 + kronecker(diag(2), w_mat/init$sigma_2)
  U_iso <- chol(Omega_psi_iso)
  mu_psi_iso <- backsolve(U_iso, forwardsolve(t(U_iso), rhs))

  z_raw <- rnorm(2 * k_free)
  W <- D_p %*% init$psi
  if (z_cand == 1) {
    psi_cand <- matrix(mu_psi_slab + backsolve(U_slab, z_raw), nrow = k_free, ncol = 2)
    W_cand <- D_p %*% psi_cand
    prior_cand <- (-0.5/init$sigma_2)*tr(t(W_cand)%*%S_mat%*%W_cand) +
      logdetK + lp_sig

    diff_cand <- as.vector(psi_cand) - as.vector(mu_psi_slab)
    posterior_cand <- sum(log(diag(U_slab))) - 0.5 * sum((U_slab %*% diff_cand)^2)

    ind_cand <- log(init$chi)
    ind_cand_bis <- log(1 - q_i)
  }
  if(z_cand == 0){
    psi_cand <- matrix(mu_psi_iso  + backsolve(U_iso,  z_raw), nrow = k_free, ncol = 2)
    W_cand <- D_p %*% psi_cand
    prior_cand <-  (-0.5 /init$sigma_2) * tr(t(psi_cand) %*% w_mat %*% psi_cand) -k_free*log(init$sigma_2) + sum(log(w)) + lp_sig 
    diff_cand <- as.vector(psi_cand) - as.vector(mu_psi_iso)
    posterior_cand <- sum(log(diag(U_iso))) - 0.5 * sum((U_iso %*% diff_cand)^2)
    ind_cand <- log(1 - init$chi)
    ind_cand_bis <- log(q_i)
  }

  if (init$z == 0) {
    diff_cand <- as.vector(init$psi) - as.vector(mu_psi_iso)
    posterior_curr <- sum(log(diag(U_iso))) - 0.5 * sum((U_iso %*% diff_cand)^2)

    prior_curr <- (-0.5 /init$sigma_2) * tr(t(init$psi) %*% w_mat %*% init$psi) -k_free*log(init$sigma_2) + lp_sig + sum(log(w))
    ind_curr <- log(1 - init$chi)
    ind_curr_bis <- log(q_i)
    }

  if (init$z == 1){
    diff_cand <- as.vector(init$psi) - as.vector(mu_psi_slab)
    posterior_curr <- sum(log(diag(U_slab))) - 0.5 * sum((U_slab %*% diff_cand)^2)
    prior_curr <-  -0.5*(1/init$sigma_2)*tr(t(W)%*%S_mat%*%W) +
      logdetK + lp_sig
    ind_curr <- log(init$chi)
    ind_curr_bis <- log(1 - q_i)
    }
  log_lik_cand <- compute_val(W_cand, init, hyper, X)
  log_lik_curr <- compute_val(W, init, hyper, X)
  #prior_curr <- hyper$traces

  log_alpha <- min((log_lik_cand - log_lik_curr) + (prior_cand - prior_curr) + (posterior_curr - posterior_cand) + (ind_cand - ind_curr ) +
    (ind_curr_bis - ind_cand_bis), 0)
  #print(log_alpha)
  log_u <- log(runif(n = 1))
  if (log_u <= log_alpha) {
    init$z <- z_cand
    init$psi <- psi_cand
    init$W <- W_cand
    hyper$traces <- prior_cand
    for (i in 1:n) {
      R_i <- init$R[, , i]
      alpha_i <- init$alphas[i]
      eta_mat <- matrix(init$eta[i, ], nrow = k_l, ncol = 2, byrow = TRUE)
      init$mean_i[, , i] <- alpha_i * (init$mean + S_mat %*% init$W) %*% R_i + eta_mat
      res <- backsolve(init$chol_c, X[, , i] - init$mean_i[, , i], transpose = TRUE)
      init$residual[, , i] <- t(res) %*% res
    }
    hyper$log_lik <- log_density(init)
    #hyper$D_p <- D_p

  }
  #print(init$z)
  init$r <- r_curr
  hyper$D_p <- D_p
  #print(init$z)

}

