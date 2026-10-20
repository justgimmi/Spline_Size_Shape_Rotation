

# sample_psi <- function(init, hyper, X) {
#   n   <- init$n # sample size 
#   k_l <- init$k_l # number of landmarks 
#   r_curr <- init$r # current value of the radius
#   D_p <- hyper$D_p #current D_p
#   k_free  <- hyper$k_free # number of free directions for the warping
#   S_mat   <- hyper$S_mat
#   S_D_p <- S_mat %*% D_p
#   K       <- crossprod(D_p, S_D_p)          # prior precision of each column of psi_i (x sigma2)
#   Lk <- chol(K)
#   Kinv    <- chol2inv(Lk)
#   Sigma_inv <- init$Sigma_inv
#   sigma_2   <- init$sigma_2
# 
#   # 1. G = D' S C^-1 S D  (common to all subjects)
#   G    <- t(S_D_p)%*% S_D_p
# 
#   # 2. Precision of vec(psi_i): the same for every subject (alpha_i and R_i cancel)
#   Omega_psi <- kronecker(Sigma_inv, G) + kronecker(diag(2), K / sigma_2)
#   U_omega   <- chol(Omega_psi)
# 
#   # 3. One draw per subject
#   for (i in 1:n) {
#     R_i     <- init$R[, , i]
#     alpha_i <- init$alphas[i]
#     eta_mat <- init$eta_array[,,i]
# 
#     Z_i <- X[, , i] - alpha_i * (init$mean %*% R_i) - eta_mat
#     M_i <- (1 / alpha_i) * (t(S_D_p)  %*% Z_i %*% t(R_i))
#     rhs <- as.vector(M_i %*% Sigma_inv)                          # = vec(M_i Sigma^-1)
# 
#     mu_psi <- backsolve(U_omega, forwardsolve(t(U_omega), rhs))
#     vec_psi_sample <- mu_psi + backsolve(U_omega, rnorm(2 * k_free))
# 
#     init$psi[, , i] <- matrix(vec_psi_sample, nrow = k_free, ncol = 2)
#     init$W[, , i]   <- D_p %*% init$psi[, , i]
# 
#     init$mean_i[, , i] <- alpha_i * (init$mean + S_mat %*% init$W[, , i]) %*% R_i + eta_mat
#     res <- X[, , i] - init$mean_i[, , i]
#     init$residual[, , i] <- crossprod(res)
#   }
# 
#   # 4. Prior term 
#   hyper$V <- median(diag(S_D_p %*% Kinv %*% t(S_D_p)))
#   # V <- 1
#   #lp_sig  <- -log(V) - log(1 + (sigma_2 / (V * hyper$a_sigma))^2)
#   lp_sig  <- hyper$a_s2 * log(hyper$V) - hyper$b_s2 * hyper$V / init$sigma_2
#   logdetK <- n*(chol_det(Lk) - k_free*log(init$sigma_2))
#   hyper$traces <- -0.5 / sigma_2 * trace_function(init$W, S_mat) + logdetK + lp_sig
# 
#   hyper$D_p <- D_p
#   init$r   <- r_curr
#   hyper$log_lik  <- log_density(init)
#   hyper$SW <- tensor(S_mat, init$W, 2, 1)
# }


sample_psi <- function(init, hyper, X) {
  n      <- init$n
  k_l    <- init$k_l
  D_p    <- hyper$D_p
  k_free <- hyper$k_free
  S_mat  <- hyper$S_mat
  mu     <- init$mean
  sigma_2   <- init$sigma_2
  Sigma_inv <- init$Sigma_inv
  
  S_D_p <- S_mat %*% D_p
  K     <- crossprod(D_p, S_D_p)
  Lk    <- chol(K)
  Kinv  <- chol2inv(Lk)
  G     <- crossprod(S_D_p)
  
  Omega_psi <- kronecker(Sigma_inv, G) + kronecker(diag(2), K / sigma_2)
  U_omega   <- chol(Omega_psi)
  
  # per-unit quantities as k x n matrices (rows = landmarks, columns = units)
  alphas <- as.numeric(init$alphas)
  cs <- matrix(init$R[1, 1, ], k_l, n, byrow = TRUE)     # cos(lambda_i)
  sn <- matrix(init$R[2, 1, ], k_l, n, byrow = TRUE)     # sin(lambda_i)
  al <- matrix(alphas,         k_l, n, byrow = TRUE)
  eta1 <- matrix(init$eta_array[, 1, ], k_l, n)
  eta2 <- matrix(init$eta_array[, 2, ], k_l, n)
  X1 <- matrix(X[, 1, ], k_l, n);  X2 <- matrix(X[, 2, ], k_l, n)
  
  # F_i - mean, with F_i = (X_i - eta_i) R_i' / alpha_i   (so M_i = (S D_p)' (F_i - mean))
  E1 <- X1 - eta1; E2 <- X2 - eta2
  G1 <- (E1 * cs - E2 * sn) / al - mu[, 1]
  G2 <- (E1 * sn + E2 * cs) / al - mu[, 2]
  M1 <- crossprod(S_D_p, G1)                             # k_free x n
  M2 <- crossprod(S_D_p, G2)
  
  # rhs_i = vec(M_i Sigma^-1), stacked as (2 k_free) x n
  rhs <- rbind(M1 * Sigma_inv[1, 1] + M2 * Sigma_inv[2, 1],
               M1 * Sigma_inv[1, 2] + M2 * Sigma_inv[2, 2])
  
  mu_psi <- backsolve(U_omega, backsolve(U_omega, rhs, transpose = TRUE))
  draw   <- mu_psi + backsolve(U_omega, matrix(rnorm(2 * k_free * n), 2 * k_free, n))
  
  psi1 <- draw[1:k_free, , drop = FALSE]
  psi2 <- draw[(k_free + 1):(2 * k_free), , drop = FALSE]
  init$psi[, 1, ] <- psi1
  init$psi[, 2, ] <- psi2
  
  # W_i = D_p psi_i,  S W_i = (S D_p) psi_i
  W1  <- D_p %*% psi1
  W2  <- D_p %*% psi2
  SW1 <- S_D_p %*% psi1
  SW2 <- S_D_p %*% psi2
  init$W[, 1, ] <- W1;   init$W[, 2, ] <- W2
  
  # mean_i = alpha_i (mean + S W_i) R_i + eta_i
  A1 <- mu[, 1] + SW1;   A2 <- mu[, 2] + SW2
  Mi1 <- al * ( A1 * cs + A2 * sn) + eta1
  Mi2 <- al * (-A1 * sn + A2 * cs) + eta2
  init$mean_i[, 1, ] <- Mi1
  init$mean_i[, 2, ] <- Mi2
  
  # residual_i = (X_i - mean_i)' (X_i - mean_i)
  Y1 <- X1 - Mi1; Y2 <- X2 - Mi2
  r12 <- colSums(Y1 * Y2)
  init$residual[1, 1, ] <- colSums(Y1 * Y1)
  init$residual[1, 2, ] <- r12
  init$residual[2, 1, ] <- r12
  init$residual[2, 2, ] <- colSums(Y2 * Y2)
  
  # prior terms
  hyper$V <- median(rowSums((S_D_p %*% Kinv) * S_D_p))
  lp_sig  <- hyper$a_s2 * log(hyper$V) - hyper$b_s2 * hyper$V / sigma_2
  logdetK <- n * (chol_det(Lk) - k_free * log(sigma_2))
  tr_W    <- sum(W1 * SW1) + sum(W2 * SW2)               # sum_i tr(W_i' S W_i)
  hyper$traces <- -0.5 / sigma_2 * tr_W + logdetK + lp_sig
  
  hyper$D_p <- D_p
  hyper$SW  <- array(c(SW1, SW2), c(k_l, n, 2)) |> aperm(c(1, 3, 2))   # k x 2 x n
  hyper$log_lik <- log_density(init)
}
