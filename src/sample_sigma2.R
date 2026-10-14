# sample_sigma2 <- function(init, hyper, burn){
#   # This function is used for sampling the tempering parameters of the gibbs distribution
#   # We avoided the Inversegamma distribution because we want to have an higher mass close to 0.
#   # In particular, we want to use a half-cauchy ditribution with location parameter equal to 0
#   # and a small scale parameter
#   # suppose sigma_2 is the variance
#   ratio_curr <- init$psi / sqrt(init$sigma_2) # non central parametriziation 
#   
#   log_sigma_curr <- log(init$sigma_2) #current log std
#   log_sigma_proposal <- log_sigma_curr + rnorm(n = 1, sd = sqrt(hyper$lambda_sigma * hyper$Sigma_sigma)) # std proposed in log
#   sigma_proposal <- exp(log_sigma_proposal) # std proposed
#   sigma_curr <- init$sigma_2 # curr std
# 
#   trace_cand <- -0.5 * 1/(sigma_proposal) * tr(t(init$W) %*% hyper$S_mat %*% init$W)
#   trace_curr <- -0.5 * 1/(sigma_curr) * tr(t(init$W) %*% hyper$S_mat %*% init$W)
# 
#   det_cand <- -log_sigma_proposal*(init$k_l - 3)
#   det_curr <- -log_sigma_curr*(init$k_l - 3)
# 
#   prior_cand <- -log(1 + (sigma_proposal/hyper$a_sigma)^2) + log_sigma_proposal
#   prior_curr <- -log(1 + (sigma_curr/hyper$a_sigma)^2) + log_sigma_curr
# 
#   # prior_cand <- -sigma_proposal*hyper$a_sigma + log_sigma_proposal
#   # prior_curr <- -sigma_curr*hyper$a_sigma + log_sigma_curr
# 
#   log_cand <- trace_cand + det_cand + prior_cand
# 
#   log_curr <- trace_curr + det_curr + prior_curr
# 
#   log_alpha <- min(0, log_cand - log_curr)
#   thresh <- log(runif(n = 1))
# 
#   if (thresh <= log_alpha) {
#     init$sigma_2 <- sigma_proposal
#     hyper$traces <- -0.5 * 1/((init$sigma_2)) * tr(t(init$W) %*% hyper$S_mat %*% init$W) +
#       determinant(t(hyper$D_p) %*% hyper$S_mat %*%hyper$D_p, logarithm = TRUE)$modulus
#     hyper$sigma_acc <-  hyper$sigma_acc + 1
#   }
# 
#   if (burn == TRUE) {
#     hyper$gamma_sigma <- min(0.01, hyper$t^(-0.5))
# 
#     hyper$lambda_sigma <- hyper$lambda_sigma*exp(hyper$gamma_sigma *(exp(log_alpha) - hyper$a_opt))
#     hyper$Sigma_sigma <- hyper$Sigma_sigma+ hyper$gamma_sigma *((log(init$sigma_2) - hyper$mu_sigma)^2 - hyper$Sigma_sigma)
#     hyper$mu_sigma <- hyper$mu_sigma + hyper$gamma_sigma *(log(init$sigma_2)- hyper$mu_sigma)
# 
#   }
# 
# 
# }
# 
# sample_sigma2 <- function(init, hyper, X, burn = TRUE){
#   n <- init$n
#   k_l <- init$k_l
#   k_free <- ncol(hyper$D_p)
#   K <- t(hyper$D_p) %*% hyper$S_mat %*% hyper$D_p
#   Kinv <- chol2inv(chol(K))
#   w    <- 1 / diag(Kinv)
#   w_mat <- diag(w, nrow = k_free, ncol = k_free)
#   ## ---- STEP A (CP): sigma_2 | psi — solo prior, nessuna verosimiglianza ----
#   quad <-  tr(t(init$psi[,,1]) %*% K %*% init$psi[,,1])
# 
#   V <- median(diag(hyper$S_mat %*% hyper$D_p %*% Kinv %*% t(hyper$D_p) %*% t(hyper$S_mat)))
#   log_sigma_curr <- log(init$sigma_2)
#   log_sigma_prop <- log_sigma_curr + rnorm(1, sd = sqrt(hyper$lambda_sigma_cp * hyper$Sigma_sigma_cp))
#   sigma_prop <- exp(log_sigma_prop)
# 
#   trace_cand <- -0.5 * quad / sigma_prop
#   trace_curr <- -0.5 * quad / init$sigma_2
#   det_cand   <- -log_sigma_prop * (k_free)
#   det_curr   <- -log_sigma_curr * (k_free)
#   prior_cand <- -log(1 + (sigma_prop/(V*hyper$a_sigma))^2) + log_sigma_prop
#   prior_curr <- -log(1 + (init$sigma_2/(V*hyper$a_sigma))^2) + log_sigma_curr
# 
#   log_alpha_A <- min(0, (trace_cand + det_cand + prior_cand) -
#                        (trace_curr + det_curr + prior_curr))
#   if (log(runif(1)) <= log_alpha_A) {
#     init$sigma_2 <- sigma_prop
#     hyper$sigma_acc_cp <- hyper$sigma_acc_cp + 1
#   }
#   if (burn) {
#     g <- min(0.01, hyper$t^(-0.5))
#     hyper$lambda_sigma_cp <- hyper$lambda_sigma_cp * exp(g * (exp(log_alpha_A) - hyper$a_opt))
#     hyper$Sigma_sigma_cp  <- hyper$Sigma_sigma_cp + g * ((log(init$sigma_2) - hyper$mu_sigma_cp)^2 - hyper$Sigma_sigma_cp)
#     hyper$mu_sigma_cp     <- hyper$mu_sigma_cp + g * (log(init$sigma_2) - hyper$mu_sigma_cp)
#   }
# 
#   # # #epsilon = psi / sqrt(sigma_2) ----
#   # eps <- init$psi / sqrt(init$sigma_2)
#   # #
#   # # ## ---- STEP B (NCP): sigma_2 | epsilon — usa la verosimiglianza dei dati ----
#   # log_sigma_curr2 <- log(init$sigma_2/V)
#   # log_sigma_prop2 <- log_sigma_curr2 + rnorm(1, sd = sqrt(hyper$lambda_sigma_ncp * hyper$Sigma_sigma_ncp))
#   # sigma_prop2 <- exp(log_sigma_prop2) * V
#   # #
#   # psi_cand <- eps * sqrt(sigma_prop2)
#   # W_cand <- hyper$D_p %*% psi_cand
#   #
#   # mean_i_cand <- array(0, dim = c(k_l, 2, n))
#   # residual_cand <- array(0, dim = c(2, 2, n))
#   # val_cand <- 0
#   # for (i in 1:n) {
#   #   R_i <- init$R[,,i]; alpha_i <- init$alphas[i]
#   #   eta_mat <- matrix(init$eta[i,], nrow = k_l, ncol = 2, byrow = TRUE)
#   #   mean_i_cand[,,i] <- alpha_i * (init$mean + hyper$S_mat %*% W_cand) %*% R_i + eta_mat
#   #   res <- backsolve(init$chol_c, X[,,i] - mean_i_cand[,,i], transpose = TRUE)
#   #   residual_cand[,,i] <- t(res) %*% res
#   #   val_cand <- val_cand + sum(init$Q_R[,,i] * residual_cand[,,i])
#   # }
#   # val_curr <- sum(init$Q_R * init$residual)
#   #
#   # prior_cand2 <- -log(1 + (sigma_prop2/(hyper$a_sigma))^2) + log_sigma_prop2
#   # prior_curr2 <- -log(1 + (init$sigma_2/(hyper$a_sigma))^2) + log_sigma_curr2
#   #
#   # log_alpha_B <- min(0, (-0.5*val_cand + prior_cand2) - (-0.5*val_curr + prior_curr2))
#   # if (log(runif(1)) <= log_alpha_B) {
#   #   init$sigma_2 <- sigma_prop2
#   #   init$psi <- psi_cand
#   #   init$W <- W_cand
#   #   init$mean_i <- mean_i_cand
#   #   init$residual <- residual_cand
#   #   hyper$sigma_acc_ncp <- hyper$sigma_acc_ncp + 1
#   # }
#   # if (burn) {
#   #   g <- min(0.01, hyper$t^(-0.5))
#   #   hyper$lambda_sigma_ncp <- hyper$lambda_sigma_ncp * exp(g * (exp(log_alpha_B) - hyper$a_opt))
#   #   hyper$Sigma_sigma_ncp  <- hyper$Sigma_sigma_ncp + g * ((log(init$sigma_2/V) - hyper$mu_sigma_ncp)^2 - hyper$Sigma_sigma_ncp)
#   #   hyper$mu_sigma_ncp     <- hyper$mu_sigma_ncp + g * (log(init$sigma_2/V) - hyper$mu_sigma_ncp)
#   # }
# 
#   lp_sig <- -log(V) - log(1 + (init$sigma_2/(V*hyper$a_sigma))^2)
#   logdetK <- as.numeric(determinant(K/init$sigma_2, logarithm = TRUE)$modulus)
# 
#   hyper$traces <- -0.5/init$sigma_2 * trace_function_psi(init$W, hyper$S_mat)+ logdetK  + lp_sig
#   hyper$log_lik <- log_density(init)
# }

# 
sample_sigma2 <- function(init, hyper, X, burn = TRUE, use_ncp = TRUE){
  n      <- init$n
  k_l    <- init$k_l
  D_p    <- hyper$D_p
  S_mat  <- hyper$S_mat
  k_free <- ncol(D_p)
  K      <- crossprod(D_p, S_mat %*% D_p)
  Kinv   <- chol2inv(chol(K))
  V      <- median(diag(S_mat %*% D_p %*% Kinv %*% t(D_p) %*% t(S_mat)))
  # V <- 1
  ## ---- STEP A (CP): sigma_2 | psi_1..n — prior only ----
  quad <- trace_function_psi_K(init$psi, K)
  log_sigma_curr <- log(init$sigma_2)
  log_sigma_prop <- log_sigma_curr + rnorm(1, sd = sqrt(hyper$lambda_sigma_cp * hyper$Sigma_sigma_cp))
  sigma_prop     <- exp(log_sigma_prop)

  lpA <- function(s2, ls2){
    -0.5 * quad / s2 - n * k_free * ls2 -
    log(1 + (s2 / (V * hyper$a_sigma))^2) + ls2

    # -0.5 * quad / s2 - k_free * ls2 -
    #   log(1 + (s2 / (V * hyper$a_sigma))^2) + ls2
  }
  log_alpha_A <- min(0, lpA(sigma_prop, log_sigma_prop) - lpA(init$sigma_2, log_sigma_curr))

  if (log(runif(1)) <= log_alpha_A) {
    init$sigma_2 <- sigma_prop
    hyper$sigma_acc_cp <- hyper$sigma_acc_cp + 1
  }
  if (burn) {
    g <- min(0.01, hyper$t^(-0.5))
    hyper$lambda_sigma_cp <- hyper$lambda_sigma_cp * exp(g * (exp(log_alpha_A) - hyper$a_opt))
    hyper$Sigma_sigma_cp  <- hyper$Sigma_sigma_cp + g * ((log(init$sigma_2) - hyper$mu_sigma_cp)^2 - hyper$Sigma_sigma_cp)
    hyper$mu_sigma_cp     <- hyper$mu_sigma_cp + g * (log(init$sigma_2) - hyper$mu_sigma_cp)
  }

  ## ---- STEP B (NCP): sigma_2 | eps, eps_i = psi_i / sigma — data likelihood ----
  if (use_ncp) {
    eps <- init$psi / sqrt(init$sigma_2)
    log_sigma_curr2 <- log(init$sigma_2)
    log_sigma_prop2 <- log_sigma_curr2 + rnorm(1, sd = sqrt(hyper$lambda_sigma_ncp * hyper$Sigma_sigma_ncp))
    sigma_prop2     <- exp(log_sigma_prop2)

    psi_cand      <- eps * sqrt(sigma_prop2)
    W_cand        <- array(NA, dim(init$W))
    mean_i_cand   <- array(NA, c(k_l, 2, n))
    residual_cand <- array(NA, c(2, 2, n))
    val_cand <- 0
    for (i in 1:n) {
      W_cand[,,i] <- D_p %*% psi_cand[,,i]
      eta_mat <- matrix(init$eta[i,], k_l, 2, byrow = TRUE)
      mean_i_cand[,,i] <- init$alphas[i] * (init$mean + S_mat %*% W_cand[,,i]) %*% init$R[,,i] + eta_mat
      res <- backsolve(init$chol_c, X[,,i] - mean_i_cand[,,i], transpose = TRUE)
      residual_cand[,,i] <- crossprod(res)
      val_cand <- val_cand + sum(init$Q_R[,,i] * residual_cand[,,i])
    }
    val_curr <- sum(init$Q_R * init$residual)

    lpB <- function(s2, ls2){
      -log(1 + (s2 / (V * hyper$a_sigma))^2) + ls2
    }
    log_alpha_B <- min(0, (-0.5 * val_cand + lpB(sigma_prop2, log_sigma_prop2)) -
                         (-0.5 * val_curr + lpB(init$sigma_2, log_sigma_curr2)))
    if (log(runif(1)) <= log_alpha_B) {
      init$sigma_2  <- sigma_prop2
      init$psi      <- psi_cand
      init$W        <- W_cand
      init$mean_i   <- mean_i_cand
      init$residual <- residual_cand
      hyper$sigma_acc_ncp <- hyper$sigma_acc_ncp + 1
    }
    if (burn) {
      g <- min(0.01, hyper$t^(-0.5))
      hyper$lambda_sigma_ncp <- hyper$lambda_sigma_ncp * exp(g * (exp(log_alpha_B) - hyper$a_opt))
      hyper$Sigma_sigma_ncp  <- hyper$Sigma_sigma_ncp + g * ((log(init$sigma_2) - hyper$mu_sigma_ncp)^2 - hyper$Sigma_sigma_ncp)
      hyper$mu_sigma_ncp     <- hyper$mu_sigma_ncp + g * (log(init$sigma_2) - hyper$mu_sigma_ncp)
    }
  }

  ## ---- sync: same formula as val_3 in eval_log_lik_theta_joint / eval_log_lik_beta ----
  s2      <- init$sigma_2
  lp_sig  <- -log(V) - log(1 + (s2 / (V * hyper$a_sigma))^2)
  logdetK <- as.numeric(determinant(K / s2, logarithm = TRUE)$modulus)
  hyper$traces  <- -0.5 / s2 * trace_function_psi(init$W, S_mat) + n * logdetK + lp_sig
  # hyper$traces  <- -0.5 / s2 * trace_function_psi(init$W, S_mat) + logdetK + lp_sig
  hyper$log_lik <- log_density(init)
}

#
sample_sigma2_ig <- function(init, hyper, X){
  n <- init$n
  D_p <- hyper$D_p
  S_mat <- hyper$S_mat
  k_free <- ncol(D_p)
  K <- crossprod(D_p, S_mat %*% D_p)
  Lk <- chol(K)
  Kinv <- chol2inv(Lk)
  V <- median(diag(S_mat %*% D_p %*% Kinv %*% t(D_p) %*% t(S_mat)))
  quad <- trace_function_psi_K(init$psi, K)
  init$sigma_2 <- 1 / rgamma(1, shape = hyper$a_s2 + n * k_free,
                             rate  = hyper$b_s2 * V + 0.5 * quad)
  s2 <- init$sigma_2
  lp_sig  <- hyper$a_s2 * log(V) - hyper$b_s2 * V / s2        # log p(sigma2 | V), theta/gamma-dependent part
  logdetK <- as.numeric(determinant(K / s2, logarithm = TRUE)$modulus)
  hyper$traces  <- -0.5 / s2 * trace_function_psi(init$W, S_mat) + n * logdetK + lp_sig
  hyper$log_lik <- log_density(init)
}

