
# 
# eval_log_lik_alpha <- function(alpha_candidate, hyper, init, k_l, X, i) {
#   # This function is useful to evaluate the log_posterior distribution for the MH-ratio
#   # alpha_candidate --> log_proposal of alpha
# 
#   alphas_curr <- exp(alpha_candidate) # candidate value for alpha
#   n <- init$n # sample size 
#   R_i <- init$R[,,i] # Rotation matrix 
#   eta_matrix <- init$eta_array[,,i] # eta matrix 
#   mean_i <- alphas_curr * (init$mean + hyper$SW[,,i]) %*% R_i + eta_matrix # mean for i-th unit
#   residual_part <- X[,,i] - mean_i # residual computation
#   residual_cand <- crossprod(residual_part)
#   Q_R <- (1 / (alphas_curr^2)) * t(R_i) %*% init$Sigma_inv %*% R_i
#   quad_i <- residual_cand[1, 1] * Q_R[1, 1] + 2 * residual_cand[1, 2] * Q_R[1, 2] + 
#     residual_cand[2, 2] * Q_R[2, 2]
#   # log_lik_total <- - 0.5 * quad_i + (-2*k_l+ hyper$a)*alpha_candidate - hyper$b*alphas_curr
#   log_lik_total <- - 0.5 * quad_i + (-2*k_l)*alpha_candidate
# 
#   return(list(log_lik = log_lik_total,
#               residual = residual_cand,
#               Q_R = Q_R,
#               mean_i = mean_i))
# }
# 
# 
# 
# sample_alpha <- function(init, hyper, X, k_l, burn = TRUE){# this function has been checked
#   n <- init$n
#   #hyper$t <- hyper$t + 1
#   for (i in 1:n) {
#     alpha_candidate <- log(init$alphas[i]) + rnorm(n = 1, sd = sqrt(hyper$lambda_a[i] * hyper$Sigma_a[i]))
#     res_cand <- eval_log_lik_alpha(alpha_candidate, hyper, init, k_l, X, i)
#     u <- runif(n = 1)
#     alpha <- res_cand$log_lik - hyper$log_lik_a[i]
#     hyper$a_aver_a[i] <- hyper$a_aver_a[i] + min(exp(alpha), 1)
#     
#     if (log(u) <= alpha) {
#       init$alphas[i] <- exp(alpha_candidate)
#       init$Q_R[,,i] <- res_cand$Q_R
#       init$mean_i[,,i] <- res_cand$mean_i
#       init$residual[,,i] <- res_cand$residual
#       hyper$log_lik_a[i] <- res_cand$log_lik
#       hyper$alpha_acc[i] <-  hyper$alpha_acc[i] + 1
#       
#     }
#     
#     if (burn == TRUE) {
#       hyper$gamma_a <- min(0.01, hyper$t^(-0.5))
#       #hyper$lambda_a[i] <- hyper$lambda_a[i]*exp(hyper$gamma_a *(hyper$a_aver_a[i]/hyper$t - hyper$a_opt))
#       hyper$lambda_a[i] <- hyper$lambda_a[i]*exp(hyper$gamma_a *(min(exp(alpha), 1) - hyper$a_opt))
#       hyper$Sigma_a[i] <- hyper$Sigma_a[i] + hyper$gamma_a *((log(init$alphas[i]) - hyper$mu_a[i])^2 - hyper$Sigma_a[i])
#       hyper$mu_a[i] <- hyper$mu_a[i] + hyper$gamma_a *(log(init$alphas[i])- hyper$mu_a[i])
#     }
#     
#   }
#   hyper$log_lik <- log_density(init)
# 
# }
sample_alpha <- function(init, hyper, X, k_l, burn = TRUE) {
  n  <- init$n
  Si <- init$Sigma_inv
  
  # per-unit rotations and quantities as k x n matrices
  cs <- init$R[1, 1, ]; sn <- init$R[2, 1, ]
  csm <- matrix(cs, k_l, n, byrow = TRUE); snm <- matrix(sn, k_l, n, byrow = TRUE)
  E1 <- matrix(X[, 1, ] - init$eta_array[, 1, ], k_l, n)
  E2 <- matrix(X[, 2, ] - init$eta_array[, 2, ], k_l, n)
  A1 <- init$mean[, 1] + matrix(hyper$SW[, 1, ], k_l, n)
  A2 <- init$mean[, 2] + matrix(hyper$SW[, 2, ], k_l, n)
  B1 <- A1 * csm + A2 * snm
  B2 <- -A1 * snm + A2 * csm
  
  # Sigma_tilde_i = R_i' Sigma^-1 R_i  (entries t11, t12, t22)
  t11 <- cs^2 * Si[1,1] + 2*cs*sn*Si[1,2] + sn^2 * Si[2,2]
  t22 <- sn^2 * Si[1,1] - 2*cs*sn*Si[1,2] + cs^2 * Si[2,2]
  t12 <- cs*sn*(Si[2,2] - Si[1,1]) + (cs^2 - sn^2) * Si[1,2]
  
  a <- t11 * colSums(E1*E1) + 2*t12 * colSums(E1*E2) + t22 * colSums(E2*E2)
  b <- t11 * colSums(E1*B1) + t12 * (colSums(E1*B2) + colSums(E2*B1)) + t22 * colSums(E2*B2)
  c <- t11 * colSums(B1*B1) + 2*t12 * colSums(B1*B2) + t22 * colSums(B2*B2)
  ll <- function(al) -0.5 * (a / al^2 - 2 * b / al + c) - 2 * k_l * log(al)
  
  alphas   <- as.numeric(init$alphas)
  lam_a    <- as.numeric(hyper$lambda_a)
  Sig_a    <- as.numeric(hyper$Sigma_a)
  mu_a     <- as.numeric(hyper$mu_a)
  
  prop     <- log(alphas) + rnorm(n, sd = sqrt(lam_a * Sig_a))
  ll_cand  <- ll(exp(prop))
  ll_curr  <- ll(alphas)
  lr       <- ll_cand - ll_curr                       # "alpha" in your code
  acc      <- log(runif(n)) <= lr
  
  hyper$a_aver_a   <- hyper$a_aver_a + pmin(exp(lr), 1)
  hyper$alpha_acc  <- hyper$alpha_acc + acc
  alphas[acc]      <- exp(prop[acc])
  init$alphas      <- alphas
  hyper$log_lik_a  <- ll(alphas)                      # kept for compatibility
  
  if (burn) {
    g <- min(0.01, hyper$t^(-0.5)); hyper$gamma_a <- g
    la <- log(alphas)
    hyper$lambda_a <- lam_a * exp(g * (pmin(exp(lr), 1) - hyper$a_opt))
    hyper$Sigma_a  <- Sig_a + g * ((la - mu_a)^2 - Sig_a)
    hyper$mu_a     <- mu_a + g * (la - mu_a)
  }
  
  # resync derived arrays with the new alphas (identical to old values for rejected units)
  al  <- matrix(alphas, k_l, n, byrow = TRUE)
  M1  <- al * B1 + matrix(init$eta_array[, 1, ], k_l, n)
  M2  <- al * B2 + matrix(init$eta_array[, 2, ], k_l, n)
  init$mean_i[, 1, ] <- M1; init$mean_i[, 2, ] <- M2
  Y1 <- matrix(X[, 1, ], k_l, n) - M1; Y2 <- matrix(X[, 2, ], k_l, n) - M2
  r12 <- colSums(Y1 * Y2)
  init$residual[1, 1, ] <- colSums(Y1 * Y1); init$residual[2, 2, ] <- colSums(Y2 * Y2)
  init$residual[1, 2, ] <- r12;              init$residual[2, 1, ] <- r12
  a2 <- alphas^2
  init$Q_R[1, 1, ] <- t11 / a2; init$Q_R[2, 2, ] <- t22 / a2
  init$Q_R[1, 2, ] <- t12 / a2; init$Q_R[2, 1, ] <- t12 / a2
  
  hyper$log_lik <- log_density(init)
}
