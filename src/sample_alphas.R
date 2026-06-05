# eval_log_lik_alpha <- function(alpha_candidate, hyper, init, k_l, X) {
#   # This function is useful to evaluate the log_posterior distribution 
#   # alpha_candidate --> log_proposal of alpha
# 
#   alphas_curr <- exp(alpha_candidate)
#   n <- init$n
#   
#   log_lik_total <- 0
#   residual_new <- array(NA, dim = c(2, 2, n))
#   mean_i_new <- array(NA, dim = c(k_l, 2, n))
#   Q_R_new <- array(NA, dim = c(2, 2, n))
#   
#   for (i in 1:n) {
#     R_i <- init$R[,,i]
#     eta_matrix <- matrix(init$eta[i,], nrow = k_l, ncol = 2, byrow = TRUE)
#     mean_i <- alphas_curr[i] * (init$mean %*% R_i) + eta_matrix
#     residual <- X[,,i] - mean_i
#     mean_i_new[,,i] <- mean_i
#     residual_new[,,i] <- t(residual)%*%residual
#     Q_R <- (1 / (alphas_curr[i]^2)) * t(R_i) %*% init$Sigma_inv %*% R_i
#     Q_R_new[,,i] <- Q_R
#     quad_i <- sum(residual[, 1]^2) * Q_R[1, 1] + 
#       2 * sum(residual[, 1] * residual[, 2]) * Q_R[1, 2] + 
#       sum(residual[, 2]^2) * Q_R[2, 2]
#     
#     log_lik_total <- log_lik_total - 0.5 * quad_i
#   }
#   
#   log_lik_total <- log_lik_total + sum((-2*k_l+ hyper$a)*alpha_candidate - hyper$b*alphas_curr) 
#   
#   return(list(log_lik = log_lik_total, 
#               residual = residual_new, 
#               Q_R = Q_R_new, 
#               mean_i = mean_i_new))
# }
# 
# 
# 
# sample_alpha <- function(init, hyper, X, k_l, burn = TRUE){
#   n <- init$n
#   A <- hyper$L_a
#   alpha_candidate <- log(init$alphas) + crossprod(A, rnorm(n))
#   res_cand <- eval_log_lik_alpha(alpha_candidate, hyper, init, k_l, X)
#   u <- runif(n = 1)
#   alpha <- res_cand$log_lik - hyper$log_lik_a
#   hyper$a_aver_a <- hyper$a_aver_a + min(exp(alpha), 1)
#   if (log(u) <=  min(alpha, 0)) {
#     init$alphas <- exp(alpha_candidate)
#     init$Q_R <- res_cand$Q_R
#     init$mean_i <- res_cand$mean_i
#     init$residual <- res_cand$residual
#     hyper$log_lik_a <- res_cand$log_lik
#     hyper$log_lik <- log_density(init)
#     
#     if (burn == TRUE) {
#       hyper$t <- hyper$t + 1
#       hyper$gamma_a <- min(0.01, hyper$t^(-0.5))
#       hyper$lambda_a <- hyper$lambda_a*exp(hyper$gamma_a *(hyper$a_aver_a/hyper$t - hyper$a_opt))
#       hyper$Sigma_a <- hyper$Sigma_a + hyper$gamma_a *((log(init$alphas) - hyper$mu_a)%*%t(log(init$alphas) - hyper$mu_a) - hyper$Sigma_a) 
#       hyper$L_a <- chol(hyper$lambda_a*hyper$Sigma_a + diag(1e-6, n))
#       hyper$mu_a <- hyper$mu_a + hyper$gamma_a *(log(init$alphas)- hyper$mu_a)
#     }
#     
#   }
#   
#   else{
#     if (burn == TRUE) {
#       hyper$t <- hyper$t + 1
#       hyper$gamma_a <- min(0.01, hyper$t^(-0.5))
#       hyper$lambda_a <- hyper$lambda_a*exp(hyper$gamma_a *(hyper$a_aver_a/hyper$t - hyper$a_opt))
#       hyper$Sigma_a <- hyper$Sigma_a + hyper$gamma_a *((log(init$alphas) - hyper$mu_a)%*%t(log(init$alphas)- hyper$mu_a) - hyper$Sigma_a) 
#       hyper$L_a <- chol(hyper$lambda_a*hyper$Sigma_a + diag(1e-6, n))
#       hyper$mu_a <- hyper$mu_a + hyper$gamma_a *(log(init$alphas) - hyper$mu_a)
#     }
#   }
#   app <- list()
#   app$init <- init
#   app$hyper <- hyper 
#   return(app)
#   
# }





eval_log_lik_alpha <- function(alpha_candidate, hyper, init, k_l, X, i) {
  # This function is useful to evaluate the log_posterior distribution
  # alpha_candidate --> log_proposal of alpha

  alphas_curr <- exp(alpha_candidate)
  n <- init$n
  R_i <- init$R[,,i]
  eta_matrix <- matrix(init$eta[i,], nrow = k_l, ncol = 2, byrow = TRUE)
  mean_i <- alphas_curr * (init$mean %*% R_i) + eta_matrix
  residual <- X[,,i] - mean_i
  Q_R <- (1 / (alphas_curr^2)) * t(R_i) %*% init$Sigma_inv %*% R_i
  quad_i <- sum(residual[, 1]^2) * Q_R[1, 1] +
      2 * sum(residual[, 1] * residual[, 2]) * Q_R[1, 2] +
      sum(residual[, 2]^2) * Q_R[2, 2]

  log_lik_total <- - 0.5 * quad_i + (-2*k_l+ hyper$a)*alpha_candidate - hyper$b*alphas_curr


  return(list(log_lik = log_lik_total,
              residual = t(residual)%*%residual,
              Q_R = Q_R,
              mean_i = mean_i))
}



sample_alpha <- function(init, hyper, X, k_l, burn = TRUE){
  n <- init$n
  hyper$t <- hyper$t + 1
  for (i in 1:n) {
    alpha_candidate <- log(init$alphas[i]) + rnorm(n = 1, sd = sqrt(hyper$lambda_a[i]*hyper$Sigma_a[i]))
    res_cand <- eval_log_lik_alpha(alpha_candidate, hyper, init, k_l, X, i)
    u <- runif(n = 1)
    alpha <- res_cand$log_lik - hyper$log_lik_a[i]
    hyper$a_aver_a[i] <- hyper$a_aver_a[i] + min(exp(alpha), 1)
    
    if (log(u) <= alpha) {
      init$alphas[i] <- exp(alpha_candidate)
      init$Q_R[,,i] <- res_cand$Q_R
      init$mean_i[,,i] <- res_cand$mean_i
      init$residual[,,i] <- res_cand$residual
      hyper$log_lik_a[i] <- res_cand$log_lik
      
    }
    
    if (burn == TRUE) {
      hyper$gamma_a <- min(0.01, hyper$t^(-0.5))
      #hyper$lambda_a[i] <- hyper$lambda_a[i]*exp(hyper$gamma_a *(hyper$a_aver_a[i]/hyper$t - hyper$a_opt))
      hyper$lambda_a[i] <- hyper$lambda_a[i]*exp(hyper$gamma_a *(min(exp(alpha), 1) - hyper$a_opt))
      hyper$Sigma_a[i] <- hyper$Sigma_a[i] + hyper$gamma_a *((log(init$alphas[i]) - hyper$mu_a[i])^2 - hyper$Sigma_a[i])
      hyper$mu_a[i] <- hyper$mu_a[i] + hyper$gamma_a *(log(init$alphas[i])- hyper$mu_a[i])
    }
    
  }
  hyper$log_lik <- log_density(init)
  # app <- list()
  # app$init <- init
  # app$hyper <- hyper
  #return(app)

}
