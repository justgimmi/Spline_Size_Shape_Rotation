

eval_log_lik_alpha <- function(alpha_candidate, hyper, init, k_l, X, i) {
  # This function is useful to evaluate the log_posterior distribution
  # alpha_candidate --> log_proposal of alpha

  alphas_curr <- exp(alpha_candidate)
  n <- init$n
  R_i <- init$R[,,i]
  eta_matrix <- matrix(init$eta[i,], nrow = k_l, ncol = 2, byrow = TRUE)
  mean_i <- alphas_curr * (init$mean %*% R_i) + eta_matrix
  
  residual_part <- backsolve(init$chol_c, X[,,i] - mean_i, transpose = TRUE)
  residual_cand <- t(residual_part)%*%residual_part
  Q_R <- (1 / (alphas_curr^2)) * t(R_i) %*% init$Sigma_inv %*% R_i
  # quad_i <- sum(residual[, 1]^2) * Q_R[1, 1] +
  #     2 * sum(residual[, 1] * residual[, 2]) * Q_R[1, 2] +
  #     sum(residual[, 2]^2) * Q_R[2, 2]
  quad_i <- residual_cand[1, 1] * Q_R[1, 1] + 2 * residual_cand[1, 2] * Q_R[1, 2] + 
    residual_cand[2, 2] * Q_R[2, 2]
  log_lik_total <- - 0.5 * quad_i + (-2*k_l+ hyper$a)*alpha_candidate - hyper$b*alphas_curr


  return(list(log_lik = log_lik_total,
              residual = residual_cand,
              Q_R = Q_R,
              mean_i = mean_i))
}



sample_alpha <- function(init, hyper, X, k_l, burn = TRUE){
  n <- init$n
  #hyper$t <- hyper$t + 1
  for (i in 1:n) {
    alpha_candidate <- log(init$alphas[i]) + rnorm(n = 1, sd = sqrt(hyper$lambda_a[i] * hyper$Sigma_a[i]))
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
      hyper$alpha_acc[i] <-  hyper$alpha_acc[i] + 1
      
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

}
