sample_phi <- function(init, hyper, X, burn = TRUE){
  phi_curr <- init$phi
  phi_curr_mod <- (phi_curr - hyper$a_phi)/(hyper$b_phi - hyper$a_phi)
  phi_curr_trans <-  g_phi_star(phi_curr, hyper$a_phi, hyper$b_phi)
  phi_candidate <- phi_curr_trans +
    rnorm(n = 1, sd = sqrt(hyper$lambda_phi * hyper$Sigma_phi))

  phi_candidate_star <- g_phi(phi_candidate, hyper$a_phi, hyper$b_phi)
  
  C_cand <- exp(-init$mat_dist/phi_candidate_star) # landmark covariance matrix
  chol_cand <- chol((C_cand))
  
  log_curr <- -init$n*log(det(init$C)) + phi_curr_trans - 2*log(1 + exp(phi_curr_trans)) -0.5*sum(init$Q_R[1,1,] * init$residual[1,1,] +
                                               init$Q_R[1,2,] * init$residual[2,1,] +
                                               init$Q_R[2,1,] * init$residual[1,2,] +
                                               init$Q_R[2,2,] * init$residual[2,2,]) 
  residual_array <- array(0, dim = c(2, 2, init$n))
  val <- 0
  for (i in 1:init$n) {
    mean_i_cand <- init$mean_i[,,i]
    res_cand <- X[,,i] - mean_i_cand
    res_cand <- backsolve(chol_cand,res_cand,  transpose = TRUE)
    res_mat <- t(res_cand) %*% res_cand
    
    residual_array[,,i] <- res_mat
    val <- val + sum(init$Q_R[,,i] * res_mat)
  }
  
  
  log_star <- -init$n*log(det(C_cand)) + phi_candidate - 2*log(1 + exp(phi_candidate)) - 0.5*val
  u <- runif(n = 1)
  alpha <- min(log_star - log_curr, 0)
  if (log(u) < alpha) {
    init$phi <- phi_candidate_star
    init$residual <- residual_array
    hyper$log_lik <- log_star - phi_candidate + 2*log(1 + exp(phi_candidate)) -2*init$k_l*sum(log(init$alphas)) - ((init$n*init$k_l)/2)*log(det(init$Sigma))
    init$C <- C_cand 
    init$chol_c <- chol_cand
    hyper$phi_acc <- hyper$phi_acc + 1
  }
  
  if (burn == TRUE) {
    hyper$gamma_phi <- min(0.01, hyper$t^(-0.5))
    hyper$lambda_phi <- hyper$lambda_phi * exp(hyper$gamma_phi * (exp(alpha) - hyper$a_opt))
    hyper$Sigma_phi  <- hyper$Sigma_phi + hyper$gamma_phi *
      ((g_phi_star(init$phi, hyper$a_phi, hyper$b_phi) - hyper$mu_phi)^2 - hyper$Sigma_phi)
    hyper$mu_phi     <- hyper$mu_phi + hyper$gamma_phi *
      (g_phi_star(init$phi, hyper$a_phi, hyper$b_phi) - hyper$mu_phi)
  }
  
  
  
}