# Sample the angle parameters -----


sample_lambda_eta <- function(init, hyper, index, X, burn = TRUE){
  # the idea of this function is to sample lambda and eta all together
  # we sample lambda from a normal distribution (recall that sin and cosine are periodic functions)
  # while we sample eta from the full conditional.
  # After a bit of math, it is possible to observe that the Metropolis-Hastings ratio
  # depend on the ratio between \pi(\lambda|... except eta) computed at the proposed and the current value
  # this implies that theoretically this is a partially collapsed metropolis step within a gibbs sampler
  
  k_l <- init$k_l # number of landmarks
  alpha_i <- init$alphas[index] # size per unit i
  lambda_curr <- init$lambdas[index] # current lambda i
  eta_curr <- init$eta[index, ] # current eta i
  ones <- matrix(1, nrow = k_l, ncol = 1)
  lambda_prop <- lambda_curr + rnorm(n = 1, mean = 0, sd = sqrt(hyper$lambda_lambda[index] *  hyper$Sigma_lambda[index])) # proposed lambda
  #lambda_star <- runif(n = 1, 0, 2*pi)
  lambda_star <- lambda_prop
  R_mat_star <- Rmat(lambda_star)
  Q_R_star <- (1/alpha_i^2) * t(R_mat_star) %*% init$Sigma_inv %*% R_mat_star # proposed inverse
  c_sum <- k_l # 1_k^T %*% 1_k
  
  mean_res_star <- alpha_i * (init$mean + hyper$SW[,,index]) %*% R_mat_star # proposed mean configuration
  tilde_E_star <- X[,,index] - mean_res_star
  
  quad_star  <- sum((t(tilde_E_star)%*% tilde_E_star) * Q_R_star)
  Sigma_star <- solve(c_sum * Q_R_star + hyper$Sigma_eta_inv)
  mu_eta_star <- t(as.matrix(Sigma_star %*% Q_R_star %*% (t(tilde_E_star) %*% ones)))
  quad_eta_star <- as.numeric(
    mu_eta_star %*%
      solve(Sigma_star) %*%
      t(mu_eta_star)
  )
  log_prop <- 0.5*log((det(Sigma_star))) - 0.5*quad_star + 0.5*quad_eta_star
  Sigma_curr <- solve(c_sum * init$Q_R[,,index] + hyper$Sigma_eta_inv)

  
  #Sigma_curr <- solve(init$k_l*init$Q_R[,,index] + hyper$Sigma_eta_inv)
  eta_matrix <- matrix(init$eta[index,], nrow = length(init$thetas), ncol = 2, byrow = T)
  tilde_E_curr <- X[,,index] - init$mean_i[,,index] + eta_matrix
  quad_curr  <- sum((t(tilde_E_curr)  %*% tilde_E_curr) * init$Q_R[,,index])
  #quad_curr <- sum((tilde_E_curr %*% init$Q_R[,,index]) * tilde_E_curr)
  mu_eta_curr <- t(as.matrix(Sigma_curr %*% init$Q_R[,,index] %*% (t(tilde_E_curr) %*% ones)))
  #mu_eta_curr <- as.matrix(t(Sigma_curr %*%init$Q_R[,,index] %*%t(tilde_E_curr) %*% matrix(1, nrow = k_l)))
  quad_eta_curr <- as.numeric(
    mu_eta_curr %*%
      solve(Sigma_curr) %*%
      t(mu_eta_curr)
  )
  log_curr <- 0.5*log((det(Sigma_curr))) - 0.5*quad_curr + 0.5*quad_eta_curr

  log_alpha <- min(log_prop - log_curr, 0)
  
  if (log(runif(1)) < log_alpha) {
    init$lambdas[index] <- lambda_star
    init$Q_R[,,index] <- Q_R_star
    init$eta[index, ] <- mvrnorm(n = 1, mu = mu_eta_star, Sigma = Sigma_star)
    init$eta_array[,,index] <- matrix(init$eta[index, ], nrow = k_l, ncol = 2, byrow = TRUE)
    init$R[,,index] <- R_mat_star
    init$mean_i[,,index] <- mean_res_star + matrix(init$eta[index, ],
                                                   nrow = k_l, ncol = 2, byrow = TRUE)
    
    E_final <- tilde_E_star - init$eta_array[,,index]
    # res_final <- backsolve(init$chol_c, E_final, transpose = TRUE)   # <- whitening added
    init$residual[,,index] <- crossprod(E_final)
    
    hyper$lambda_acc[index] <- hyper$lambda_acc[index] + 1
  }
  
  if (burn == TRUE) {
    
    hyper$gamma_lambda <- min(0.01, hyper$t^(-0.5))
    hyper$lambda_lambda[index] <- hyper$lambda_lambda[index] * exp(hyper$gamma_lambda * (min(exp(log_alpha), 1) - hyper$lambda_opt))
    hyper$Sigma_lambda[index] <- hyper$Sigma_lambda[index] + hyper$gamma_lambda *((init$lambdas[index] - hyper$mu_lambda[index])^2 - hyper$Sigma_lambda[index])
    hyper$mu_lambda[index] <- hyper$mu_lambda[index] + hyper$gamma_lambda *(init$lambdas[index] - hyper$mu_lambda[index])
  }
  
  
}


# sample_lambda_eta <- function(init, hyper, index, X, burn = TRUE){
#   # Collapsed Metropolis-Hastings: eta_i is analytically integrated out of
#   # the acceptance ratio for lambda_i, then drawn from its exact conditional
#   # once lambda_i is accepted.
#   
#   k_l <- init$k_l
#   alpha_i <- init$alphas[index]
#   lambda_curr <- init$lambdas[index]
#   eta_curr <- init$eta[index, ]
#   ones <- matrix(1, nrow = k_l, ncol = 1)
#   
#   C_chol <- init$chol_c          # k_l x k_l, inverse of the landmark-correlation chol matrix
#   s_part <- backsolve(C_chol,ones, transpose = TRUE)              # scalar: 1_k' Cinv 1_k
#   s <- as.numeric(t(s_part)%*%s_part)
#   #Cinv_1 <- rowSums(Cinv)    # Cinv %*% 1_k, precomputed once
#   
#   lambda_prop <- lambda_curr + rnorm(1, mean = 0, sd = sqrt(hyper$lambda_lambda[index] * hyper$Sigma_lambda[index]))
#   lambda_star <- lambda_prop
#   R_mat_star <- matrix(c(cos(lambda_star), -sin(lambda_star),
#                          sin(lambda_star),  cos(lambda_star)), byrow = TRUE, nrow = 2)
#   Q_R_star <- (1/alpha_i^2) * t(R_mat_star) %*% init$Sigma_inv %*% R_mat_star
#   mean_res_star <- alpha_i * (init$mean + hyper$S_mat%*%init$W) %*% R_mat_star
#   tilde_E_star <- X[,,index] - mean_res_star     # E_i(lambda_star)
#   res_star <- backsolve(C_chol, tilde_E_star, transpose = TRUE)
#   res_star_bis <- t(res_star)%*%res_star
#   # const_i = tr[Q_i * E_i' Cinv E_i]
#   quad_star <- res_star_bis[1, 1] * Q_R_star[1, 1] + 2 * res_star_bis[1, 2] * Q_R_star[1, 2] + 
#     res_star_bis[2, 2] * Q_R_star[2, 2]
#   
#   # posterior precision/covariance for eta_i, using s = 1_k' Cinv 1_k (NOT k_l)
#   Sigma_star <- solve(s * Q_R_star + hyper$Sigma_eta_inv)
#   
#   # b_i = E_i' Cinv 1_k ; mu_eta = Sigma_post %*% Q_i %*% b_i
#   b_star  <- t(res_star) %*% s_part
#   #b_star <- v %*% ones
#   mu_eta_star <- as.numeric(Sigma_star %*% Q_R_star %*% b_star)
#   
#   quad_eta_star <- as.numeric(t(mu_eta_star) %*% solve(Sigma_star) %*% mu_eta_star)
#   
#   log_prop <- 0.5*log(det(Sigma_star)) - 0.5*quad_star + 0.5*quad_eta_star
#   
#   # ---- current state ----
#   Q_R_curr <- init$Q_R[,,index]
#   eta_matrix <- matrix(eta_curr, nrow = k_l, ncol = 2, byrow = TRUE)
#   tilde_E_curr <- X[,,index] - init$mean_i[,,index] + eta_matrix   # E_i(lambda_curr)
#   res_curr <- backsolve(C_chol, tilde_E_curr, transpose = TRUE)
#   res_curr_bis <- t(res_curr)%*%res_curr
#   
#   quad_curr <- res_curr_bis[1, 1] * Q_R_curr[1, 1] + 2 * res_curr_bis[1, 2] * Q_R_curr[1, 2] + 
#     res_curr_bis[2, 2] * Q_R_curr[2, 2]
#   Sigma_curr <- solve(s * Q_R_curr + hyper$Sigma_eta_inv)
#   b_curr  <- t(res_curr) %*% s_part
#   mu_eta_curr <- as.numeric(Sigma_curr %*% Q_R_curr %*% b_curr)
#   quad_eta_curr <- as.numeric(t(mu_eta_curr) %*% solve(Sigma_curr) %*% mu_eta_curr)
#   
#   log_curr <- 0.5*log(det(Sigma_curr)) - 0.5*quad_curr + 0.5*quad_eta_curr
#   
#   log_alpha <- min(log_prop - log_curr, 0)
#   
#   if (log(runif(1)) < log_alpha) {
#     init$lambdas[index] <- lambda_star
#     init$Q_R[,,index] <- Q_R_star
#     init$eta[index, ] <- MASS::mvrnorm(n = 1, mu = mu_eta_star, Sigma = Sigma_star)
#     init$R[,,index] <- R_mat_star
#     
#     eta_new_matrix <- matrix(init$eta[index, ], nrow = k_l, ncol = 2, byrow = TRUE)
#     init$mean_i[,,index] <- mean_res_star + eta_new_matrix
#     resid_star <- tilde_E_star - eta_new_matrix
#     res_stars <- backsolve(C_chol, resid_star, transpose = TRUE)
#     init$residual[,,index] <- t(res_stars)%*%res_stars  # keep consistent with C-weighted residual, if that's what downstream expects
#     
#     hyper$lambda_acc[index] <- hyper$lambda_acc[index] + 1
#   }
#   
#   if (burn == TRUE) {
#     hyper$gamma_lambda <- min(0.01, hyper$t^(-0.5))
#     hyper$lambda_lambda[index] <- hyper$lambda_lambda[index] * exp(hyper$gamma_lambda * (min(exp(log_alpha), 1) - hyper$lambda_opt))
#     hyper$Sigma_lambda[index] <- hyper$Sigma_lambda[index] + hyper$gamma_lambda * ((init$lambdas[index] - hyper$mu_lambda[index])^2 - hyper$Sigma_lambda[index])
#     hyper$mu_lambda[index] <- hyper$mu_lambda[index] + hyper$gamma_lambda * (init$lambdas[index] - hyper$mu_lambda[index])
#   }
#   
#   return(list(init = init, hyper = hyper))
# }
