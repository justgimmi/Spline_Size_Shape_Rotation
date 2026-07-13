# Sample the angle parameters -----

# eval_log_lik_i <- function(lambda_candidate, hyper, init, i, X) {
#   # this function is just a wrap up to compute the contribution of the new element inside the log-lik
#   lambda_mod <- lambda_candidate %% (2*pi)
#   R <- matrix(c(cos(lambda_mod), -sin(lambda_mod),
#                 sin(lambda_mod), cos(lambda_mod)),
#               byrow = T, nrow = 2, ncol = 2)
#   eta_matrix <- matrix(init$eta[i,], nrow = length(init$thetas), ncol = 2, byrow = T)
#   mean_i <- init$alphas[i]*init$mean%*%R+ eta_matrix
#   Q_R <- (1/(init$alphas[i]^2)) *t(R)%*%init$Sigma_inv%*%R
#   residual <- X[,,i] - mean_i
# 
#   quad_new <- sum(residual[,1]^2) * Q_R[1,1] +
#     2 * sum(residual[,1] * residual[,2]) * Q_R[1,2] +
#     sum(residual[,2]^2) * Q_R[2,2]
# 
# 
#   return(list(log_lik = -0.5 * sum(quad_new), residual  = residual, Q_R = Q_R,
#               mean_i = mean_i, R = R))
# }
# 
# 
# 
# sample_lambda <- function(init, hyper, X, k_l){
#   # This function build the slice sampler for the angle parameters
#   n <- init$n # number of samples
#   k <- k_l # number of landmarks
#   # In this way, for each unit, we are decomposing the R^TSigma^-1R in three different elements
#   m <- hyper$m
#   w <- hyper$width_lambda
#   for (i in 1:n) {
# 
#     Q_R_11 <- init$Q_R[1,1,i]
#     Q_R_12 <- init$Q_R[1,2,i]
#     Q_R_22 <- init$Q_R[2,2,i]
# 
#     e <- rexp(1)
#     v <- runif(1)
#     left <- init$lambdas[i] - v * w
#     right <- left + w
#     bool <- TRUE
# 
#     v_old_x <- X[, 1, i] - init$mean_i[, 1, i] # compute the old value for the residuals x of unit i
#     v_old_y <- X[, 2, i] - init$mean_i[, 2, i] # compute the old value for the residuals y of unit i
#     quad_old <- sum(v_old_x^2) * Q_R_11 +
#       2 * sum(v_old_x * v_old_y) * Q_R_12 +
#       sum(v_old_y^2) * Q_R_22 # this is the contribution of the i-th landmark on the log-likelihood
#     # it is very easy to evaluate
# 
#     log_lik_i_old <- -0.5 * sum(quad_old)
#     log_thresh <- log_lik_i_old - e
#     u <- runif(1)
#     J <- floor(m*u)
#     up <- (m-1) - J
# 
#     while ((right - left) < 2*pi && eval_log_lik_i(left, hyper, init, i, X)$log_lik >= log_thresh && J > 0) {
#       left <- left - w
#       J <- J - 1
#     }
#     while ((right - left) < 2*pi && eval_log_lik_i(right, hyper, init,i, X)$log_lik >= log_thresh && up > 0) {
# 
#       right <- right + w
#       up <- up - 1
#     }
# 
#     if ((right - left) >= 2*pi) {
#       right <- left + 2*pi
#     }
# 
# 
#     while (bool == TRUE) {
#       lambda_raw <- left + runif(n = 1)*(right - left)
# 
#       res_cand <- eval_log_lik_i(lambda_raw, hyper, init, i, X)
# 
#       log_lambda_prop <- res_cand$log_lik
# 
#       if (log_lambda_prop >= log_thresh) {
#         init$lambdas[i] <- lambda_raw  %% (2*pi)
#         init$mean_i[, , i] <- res_cand$mean_i
#         init$R[,,i] <- res_cand$R
#         init$residual[, , i] <- t(res_cand$residual)%*%res_cand$residual
#         init$Q_R[,,i] <- res_cand$Q_R
#         hyper$log_lik <- hyper$log_lik +  res_cand$log_lik - log_lik_i_old
#         bool <- FALSE
#       } else {
#         if (lambda_raw < init$lambdas[i]) {
#           left <- lambda_raw
#         } else {
#           right <- lambda_raw
#         }
#       }
#     }
#   }
# }


sample_lambda_eta <- function(init, hyper, index, X, burn = TRUE){
  # the idea of this function is to sample lambda and eta all together
  # we sample lambda from either a uniform or a von-mises distribution 
  # while we sample eta from the full conditional.
  # After a bit of math, it is possible to observe that the Metropolis-Hastings ratio
  # depend on the ratio between \pi(\lambda|... except eta) computed at the proposed and the current value
  
  k_l <- init$k_l # number of landmarks
  alpha_i <- init$alphas[index] # size per unit i
  lambda_curr <- init$lambdas[index] # current lambda i
  eta_curr <- init$eta[index, ] # current eta i
  lambda_prop <- lambda_curr + rnorm(n = 1, mean = 0, sd = sqrt(hyper$lambda_lambda[index] *  hyper$Sigma_lambda[index])) # proposed lambda
  #lambda_star <- runif(n = 1, 0, 2*pi)
  lambda_star <- lambda_prop
  R_mat_star <- matrix(c(cos(lambda_star), -sin(lambda_star),
                    sin(lambda_star),  cos(lambda_star)), byrow = TRUE, nrow = 2) # proposed Rotation
  Q_R_star <- (1/alpha_i^2) * t(R_mat_star) %*% init$Sigma_inv %*% R_mat_star # proposed inverse 
  mean_res_star <- alpha_i * (init$mean %*% R_mat_star) # proposed mean configuration
  tilde_E_star <- X[,,index] - mean_res_star
  quad_star <- sum((tilde_E_star %*% Q_R_star) * tilde_E_star)
  Sigma_star <- solve(k_l*Q_R_star + hyper$Sigma_eta_inv)
  mu_eta_star <- as.matrix(t(Sigma_star %*%Q_R_star %*%t(tilde_E_star) %*% matrix(1, nrow = k_l)))
  quad_eta_star <- as.numeric(
    mu_eta_star %*%
      solve(Sigma_star) %*%
      t(mu_eta_star)
  )
  log_prop <- 0.5*log((det(Sigma_star))) - 0.5*quad_star + 0.5*quad_eta_star
  
  Sigma_curr <- solve(init$k_l*init$Q_R[,,index] + hyper$Sigma_eta_inv)
  eta_matrix <- matrix(init$eta[index,], nrow = length(init$thetas), ncol = 2, byrow = T)
  tilde_E_curr <- X[,,index] - init$mean_i[,,index] + eta_matrix
  quad_curr <- sum((tilde_E_curr %*% init$Q_R[,,index]) * tilde_E_curr)
  mu_eta_curr <- as.matrix(t(Sigma_curr %*%init$Q_R[,,index] %*%t(tilde_E_curr) %*% matrix(1, nrow = k_l)))
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
    #mu_eta <- t(Sigma_star %*% init$Q_R[,,index]%*% t(tilde_E_star)%*%matrix(1, nrow = init$k_l))
    
    init$eta[index, ] <- mvrnorm(n = 1, mu = mu_eta_star, Sigma = Sigma_star)
    init$R[,,index] <- R_mat_star
    init$mean_i[,,index] <- mean_res_star + matrix(init$eta[index, ], 
                                                   nrow = k_l, ncol = 2, byrow = TRUE)
    
    init$residual[,,index] <- t(tilde_E_star - matrix(init$eta[index, ], 
                                                      nrow = k_l, ncol = 2, byrow = TRUE)) %*% (tilde_E_star - matrix(init$eta[index, ], 
                                                                                                                   nrow = k_l, ncol = 2, byrow = TRUE))
    #print("accepted") 
    hyper$lambda_acc[index] <- hyper$lambda_acc[index] + 1
  }
  
  if (burn == TRUE) {
    
    hyper$gamma_lambda <- min(0.01, hyper$t^(-0.5))
    hyper$lambda_lambda[index] <- hyper$lambda_lambda[index] * exp(hyper$gamma_lambda * (min(exp(log_alpha), 1) - hyper$lambda_opt))
    # hyper$lambda_lambda[index] <- hyper$lambda_lambda[index]*exp(hyper$gamma_lambda *(min(exp(log_alpha), 1) - hyper$lambda_opt))
    hyper$Sigma_lambda[index] <- hyper$Sigma_lambda[index] + hyper$gamma_lambda *((init$lambdas[index] - hyper$mu_lambda[index])^2 - hyper$Sigma_lambda[index])
    hyper$mu_lambda[index] <- hyper$mu_lambda[index] + hyper$gamma_lambda *(init$lambdas[index] - hyper$mu_lambda[index])
    #print(min(exp(log_alpha), 1) )
  }
  
  
}
