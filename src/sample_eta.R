# Gibbs step for the eta parameters -----

sample_eta <- function(init, hyper, index, X){
  hyper$Sigma_eta_new <- solve(init$k_l*init$Q_R[,,index] + hyper$Sigma_eta_inv)
  eta_matrix <- matrix(init$eta[index,], nrow = length(init$thetas), ncol = 2, byrow = T)
  hyper$mu_eta[index, ] <- t(hyper$Sigma_eta_new %*% init$Q_R[,,index]%*% t(X[,,index] - init$mean_i[,,index] + eta_matrix)%*%matrix(1, nrow = init$k_l))
  init$eta[index, ] <- mvrnorm(n = 1, mu = hyper$mu_eta[index, ], Sigma = hyper$Sigma_eta_new)
  
  
  eta_matrix <- matrix(init$eta[index,], nrow = length(init$thetas), ncol = 2, byrow = T)
  init$mean_i[,,index] <- init$alphas[index]*init$mean%*%init$R[,,index] + eta_matrix # compute the mean configuration for every unit
  init$residual[,,index] <- t(X[,,index] - init$mean_i[,,index])%*%(X[,,index] - init$mean_i[,,index])
}

