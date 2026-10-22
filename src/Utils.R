require(here)
source(file = here("src/Packages.R"))

# Basis Construction ------
# Here we define a function to define the cubic B-splines 

Basis_Construction <- function(thetas, L, degree = 3, derivatives = FALSE){
  # I do not like to use already implemented function because they are based
  # on the range of the vector to define the knots so every iteration of the MCMC
  # we could have different basis following this approach

  # thetas: --> vector of theta values where to evaluate the basis
  # L: --> how many internal knots do we want?
  # degree: --> spline degree

  # The idea of this function is the following. Given the fact that I want the starting value in 0
  # and the end point in 2 \pi, I first define equispaced knots on this interval. Then, we have to
  # add padding on the left and on the right. A conservative choice would be to add degree times 0 before
  # and degree times 2\pi after. Talking with Johannes, I understood that the best approach to have
  # meaningfull penalties is to have equispaced knots also on the left and right part of the interval
  delta <- (2 * pi) / L # equispaced shift

  knots_base <- seq(0, 2 * pi, by = delta)[-c(1, L+1)] # base knots

  total_knots <- c(seq(0 - degree*delta, 0, by = delta), knots_base, seq(2*pi , 2*pi + degree*delta,
                                                                         by = delta))
  # Basis <- splineDesign(total_knots, thetas, ord = degree + 1, outer.ok = FALSE,
  #                       derivs = F) # Basis construction
  
  if (derivatives == FALSE) {
    Basis <- splineDesign(total_knots, thetas, ord = degree + 1, outer.ok = FALSE,
                          derivs = F) # Basis construction
  }
  
  else{
    Basis <- splineDesign(total_knots, thetas, ord = degree + 1, outer.ok = FALSE,
                          derivs = derivatives) # Basis construction
    
  }
  # At the End we obtain a matrix of size length(theta) x (L + degree)
  return(Basis)
}
# 

# Useful function for K1 and K2 -----
K1_construction <- function(J){ # smoothness constraint
  # J: --> number of basis 
  I <- diag(1, J) # define identity matrix
  D2  <-  matrix(0, J-2, J) # define empty matrix
  for (i in 1:nrow(D2)) {
    D2[i, i:(i+2)] <- c(-1, 2, -1)
  }
  return(t(D2) %*% D2)
  
}

K2_construction <- function(J){ # symmetry constraint
  # J: --> number of basis 
  I <- diag(1, J) # define identity matrix
  R_J <- matrix(0, J, J) # reverse identity matrix 
  for (i in 1:nrow(R_J)) {
    R_J[i, J - (i-1)] <- 1
  } 
  return(I - R_J)
}

# Simulation Function ----

# First of all we need a function to sample from a rank deficient normal distribution 
rmvnorm_rd <- function(n, mu, Precision, tol) {
  
  # n: --> how many samples?
  # mu: --> mean of the normal
  # Precison: --> Precision matrix
  # tol: --> degree of tolerance 
  
  eig <- eigen(Precision, symmetric = TRUE)
  
  keep <- eig$values > tol # define which col to take
  
  U <- eig$vectors[, keep, drop = FALSE] # consider just the first keep columns
  lambda <- eig$values[keep] # consider lambda
  
  r <- length(lambda) # rank of the Precision matrix
  
  Z <- matrix(rnorm(n * r), nrow = n) # sample from the normal (0, I)
  
  X <- t(U%*%diag(1 / sqrt(lambda), r) %*%t(Z)) # reproject back
  X <- sweep(X, 2, mu, "+") # add the mean
  
  return(as.vector(X))
}



in_model_sample <- function(n = 1, K_l,  thetas = NA, n_int_knots, degree = 3, tau = c(0.1, 0.1),sigma_2 = 0.01,
                            beta_values = NA, lambda = NA, alphas = NA, eta = NA, 
                            Sigma_e = NA, indep_error = TRUE){ # this functions can be used to sample directly from the model
  # K_l --> lanmdarks number 
  # thetas --> vector of angles
  # n_int_knots --> number of internal_knots
  # degree --> degree of the polynomial
  # tau --> strength of the penalty
  # thetas block ----
  if(length(thetas) == 1){
    thetas <- runif(K_l) *2*pi # sample angles if not given
  }
  thetas[thetas < 0]  <- thetas[thetas < 0] + 2*pi # be aware of the domain [0, 2pi]
  thetas <- sort(thetas) # sort the angles
  
  # basis block -----
  n_basis     <- n_int_knots + degree # final number of basis
  B_sim <- Basis_Construction(thetas, L = n_int_knots, degree = degree) # Basis lenght(thetas) X n_basis
  K1 <- K1_construction(n_basis) # smoothness
  K2 <- K2_construction(n_basis) # symmetry 
  K <- K1 + K2
  P <- (1/(tau[1]))*K1 + (1/(tau[2]))*K2  # Precision Matrix
  tol <- 1e-6
  L <- qr(K)$rank
  eig <- eigen(K, symmetric = TRUE)
  ord <- order(eig$values, decreasing = TRUE)
  eig_vals <- eig$values[ord]
  eig_vecs <- eig$vectors[, ord, drop = FALSE]
  keep <- eig_vals > tol
  #hyper$L <- sum(keep) # approximation for the normal distribution ---> non null eigenvalues
  Eigen_matrix <- diag(1/sqrt(eig_vals[keep]), nrow = L) # using for sample from the prior
  Eigen_vector <- eig_vecs[, keep, drop = FALSE] # Eigen vectors with non null eigen values
  Eigen_vector_null <- t(eig_vecs[, !keep, drop = FALSE]) # eigen vectors eith null eigen values
  R <- Basis_Construction(0, n_int_knots, degree) - Basis_Construction(2*pi, n_int_knots, degree) # (B(0) - B(2pi))
  R1 <- Basis_Construction(0, n_int_knots, degree, derivatives = 1) - Basis_Construction(2*pi, n_int_knots, degree, derivatives = 1) # (B'(0) - B'(2pi))
  R2 <- Basis_Construction(0, n_int_knots, degree, derivatives = 2) - Basis_Construction(2*pi, n_int_knots, degree, derivatives = 2)
  A <- rbind(Eigen_vector_null, R, R1, R2)
  # A <- rbind(Eigen_vector_null, R)
  aux_eig <- eigen(t(A) %*% A)
  A_bar <- t(aux_eig$vectors[,aux_eig$values < tol])
  union_matrix <- rbind(A, A_bar)
  inv_union_matrix <- solve(union_matrix)
  C_bar <- inv_union_matrix[, (nrow(A)+1):n_basis ]
  K_gamma1 <- t(C_bar) %*% K1 %*% C_bar
  K_gamma2 <- t(C_bar) %*% K2 %*% C_bar
  L <- ncol(C_bar)
  Sigma_gamma_inv <-  (1/tau[1])*K_gamma1 + 
    (1/tau[2])*K_gamma2
  U_lambda <- t(chol(Sigma_gamma_inv))
  gamma_samp <- solve(t(U_lambda), rnorm(L))
  beta_values <- C_bar %*%  gamma_samp
  
  
  # if (length(beta_values) == 1 ) { # the user can supply the values of beta
  #   beta_values <- rmvnorm_rd(n = 1, mu = rep(0, n_basis), P, tol = 1e-7)
  #   
  # }
   # Mean block ----
  r <- exp(B_sim %*% as.vector(beta_values)) # radius
  mu_mean_x <- r*cos(thetas)
  mu_mean_y <- r*sin(thetas)
  mu_mean <- cbind(mu_mean_x, mu_mean_y) # mean configuration 
  k_free <- K_l - 3
  S_mat <- S_mat_constructor(mu_mean)
  D_p <- D_p_constructor(mu_mean, K_l)
  V <- median(diag(D_p %*% S_mat %*% solve(t(D_p) %*% S_mat %*% D_p) %*% S_mat %*% t(D_p)))
  psi <- array(0, dim = c(k_free, 2, n))
  if (sigma_2 == 0 && indep_error == TRUE) {
    Precision <- diag(1, nrow = k_free)
    U_prec <- chol(Precision)
  }
  
  else if(sigma_2 != 0 && indep_error == TRUE){
    Precision <- diag(1/sigma_2, nrow = k_free)
    U_prec <- chol(Precision)
    k_free <- ncol(D_p)  # k - 5
    for (i in 1:n) {
      Z <- matrix(rnorm(k_free * 2), nrow = k_free, ncol = 2)
      psi[,,i] <- backsolve(U_prec, Z)
    }
  }
  else{
    Precision <- (1 / sigma_2) * t(D_p) %*% S_mat %*% D_p
    U_prec <- chol(Precision)
    
    for (i in 1:n) {
      Z <- matrix(rnorm(k_free * 2), nrow = k_free, ncol = 2)
      psi[,,i] <- backsolve(U_prec, Z)
    }
  }
  W <- tensor(D_p, psi, 2, 1)

  R <- array(NA, dim = c(2, 2, n))
  mu_i <- array(NA, dim = c(K_l, 2, n))
  if (length(lambda) == 1) {
    lambda <- runif(n = n)*2*pi # generate rotation angles
  }
  
  if (length(alphas) == 1) {
    alphas <- rgamma(n = n, shape = 3, rate = 0.5) # sample the size effect
    
  }
  
  if (length(eta) == 1) {
    eta <- matrix(rnorm(n = n*2, sd = 8), ncol = 2) # sample the translation effect 
    
  }
  
  for (i in 1:n) { # compute rotation matrix 
    R[,,i] <- matrix(c(cos(lambda[i]), -sin(lambda[i]), sin(lambda[i]), cos(lambda[i])), 
                     byrow = T, nrow = 2, ncol = 2)
    eta_matrix <- matrix(eta[i,], nrow = K_l, ncol = 2, byrow = T)
    mu_i[,,i] <- alphas[i]*(mu_mean + S_mat %*% W[,,i]) %*%R[,,i] + eta_matrix # compute the mean configuration for every unit
  }
  
  s_x2 <- 9 * 1e-4
  s_y2 <- 9 * 1e-4
  rho <- 0.2 
  s_xy <- rho * sqrt(s_x2 * s_y2)
  S <- matrix(c(s_x2, s_xy,
                s_xy, s_y2), nrow = 2, ncol = 2)
  Sigma_e <- riwish(4, S)
  #Sigma_e <- diag(2)*0.00035
  X <- array(NA, dim = c(K_l, 2, n))
  for (i in 1:n) {
    # to sample we use the fact that X \sim MN(mu_{i}, I_k, Sigma)
    Sigmas <- alphas[i]^2 *t(R[,,i])%*%Sigma_e%*%R[,,i]
    I <- diag(K_l)
    X[,,i] <- matrix(rnorm(n = K_l * 2), nrow = K_l, ncol = 2)%*%chol(Sigmas) + mu_i[,,i]
    #X[,,i] <- chol(I)%*%matrix(rnorm(n = K_l * 2), nrow = K_l, ncol = 2)%*%t(chol(Sigma_e)) +  mu_i[,,i]
    }
  
  
  
  
  final_param <- list()
  final_param$mu <- mu_mean
  final_param$theta <- thetas
  final_param$r <- r
  final_param$tau <- tau
  final_param$degree <- degree
  final_param$K <- K
  final_param$n_int_knots <- n_int_knots
  final_param$R <- R
  final_param$lambda <- lambda
  final_param$alphas <- alphas
  final_param$eta <- eta
  final_param$mu_i <- mu_i
  final_param$Sigma_e <- Sigma_e
  final_param$X <- X
  final_param$betas <- beta_values
  final_param$psi <- psi
  final_param$S_mat <- S_mat
  final_param$W <- W
  final_param$D_p <- D_p
  final_param$sigma_2 <- sigma_2
  final_param$sigma_2_epur <- sigma_2/V
  return(final_param)

}



sample_from_mean <- function(mean_shape, n, sigma_2, plot = TRUE){ # this functions can be used to sample directly from the model
  # mean_shape: star-shaped-mean given by the user
  # n: sample size
  k <- dim(mean_shape)[1]
  
  k_free <- k - 3 # free values 
  S_mat <- S_mat_constructor(mean_shape) # k x k matrix 
  D_p <- D_p_constructor(mean_shape, k) # k x k-3 matrix 
  V <- median(diag(S_mat %*% D_p %*% solve(t(D_p) %*% S_mat %*% D_p) %*% t(D_p) %*% S_mat)) # define V
  psi <- array(0, dim = c(k_free, 2, n)) # define where to save psi for each unit
  
  if (sigma_2 == 0) {
    Precision <- diag(1, nrow = k_free, ncol = k_free)
  }
  else{
    Precision <- (1 / sigma_2) * t(D_p) %*% S_mat %*% D_p # precision of \psi_i
  }
  U_prec <- chol(Precision) # upper triangular matrix
  
  for (i in 1:n) {
    Z <- matrix(rnorm(k_free * 2), nrow = k_free, ncol = 2)
    psi[,,i] <- backsolve(U_prec, Z)
  }
  W <- tensor(D_p, psi, 2, 1) # define W_i for each unit
  R <- array(NA, dim = c(2, 2, n)) 
  mu_i <- array(NA, dim = c(k, 2, n))
  lambda <- runif(n = n)*2*pi #generate rotation angles completely at random
  alphas <- rgamma(n = n, shape = 3, rate = 0.5) # sample the size effect
  eta <- matrix(rnorm(n = n*2, sd = 2), ncol = 2) # sample the translation effect 
  if (sigma_2 != 0) {
    for (i in 1:n) { # compute rotation matrix 
      R[,,i] <- Rmat(lambda[i])
      eta_matrix <- matrix(eta[i,], nrow = k, ncol = 2, byrow = T)
      mu_i[,,i] <- alphas[i]*(mean_shape + S_mat %*% W[,,i]) %*%R[,,i] + eta_matrix # compute the mean configuration for every unit
    }
    
  }
  
  else{
    for (i in 1:n) { # compute rotation matrix 
      R[,,i] <- Rmat(lambda[i])
      eta_matrix <- matrix(eta[i,], nrow = k, ncol = 2, byrow = T)
      mu_i[,,i] <- alphas[i]*mean_shape%*%R[,,i] + eta_matrix # compute the mean configuration for every unit
    }
    
  }
  s_x2 <- 2 * 1e-5
  s_y2 <- 2 * 1e-5
  rho <- 0.2 
  s_xy <- rho * sqrt(s_x2 * s_y2)
  S <- matrix(c(s_x2, s_xy,
                s_xy, s_y2), nrow = 2, ncol = 2)
  Sigma_e <- riwish(4, S)
  X <- array(NA, dim = c(k, 2, n))
  for (i in 1:n) {
    # to sample we use the fact that X \sim MN(mu_{i}, I_k, Sigma)
    Sigmas <- alphas[i]^2 *t(R[,,i])%*%Sigma_e%*%R[,,i]
    X[,,i] <- matrix(rnorm(n = k * 2), nrow = k, ncol = 2)%*%chol(Sigmas) + mu_i[,,i]
  }
  
  if (plot) {
    rot <- mean_shape
    mu_i <- mu_i[,,1]
    X_i <- X[,,1]
    par(mfrow = c(1, 3))
    plot(rot, asp = 1, pch = 21, bg = "darkgreen", main = "Latent mu")
    lines(c(rot[,1], rot[1,1]), c(rot[,2], rot[1,2]), col = "forestgreen", lwd = 2)
    
    plot(mu_i, asp = 1, pch = 21, bg = "darkred", main = "mu for the first unit")
    lines(c(mu_i[,1], mu_i[1,1]), c(mu_i[,2], mu_i[1,2]), col = "red", lwd = 2)
    
    
    plot(X_i, asp = 1, pch = 21, bg = "darkred", main = " X_1 unit")
    lines(c(X_i[,1], X_i[1,1]), c(X_i[,2], X_i[1,2]), col = "red", lwd = 2)
    
  }
  
  
  final_param <- list()
  final_param$mu <- mean_shape
  final_param$lambda <- lambda
  final_param$alphas <- alphas
  final_param$eta <- eta
  final_param$mu_i <- mu_i
  final_param$Sigma_e <- Sigma_e
  final_param$X <- X
  final_param$psi <- psi
  final_param$S_mat <- S_mat
  final_param$W <- W
  final_param$D_p <- D_p
  final_param$sigma_2 <- sigma_2
  final_param$sigma_2_epur <- sigma_2/V
  return(final_param)
  
}

# MCMC Util functions ----
Rmat <- function(angle) { # rotation matrix 
  matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)), nrow = 2, ncol = 2, byrow = TRUE)
}


# S_mat_constructor <- function(mu){
#   # mu: value of the latent mean K \times K
#   # this function build a conditional positive definite matrix 
#   mat <- dist(mu)
#   mat[mat > 0] <- mat[mat > 0]^2 * log(mat[mat > 0])
#   return(as.matrix(mat))
# }

S_mat_constructor <- function(mu){
  # mu: value of the latent mean K \times K
  # this function build a conditional positive definite matrix 
  dx <- outer(mu[, 1], mu[, 1], "-")
  dy <- outer(mu[, 2], mu[, 2], "-")
  d2 <- dx * dx + dy * dy
  0.5 * d2 * log(d2 + (d2 == 0))
}

# D_p_constructor <- function(mu){
#   k_l <- nrow(mu)
#   ones <- matrix(1, nrow = 1, ncol = k_l)
#   B <- rbind(ones, t(mu))
#   mat <- t(B) %*% B
#   eigs <- eigen(mat, symmetric = TRUE)
#   # take the k_l - 3 smallest eigenvalues, by rank, not by absolute threshold
#   ord <- order(eigs$values)                # ascending
#   eg_vec <- eigs$vectors[, ord[1:(k_l - 3)], drop = FALSE]
#   return(eg_vec)
# }

# D_p_raw <- function(mu) {
#   k_l <- nrow(mu)
#   ones <- matrix(1, nrow = 1, ncol = k_l)
#   B <- rbind(ones, t(mu))
#   
#   mat <- t(B) %*% B
#   eigs <- eigen(mat, symmetric = TRUE)
#   eg_vec <- eigs$vectors[,eigs$values< 1e-10, drop = FALSE] 
#   
#   matts <- rbind(B , t(eg_vec))
#   inv <- solve(matts)
#   inv
#   
#   C_bar_w <- inv[, (nrow(B) + 1):k_l]
#   return(C_bar_w)
#   
# }
# D_p_raw <- function(mu, k) {
#     A <- rbind(1, t(mu))          # affine + zero mean displacement (3 x k)
#     sv <- svd(A, nu = 0, nv = k)
#     sv$v[, 4:k, drop = FALSE]                          # k - 3 columns
# 
# }


D_p_constructor <- function(mu, k) {
  # this function is used to define the basis function for the warping parameters
  # we know that the svd could have a flipping sign or whatever but we align everything a posteriori
  A <- rbind(1, t(mu))          # zero mean displacement (3 x k)
  sv <- svd(A, nu = 0, nv = k)  # svd 
  sv$v[, 4:k, drop = FALSE]   # free dimensions
  #qr.Q(qr(t(A)), complete = TRUE)[, 4:k, drop = FALSE]
}
# covariance_mat <- function(mat, thetas){
#   # thetas: vector of angles
#   # mat: K x K covariance matrix
#   # phi: scale parameter of the exponential covariance function
#   k_l <- length(thetas)
#   for (i in 1:k_l) {
#     for (j in 1:k_l) {
#       #d <- min(2*pi - abs(thetas[i] - thetas[j]), abs(thetas[i] - thetas[j]))
#       mat[i, j] <- min(2*pi - abs(thetas[i] - thetas[j]), abs(thetas[i] - thetas[j]))
#       mat[j, i] <- mat[i, j]
#     }
#   }
#   return(mat)
# }

# covariance_mat <- function(mat, mu){
#   # thetas: vector of angles
#   # mat: K x K covariance matrix
#   # phi: scale parameter of the exponential covariance function
#   # k_l <- nrow(mu)
#   # for (i in 1:k_l) {
#   #   for (j in 1:k_l) {
#   #     #d <- min(2*pi - abs(thetas[i] - thetas[j]), abs(thetas[i] - thetas[j]))
#   #     mat[i, j] <- min(2*pi - abs(thetas[i] - thetas[j]), abs(thetas[i] - thetas[j]))
#   #     mat[j, i] <- mat[i, j]
#   #   }
#   # }
#   D  <- as.matrix(dist(mu))
#   k <- nrow(mu)
#   # D <- D / mean(D[cbind(1:k, c(2:k, 1))])
#   return(D)
# }



init_param <- function(tau, lambdas, thetas, betas, Sigma, eta, alphas, gammas, phi, n){
  # n --> number of units 
  param <- list()
  param$n <- n # number of units
  param$k_l <- length(thetas) # number of landmarks 
  param$tau <- tau # 1x2 vector of variances 
  param$lambdas <- lambdas # 1 x k vector of rotation parameters
  param$thetas <- thetas # 1 x k vector of angles
  param$gammas <- gammas # 1 x (k-3) vector of free spline coefficients
  param$gaps <- c(thetas[1], diff(thetas), 2*pi- thetas[param$k_l])/(2*pi) # vector of k + 1 gaps
  param$omega_theta <- log(param$gaps) - log(param$gaps[param$k_l + 1]) # vector of k + 1 omegas
  param$betas <- as.matrix(betas) # 1 x k vector of spline coefficients
  param$Sigma <- Sigma # 2 x 2 variance-covariance matrix
  param$eta <- eta # n x 2 matrix of locations parameters
  param$eta_array <- array(0, dim = c(length(thetas), 2, n))  #k x 2 x n array of shifting locations
  param$alphas <- alphas # 1 x n vector of size parameters 
  
  param$r <- numeric(length(thetas)) # deterministic vector of radius
  param$mean <- matrix(0, nrow = length(thetas), ncol = 2) # deterministic matrix of star shape mean
  param$mean_i <- array(0, dim = c(length(thetas), 2, n))  # deterministic array of means
  param$R <- array(NA, dim = c(2, 2, n)) #deterministic array of rotation matrices
  
  param$Sigma_inv <- solve(param$Sigma) # deterministic inverse 
  param$Q_R <- array(NA, dim = c(2, 2, n)) # deterministic  R^T \Sigma^-1 R
  param$residual <-  array(NA, dim = c(2, 2, n)) #deterministic (X_i - \mu_i)^T(X_i - \mu_i) 
  param$sigma_2 <- 1e-4 # temperature of the warping 
  param$psi <- array(NA, dim = c(length(thetas) - 3, 2, n)) # k_free x 2 x n array of free warping
  param$psi_flatten <- matrix(param$psi, nrow = length(thetas) - 3) # flatten version of psi
  param$W <- array(NA, dim = c(length(thetas), 2, n)) # k x 2 x n array of constrained warping
  return(param)
  
}

hyperparameters <- function(a_tau, b_tau, n_basis, degree, nu, psi, n, a, b, Sigma_eta,
                            tol = 1e-10,a_phi,b_phi,X){
  
  hyper <- list()
  # tau block -----
  hyper$a_tau <- a_tau
  hyper$b_tau <- b_tau
  hyper$a_new <- 0
  hyper$b_new <- 0
  
  # Basis block -----
  hyper$n_basis <- n_basis
  hyper$degree <- degree
  J <- n_basis + degree # find the number of basis
  #J <- n_basis
  
  # tau1 and tau2 block
  hyper$K1 <- K1_construction(J) # construct K1
  hyper$K2 <- K2_construction(J) #construct K2
  hyper$K <- (hyper$K1 + hyper$K2)
  hyper$L <- qr(hyper$K)$rank
  eig <- eigen(hyper$K, symmetric = TRUE)
  ord <- order(eig$values, decreasing = TRUE)
  eig_vals <- eig$values[ord]
  eig_vecs <- eig$vectors[, ord, drop = FALSE]
  keep <- eig_vals > tol
  #hyper$L <- sum(keep) # approximation for the normal distribution ---> non null eigenvalues
  hyper$Eigen_matrix <- diag(1/sqrt(eig_vals[keep]), nrow = hyper$L) # using for sample from the prior
  hyper$Eigen_vector <- eig_vecs[, keep, drop = FALSE] # Eigen vectors with non null eigen values
  hyper$Eigen_vector_null <- t(eig_vecs[, !keep, drop = FALSE]) # eigen vectors eith null eigen values
  hyper$R <- Basis_Construction(0, n_basis, degree) - Basis_Construction(2*pi, n_basis, degree)
  hyper$R1 <- Basis_Construction(0, n_basis, degree, derivatives = 1) - Basis_Construction(2*pi, n_basis, degree, derivatives = 1) # (B'(0) - B'(2pi))
  hyper$R2 <- Basis_Construction(0, n_basis, degree, derivatives = 2) - Basis_Construction(2*pi, n_basis, degree, derivatives = 2)
  #hyper$A <- rbind(hyper$Eigen_vector_null, hyper$R)
  hyper$A <- rbind(hyper$Eigen_vector_null, hyper$R, hyper$R1,hyper$R2)
  aux_eig <- eigen(t(hyper$A) %*% hyper$A)
  hyper$A_bar <- t(aux_eig$vectors[,aux_eig$values < tol])
  union_matrix <- rbind(hyper$A, hyper$A_bar)
  inv_union_matrix <- solve(union_matrix)
  hyper$C_bar <- inv_union_matrix[, (nrow(hyper$A)+1):J ]
  hyper$L <- ncol(hyper$C_bar)
  #hyper$C_bar <- hyper$Eigen_vector %*% hyper$Eigen_matrix # \beta = C_bar \gamma
  #hyper$A_bar <- diag(sqrt(eig_vals[keep]), nrow = hyper$L) %*% t(hyper$Eigen_vector) # \gamma = A_bar \beta
  hyper$U_lambda <- t(chol(t(hyper$C_bar) %*% hyper$K %*% hyper$C_bar))
  hyper$K_gamma1 <- t(hyper$C_bar) %*% hyper$K1 %*% hyper$C_bar
  hyper$K_gamma2 <- t(hyper$C_bar) %*% hyper$K2 %*% hyper$C_bar
  hyper$lambda_tau <- (2.38^2)/2
  hyper$mu_tau <- matrix(0, nrow = 2)
  hyper$Sigma_tau <- diag(0.1, nrow = 2)
  hyper$a_opts <- 0.234
  hyper$tau_acc <- 0
  
  # Theta block -----

  hyper$log_lik <- 0
  k_l <- dim(X)[1]
  
  # Sigma block ----
  hyper$nu <- nu
  hyper$psi <- psi
  hyper$nu_post <- 0
  hyper$psi_post <- psi
  
  # Alphas block -----
  hyper$lambda_a <- matrix(2.38^2, nrow = n)
  hyper$mu_a <- matrix(0, nrow = n)
  hyper$Sigma_a <- matrix(0.1, nrow = n)
  #hyper$L_a <- chol(hyper$lambda_a*hyper$Sigma_a)
  hyper$alpha_acc <- matrix(0, nrow = n)
  hyper$gamma_a <- 1
  hyper$a_opt <- 0.44
  hyper$a_aver_a <- numeric(n)
  hyper$t <- 0
  hyper$a <- a
  hyper$b <- b
  hyper$log_lik_a <- numeric(n)
  
  
  # Eta block ----
  hyper$Sigma_eta_inv <- solve(Sigma_eta)
  hyper$mu_eta <- matrix(0, nrow = n, ncol = 2)
  hyper$Sigma_eta_new <- Sigma_eta
  hyper$mu_new <- matrix(0, nrow = n, ncol = 2)
  
  hyper$lambda_lambda <- matrix(2.38^2, nrow = n)
  hyper$mu_lambda <- matrix(0, nrow = n)
  hyper$Sigma_lambda <- matrix(0.1, nrow = n)
  hyper$gamma_lambda <- 1
  hyper$lambda_opt <- 0.44
  hyper$lambda_acc <- matrix(0, nrow = n)
  
  # phi block
  hyper$lambda_phi <- 2.38^2
  hyper$mu_phi <- 0
  hyper$Sigma_phi <- 0.1
  #hyper$L_a <- chol(hyper$lambda_a*hyper$Sigma_a)
  hyper$alpha_phi <- 0
  hyper$gamma_phi <- 1
  hyper$a_opt <- 0.44
  hyper$a_phi <- a_phi
  hyper$b_phi <- b_phi
  hyper$phi_acc <- 0
  
  
  hyper$S_mat <- matrix(0, nrow = k_l, ncol = k_l)
  hyper$a_sigma <- 1
  hyper$lambda_sigma_cp <- 2.38^2
  hyper$Sigma_sigma_cp  <- 0.1
  hyper$mu_sigma_cp     <- 0
  hyper$sigma_acc_cp    <- 0
  
  hyper$lambda_sigma_ncp <- 2.38^2
  hyper$Sigma_sigma_ncp  <- 0.1
  hyper$mu_sigma_ncp     <- 0
  hyper$sigma_acc_ncp    <- 0
  
  
  hyper$a_chi <- 2
  hyper$b_chi <- 2
  hyper$a_s2 <- 1
  hyper$b_s2 <- 1
  return(hyper)
}


g_phi <- function(phi_star, a_phi, b_phi){
  phi <- (exp(phi_star)/(1 + exp(phi_star)))*(b_phi - a_phi) + a_phi
  return(phi)
  
}

# project_warp <- function(S_mat, W, thetas) {
#   #t_mat <- cbind(-sin(thetas), cos(thetas))
#   S_mat %*% W
#   #s <- rowSums(Delta * t_mat)
#   #s * t_mat   # k x 2, row h = s_h * t_h
# }

g_phi_star <- function(phi, a_phi, b_phi){
  app <- (phi - a_phi)/(b_phi - a_phi)
  phi_star <- log(app/(1 - app))
  return(phi_star)
  
}

create_output <- function(mcmc_iter, k_l, n, L){
  output <- list()
  output$tau <- matrix(NA, nrow = mcmc_iter, ncol = 2)
  output$theta <- matrix(NA, nrow = mcmc_iter, ncol = k_l)
  output$lambdas <- matrix(NA, nrow = mcmc_iter, ncol = n)
  output$Sigma <- array(NA, dim = c(mcmc_iter, 2, 2))
  output$alphas <- matrix(NA, nrow = mcmc_iter, ncol = n)
  output$betas <- matrix(NA, nrow = mcmc_iter, ncol = L)
  output$gammas <- matrix(NA, nrow = mcmc_iter, ncol = L - 4)
  output$eta <- array(NA, c(mcmc_iter, n, 2))
  output$mean_i <- list()
  output$pred_i <- list()
  output$mean <- array(NA, c(mcmc_iter, k_l, 2))
  output$phi <- matrix(NA, nrow = mcmc_iter, ncol =1)
  output$psi <- list()
  output$W <- list()
  output$S <- list()
  output$SW <- list()
  # output$psi <- array(NA, dim = c(mcmc_iter, k_l - 5, 2))
  output$sigma_2 <- matrix(NA, nrow = mcmc_iter, ncol = 1)
  output$sigma_2_epurated <- matrix(NA, nrow = mcmc_iter, ncol = 1)
  output$dis <- matrix(NA, nrow = mcmc_iter, ncol = n)
  #output$z <- matrix(NA, nrow = mcmc_iter, ncol = 1)
  return(output)
}



compute_proc_dist <- function(init, X, X_pred){
  n <- dim(X)[3]
  k_l <- dim(X)[1]
  p <- 2
  Sigma_e <- init$Sigma
  dis <- numeric(n)
  for (i in 1:n) {
    M_subject <- init$mean_i[,,i]
    
    # 2) Affine component
    alpha_i   <- init$alphas[i]
    lambda_i  <- init$lambdas[i]
    
    # 3) Rotation matrix 
    R_lambda <- Rmat(lambda_i)
    
    # 4) Spatial covariance
    Sigma_spatial_i <- (alpha_i^2) * (t(R_lambda) %*% Sigma_e %*% R_lambda)
    L_Sigma_spatial <- t(chol(Sigma_spatial_i)) 
    Z <- matrix(rnorm(k_l * p), nrow = k_l, ncol = p)
    # E_i <- chol(C) %*% Z %*% t(L_Sigma_spatial)
    E_i <- Z %*% t(L_Sigma_spatial)
    X_pred[,,i]<- M_subject + E_i
    dis[i] <- sin(riemdist(X[,,i],X_pred[,,i]))
  }
  return(list(dis = dis, X_pred = X_pred))
}


mean_constructor <- function(n_basis, degree, init, hyper, X){
  # n_basis --> number of internal knots
  # degree --> degree of the polynomial
  
  n <- init$n
  B_sim <- Basis_Construction(init$thetas, n_basis, degree) # build the basis with the current value of theta
  init$r <- exp(B_sim %*% as.vector(init$betas))# compute the current value of the radius
  mu_mean_x <- init$r*cos(init$thetas)
  mu_mean_y <- init$r*sin(init$thetas)
  mean <- cbind(mu_mean_x, mu_mean_y)
  init$mean <- mean # average configuration 
  
  
  hyper$S_mat <- S_mat_constructor(mean)
  hyper$D_p <- D_p_constructor(mean, init$k_l)
  hyper$k_free <- ncol(hyper$D_p)
  Precision <- t(hyper$D_p) %*% hyper$S_mat %*% hyper$D_p
  K_inv <- chol2inv(chol(Precision)) # variance
  k_free <- ncol(hyper$D_p)
  # V <- diag(1/diag(K_inv), nrow = k_free, ncol =k_free) #precision 
  # V <- (1/init$sigma_2) * V
  U_prec <- chol(Precision/init$sigma_2)
    # k - 
  for (i in 1:n) { # compute rotation matrix 
    z_raw <- matrix(rnorm(2*k_free), k_free, 2)
    init$psi[,,i] <-  matrix(backsolve(U_prec, z_raw), nrow = k_free, ncol = 2)
    # init$psi[,,i] <- app
    #init$psi[,,i] <- matrix(0, nrow = k_free, ncol = 2 )
    
    init$W[,,i] <- hyper$D_p %*% init$psi[,,i]
    init$R[,,i] <- matrix(c(cos(init$lambdas[i]), -sin(init$lambdas[i]), 
                                  sin(init$lambdas[i]), cos(init$lambdas[i])), 
                                byrow = T, nrow = 2, ncol = 2)
    init$eta_array[,,i] <- matrix(init$eta[i,], nrow = length(init$thetas), ncol = 2, byrow = T)
    init$mean_i[,,i] <- init$alphas[i]*(init$mean + hyper$S_mat%*%init$W[,,i])%*%init$R[,,i] + init$eta_array[,,i] # compute the mean configuration for every unit
    #init_param$mean_i[,,i] <- init_param$alphas[i]*(init_param$mean + project_warp(hyper$S_mat, init_param$W, init_param$thetas))%*%init_param$R[,,i] +
    #  eta_matrix # compute the mean configuration for every unit
    
    init$Q_R[,,i] <- (1/(init$alphas[i]^2)) *t(init$R[,,i])%*%init$Sigma_inv%*%init$R[,,i] 
    
    res <- X[,,i] - init$mean_i[,,i]
    #init_param$residual[,,i] <- t(X[,,i] - init_param$mean_i[,,i])%*%init_param$C %*% (X[,,i] - init_param$mean_i[,,i])
    init$residual[,,i] <- t(res)%*%res
  }
  hyper$mean_i_array <- init$mean_i
  hyper$residual_array <- init$residual
  K <- t(hyper$D_p) %*% hyper$S_mat%*% hyper$D_p
  Kinv <- chol2inv(chol(K))
  w    <- 1 / diag(Kinv)
  w_mat <- diag(w, nrow = k_free)
  V <-  median(diag(hyper$S_mat %*% hyper$D_p %*% solve(K) %*% t(hyper$D_p) %*% t(hyper$S_mat)))
  hyper$V <- V
  # V <- 1
  # lp_sig <- -log(V) - log(1 + (init$sigma_2/(V*hyper$a_sigma))^2)
  lp_sig  <- hyper$a_s2 * log(V) - hyper$b_s2 * V / init$sigma_2
  logdetK <- as.numeric(determinant(K, logarithm = TRUE)$modulus)
  hyper$traces <- -0.5 * (1/init$sigma_2) *  trace_function(init$W, hyper$S_mat) +
      n*logdetK - n*k_free*log(init$sigma_2) + lp_sig
  #hyper$D_ref <- D_p_raw(mean, init$k_l)
  hyper$SW <- tensor(hyper$S_mat, init$W, 2, 1)# candidate value of SW
  # return(init_param)
}


# candidates <- function(init, hyper){ # this is the most useful object in this entire code
#   # the main problem of using elliptical sampling is that you do not know how many times you
#   # are going to call the density function and it is possible that in this way you create quite
#   # a lot of objects
#   cand <- list()
#   cand$mean_i_array <- init$mean_i
#   cand$residual_array <- init$residual
#   cand$theta_cand <- init$thetas
#   
#   
#   k    <- nrow(mean_mat)
#   S    <- S_mat_constructor(mean_mat)
#   D    <- D_p_constructor(mean_mat, k)
#   SD   <- S %*% D
#   K    <- crossprod(D, SD)
#   Lk   <- chol(K)
#   Kinv <- chol2inv(Lk)
# }


trace_function <- function(W, S_mat) {
  # this function is very useful to compute trace functions and it is surprisingly fast
  # it is used either to compute the \sum_{i = 1}^{n} tr(W_{i}^T S_mat W_i)
  # but also to compute the other versione \sum_{i = 1}^{n} tr(\psi_{i}^T K \psi_o)
  d <- dim(W) # k \times 2 \times n
  W_cat <- matrix(W, nrow = d[1]) # flatten the array in one direction to obtain k \times 2n
  sum(W_cat * (S_mat %*% W_cat)) # compute the wanted quantity
}

chol_det <- function(Lk){
  # given a upper/lower triangular matrix obtained via cholesky, it computes the corresponding log_determinad
  2*sum(log(diag(Lk)))
}



log_density <- function(init) {
  Sg  <- init$Sigma
  ldS <- log(Sg[1, 1] * Sg[2, 2] - Sg[1, 2] * Sg[2, 1])        # log det of a 2x2
  val   <- sum(init$Q_R * init$residual)                        # = sum_i tr(Q_i res_i), res_i symmetric
  val_2 <- -2 * init$k_l * sum(log(init$alphas)) - 0.5 * init$n * init$k_l * ldS
  -0.5 * val + val_2
}


extract_mean_i <- function(output_mcmc, subject_idx){
  n_saved_iter <- dim(output_mcmc$mean)[1] # number of iterations
  k_l <- dim(output_mcmc$mean)[2] # number of landmarks
  p <- 2     
  mean_i <- array(0, dim = c(k_l, p, n_saved_iter))
  pred_i <- array(0, dim = c(k_l, p, n_saved_iter))
  for (i in 1:n_saved_iter) {
    mean_i[,,i] <- output_mcmc$mean_i[[i]][, , subject_idx]
    pred_i[,,i] <- output_mcmc$pred_i[[i]][, , subject_idx]
  }
  return(list(mean_i = mean_i, pred_i = pred_i))
}

sample_posterior_predictive <- function(output_mcmc, subject_idx) { 
  # This function is quite useful because it allow a direct sampling from the posterior predictive distrubution 
  n_saved_iter <- dim(output_mcmc$mean)[1] # number of MCMC iterations 
  k_l <- dim(output_mcmc$mean)[2] # number of landmarks
  p <- 2                           # dimension
  X_pred <- array(0, dim = c(k_l, p, n_saved_iter)) # final object 
  for (i in 1:n_saved_iter) {
    # 1) Mean per subject
    M_subject <- output_mcmc$mean_i[[i]][, , subject_idx]
    
    # 2) Affine component
    alpha_i   <- output_mcmc$alphas[i, subject_idx]
    lambda_i  <- output_mcmc$lambdas[i, subject_idx]
    Sigma_e   <- output_mcmc$Sigma[i, , ]
    
    # 3) Rotation matrix 
    R_lambda <- Rmat(lambda_i)
    
    # 4) Spatial covariance
    Sigma_spatial_i <- (alpha_i^2) * (t(R_lambda) %*% Sigma_e %*% R_lambda)
    L_Sigma_spatial <- t(chol(Sigma_spatial_i)) 
    # mat_dist <- covariance_mat(mat_dist, output_mcmc$theta[i,])
    # C <- exp(-mat_dist/output_mcmc$phi[i])
    # 5) Sampling scheme 
    
    Z <- matrix(rnorm(k_l * p), nrow = k_l, ncol = p)
    # E_i <- chol(C) %*% Z %*% t(L_Sigma_spatial)
    E_i <- Z %*% t(L_Sigma_spatial)
    X_pred[, , i] <- M_subject + E_i
  }
  
  return(X_pred)
}

## enforce a minimum cyclic gap, keeping order and total 2*pi
spread_angles <- function(theta, g_min){
  k <- length(theta)
  if (k * g_min >= 2 * pi) stop("min gap too large for k landmarks")
  g <- diff(c(theta, theta[1] + 2 * pi))            # k cyclic gaps (g[k] = wrap-around gap)
  for (it in 1:100) {
    small <- g < g_min
    if (!any(small)) break
    deficit <- sum(g_min - g[small])
    excess  <- g[!small] - g_min
    g[small]  <- g_min
    g[!small] <- g[!small] - deficit * excess / sum(excess)   # take it from the large gaps
  }
  th <- theta[1] + c(0, cumsum(g[-k]))
  s  <- (th[1] + th[k] - 2 * pi) / 2                 # shift so that 0 < th[1] and th[k] < 2*pi
  th - s
}

init_mu_star <- function(X,
                         n_basis = 12, degree = 3,
                         nu      = 4,
                         a_tau   = 2, b_tau = 0.5,
                         start   = c("gpa", "circle"),
                         ridge   = 1e-1,
                         tau0    = c(0.1, 0.1), phi0 = NULL,
                         Sigma_eta = diag(1, 2),
                         min_gap_frac = 0.05,          ## NEW: min gap = frac * (2*pi/k)
                         plot    = TRUE) {
  start <- match.arg(start)
  if (nu <= 3) stop("nu must be > 3 for E[Sigma] to exist")
  k <- dim(X)[1]; n <- dim(X)[3]
  Rot <- function(l) matrix(c(cos(l), -sin(l), sin(l), cos(l)), 2, byrow = TRUE)
  ang <- function(M) atan2(M[, 2], M[, 1]) %% (2 * pi)
  
  ## 1. order landmarks by the GPA angles, rerun GPA, fix a possible cyclic shift
  mean_shape <- shapes::procGPA(X, scale = TRUE)$mshape
  theta_gpa <- atan2(mean_shape[,2], mean_shape[,1])
  theta_gpa <- ifelse(theta_gpa < 0, theta_gpa + 2*pi, theta_gpa)
  orders <- order(theta_gpa)
  X <- X[orders, , ]
  mean_shape <- shapes::procGPA(X, scale = TRUE)$mshape
  theta_gpa  <- atan2(mean_shape[,2], mean_shape[,1])
  theta_gpa <- ifelse(theta_gpa < 0, theta_gpa + 2*pi, theta_gpa)

  ## 1b. minimum gap between consecutive angles
  g_min <- min_gap_frac * 2 * pi / k
  theta <- spread_angles(theta_gpa, g_min)
  moved <- max(abs(theta - theta_gpa))
  
  ## 2. one single hyper (C_bar does not depend on psi); phi bounds from the spread angles
  a_phi <- 0.01                      # now >= g_min / 3
  b_phi <- 5
  build_hyper <- function(psi)
    hyperparameters(a_tau, b_tau, n_basis, degree, nu = nu, psi = psi, n = n,
                    a = 0.1, b = 0.1, Sigma_eta = Sigma_eta,
                    a_phi = a_phi, b_phi = b_phi, X = X)
  hyper <- build_hyper(diag(1, 2))
  
  ## 3. mu*: fit log r at the GPA angles, evaluate the curve at the spread angles
  B_fit <- Basis_Construction(theta_gpa, n_basis, degree)   ## NEW: fit where the radii were observed
  # B     <- Basis_Construction(theta,     n_basis, degree)   ##      evaluate at the angles used by mu*
  if (start == "gpa") {
    log_r <- log(sqrt(rowSums(mean_shape^2)))
    log_r <- log_r - mean(log_r)
    beta0 <- solve(crossprod(B_fit) + ridge * diag(ncol(B_fit)), crossprod(B_fit, log_r))
    gamma <- solve(crossprod(hyper$C_bar), crossprod(hyper$C_bar, beta0))
  } else {
    gamma <- rep(0, ncol(hyper$C_bar))
  }
  beta <- as.vector(hyper$C_bar %*% gamma)
  r    <- as.numeric(exp(B_fit %*% beta))
  mu   <- cbind(r * cos(theta_gpa), r * sin(theta_gpa))
  mu_c <- scale(mu, scale = FALSE)
  cm   <- attr(mu_c, "scaled:center")
  hyper$b_phi <- max(dist(mu_c))/3
  hyper$a_phi <- min(dist(mu_c)[dist(mu_c) != 0])/3
  phi0 <- (hyper$a_phi + hyper$b_phi) / 2
  ## 4. similarity per subject
  alpha <- lambda <- numeric(n); eta <- matrix(NA, n, 2); refl <- logical(n)
  for (i in 1:n) {
    Xc <- scale(X[, , i], scale = FALSE)
    s  <- svd(crossprod(Xc, mu_c))
    Ri <- s$v %*% t(s$u)
    if (det(Ri) < 0) { refl[i] <- TRUE; s$u[, 2] <- -s$u[, 2]; Ri <- s$v %*% t(s$u) }
    lambda[i] <- atan2(Ri[2, 1], Ri[1, 1]) %% (2 * pi)
    alpha[i]  <- sum(Xc * (mu_c %*% Ri)) / sum(mu_c^2)
    eta[i, ]  <- colMeans(X[, , i]) - alpha[i] * as.vector(cm %*% Ri)
  }
  if (any(refl))
    warning(sum(refl), " subject(s) are mirror images of the mean: ",
            paste(which(refl), collapse = ", "))
  
  ## 5. S_hat in the mu* frame
  S_hat <- matrix(0, 2, 2)
  for (i in 1:n) {
    Ri <- Rot(lambda[i])
    E  <- X[, , i] - (alpha[i] * (mu %*% Ri) + matrix(eta[i, ], k, 2, byrow = TRUE))
    S_hat <- S_hat + crossprod((E %*% t(Ri)) / alpha[i])
  }
                     # E[Sigma] = Sigma_target
  S_hat <- S_hat / (n * k)
  #hyper <- build_hyper((100 - 3) * S_hat)
  hyper$min_gap <- g_min                              ## NEW: used by the theta sampler
  #hyper$psi <- S_hat
  init <- init_param(tau0, lambda, theta_gpa, beta, S_hat, eta, alpha, gamma,
                     phi = phi0, n = n)
  

  rms <- sapply(1:n, function(i)
    sqrt(mean((X[, , i] - (alpha[i] * (mu %*% Rot(lambda[i])) +
                             matrix(eta[i, ], k, 2, byrow = TRUE)))^2)))
  
  ## 7. diagnostic plot
  if (plot) {
    par(mfrow = c(1,1))
    th_grid <- seq(0, 2 * pi, length.out = 400)
    r_grid  <- as.numeric(exp(Basis_Construction(th_grid, n_basis, degree) %*% beta))
    curve   <- cbind(r_grid * cos(th_grid), r_grid * sin(th_grid))
    aligned <- lapply(1:n, function(i)
      ((X[, , i] - matrix(eta[i, ], k, 2, byrow = TRUE)) %*% t(Rot(lambda[i]))) / alpha[i])
    rng <- range(unlist(aligned), curve)
    plot(NA, xlim = rng, ylim = rng, asp = 1, xlab = "", ylab = "",
         main = paste0("Initial mu* (", start, ")"))
    for (A in aligned) points(A, col = adjustcolor("grey40", 0.3), pch = 16, cex = 0.6)
    lines(curve, col = "red", lwd = 2)
    points(mu, col = "red", pch = 19)
    text(mu, labels = 1:k, pos = 3, cex = 0.7, col = "red")
  }
  
  list(init = init, hyper = hyper, X = X, theta = theta_gpa, theta_gpa = theta_gpa,
       mu = mu, S_hat = S_hat, rms = rms, reflected = which(refl))
}
gg_mcmc_diagnostics <- function(data, param_name = "NA", real_values = "NA" ,matrix_param = FALSE) {
  # This function is very useful to build fast plots for our model
  
  if (is.null(dim(data))) {
    data <- matrix(data, ncol = 1)
    colnames(data) <- param_name
  } else if (is.null(colnames(data))) {
    if (matrix_param == TRUE) {
      colnames(data) <- paste(param_name, 1:ncol(data))
    }
    
    else{
      
      
    }
  }
  
  plot_list <- list()
  n_iter <- nrow(data)
  i <- 1
  for (col in colnames(data)) {
    chain <- data[, col]
    
    ess_val <- round(LaplacesDemon::ESS(chain), 1)
    quantiles <- quantile(chain, probs = c(0.025, 0.975))
    mean_val <- mean(chain)
    
    df <- data.frame(
      Iteration = 1:n_iter,
      Value = chain
    )
    p <- ggplot(df, aes(x = Iteration, y = Value)) +
      annotate("rect", xmin = -Inf, xmax = Inf, 
               ymin = quantiles[1], ymax = quantiles[2], 
               fill = "#3182bd", alpha = 0.15) +
      geom_line(color = "gray25", linewidth = 0.4) +
      geom_hline(yintercept = mean_val, color = "#e41a1c", 
                 linetype = "dashed", linewidth = 0.8) +
      geom_hline(yintercept = quantiles[1], color = "#3182bd", linetype = "dotted") +
      geom_hline(yintercept = quantiles[2], color = "#3182bd", linetype = "dotted") +
      annotate("label", x = Inf, y = Inf, 
               label = paste("ESS:", ess_val), 
               hjust = 1.1, vjust = 1.1, 
               fill = "white", alpha = 0.85, 
               fontface = "bold", size = 3.5, color = "gray10") +
      labs(
        title = paste("Traceplot of", col),
        subtitle = paste0("95% Credibility Interval: [", 
                          round(quantiles[1], 4), ", ", 
                          round(quantiles[2], 4), "]"),
        x = "Iteration",
        y = "Value"
      ) +
      theme_minimal(base_size = 11) +
      theme(
        plot.title = element_text(face = "bold", color = "gray10"),
        plot.subtitle = element_text(color = "gray40", size = 9),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "gray92")
      )
    if (all(!is.na(real_values)) ) {
      p <- p + geom_hline(yintercept = real_values[i], color = "cyan", 
                          linetype = "dashed", linewidth = 0.8)
    }
    plot_list[[col]] <- p
    i <- i + 1
    }
  if (length(plot_list) == 1) {
    return(plot_list[[1]])
  } else {
    return(plot_list)
  }
}
show <- function(p) {
  if (inherits(p, "ggplot")) print(p) else for (q in p) print(q)
}



visualize_results <- function(output, hypers, X, file_name, sam = NULL){
  options(error = NULL)
  if (!is.null(dev.list())) { # controllino anti-bestemmie
    dev.off()
  }
  pdf(paste(file_name, sep = "_predictions", ".pdf"), width = 20, height = 5.5) # open pdf
  par(mfrow = c(1, 4)) # define the margins for the first window
  
  n_iter <- dim(output$mean)[1] # number of iterations
  k_l <- dim(output$mean)[2] # number of landmarks
  n <- dim(X)[3] # sample size
  mean_samp <- output$mean # common star shaped mean  [iter, k, 2]
  mean_low <- apply(mean_samp, c(2,3), quantile, probs = 0.025) # compute quantiles
  mean_high <- apply(mean_samp, c(2,3), quantile, probs = 0.975)
  mean_samp_mean <- apply(mean_samp, c(2,3), mean) # compute average value
  for (i in 1:n) {
    X_i <- X[,,i] # unit i # 
    
    plot(mean_samp_mean, asp = 1, pch = 21, bg = "darkblue",
         main = "Estimated Mean + CI", xlim = range(mean_samp_mean[,1]) + c(-0.1, 0.1),
         ylim = range(mean_samp_mean[,2]) + c(-0.1, 0.1))
    
    lines(c(mean_samp_mean[,1], mean_samp_mean[1,1]),
          c(mean_samp_mean[,2], mean_samp_mean[1,2]),
          col = "blue", lwd = 2)
    
    for (k in 1:k_l) {
      
      polygon(
        x = c(mean_low[k,1], mean_high[k,1],
              mean_high[k,1], mean_low[k,1]),
        y = c(mean_low[k,2], mean_low[k,2],
              mean_high[k,2], mean_high[k,2]),
        border = NA,
        col = "red"
      )
    }
    
    if (!is.null(sam)) {
      points(sam$m, col = "yellow")
      
    }
    # 
    plot(X_i, asp = 1, pch = 21, bg = "darkred",
         main = "i-th unit", xlim = range(X_i[,1]) + c(-0.01, 0.01),
         ylim = range(X_i[,2]) + c(-0.01, 0.01) )
    
    lines(c(X_i[,1], X_i[1,1]),
          c(X_i[,2], X_i[1,2]),
          col = "red", lwd = 2)
    
    app <- extract_mean_i(output, i)
    samp_i <- app$mean_i # mean of unit i
    mean_samp_i <- apply(samp_i, c(1,2), mean)
    low_i <- apply(samp_i, c(1,2), quantile, 0.025)
    high_i <- apply(samp_i, c(1,2), quantile, 0.975)
    plot(mean_samp_i, asp = 1, pch = 21, bg = "black",
         main = "Estimated config + CI",  xlim = range(mean_samp_i[,1]) + c(-0.01, 0.01),
         ylim = range(mean_samp_i[,2]) + c(-0.01, 0.01))
    # 
    lines(c(mean_samp_i[,1], mean_samp_i[1,1]),
          c(mean_samp_i[,2], mean_samp_i[1,2]),
          col = "brown", lwd = 2)
    for (k in 1:k_l) {
      
      polygon(
        x = c(low_i[k,1], high_i[k,1],
              high_i[k,1], low_i[k,1]),
        y = c(low_i[k,2], low_i[k,2],
              high_i[k,2], high_i[k,2]),
        col = rgb(0,0,1,0.40),
        border = NA
      )
      
    }
    posterior_pred <- app$pred_i # posterior predictive distribution
    post_samp_i <- apply(posterior_pred, c(1,2), mean)
    low_post_i <- apply(posterior_pred, c(1,2), quantile, 0.025)
    high_post_i <- apply(posterior_pred, c(1,2), quantile, 0.975)
    plot(post_samp_i, asp = 1, pch = 21, bg = "black",
         main = "Posterior Predictive", xlim = range(post_samp_i[,1]) + c(-0.01, 0.01),
         ylim = range(post_samp_i[,2]) + c(-0.01, 0.01))
    lines(c(post_samp_i[,1], post_samp_i[1,1]),
          c(post_samp_i[,2], post_samp_i[1,2]),
          col = "brown", lwd = 2)
    for (k in 1:k_l) {
      
      polygon(
        x = c(low_post_i[k,1], high_post_i[k,1],
              high_post_i[k,1], low_post_i[k,1]),
        y = c(low_post_i[k,2], low_post_i[k,2],
              high_post_i[k,2], high_post_i[k,2]),
        col = rgb(0,0,1,0.40),
        border = NA
      )
    }
    points(X_i, pch = 16, col = "red")
    
  }
  dev.off()
  pdf(paste(file_name, sep = "_diagnostics", ".pdf"), width = 20, height = 5.5) # open pdf
  if (is.null(sam)) {
    show(gg_mcmc_diagnostics(output$tau[,1], param_name = "tau_1", real_values = NA))
    show(gg_mcmc_diagnostics(output$tau[,2], param_name = "tau_2", real_values = NA))
    show(gg_mcmc_diagnostics(output$alphas/output$alphas[,1], param_name = "alpha ratios others vs 1", real_values = NA, TRUE))
    show(gg_mcmc_diagnostics((output$lambdas - output$lambdas[,1])%%(2*pi), param_name = "lambda diff others vs 1", real_values = NA, TRUE))
    show(gg_mcmc_diagnostics(output$sigma_2[,], param_name = "sigma_2", real_values = NA))
    show(gg_mcmc_diagnostics(output$sigma_2_epurated[,], param_name = "sigma_2 epur", real_values = NA))
  }
  
  else{
    show(gg_mcmc_diagnostics(output$tau[,1], param_name = "tau_1", real_values = sam$tau[1]))
    show(gg_mcmc_diagnostics(output$tau[,2], param_name = "tau_2", real_values = sam$tau[2]))
    show(gg_mcmc_diagnostics(output$alphas/output$alphas[,1], param_name = "alpha ratios others vs 1", real_values = sam$alphas/sam$alphas[1], TRUE))
    show(gg_mcmc_diagnostics((output$lambdas - output$lambdas[,1])%%(2*pi), param_name = "lambda diff others vs 1", real_values = (sam$lambda - sam$lambda[1]) %%(2*pi), TRUE))
    V <- median(diag(sam$S_mat %*% sam$D_p %*% solve(t(sam$D_p)%*%sam$S_mat %*% sam$D_p) %*% t(sam$D_p) %*% sam$S_mat))
    show(gg_mcmc_diagnostics(output$sigma_2[,], param_name = "sigma_2", real_values = sam$sigma_2))
    show(gg_mcmc_diagnostics(output$sigma_2_epurated[,], param_name = "sigma_2_epur", real_values = sam$sigma_2/V))
    # show(gg_mcmc_diagnostics(output$eta[,,1] - output$eta[,1,1] , param_name = "eta diff others vs 1", real_values = sam$eta[,1] - sam$eta[1,1], TRUE))
    # show(gg_mcmc_diagnostics(output$eta[,,2] - output$eta[,1,2] , param_name = "eta diff others vs 1", real_values = sam$eta[,2] - sam$eta[1,2], TRUE))
  }
  dev.off()
  
  
  X0 <- X[, , 2]
  pair <- c(2, 14)
  distt <- 0
  for (i in 1:nrow(X0)) {
    for (j in 1:nrow(X0)) {
      
      dis <- sum((X0[i,] - X0[j,])^2)
      
      if (dis > distt) {
        pair <- c(i, j)
        distt <- dis
      }
    }
  }
  # pair # those are the two landmarks w.r.t we are going to identify stuff
  nsamp <- n_iter
  k <- k_l
  mean_orig <- array(NA, c(k, 2, nsamp))
  mean_noloc <- array(NA, c(k, 2, nsamp))
  mean_noloc_i <- array(NA, c(k, 2, nsamp))
  mean_norot <- array(NA, c(k, 2, nsamp))
  mean_norot_i <- array(NA, c(k, 2, nsamp))
  mean_nosize <- array(NA, c(k, 2, nsamp))
  mean_nosize_i <- array(NA, c(k, 2, nsamp))
  sigma_orig <- array(NA, c(2, 2, nsamp))
  sigma_noloc <- array(NA, c(2, 2, nsamp))
  sigma_norot <- array(NA, c(2, 2, nsamp))
  sigma_nosize <- array(NA, c(2, 2, nsamp))
  warp_mean <- array(NA, c(k, 2, n))
  pdf(paste(file_name, sep = "_mean_ident", ".pdf"), width = 15, height = 8.5) # open pdf
  
  for (iobs in 1:n) {
    
    # 1. Feedback reale sulla console R
    cat(sprintf("\rElaborazione Unità: %3d", iobs))
    flush.console()
    
    krot <- pair[2]
    for (isamp in 1:nsamp){
      ## Creo oggetti originali
      rrr <- Rmat(output$lambdas[isamp, iobs])
      sigma_orig[, , isamp] <- t(rrr) %*% output$Sigma[isamp, , ] %*% rrr * output$alphas[isamp, iobs]^2
      mean_orig[, , isamp] <- output$mean[isamp, ,] # mean
      S <- S_mat_constructor(mean_orig[, , isamp])
      warping <- S%*%output$W[[isamp]][,,iobs]
      # Metto una coordinata a zero (Rimuovo traslazione)
      mean_noloc[, , isamp] <- mean_orig[, , isamp] - matrix(mean_orig[pair[1], , isamp], nrow = k, ncol = 2, byrow = T)
      warping <- warping - matrix(mean_orig[pair[1], , isamp], nrow = k, ncol = 2, byrow = T)
      mean_noloc_i[, , isamp] <- mean_noloc[, , isamp] + warping
      sigma_noloc[, , isamp] <- sigma_orig[, , isamp]
      
      # Tolgo la rotazione
      angle <- atan2(mean_noloc[krot, 2, isamp], mean_noloc[krot, 1, isamp])
      mean_norot[, , isamp] <- mean_noloc[, , isamp] %*% Rmat(angle)
      warping <- warping %*% Rmat(angle)
      mean_norot_i[, , isamp] <- mean_norot[, , isamp] +  warping
      sigma_norot[, , isamp] <- t(Rmat(angle)) %*% sigma_noloc[, , isamp] %*% Rmat(angle)
      
      # Stabilizzo elevando la media campionaria
      #mean_norot[, , isamp] <- mean_norot[, , isamp] - matrix(colMeans(mean_norot[, , isamp]), nrow = k, ncol = 2, byrow = T)
      #mean_norot_i[, , isamp] <- mean_norot_i[, , isamp] - matrix(colMeans(mean_norot_i[, , isamp]), nrow = k, ncol = 2, byrow = T)
      # normx <- norm(mean_norot[, , isamp])
      normx <- norm(sweep(mean_orig[, , isamp], 2, colMeans(mean_orig[, , isamp])), "F")
      
      mean_nosize[, , isamp] <- mean_norot[, , isamp] / normx
      warping <- warping/normx
      mean_nosize_i[, , isamp] <- mean_nosize[, , isamp] + warping
      sigma_nosize[, , isamp] <- sigma_norot[, , isamp] / normx^2
    }
    
    # Calcoliamo la forma media a posteriori per disegnare una linea guida pulita
    mean_shape_posterior <- apply(mean_nosize, c(1, 2), mean)
    mean_shape_posterior_i <- apply(mean_nosize_i, c(1, 2), mean)
    warp_mean[, , iobs] <- mean_shape_posterior_i
    # ====================================================================
    # PAGINA 1: Il "Cloud" dei Landmark allineati + Forma Media
    # ====================================================================
    par(mfrow = c(1, 1), mar = c(5, 5, 4, 2))
    
    plot(0, 0, type = "n", 
         xlim = range(mean_nosize[, 1, ]), 
         ylim = range(mean_nosize[, 2, ]),
         asp = 1, # FONDAMENTALE: mantiene le proporzioni geometriche reali
         xlab = "X", ylab = "Y",
         main = paste("Landmark Allineati Unità", iobs))
    
    # Disegniamo i punti di ogni iterazione (usando una trasparenza per vedere la densità)
    colors_palette <- rainbow(k)
    for (iii in 1:k) {
      points(t(mean_nosize[iii, , ]), col = adjustcolor(colors_palette[iii], alpha.f = 0.4), cex = 1, pch = 20)
    }
    
    # Disegniamo la sagoma della forma media a posteriori che unisce i punti
    lines(c(mean_shape_posterior[,1], mean_shape_posterior[1,1]), 
          c(mean_shape_posterior[,2], mean_shape_posterior[1,2]), 
          col = "black", lwd = 2, lty = 1)
    
    
    # Evidenziamo i landmark medi stimati
    points(mean_shape_posterior, col = "black", bg = colors_palette, pch = 21, cex = 1.4, lwd = 1.5)
    # ====================================================================
    # PAGINE 2+: I Traceplot dei Landmark (X e Y)
    # ====================================================================
    par(mfrow = c(3, 2), mar = c(4, 4, 2.5, 1)) # Ripristiniamo la griglia ordinata 3x2
    
    for (iii in 1:k) {
      for (ip in 1:2) {
        coord_label <- ifelse(ip == 1, "Coord X", "Coord Y")
        
        plot(mean_nosize_i[iii, ip, ], type = "l", col = "steelblue",
             xlab = "Iterazioni MCMC", ylab = coord_label,
             main = paste("Media", iii, "-", coord_label))
        print(ESS(mean_nosize_i[iii, ip, ]))
        
        # Linea rossa tratteggiata sul valore medio stimato per controllare la convergenza
        abline(h = mean_shape_posterior_i[iii, ip], col = "red", lty = 2, lwd = 1.5)
      }
    }
    
    par(mfrow = c(2, 2), mar = c(4, 4, 2.5, 1))
    
    sigma_labels <- matrix(c("Var(X)", "Cov(X,Y)", "Cov(Y,X)", "Var(Y)"), 2, 2)
    for (i in 1:2) {
      for (j in 1:2) {
        plot(sigma_nosize[i, j, ], type = "l", col = "darkorange",
             xlab = "Iterazioni MCMC", ylab = "Valore",
             main = paste("Sigma -", sigma_labels[i, j]))
        
        # Linea di riferimento sul valore medio di Sigma
        abline(h = mean(sigma_nosize[i, j, ]), col = "red", lty = 2, lwd = 1.5)
      }
    }
  }
  
  dev.off()
  
  pdf(paste(file_name, sep = "_funny", ".pdf"), width = 10, height = 8.5) # open pdf
  
  amp  <- 1       
  cols <- rainbow(k)
  
  for (iobs in 1:n) {
    m0 <- mean_shape_posterior       # mu*
    m1 <- warp_mean[,,iobs]     # mu* + S W_i
    rng  <- range(c(mean_nosize[,1,], mean_nosize_i[,1,]))
    rngy <- range(c(mean_nosize[,2,], mean_nosize_i[,2,]))
    
    par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
    
    # (a) forme medie posteriori e spostamento
    plot(0, 0, type = "n", xlim = rng, ylim = rngy, asp = 1, xlab = "X", ylab = "Y",
         main = paste0("Unità ", iobs, ": mu* (tratteggio) vs mu* + S W (rosso)"))
    lines(rbind(m0, m0[1, ]), lty = 2, lwd = 2)
    lines(rbind(m1, m1[1, ]), col = "red", lwd = 2)
    arrows(m0[,1], m0[,2], m0[,1] + amp*(m1[,1]-m0[,1]), m0[,2] + amp*(m1[,2]-m0[,2]),
           length = 0.06, col = "blue", lwd = 1.5)
    points(m0, pch = 21, bg = "white", cex = 1.2)
    points(m1, pch = 19, col = "red", cex = 0.9)
    text(m0, labels = 1:k, pos = 3, cex = 0.7)
    X_i <- X[,,iobs]
    plot(X_i, asp = 1, pch = 21, bg = "darkred",
         main = "i-th unit", xlim = range(X_i[,1]) + c(-0.001, 0.001),
         ylim = range(X_i[,2]) + c(-0.001, 0.001) )
    
    lines(c(X_i[,1], X_i[1,1]),
          c(X_i[,2], X_i[1,2]),
          col = "red", lwd = 2)
  }
  dev.off() 
  
  pdf(paste(file_name, sep = "_continous_shape", ".pdf"), width = 8, height = 5.5)
  par(mfrow = c(1,1))
  th_grid <- seq(0, 2 * pi, length.out = 400)
  r_grid <- matrix(NA, nrow = nsamp, ncol = 400)
  curved <- array(NA, dim = c(400, 2, nsamp))
  for (s in 1:nsamp) {
    r_grid[s, ] <- as.numeric(exp(Basis_Construction(th_grid, hypers$n_basis, hypers$degree) %*% output$betas[s,]))
    curved[,,s] <- cbind(r_grid[s, ] * cos(th_grid), r_grid[s, ] * sin(th_grid))
  }
  
  mean_low <- apply(curved, c(1,2), quantile, probs = 0.025) # compute quantiles
  mean_high <- apply(curved, c(1,2), quantile, probs = 0.975)
  mean_samp_mean <- apply(curved, c(1,2), mean) # compute average value
  
  
  plot(mean_samp_mean, asp = 1, pch = 21, bg = "black",
       main = "Continous shape", xlim = range(mean_samp_mean[,1]) + c(-0.001, 0.001),
       ylim = range(mean_samp_mean[,2]) + c(-0.001, 0.001))
  lines(c(mean_samp_mean[,1], mean_samp_mean[1,1]),
        c(mean_samp_mean[,2], mean_samp_mean[1,2]),
        col = "brown", lwd = 2)
  for (k in 1:400) {
    
    polygon(
      x = c(mean_low[k,1], mean_high[k,1],
            mean_high[k,1], mean_low[k,1]),
      y = c(mean_low[k,2], mean_low[k,2],
            mean_high[k,2], mean_high[k,2]),
      col = rgb(0,0,1,0.40),
      border = NA
    )
  }
  dev.off()
  
}

plot_status <- function(init_env, hyper_env, X, file, it, i_unit = 1) {
  # tmp <- paste0(file, ".tmp")
  closed <- function(M) rbind(M, M[1, ])
  # pdf(tmp, width = 13, height = 4.8)
  # ok <- FALSE
  # on.exit({ dev.off(); if (ok) file.rename(tmp, file) }, add = TRUE)
  par(mfrow = c(1, 3), mar = c(3, 3, 3.5, 1))
  
  # 1. mu* e mu* + S W
  mu_star <- init_env$mean
  mu_warp <- mu_star + hyper_env$S_mat %*% init_env$W[,,i_unit]
  plot(closed(mu_star), type = "o", pch = 19, col = "darkblue", asp = 1,
       xlim = range(mu_star[, 1], mu_warp[, 1]), ylim = range(mu_star[, 2], mu_warp[, 2]),
       xlab = "", ylab = "",
       main = sprintf("mu*  | iter %d\nsigma2 = %.3g  tau = (%.3g, %.3g)",
                      it, init_env$sigma_2, init_env$tau[1], init_env$tau[2]))
  lines(closed(mu_warp), lty = 2, lwd = 2, col = "red")
  text(mu_star, labels = seq_len(nrow(mu_star)), pos = 3, cex = 0.6)
  legend("topright", c("mu*", "mu* + SW"), col = c("darkblue", "red"),
         lty = c(1, 2), bty = "n", cex = 0.8)
  
  # 2. mu_i
  mu_i <- init_env$mean_i[, , i_unit]
  plot(closed(mu_i), type = "o", pch = 19, col = "brown", asp = 1,
       xlab = "", ylab = "", main = sprintf("mu_%d", i_unit))
  
  # 3. X_i vs mu_i
  Xi <- X[, , i_unit]
  plot(closed(Xi), type = "o", pch = 21, bg = "black", asp = 1,
       xlim = range(Xi[, 1], mu_i[, 1]), ylim = range(Xi[, 2], mu_i[, 2]),
       xlab = "", ylab = "", main = sprintf("X_%d (nero) vs mu_%d (rosso)", i_unit, i_unit))
  lines(closed(mu_i), col = "red", lwd = 2)
  
  # ok <- TRUE
}

##### Rccp for some functions #####
Rcpp::cppFunction(depends = "RcppArmadillo", code = '
Rcpp::List units_loglik(const arma::mat& mean, const arma::mat& Dp,
                        const arma::mat& SDp, const arma::cube& psi,
                        const arma::cube& R, const arma::cube& eta,
                        const arma::cube& X, const arma::cube& Q,
                        const arma::vec& alphas, bool keep = false) {
  const arma::uword k = X.n_rows, n = X.n_slices;
  double val = 0.0, trace = 0.0;
  arma::cube mean_i, resid, SW, W;
  if (keep) {
    mean_i.set_size(k, 2, n); resid.set_size(2, 2, n);
    SW.set_size(k, 2, n);     W.set_size(k, 2, n);
  }
  for (arma::uword i = 0; i < n; ++i) {
    arma::mat wi  = Dp  * psi.slice(i);           // W_i
    arma::mat swi = SDp * psi.slice(i);           // S W_i
    arma::mat M   = alphas(i) * (mean + swi) * R.slice(i) + eta.slice(i);
    arma::mat Y   = X.slice(i) - M;
    arma::mat r   = Y.t() * Y;
    val   += arma::accu(Q.slice(i) % r);          // sum(Q * residual)
    trace += arma::accu(wi % swi);                // tr(W_i^T S W_i)
    if (keep) { W.slice(i) = wi; SW.slice(i) = swi;
                mean_i.slice(i) = M; resid.slice(i) = r; }
  }
  if (!keep) return Rcpp::List::create(Rcpp::Named("val") = val,
                                       Rcpp::Named("trace") = trace);
  return Rcpp::List::create(Rcpp::Named("val") = val, Rcpp::Named("trace") = trace,
                            Rcpp::Named("W") = W, Rcpp::Named("SW") = SW,
                            Rcpp::Named("mean_i") = mean_i,
                            Rcpp::Named("residual") = resid);
}')
