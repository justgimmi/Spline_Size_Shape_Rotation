source(file = file.path(getwd(), "src/Packages.R"))

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
  Basis <- splineDesign(total_knots, thetas, ord = degree + 1, outer.ok = FALSE,
                        derivs = F) # Basis construction
  
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



in_model_sample <- function(n = 1, K_l,  thetas = NA, n_int_knots, degree = 3, tau = c(0.1, 0.1),phi = 0.3, sigma_2 = 0.01,
                            beta_values = NA, lambda = NA, alphas = NA, eta = NA, 
                            Sigma_e = NA, indep_error = TRUE){
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
  mat_dist <- diag(0, nrow = length(thetas), ncol = length(thetas)) 
  mat_dist <- covariance_mat(mat_dist, thetas)
  C <- exp(-mat_dist/phi) 
  
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
  #R1 <- Basis_Construction(0, n_int_knots, degree, derivatives = 1) - Basis_Construction(2*pi, n_int_knots, degree, derivatives = 1) # (B'(0) - B'(2pi))
  #R2 <- Basis_Construction(0, n_int_knots, degree, derivatives = 2) - Basis_Construction(2*pi, n_int_knots, degree, derivatives = 2)
  #A <- rbind(Eigen_vector_null, R, R1, R2)
  A <- rbind(Eigen_vector_null, R)
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
  k_free <- K_l - 5
  S_mat <- matrix(0, nrow = K_l, ncol = K_l)
  S_mat <- S_mat_constructor(S_mat, mu_mean)
  # init <- list()
  # hyper <- list()
  # init$mean <- mu_mean
  # init$betas <- beta_values
  # init$r <- r
  # init$thetas <- thetas
  # hyper$n_basis <- n_int_knots
  # hyper$degree <- degree
  
  
  D_p <- D_p_raw(mu_mean, thetas, beta_values, r, n_int_knots, degree)
  if (sigma_2 == 0 && indep_error == TRUE) {
    Precision <- diag(1, nrow = k_free)
    U_prec <- chol(Precision)
    k_free <- ncol(D_p)  # k - 5
    psi <- matrix(0, nrow = k_free, ncol = 2)
    W <- D_p %*% psi
  }
  
  else if(sigma_2 != 0 && indep_error == TRUE){
    Precision <- diag(1/sigma_2, nrow = k_free)
    U_prec <- chol(Precision)
    k_free <- ncol(D_p)  # k - 5
    
    Z <- matrix(rnorm(k_free * 2), nrow = k_free, ncol = 2)
    psi <- backsolve(U_prec, Z)
    W <- D_p %*% psi
  }
  else{
    Precision <- (1 / sigma_2) * t(D_p) %*% S_mat %*% D_p
    U_prec <- chol(Precision)
    k_free <- ncol(D_p)  # k - 3
    
    Z <- matrix(rnorm(k_free * 2), nrow = k_free, ncol = 2)
    psi <- backsolve(U_prec, Z)
    W <- D_p %*% psi
  }
  U_prec <- chol(Precision)
  k_free <- ncol(D_p)  # k - 3

  # Z <- matrix(rnorm(k_free * 2), nrow = k_free, ncol = 2)
  # psi <- backsolve(U_prec, Z)
  # W <- D_p %*% psi
  R <- array(NA, dim = c(2, 2, n))
  mu_i <- array(NA, dim = c(K_l, 2, n))
  if (length(lambda) == 1) {
    lambda <- runif(n = n)*2*pi # generate rotation angles
  }
  
  if (length(alphas) == 1) {
    alphas <- rgamma(n = n, shape = 3, rate = 0.5) # sample the size effect
    
  }
  
  if (length(eta) == 1) {
    eta <- matrix(rnorm(n = n*2, sd = 4.5), ncol = 2) # sample the translation effect 
    
  }
  
  for (i in 1:n) { # compute rotation matrix 
    R[,,i] <- matrix(c(cos(lambda[i]), -sin(lambda[i]), sin(lambda[i]), cos(lambda[i])), 
                     byrow = T, nrow = 2, ncol = 2)
    eta_matrix <- matrix(eta[i,], nrow = K_l, ncol = 2, byrow = T)
    mu_i[,,i] <- alphas[i]*(mu_mean + S_mat %*% W) %*%R[,,i] + eta_matrix # compute the mean configuration for every unit
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
    X[,,i] <- t(chol(C))%*%matrix(rnorm(n = K_l * 2), nrow = K_l, ncol = 2)%*%chol(Sigmas) + mu_i[,,i]
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
  final_param$phi <- phi
  final_param$psi <- psi
  final_param$S_mat <- S_mat
  final_param$W <- W
  final_param$D_p <- D_p
  
  return(final_param)

}


# MCMC Util functions ----
Rmat <- function(angle) {
  matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)), nrow = 2, ncol = 2)
}


S_mat_constructor <- function(mat, mu){
  # mat: K \times K conditional positive definite matrix
  # mu: value of the latent mean K \times K
  k_l <- nrow(mu) # number of landmarks
  for (i in 1:k_l) {
    for (j in 1:k_l) {
      normm <- sum((mu[i, ] - mu[j, ])^2)   # = r^2
      mat[i, j] <- ifelse(normm > 0, normm * log(sqrt(normm)), 0)  # r^2 log(r)
    }
  }
  return(mat)
}

# 
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
D_p_raw <- function(mu, thetas, betas, r, n_basis, degree) {
  k_l <- nrow(mu)
  B_deriv <- Basis_Construction(thetas, n_basis, degree, derivatives = 1)
  deriv_x <- r * as.numeric(B_deriv %*% betas) * cos(thetas) - mu[,2]
  deriv_y <- r * as.numeric(B_deriv %*% betas) * sin(thetas) + mu[,1]
  B <- rbind(matrix(1, 1, k_l), t(mu), t(deriv_x), t(deriv_y))
  sv <- svd(B, nu = 0, nv = k_l)
  sv$v[, (nrow(B) + 1):k_l, drop = FALSE]
}

D_p_constructor <- function(mu, thetas, betas, r, n_basis, degree, D_p_prev = NULL) {
  D_raw <- D_p_raw(mu, thetas, betas, r, n_basis, degree)
  # if (is.null(D_p_prev)) return(D_raw)
  # M  <- crossprod(D_p_prev, D_raw)
  # sm <- svd(M)
  # Rot <- sm$v %*% t(sm$u)
  # D_raw %*% Rot
}
covariance_mat <- function(mat, thetas){
  # thetas: vector of angles
  # mat: K x K covariance matrix
  # phi: scale parameter of the exponential covariance function
  k_l <- length(thetas)
  for (i in 1:k_l) {
    for (j in 1:k_l) {
      #d <- min(2*pi - abs(thetas[i] - thetas[j]), abs(thetas[i] - thetas[j]))
      mat[i, j] <- min(2*pi - abs(thetas[i] - thetas[j]), abs(thetas[i] - thetas[j]))
      mat[j, i] <- mat[i, j]
    }
  }
  return(mat)
}

init_param <- function(tau, lambdas, thetas, betas, Sigma, eta, alphas, gammas, phi, n){
  # n --> number of units 
  param <- list()
  param$n <- n
  param$k_l <- length(thetas)
  param$tau <- tau
  param$lambdas <- lambdas
  param$thetas <- thetas
  param$gammas <- gammas
  param$gaps <- c(thetas[1], diff(thetas), 2*pi- thetas[param$k_l])/(2*pi)
  #param$omega_theta <- log(param$gaps)
  param$omega_theta <- log(param$gaps) - log(param$gaps[param$k_l + 1])
  param$thetas_sorted <- sort(thetas)
  param$betas <- as.matrix(betas)
  param$Sigma <- Sigma 
  param$eta <- eta
  param$alphas <- alphas
  param$r <- numeric(length(thetas)) # deterministic
  param$mean <- matrix(0, nrow = length(thetas), ncol = 2) #deterministic
  param$mean_i <- array(0, dim = c(length(thetas), 2, n))  # deterministic 
  param$R <- array(NA, dim = c(2, 2, n)) #deterministic
  
  param$Sigma_inv <- solve(param$Sigma) # deterministic
  param$Q_R <- array(NA, dim = c(2, 2, n)) # deterministic 
  param$residual <-  array(NA, dim = c(2, 2, n)) #deterministic 
  param$mat_dist <- diag(0, nrow = length(thetas), ncol = length(thetas)) 
  param$mat_dist <- covariance_mat(param$mat_dist, thetas)
  param$C <- exp(-param$mat_dist/phi) # landmark covariance matrix
  param$chol_c <- chol((param$C)) # cholesky decompoition of the covariance matrix -->
  # R is dumb so it is upper-triangular
  param$phi <- phi
  
  #param$psi <- matrix(rnorm(2*length(thetas)), nrow = length(thetas), ncol = 2)
  param$sigma_2 <- 1e-4
  param$z <- 0
  param$chi <- 0.1
  return(param)
  
}

hyperparameters <- function(a_tau, b_tau, n_basis, degree,
                            width_theta = pi/4, m = 8, nu, psi, n, a, b, Sigma_eta = diag(1000, nrow = 2, ncol = 2),
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
  hyper$Sigma_tau <- diag(1e-4, nrow = 2)
  hyper$a_opts <- 0.234
  hyper$tau_acc <- 0
  
  # Theta block -----
  hyper$width_theta <- width_theta
  hyper$m <- m # it is used to have a slice of size w*m
  # Lambda block ----
  hyper$width_lambda <- width_theta
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
  hyper$Sigma_a <- matrix(1e-4, nrow = n)
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
  hyper$Sigma_lambda <- matrix(1e-4, nrow = n)
  hyper$gamma_lambda <- 1
  hyper$lambda_opt <- 0.44
  hyper$lambda_acc <- matrix(0, nrow = n)
  
  # phi block
  hyper$lambda_phi <- 2.38^2
  hyper$mu_phi <- 0
  hyper$Sigma_phi <- 1e-4
  #hyper$L_a <- chol(hyper$lambda_a*hyper$Sigma_a)
  hyper$alpha_phi <- 0
  hyper$gamma_phi <- 1
  hyper$a_opt <- 0.44
  hyper$a_phi <- a_phi
  hyper$b_phi <- b_phi
  hyper$phi_acc <- 0
  
  
  hyper$S_mat <- matrix(0, nrow = k_l, ncol = k_l)
  hyper$a_sigma <- 1e-4
  hyper$lambda_sigma_cp <- 2.38^2
  hyper$Sigma_sigma_cp  <- 1e-8
  hyper$mu_sigma_cp     <- 0
  hyper$sigma_acc_cp    <- 0
  
  hyper$lambda_sigma_ncp <- 2.38^2
  hyper$Sigma_sigma_ncp  <- 1e-8
  hyper$mu_sigma_ncp     <- 0
  hyper$sigma_acc_ncp    <- 0
  
  
  hyper$a_chi <- 2
  hyper$b_chi <- 40
  
  return(hyper)
}


g_phi <- function(phi_star, a_phi, b_phi){
  phi <- (exp(phi_star)/(1 + exp(phi_star)))*(b_phi - a_phi) + a_phi
  return(phi)
  
}

project_warp <- function(S_mat, W, thetas) {
  #t_mat <- cbind(-sin(thetas), cos(thetas))
  S_mat %*% W
  #s <- rowSums(Delta * t_mat)
  #s * t_mat   # k x 2, row h = s_h * t_h
}

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
  output$mean <- array(NA, c(mcmc_iter, k_l, 2))
  output$phi <- matrix(NA, nrow = mcmc_iter, ncol =1)
  output$psi <- array(NA, dim = c(mcmc_iter, k_l - 5, 2))
  output$sigma_2 <- matrix(NA, nrow = mcmc_iter, ncol = 1)
  output$z <- matrix(NA, nrow = mcmc_iter, ncol = 1)
  return(output)
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
  init$mean <- cbind(mean) # average configuration 
  hyper$S_mat <- S_mat_constructor(hyper$S_mat, mean)
  hyper$D_p <- D_p_constructor(mean, init$thetas, init$betas,init$r, n_basis, degree = degree)
  init$z = 1
  Precision <- t(hyper$D_p) %*% hyper$S_mat %*% hyper$D_p
  K_inv <- chol2inv(chol(Precision)) # variance
  k_free <- ncol(hyper$D_p)
  V <- diag(1/diag(K_inv), nrow = k_free, ncol =k_free) #precision 
  V <- (1/init$sigma_2) * V
  if (init$z == 1) {
    U_prec <- chol(Precision/init$sigma_2)
  }
  else{
    U_prec <- chol(V)
  }
    # k - 3
  z_raw <- rnorm(2 * k_free)
  init$psi <-  matrix(backsolve(U_prec, z_raw), nrow = k_free, ncol = 2)
  
  #init$psi <- matrix(0, nrow = k_free, ncol = 2 )
  
  init$W <- hyper$D_p %*% init$psi
  
  for (i in 1:n) { # compute rotation matrix 
    init$R[,,i] <- matrix(c(cos(init$lambdas[i]), -sin(init$lambdas[i]), 
                                  sin(init$lambdas[i]), cos(init$lambdas[i])), 
                                byrow = T, nrow = 2, ncol = 2)
    eta_matrix <- matrix(init$eta[i,], nrow = length(init$thetas), ncol = 2, byrow = T)
    init$mean_i[,,i] <- init$alphas[i]*(init$mean + hyper$S_mat%*%init$W)%*%init$R[,,i] + eta_matrix # compute the mean configuration for every unit
    #init_param$mean_i[,,i] <- init_param$alphas[i]*(init_param$mean + project_warp(hyper$S_mat, init_param$W, init_param$thetas))%*%init_param$R[,,i] +
    #  eta_matrix # compute the mean configuration for every unit
    
    init$Q_R[,,i] <- (1/(init$alphas[i]^2)) *t(init$R[,,i])%*%init$Sigma_inv%*%init$R[,,i] 
    
    res <- backsolve(init$chol_c,X[,,i] - init$mean_i[,,i],  transpose = TRUE)
    #init_param$residual[,,i] <- t(X[,,i] - init_param$mean_i[,,i])%*%init_param$C %*% (X[,,i] - init_param$mean_i[,,i])
    init$residual[,,i] <- t(res)%*%res
  }
  K <- t(hyper$D_p) %*% hyper$S_mat%*% hyper$D_p
  Kinv <- chol2inv(chol(K))
  w    <- 1 / diag(Kinv)
  w_mat <- diag(w, nrow = k_free)
  V <-  mean(diag(hyper$S_mat %*% hyper$D_p %*% solve(K) %*% t(hyper$D_p) %*% t(hyper$S_mat)))
  lp_sig <- -log(V) - log(1 + (init$sigma_2/(V*hyper$a_sigma))^2)
  logdetK <- as.numeric(determinant(K, logarithm = TRUE)$modulus)
  if (init$z == 1) {
    hyper$traces <- -0.5 * (1/init$sigma_2) * tr(t(init$W) %*% hyper$S_mat %*% init$W) +
      logdetK - k_free*log(init$sigma_2) + lp_sig
  }
  else{
    
    hyper$traces <- -0.5 * (1/init$sigma_2) * sum(init$psi * (w * init$psi)) + sum(log(w)) - k_free*log(init$sigma_2) + lp_sig
  }
  hyper$D_ref <- D_p_raw(mean, init$thetas, init$betas, init$r, n_basis, degree)
  
  # return(init_param)
}


log_density <- function(init){
  val <- sum(init$Q_R[1,1,] * init$residual[1,1,] +
               init$Q_R[1,2,] * init$residual[2,1,] +
               init$Q_R[2,1,] * init$residual[1,2,] +
               init$Q_R[2,2,] * init$residual[2,2,])
  
  val_2 <- -2*init$k_l*sum(log(init$alphas)) - ((init$n*init$k_l)/2)*log(det(init$Sigma)) - ((init$n))*log(det(init$C))
  return(-0.5 * val + val_2)
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
