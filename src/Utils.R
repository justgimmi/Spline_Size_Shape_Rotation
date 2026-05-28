source(file = file.path(getwd(), "src/Packages.R"))

# Basis Construction ------
# Here we define a function to define the cubic B-splines 

Basis_Construction <- function(thetas, L, degree = 3){
  # I do not like to use already implemented function because they are based
  # on the range of the vector to define the knots so every iteration of the MCMC
  # we could have different basis following this approach
  
  # thetas: --> vector of theta values where to evaluate the basis 
  # L: --> how many internal knots do we want?
  # degree: --> spline degree
  
  # The idea of this function is the following. Given the fact that I want the starting value in 0 
  # and the end point in 2 \pi, I first define equispaced knots on this interval. Then, we habe to
  # add padding on the left and on the right. A conservative choice would be to add degree times 0 before 
  # and degree times 2\pi after. Talking with Johannes, I understood that the best approach to have 
  # meaningfull penalties is to have equispaced knots also on the left and right part of the interval
  delta <- (2 * pi) / L # equispaced shift 
  
  knots_base <- seq(0, 2 * pi, by = delta)[-c(1, L+1)] # base knots
  
  total_knots <- c(seq(0 - degree*delta, 0, by = delta), knots_base, seq(2*pi , 2*pi + degree*delta, 
                                                                         by = delta))
  Basis <- splineDesign(total_knots, thetas, ord = degree + 1, outer.ok = FALSE,
                        derivs = F) # Basis construction
  # At the End we obtain a matrix of size length(theta) x (L + degree)
  return(Basis)
}


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



in_model_sample <- function(n = 1, K_l,  thetas = NA, n_int_knots, degree = 3, tau = 0.1, 
                            beta_values = NA, lambda = NA, alphas = NA, eta = NA, 
                            Sigma_e = NA){
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
  K <- K1 + 2*K2
  P <- (1/(tau^2))*K # Precision Matrix
  if (length(beta_values) == 1 ) { # the user can supply the values of beta
    beta_values <- rmvnorm_rd(n = 1, mu = rep(0, n_basis), P, tol = 1e-7)
    
  }
  
  # Mean block ----
  r <- exp(B_sim %*% as.vector(beta_values)) # radius
  mu_mean_x <- r*cos(thetas)
  mu_mean_y <- r*sin(thetas)
  mu_mean <- cbind(mu_mean_x, mu_mean_y) # mean configuration 
  R <- array(NA, dim = c(2, 2, n))
  mu_i <- array(NA, dim = c(K_l, 2, n))
  if (length(lambda) == 1) {
    lambda <- runif(n = n)*2*pi # generate rotation angles
  }
  
  if (length(alphas) == 1) {
    alphas <- rgamma(n = n, shape = 2, rate = 2) # sample the size effect
    
  }
  
  if (length(eta) == 1) {
    eta <- matrix(rnorm(n = n*2, sd = 4), ncol = 2) # sample the translation effect 
    
  }
  
  for (i in 1:n) { # compute rotation matrix 
    R[,,i] <- matrix(c(cos(lambda[i]), -sin(lambda[i]), sin(lambda[i]), cos(lambda[i])), 
                     byrow = T, nrow = 2, ncol = 2)
    eta_matrix <- matrix(eta[i,], nrow = K_l, ncol = 2, byrow = T)
    mu_i[,,i] <- alphas[i]*mu_mean%*%R[,,i] + eta_matrix # compute the mean configuration for every unit
  }
  
  # Variance block ---- 
  S <- matrix(c(2.8*1e-4, 0,
                0, 2.1*1e-4), 2, 2)
  
  Sigma_e <- riwish(4, S)# sample measurement error on the p x p
  #Sigma_e <- diag(2)*0.00035
  X <- array(NA, dim = c(K_l, 2, n))
  for (i in 1:n) {
    # to sample we use the fact that X \sim MN(mu_{i}, I_k, Sigma)
    Sigmas <- alphas[i]^2 *t(R[,,i])%*%Sigma_e%*%R[,,i]
    I <- diag(K_l)
    X[,,i] <- chol(I)%*%matrix(rnorm(n = K_l * 2), nrow = K_l, ncol = 2)%*%t(chol(Sigmas)) + mu_i[,,i]
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
  
  return(final_param)

}

