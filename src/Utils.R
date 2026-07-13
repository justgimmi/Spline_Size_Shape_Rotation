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
# 
# Basis_Construction <- function(thetas, L, degree = 3){
#   # thetas: vettore dei landmark angles in (0, 2*pi)
#   # L: numero di nodi equispaziati nell'intervallo circolare
#   # degree: grado della spline (3 = cubica)
#   
#   # 1. Definiamo i nodi base equispaziati su [0, 2*pi]
#   delta <- (2 * pi) / L
#   knots_base <- seq(0, 2 * pi, length.out = L + 1) # include 0 e 2*pi esattamente una volta
#   
#   # 2. Per le spline periodiche, estendiamo i nodi all'esterno "avvolgendo" il cerchio
#   # Nodi a sinistra: prendiamo gli ultimi nodi prima di 2*pi e li trasliamo a sinistra di 2*pi
#   left_knots <- knots_base[(L + 1 - degree):L] - 2 * pi
#   # Nodi a destra: prendiamo i primi nodi dopo lo 0 e li trasliamo a destra di 2*pi
#   right_knots <- knots_base[2:(degree + 1)] + 2 * pi
#   
#   # Uniamo i nodi in una griglia perfettamente equispaziata e coerente
#   total_knots <- c(left_knots, knots_base, right_knots)
#   
#   # 3. Generiamo la base standard sulla griglia estesa
#   # splineDesign richiede l'ordine (degree + 1)
#   ord <- degree + 1
#   Basis_raw <- splineDesign(knots = total_knots, x = thetas, ord = ord, outer.ok = TRUE)
#   
#   # 4. TRICK DELLA PERIODICITÀ: Uniamo le colonne corrispondenti ai nodi che si sovrappongono
#   # Le ultime 'degree' colonne si sovrappongono esattamente alle prime 'degree' colonne del cerchio
#   n_col_raw <- ncol(Basis_raw)
#   
#   Basis_periodic <- Basis_raw[, 1:(n_col_raw - degree), drop = FALSE]
#   
#   # Sommiamo l'effetto ciclico delle ultime colonne sulle prime
#   for(k in 1:degree){
#     Basis_periodic[, k] <- Basis_periodic[, k] + Basis_raw[, n_col_raw - degree + k]
#   }
#   
#   # Ora la matrice ha dimensione esattamente: length(thetas) x L
#   return(Basis_periodic)
# }

# L = 10

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
  K <- K1 + K2
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
  S <- matrix(c(3.5*1e-5, 0,
                0, 2.2*1e-5), 2, 2)
  
  Sigma_e <- riwish(4, S)# sample measurement error on the p x p
  #Sigma_e <- diag(2)*0.00035
  X <- array(NA, dim = c(K_l, 2, n))
  for (i in 1:n) {
    # to sample we use the fact that X \sim MN(mu_{i}, I_k, Sigma)
    Sigmas <- alphas[i]^2 *t(R[,,i])%*%Sigma_e%*%R[,,i]
    I <- diag(K_l)
    X[,,i] <- chol(I)%*%matrix(rnorm(n = K_l * 2), nrow = K_l, ncol = 2)%*%chol(Sigmas) + mu_i[,,i]
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
  
  return(final_param)

}


# MCMC Util functions ----



init_param <- function(tau, lambdas, thetas, betas, Sigma, eta, alphas, n){
  # n --> number of units 
  param <- list()
  param$n <- n
  param$k_l <- length(thetas)
  param$tau <- tau
  param$lambdas <- lambdas
  param$thetas <- thetas
  param$gaps <- c(thetas[1], diff(thetas), 2*pi- thetas[param$k_l])/(2*pi)
  #param$omega_theta <- log(param$gaps)
  param$omega_theta <- log(param$gaps) - mean(log(param$gaps))
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
  
  return(param)
}

hyperparameters <- function(a_tau, b_tau, n_basis, degree,
                            width_theta = pi/4, m = 8, nu, psi, n, a, b, Sigma_eta = diag(1000, nrow = 2, ncol = 2),
                            tol = 1e-10, X){
  
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
  K1 <- K1_construction(J) # construct K1
  K2 <- K2_construction(J) #construct K2
  hyper$K <- (K1 + K2)
  hyper$L <- qr(hyper$K)$rank
  eig <- eigen(hyper$K, symmetric = TRUE)
  ord <- order(eig$values, decreasing = TRUE)
  eig_vals <- eig$values[ord]
  eig_vecs <- eig$vectors[, ord, drop = FALSE]
  
  
  keep <- eig_vals > tol
  hyper$L <- sum(keep) # approximation for the normal distribution
  hyper$Eigen_matrix <- diag(1/sqrt(eig_vals[keep]), nrow = hyper$L) # using for sample from the prior
  hyper$Eigen_vector <- eig_vecs[, keep, drop = FALSE]
  hyper$Eigen_vector_null <- t(eig_vecs[, !keep, drop = FALSE])
  #hyper$U_lambda <- hyper$Eigen_vector %*% hyper$Eigen_matrix
  
  hyper$C_bar <- hyper$Eigen_vector %*% hyper$Eigen_matrix # \beta = C_bar \gamma
  hyper$A_bar <- diag(sqrt(eig_vals[keep]), nrow = hyper$L) %*% t(hyper$Eigen_vector) # \gamma = A_bar \beta
  
  
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
  # hyper$lambda_a <- (2.38^2)/n
  # hyper$mu_a <- matrix(0, nrow = n)
  # hyper$Sigma_a <- diag(0.01, nrow = n, ncol = n)
  # hyper$L_a <- chol(hyper$lambda_a*hyper$Sigma_a)
  # hyper$gamma_a <- 1
  # hyper$a_opt <- 0.234
  # hyper$a_aver_a <- 0
  # hyper$t <- 0
  # hyper$a <- a
  # hyper$b <- b
  # hyper$log_lik_a <- 0
  
  hyper$lambda_a <- matrix(2.38^2, nrow = n)
  hyper$mu_a <- matrix(0, nrow = n)
  hyper$Sigma_a <- matrix(1e-8, nrow = n)
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
  hyper$Sigma_eta_new <- hyper$Sigma_eta
  hyper$mu_new <- matrix(0, nrow = n, ncol = 2)
  
  hyper$lambda_lambda <- matrix(2.38^2, nrow = n)
  hyper$mu_lambda <- matrix(0, nrow = n)
  hyper$Sigma_lambda <- matrix(1e-8, nrow = n)
  hyper$gamma_lambda <- 1
  hyper$lambda_opt <- 0.44
  hyper$lambda_acc <- matrix(0, nrow = n)
  return(hyper)
}

create_output <- function(mcmc_iter, k_l, n, L){
  output <- list()
  output$tau <- numeric(mcmc_iter)
  output$theta <- matrix(NA, nrow = mcmc_iter, ncol = k_l)
  output$lambdas <- matrix(NA, nrow = mcmc_iter, ncol = n)
  output$Sigma <- array(NA, dim = c(mcmc_iter, 2, 2))
  output$alphas <- matrix(NA, nrow = mcmc_iter, ncol = n)
  output$betas <- matrix(NA, nrow = mcmc_iter, ncol = L)
  output$eta <- array(NA, c(mcmc_iter, n, 2))
  output$mean_i <- list()
  output$mean <- array(NA, c(mcmc_iter, k_l, 2))
  return(output)
}

mean_constructor <- function(n_basis, degree, init_param, X){
  # n_basis --> number of internal knots
  # degree --> degree of the polynomial
  
  n <- init_param$n
  B_sim <- Basis_Construction(init_param$thetas, n_basis, degree) # build the basis with the current value of theta
  init_param$r <- exp(B_sim %*% as.vector(init_param$betas))# compute the current value of the radius
  mu_mean_x <- init_param$r*cos(init_param$thetas)
  mu_mean_y <- init_param$r*sin(init_param$thetas)
  init_param$mean <- cbind(mu_mean_x, mu_mean_y) # average configuration 
  
  for (i in 1:n) { # compute rotation matrix 
    init_param$R[,,i] <- matrix(c(cos(init_param$lambdas[i]), -sin(init_param$lambdas[i]), 
                                  sin(init_param$lambdas[i]), cos(init_param$lambdas[i])), 
                                byrow = T, nrow = 2, ncol = 2)
    eta_matrix <- matrix(init_param$eta[i,], nrow = length(init_param$thetas), ncol = 2, byrow = T)
    init_param$mean_i[,,i] <- init_param$alphas[i]*init_param$mean%*%init_param$R[,,i] + eta_matrix # compute the mean configuration for every unit
    init_param$Q_R[,,i] <- (1/(init_param$alphas[i]^2)) *t(init_param$R[,,i])%*%init_param$Sigma_inv%*%init_param$R[,,i] 
    init_param$residual[,,i] <- t(X[,,i] - init_param$mean_i[,,i])%*%(X[,,i] - init_param$mean_i[,,i])
  }
  # return(init_param)
}


log_density <- function(init){
  val <- sum(init$Q_R[1,1,] * init$residual[1,1,] +
               init$Q_R[1,2,] * init$residual[2,1,] +
               init$Q_R[2,1,] * init$residual[1,2,] +
               init$Q_R[2,2,] * init$residual[2,2,])
  
  val_2 <- -2*init$k_l*sum(log(init$alphas)) - ((init$n*init$k_l)/2)*log(det(init$Sigma))
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
