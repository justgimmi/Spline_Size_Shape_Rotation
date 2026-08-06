source("src/Packages.R")
source("src/Utils.R")
source("src/sample_tau.R")
source("src/sample_thetas.R")
source("src/sample_lambdas.R")
source("src/sample_Sigma.R")
source("src/sample_alphas.R")
source("src/sample_betas.R")
source("src/sample_eta.R")
source("src/sample_phi.R")
total_iter = 100000
burnin = 30000
thinning  = 10
par(mfrow = c(1, 3))
# MCMC Sampler for all the parameters -----
RB_MCMC <- function(total_iter, burnin, thinning, 
                    X, init, hyper){
  # total_iter : --> total number of iterations for the sampler
  # burnin : --> burnin for the chain
  # thinning: --> how much thinning do we apply
  # X: --> array of size K x p x n (where K is the number of lanmdarks, p the dimension and n the number of observations)
  # init_param --> starting values for the parameters (it depends on how we decide to build the order)
  # hyper --> list of all the hyperparameters 
  
  init_env <- as.environment(init)
  hyper_env <- as.environment(hyper)
  #init_env$gammas <- hyper_env$A_bar %*% init_env$betas
  mcmc_iter <- floor((total_iter - burnin)/thinning) # how many iterations do we need to save?
  total_steps <- burnin + (mcmc_iter - 1) * thinning
  #pbar <- cli_progress_bar(total = total_steps)
  k_l <- dim(X)[1] # number of landmarks
  n <- init$n
  output <- create_output(mcmc_iter, k_l, init$n, hyper$degree + hyper$n_basis) # create the output list
  # pb <- txtProgressBar(min = 1, max = mcmc_iter, style = 3, char = "*")
  pbar <- cli_progress_bar(
    name = "MCMC Sampler",
    total = total_steps,
    type = "iterator",
    format = "{cli::pb_name} {cli::pb_bar} {cli::pb_percent} | ETA: {cli::pb_eta} | Time Elapsed: {cli::pb_elapsed}"
  )
  thin <- burnin
  mean_constructor(hyper_env$n_basis, hyper_env$degree, init_env, X) # construct the mean configuration

  burn  = TRUE
  hyper_env$B_sim <- Basis_Construction(init_env$thetas, hyper_env$n_basis, hyper_env$degree)
  hyper_env$log_lik <- log_density(init_env)
  dev.off()
  par(mfrow = c(1, 3))
  for (i in 1:mcmc_iter) {
    
    for (j in 1:thin) {
      hyper_env$t <- hyper_env$t + 1
      cli_progress_update(id = pbar, inc = 1)
      sample_phi(init_env, hyper_env, X, burn = burn)
      sample_Sigma(init_env, hyper_env, k_l)
      # sample_theta_gaps(init_env, hyper_env, X)
      #sample_theta_joint(init_env, hyper_env, X)
      #sample_theta_joint(init_env, hyper_env, X, burn)
      #sample_theta_logistic(init_env, hyper_env, X, burn)
      sample_theta_logit_normal(init_env, hyper_env, X)
      sample_beta(init_env, hyper_env, X, k_l)
      for (indice in 1:n) {
        #sample_eta(init_env, hyper_env, indice, X) 
        sample_lambda_eta(init_env, hyper_env, indice, X, burn)
      }
      hyper_env$log_lik <- log_density(init_env)
      
      #sample_lambda(init_env, hyper_env, X, k_l)
      
      #hyper_env$a_new <- hyper_env$L/2 + hyper_env$a_tau # shape parameter
      hyper_env$a_new <- ncol(hyper_env$C_bar)/2 + hyper_env$a_tau
      #hyper_env$b_new <- (t(init_env$betas)%*%hyper_env$K%*%init_env$betas)/2+ hyper_env$b_tau # scale parameter
      #hyper_env$b_new <- (t(init_env$gammas)%*%init_env$gammas)/2+ hyper_env$b_tau
      hyper_env$b_new <-  (t(init_env$betas) %*% hyper_env$K %*% init_env$betas)/2 + hyper_env$b_tau
      sample_tau(init_env, a_new = hyper_env$a_new, b_new = hyper_env$b_new) # sample new tau^2

      for (y in 1:n) {
        hyper_env$log_lik_a[y] <- eval_log_lik_alpha(log(init_env$alphas[y]), hyper_env, init_env, k_l, X, y)$log_lik

      }
      sample_alpha(init_env, hyper_env, X, k_l, burn)

      
      if (i == 1 & j %% 50 == 0) {
        # cli_inform(paste0("\n[Burn-in - Iterazione globale: ", hyper$t, "]"))
        #print(hyper_env$a_aver_a/hyper_env$t)
        # print(init$Sigma)
        # print(init_env$tau)
        #print(init_env$betas)
        # print(hyper_env$Eigen_vector_null%*%init_env$betas)
        # print(init_env$Sigma)
        # print(init_env$alphas)
        print(init_env$thetas)
        #print(init_env$thetas_sorted)
        #print(init_env$lambdas)
        #print(init_env$alphas)
        print(init_env$tau)
        print(init_env$Sigma)
        #print(init_env$eta)
        print(hyper_env$lambda_acc/j)
        print(hyper_env$a_aver_a/j)
        print(hyper_env$phi_acc/j)
        print(init_env$phi)
        X_mean <- init_env$mean
        plot(X_mean, asp = 1, pch = 21, bg = "darkblue", main = "Latent mu")
        lines(c(X_mean[,1], X_mean[1,1]), c(X_mean[,2], X_mean[1,2]), col = "blue", lwd = 2)
        X_new_i <- init_env$mean_i[,,1]
        plot(X_new_i, asp = 1, pch = 21, bg = "black", main = "Latent mu unit i")
        lines(c(X_new_i[,1], X_new_i[1,1]), c(X_new_i[,2], X_new_i[1,2]), col = "brown", lwd = 2)
        
        plot(X[,,1], asp = 1, pch = 21, bg = "black", main = "Observed Unit i")
        lines(c(X[,1,1], X[1,1,1]), c(X[,2,1], X[1,2,1]), col = "brown", lwd = 2)
      }
    
      
      
    }
    if(i == 1){ # for the first iteration we are doing burnin, then we are doing thinning
    thin <- thinning
    burn <- FALSE
    }
    
    if (i %% 50 == 0) {
      # cli_inform(paste0("\n[Burn-in - Iterazione globale: ", hyper$t, "]"))
      #print(hyper_env$a_aver_a/hyper_env$t)
      # print(init$Sigma)
      # print(init_env$tau)
      # print(init_env$betas)
      # print(hyper_env$Eigen_vector_null%*%init_env$betas)
      # print(init_env$Sigma)
      # print(init_env$alphas)
      # print(init_env$thetas)
      # print(init_env$lambdas)
      # print(init_env$alphas)
      # print(init_env$tau)
      # print(init_env$Sigma)
      X_mean <- init_env$mean
      plot(X_mean, asp = 1, pch = 21, bg = "darkblue", main = "Estimated Mean")
      lines(c(X_mean[,1], X_mean[1,1]), c(X_mean[,2], X_mean[1,2]), col = "blue", lwd = 2)
      X_new_i <- init_env$mean_i[,,1]
      plot(X_new_i, asp = 1, pch = 21, bg = "black", main = "estimated Conf for unit i")
      lines(c(X_new_i[,1], X_new_i[1,1]), c(X_new_i[,2], X_new_i[1,2]), col = "brown", lwd = 2)
      
      plot(X[,,1], asp = 1, pch = 21, bg = "black", main = "estimated Conf for unit i")
      lines(c(X[,1,1], X[1,1,1]), c(X[,2,1], X[1,2,1]), col = "brown", lwd = 2)
    }
    output$tau[i] <- init_env$tau
    output$theta[i,] <- init_env$thetas
    output$lambdas[i,] <- init_env$lambdas%%(2*pi)
    output$Sigma[i,,] <- init_env$Sigma
    output$alphas[i, ] <- init_env$alphas
    output$betas[i, ] <- init_env$betas
    output$eta[i,,] <- init_env$eta
    output$mean_i[[i]] <- init_env$mean_i
    output$mean[i, ,] <- init_env$mean
    output$gammas[i, ] <- init_env$gammas
    output$phi[i,1] <- init_env$phi
    #setTxtProgressBar(pb, value = i, title = "Ziocan")
    #
    
  }
  cli_progress_done(id = pbar)
  out <- list()
  out$output <- output
  out$init <- init_env
  return(out)
  
}


X <- sam$X
# X <- X/100
X <- sam/100
#X <- sam

X_centered <- array(NA, dim = dim(X))


for (i in 1:dim(X)[3]) {
  X_centered[,,i] <- X[,,i] - matrix(colMeans(X[,,i]), nrow = dim(X)[1], ncol = 2, byrow = TRUE)
}

gpa <- procGPA(X_centered, scale = FALSE)
#?procGPA
mean_shape <- gpa$mshape
theta_init <- atan2(mean_shape[,2], mean_shape[,1])
theta_init <- ifelse(theta_init < 0, theta_init + 2*pi, theta_init)

mat_dist <- diag(0, nrow = length(theta_init), ncol = length(theta_init)) 
mat_dist <- covariance_mat(mat_dist,theta_init)
b_phi <- 2*max(mat_dist)/3
#b_phi <- pi
3/max(mat_dist)
3/min(mat_dist[mat_dist != 0])
a_phi <- min(mat_dist[mat_dist != 0])/3
hyper <- hyperparameters(1, 5e-4, 12, 3, width_theta = 1, m = 6, nu = 4,
                         psi = diag(0.1, nrow = 2),n = dim(X)[3], a = 0.01, b = 0.01, a_phi = a_phi, b_phi = b_phi, X = X)
hyper$n_basis
ord <- order(theta_init)
theta_init <- theta_init[ord]
mean_shape_ordered <- mean_shape[ord, ]


X_new <- array(NA, dim = dim(X))
for(i in 1:dim(X)[3]) {
  X_new[,,i] <- X[ord, , i] # Ordiniamo i landmark per il sarago i
}

plot(X[,,1])
for (i in 1:18) {
  text(X[i,1,1], X[i, 2, 1], i)
  
}

plot(X_new[,,1])
for (i in 1:18) {
  text(X_new[i,1,1], X_new[i, 2, 1], i)
  
}

X <- X_new
#theta_init <- sort(theta_init)

r_emp <- sqrt(mean_shape_ordered[,1]^2 + mean_shape_ordered[,2]^2)
log_r_emp <- log(r_emp)
B_init <- Basis_Construction(theta_init, 12, 3)

beta_init <- solve(t(B_init) %*% B_init + 1e-1 * diag(ncol(B_init)), t(B_init) %*% log_r_emp)

sum(beta_init)
gamma_init <- solve(
  t(hyper$C_bar) %*% hyper$C_bar
) %*%
  t(hyper$C_bar) %*%
  beta_init

#gamma_init_ridotto <- as.vector(hyper$A_bar %*% beta_init)
beta_init <- as.vector(hyper$C_bar %*% gamma_init)
sum(beta_init)
r_pure <- as.numeric(exp(B_init %*% beta_init))
mu_pure <- cbind(r_pure * cos(theta_init), r_pure * sin(theta_init))
mu_pure_centered <- scale(mu_pure, scale = FALSE)


alpha_sync  <- numeric(dim(X)[3])
lambda_sync <- numeric(dim(X)[3])
eta_sync    <- matrix(NA, nrow = dim(X)[3], ncol = 2)

for (i in 1:dim(X)[3]) {
  X_i_centered  <- scale(X[,,i], scale = FALSE)
  eta_sync[i, ] <- colMeans(X[,,i])
  
  M <- t(X_i_centered) %*% mu_pure_centered
  svd_decomp <- svd(M)
  
  # 2. Derive the optimal rotation matrix mapping directly from MU -> X_i
  # We use V %*% t(U) instead of U %*% t(V)
  R_i <- svd_decomp$v %*% t(svd_decomp$u)
  
  # 3. Extract the rotation angle from the correct matrix orientation
  lambda_sync[i] <- atan2(R_i[2,1], R_i[1,1])
  
  # 4. Update the scaling factor calculation
  # We project mu into the space of X_i. Since they now align perfectly, 
  # this dot product is guaranteed to be highly positive.
  alpha_sync[i]  <- sum(X_i_centered * (mu_pure_centered %*% R_i)) / sum(mu_pure_centered^2)
}

lambda_sync <- ifelse(lambda_sync < 0, lambda_sync + 2*pi, lambda_sync)
#lambda_sync <- ifelse(lambda_sync < 0, lambda_sync + 2*pi, lambda_sync)

init_residuals <- array(NA, dim = c(dim(X)[1], 2, dim(X)[3]))
for (i in 1:dim(X)[3]) {
  # Costruisci la rotazione iniziale del pesce i
  R_i <- matrix(c(cos(lambda_sync[i]), -sin(lambda_sync[i]),
                  sin(lambda_sync[i]),  cos(lambda_sync[i])), byrow = TRUE, nrow = 2)
  # Configurazione media stimata ruotata, scalata e traslata
  mean_i <- alpha_sync[i] * (mu_pure %*% R_i) + matrix(eta_sync[i,], nrow = dim(X)[1], ncol = 2, byrow = TRUE)
  init_residuals[,,i] <- X[,,i] - mean_i
}

# Ricava la matrice di covarianza 2x2 reale dei residui iniziali
# init_residuals is k_l x 2 x n from your GPA-based initialization
S_hat <- matrix(0, 2, 2)
for (i in 1:dim(X)[3]) {
  S_hat <- S_hat + t(init_residuals[,,i]) %*% init_residuals[,,i]
}
S_hat <- S_hat / (dim(X)[3] * dim(X)[1])   # average per-landmark-per-subject outer product

nu_prior <- 4
psi <- (nu_prior - 2 - 1) * S_hat   # so that E[Sigma_e] = S_hat
#Sigma_init <- diag(0.001, nrow = 2)
#lambda_init <- rep(0, 120)
init <- init_param(2,lambda_sync, theta_init,beta_init,
                   S_hat,
                   eta_sync,alpha_sync, gamma_init, phi = 0.3, n = dim(X)[3])
# 
# hyper <- hyperparameters(1, 5e-4, 10, 3, width_theta = 1, m = 6, nu = 4,
#                          psi = psi,n = dim(X)[3], a = 0.01, b = 0.01, a_phi = a_phi, b_phi = b_phi, X = X)
# plot(X[,,1])
# for (i in 1:25) {
#   text(X[i,1,1], X[i, 2, 1], i)
#   
# }
# init <- init_param(0.01, runif(n = 100)*2*pi, sort(runif(n = 25)*2*pi),rmvnorm_rd(n = 1, mu = 0, Precision =  (1/0.01)*hyper$K, tol = 1e-6),
#                    diag(0.01, nrow = 2),
#                    eta_sync,rgamma(n = 100, 1, 1), n <- dim(X)[3])
# init <- init_param(0.01, runif(n = 100)*2*pi, sam$theta,sam$betas, diag(0.01, nrow = 2),
#                    sam$eta,rgamma(n = 100, 1, 1), n <- dim(sam$X)[3])

tic()
par(mfrow = c(1, 3))
MCMC_samp <- RB_MCMC(total_iter = 15000, burnin = 5000, thinning  = 10, X = X, 
                     init = init, hyper = hyper)
toc()
gcinfo(FALSE)
total_iter = 100000
burnin = 30000
thinning  = 20
X <- sam$X


sample_posterior_predictive_identita <- function(output_mcmc, subject_idx) {
  
  n_saved_iter <- dim(output_mcmc$mean)[1]
  k_l <- dim(output_mcmc$mean)[2] # Numero di landmark (K)
  p <- 2                           # Dimensione spaziale (x, y)
  
  X_pred <- array(0, dim = c(k_l, p, n_saved_iter))
  
  for (i in 1:n_saved_iter) {
    # 1. Recupera la media già calcolata del soggetto i a questa iterazione
    M_subject <- output_mcmc$mean_i[[i]][, , subject_idx]
    
    # 2. Recupera i parametri di trasformazione del soggetto e la Sigma_e globale
    alpha_i   <- output_mcmc$alphas[i, subject_idx]
    lambda_i  <- output_mcmc$lambdas[i, subject_idx]
    Sigma_e   <- output_mcmc$Sigma[i, , ]
    
    # 3. Costruisci la matrice di rotazione R(\lambda_i)
    R_lambda <- matrix(c(cos(lambda_i), -sin(lambda_i),
                         sin(lambda_i),  cos(lambda_i)), 
                       nrow = 2, ncol = 2, byrow = TRUE)
    
    # 4. Calcola la covarianza spaziale specifica del soggetto (2 x 2)
    # \alpha_i^2 * t(R_lambda) %*% Sigma_e %*% R_lambda
    Sigma_spatial_i <- (alpha_i^2) * (t(R_lambda) %*% Sigma_e %*% R_lambda)
    L_Sigma_spatial <- t(chol(Sigma_spatial_i)) # Cholesky inferiore (2 x 2)
    
    # 5. Campiona la matrice di rumore bianco standard Z (K x 2)
    Z <- matrix(rnorm(k_l * p), nrow = k_l, ncol = p)
    
    # 6. Trasforma il rumore solo a destra (poiché L_C è l'identità)
    # E_i = Z %*% t(L_Sigma_spatial)
    E_i <- Z %*% t(L_Sigma_spatial)
    
    # Configurazione finale predetta per l'iterazione i
    X_pred[, , i] <- M_subject + E_i
  }
  
  return(X_pred)
}


posterior_pred <- sample_posterior_predictive_identita(MCMC_samp$output, 1)

dim(posterior_pred)
length(MCMC_samp$output$mean_i[[1]])


save(MCMC_samp,X, file = "Fish_new.RData")
load("Sim1.RData")

par(mfrow = c(1, 1))
plot((MCMC_samp$output$lambdas[,5] - MCMC_samp$output$lambdas[,1])%%(2*pi), type ="l")
library(LaplacesDemon)
ESS(MCMC_samp$output$lambdas[,2] - MCMC_samp$output$lambdas[,1])
plot((MCMC_samp$output$theta[,9] - MCMC_samp$output$theta[,1]) %%(2*pi), type = "l")
ESS((MCMC_samp$output$theta[,2] - MCMC_samp$output$theta[,1]) %%(2*pi))
plot(MCMC_samp$output$phi, type = "l")
plot(density(MCMC_samp$output$phi))

ESS(MCMC_samp$output$theta[,5] - MCMC_samp$output$theta[,2])
sam$theta[3] - sam$theta[1]
ESS(MCMC_samp$output$alphas[,2]/MCMC_samp$output$alphas[,1])
sam$lambda[2] - sam$lambda[1]
sam$alphas[2]/sam$alphas[1]
gg_mcmc_diagnostics(MCMC_samp$output$tau, param_name = "tau", real_values = NA)
gg_mcmc_diagnostics(MCMC_samp$output$theta, param_name = "theta", real_values = NA, TRUE)
gg_mcmc_diagnostics(MCMC_samp$output$lambda, param_name = "lambda", real_values = NA, TRUE)
Sigma_samp <- cbind(MCMC_samp$output$Sigma[,1,1], MCMC_samp$output$Sigma[,1,2], MCMC_samp$output$Sigma[,1,2], MCMC_samp$output$Sigma[,2,2])
gg_mcmc_diagnostics(Sigma_samp, param_name = "Sigma", real_values = c(NA), TRUE)
gg_mcmc_diagnostics(MCMC_samp$output$phi, param_name = "phi", real_values = c(NA), TRUE)
gg_mcmc_diagnostics(MCMC_samp$output$alphas, param_name = "alpha", real_values = NA, TRUE)
gg_mcmc_diagnostics(MCMC_samp$output$betas, param_name = "betas", real_values = NA, TRUE)
gg_mcmc_diagnostics(MCMC_samp$output$gammas, param_name = "betas", real_values = NA, TRUE)
gg_mcmc_diagnostics(MCMC_samp$output$eta[,,1], param_name = "etas", real_values = NA, TRUE)
par(mfrow = c(1,1))
acf(MCMC_samp$output$theta[seq(1, 4000, 50),1])

# for (i in 1:n) {
#   
#   rot <- sam$mu
#   mu_i <- sam$mu_i[,,i]
#   X_i <- sam$X[,,i]
#   X_new_i <- init_env$mean_i[,,i]
#   X_mean <- init_env$mean
# 
#   plot(rot, asp = 1, pch = 21, bg = "darkgreen", main = "Latent mu")
#   lines(c(rot[,1], rot[1,1]), c(rot[,2], rot[1,2]), col = "forestgreen", lwd = 2)
#   
#   plot(X_mean, asp = 1, pch = 21, bg = "darkblue", main = "Estimated Mean")
#   lines(c(X_mean[,1], X_mean[1,1]), c(X_mean[,2], X_mean[1,2]), col = "blue", lwd = 2)
#   
#   plot(mu_i, asp = 1, pch = 21, bg = "darkred", main = "mu for the i-th unit")
#   lines(c(mu_i[,1], mu_i[1,1]), c(mu_i[,2], mu_i[1,2]), col = "red", lwd = 2)
#   
#   plot(X_new_i, asp = 1, pch = 21, bg = "black", main = "estimated Conf for unit i")
#   lines(c(X_new_i[,1], X_new_i[1,1]), c(X_new_i[,2], X_new_i[1,2]), col = "brown", lwd = 2)
#   
#   
# }

apply(MCMC_samp$output$betas[,], MARGIN = 2, FUN = mean)
sam$betas
thetas <- output$theta[,1]
betas <- output$betas[,5]
for (i in 1:init_env$k_l) {
  for (j in 1:13) {
    thetas <- output$theta[,i]
    betas <- output$betas[,j]
    print(cor(thetas[!is.na(thetas)], betas[!is.na(betas)]))
  }
  
}

pdf("Risultati_Pesciolini_Rotation.pdf")
par(mfrow = c(1, 3))
n = 120
n_iter <- dim(MCMC_samp$output$mean)[1]
#i = 1
k_l <- dim(MCMC_samp$output$mean)[2]
par(mfrow = c(1, 3))
for (i in 1:n) {
  
  X_i <- X[,,i]
  
  X_mean <- MCMC_samp$init$mean
  
  mean_samp <- MCMC_samp$output$mean  # [iter, k, 2]
  
  mean_low <- apply(mean_samp, c(2,3), quantile, probs = 0.025)
  mean_high <- apply(mean_samp, c(2,3), quantile, probs = 0.975)
  mean_samp_mean <- apply(mean_samp, c(2,3), mean)
  
  plot(mean_samp_mean, asp = 1, pch = 21, bg = "darkblue",
       main = "Estimated Mean + CI")
  
  lines(c(mean_samp_mean[,1], mean_samp_mean[1,1]),
        c(mean_samp_mean[,2], mean_samp_mean[1,2]),
        col = "blue", lwd = 2)
  
  for (k in 1:nrow(X_mean)) {
    
    polygon(
      x = c(mean_low[k,1], mean_high[k,1],
            mean_high[k,1], mean_low[k,1]),
      y = c(mean_low[k,2], mean_low[k,2],
            mean_high[k,2], mean_high[k,2]),
      border = NA,
      col = "red"
    )
  }
  
  plot(X_i, asp = 1, pch = 21, bg = "darkred",
       main = "i-th unit")
  
  lines(c(X_i[,1], X_i[1,1]),
        c(X_i[,2], X_i[1,2]),
        col = "red", lwd = 2)
  
  samp_i <- array(NA, dim = c(n_iter, k_l, 2))
  for (j in 1:n_iter) {
    samp_i[j,,] <- MCMC_samp$output$mean_i[[j]][,,i]
    
  }
  
  
  #samp_i <- MCMC_samp$output$mean_i[[i]]
  mean_samp_i <- apply(samp_i, c(2,3), mean)
  low_i <- apply(samp_i, c(2,3), quantile, 0.025)
  high_i <- apply(samp_i, c(2,3), quantile, 0.975)
  plot(mean_samp_i, asp = 1, pch = 21, bg = "black",
       main = "Estimated config + CI")
  
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
  
  # post_samp_i <- apply(posterior_pred, c(1,2), mean)
  # low_post_i <- apply(posterior_pred, c(1,2), quantile, 0.025)
  # high_post_i <- apply(posterior_pred, c(1,2), quantile, 0.975)
  # plot(post_samp_i, asp = 1, pch = 21, bg = "black",
  #      main = "Posterior Predictive")
  # lines(c(post_samp_i[,1], post_samp_i[1,1]),
  #       c(post_samp_i[,2], post_samp_i[1,2]),
  #       col = "brown", lwd = 2)
}


gg_mcmc_diagnostics(MCMC_samp$output$tau, param_name = "tau", real_values = NA)
gg_mcmc_diagnostics(MCMC_samp$output$theta, param_name = "theta", real_values = NA, TRUE)
gg_mcmc_diagnostics(MCMC_samp$output$lambda, param_name = "lambda", real_values = NA, TRUE)
Sigma_samp <- cbind(MCMC_samp$output$Sigma[,1,1], MCMC_samp$output$Sigma[,1,2], MCMC_samp$output$Sigma[,1,2], MCMC_samp$output$Sigma[,2,2])
gg_mcmc_diagnostics(Sigma_samp, param_name = "Sigma", real_values = c(NA), TRUE)
gg_mcmc_diagnostics(MCMC_samp$output$alphas, param_name = "alpha", real_values = NA, TRUE)
gg_mcmc_diagnostics(MCMC_samp$output$betas, param_name = "betas", real_values = NA, TRUE)

save(MCMC_samp,sam, X,  file = "prova.RData")
par(mfrow = c(1, 3))
n = 100
for (i in 1:n) {

  #rot <- sam$mu
  #mu_i <- sam$mu
  X_i <- X[,,i]
  X_new_i <- init_env$mean_i[,,i]
  X_mean <- init_env$mean

  # plot(rot, asp = 1, pch = 21, bg = "darkgreen", main = "Latent mu")
  # lines(c(rot[,1], rot[1,1]), c(rot[,2], rot[1,2]), col = "forestgreen", lwd = 2)

  plot(X_mean, asp = 1, pch = 21, bg = "darkblue", main = "Estimated Mean")
  lines(c(X_mean[,1], X_mean[1,1]), c(X_mean[,2], X_mean[1,2]), col = "blue", lwd = 2)

  plot(X_i, asp = 1, pch = 21, bg = "darkred", main = "i-th unit")
  lines(c(X_i[,1], X_i[1,1]), c(X_i[,2], X_i[1,2]), col = "red", lwd = 2)

  # plot(mu_i, asp = 1, pch = 21, bg = "darkred", main = "mu for the i-th unit")
  # lines(c(mu_i[,1], mu_i[1,1]), c(mu_i[,2], mu_i[1,2]), col = "red", lwd = 2)

  plot(X_new_i, asp = 1, pch = 21, bg = "black", main = "estimated Conf for unit i")
  lines(c(X_new_i[,1], X_new_i[1,1]), c(X_new_i[,2], X_new_i[1,2]), col = "brown", lwd = 2)


}

dev.off()


cor(thetas[!is.na(thetas)], betas[!is.na(betas)])
plot(thetas, betas)

for (i in 1:init_env$k_l) {
  
  thetas <- output$theta[,i]
  betas <- output$betas[,i]
  plot(thetas, betas)
  
  
}



MCMC_samp$output$mean[1, ,]

X <- sam

output <- MCMC_samp$output

X0 <- X[, , 1]
sam$X
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
pair
load("risultati_new.RData")
output <- MCMC_samp$output
#pair <- c(1, 12)
align_prof_method <- function(output, pair, third_landmark = NULL) {
  
  n_iter <- dim(output$mean)[1]
  k <- dim(output$mean)[2]
  n <- dim(output$lambdas)[2]
  mean_aligned <- array(0, dim = c(n_iter, k, 2))
  eta_aligned <- array(0, dim = dim(output$eta))
  
  alphas_aligned  <- matrix(0, n_iter, ncol(output$alphas))
  betas_aligned   <- matrix(0, n_iter, ncol(output$betas))
  lambdas_aligned <- matrix(0, n_iter, ncol(output$lambdas))
  theta_aligned   <- matrix(0, n_iter, ncol(output$theta))
  
  Sigma_aligned <- array(0, dim(output$Sigma))
  for (i in 1:n_iter) {
    
    M_raw <- output$mean[i,,]

    # Remove Translatiomn
    
    shift <- M_raw[pair[1],]
    M <- sweep(M_raw, 2, shift, "-")
    mat <- matrix(shift, byrow = T, nrow = k, ncol = 2)
    for (j in 1:n) {
      R <- matrix(c(cos(output$lambdas[i, j]), -sin(output$lambdas[i, j]),
                    sin(output$lambdas[i, j]),  cos(output$lambdas[i, j])),
                  2,2,byrow=TRUE)
      eta_aligned[i,j,] <- output$eta[i,j,] + output$alphas[i, j] * (mat%*%R)[1,]
      
    }
    theta_aligned[i,] <- atan2(M[,2], M[,1]) %% (2*pi)
    r_new <- sqrt(M[,1]^2 + M[,2]^2)
    idx_zero <- which.min(r_new)
    theta_fit <- theta_aligned[i, -idx_zero]
    r_fit <- r_new[-idx_zero]
    #idx_order <- order(theta_fit)
    #theta_fit <- theta_fit[idx_order]
    #r_fit <- r_fit[idx_order]
    log_r_emp <- log(r_fit)
    B_init <- Basis_Construction(theta_fit, 10, 3)
    K <- hyper$K
    tau2 <- output$tau[i] 
    
    betas_aligned[i,] <- solve(t(B_init) %*% B_init + (1/sqrt(tau2)) * K) %*% t(B_init) %*% log_r_emp
    
    v <- M[pair[2],]

    phi <- atan2(v[2], v[1])%%(2*pi)

    R <- matrix(c(cos(phi), -sin(phi),
                  sin(phi),  cos(phi)),
                2,2,byrow=TRUE)

 
    M <- M %*% R

    lambdas_aligned[i,] <- (output$lambdas[i,] - phi)%%(2*pi)
    theta_aligned[i,]   <- (theta_aligned[i,] + phi) %%(2*pi)
  

    scale <- sqrt(sum(M^2))
    M <- M / scale

    alphas_aligned[i,] <- output$alphas[i,] * scale
    betas_aligned[i,]  <- betas_aligned[i,] - log(scale)


    mean_aligned[i,,] <- M
    
    Sigma_i <- (t(R) %*% output$Sigma[i,,] %*% R) / (scale^2)
    Sigma_aligned[i,,] <- Sigma_i
  }
  
  output$mean <- mean_aligned
  output$eta <- eta_aligned
  output$alphas <- alphas_aligned
  output$betas <- betas_aligned
  output$lambdas <- lambdas_aligned
  output$theta <- theta_aligned
  output$Sigma <- Sigma_aligned
  
  output
}
output = MCMC_samp$output
output_identificato <- align_prof_method(output, pair = pair)
gg_mcmc_diagnostics(output_identificato$theta[,7] - output_identificato$theta[,1], param_name = "alpha Ratio", real_values = NA, TRUE)
dim(output_identificato$mean[,,1])
output$mean[,,1]
ESS(output_identificato$mean[,,2])
prova <- array(0, dim = c(dim(output$alphas)[1],2,2))

for (i in 1:dim(output$alphas)[1]) {
  phi <- output_identificato$lambdas[i,2]
  R <- matrix(c(cos(phi), -sin(phi),
                sin(phi),  cos(phi)),
              2,2,byrow=TRUE)
  prova[i,,] <- (output_identificato$alphas[i,2]^2)*(t(R) %*%output_identificato$Sigma[i,,] %*% R)
  
}
prova[,1,1]

solve(t(R))
pdf("Risulati_Pesci_new.pdf")
k = 19
k_l = 19
prova <- array(0, dim  = c(120, k, 2))
# for (i in 1:100) {
#   eta_matrix <- matrix(sam$eta[i, ], nrow = k, ncol = 2, byrow = TRUE)
#   mu_i <- sam$X[,,i] - eta_matrix
#   
#   mu_i <- mu_i / sam$alphas[i]
#   
#   R_inv <- Rmat(sam$lambda[i])
#   mu_i <- mu_i %*% R_inv
#   
#   # 4. PLOT
#   prova[i,,] <- mu_i
#   # plot(mu_i, asp = 1, pch = 21, bg = "darkred", 
#   #      main = paste("Mu ricostruito per l'unità", i),
#   #      xlim = c(-2, 2), ylim = c(-1.5, 1.5)) 
#   # lines(c(mu_i[,1], mu_i[1,1]), c(mu_i[,2], mu_i[1,2]), col = "red", lwd = 2)
# }
par(mfrow = c(1, 1))
mean_behav <- apply(prova, MARGIN = c(2, 3), FUN = mean)
plot(mean_behav, xlim = c(-1, 1.5), pch = 16, ylim = c(-1, 1), asp = 1)
#plot()
for (i in 1:120) {
  points(prova[i,,], asp = 1, pch = 16, col = "red")
  
}
points(apply(prova, MARGIN = c(2, 3), FUN = mean), asp = 1, pch = 16,
       xlim = c(-0.5, 0.5), ylim = c(-1, 1))
output <- MCMC_samp$output
gg_mcmc_diagnostics(MCMC_samp$output$tau, param_name = "tau", real_values = NA)
#gg_mcmc_diagnostics(output$theta - output$theta[,1], param_name = "theta", real_values = NA, TRUE)
gg_mcmc_diagnostics(output$lambda[,]%%(2*pi) - output$lambda[,1]%%(2*pi)  , param_name = "lambda diff", real_values = NA, TRUE)
# Sigma_samp <- cbind(prova[,1,1], prova[,1,2],prova[,1,2], prova[,2,2])
# Sigma_samp <- cbind(output_identificato$Sigma[,1,1],
#                     output_identificato$Sigma[,1,2],output_identificato$Sigma[,1,2],
#                     output_identificato$Sigma[,2,2])
# gg_mcmc_diagnostics(Sigma_samp, param_name = "Sigma", real_values = c(NA), TRUE)
gg_mcmc_diagnostics(output$alphas/output$alphas[,1], param_name = "alpha Ratio", real_values = NA, TRUE)
#gg_mcmc_diagnostics(output_identificato$betas, param_name = "betas", real_values =NA, TRUE)
n_iter <- dim(MCMC_samp$output$alphas)[1]
par(mfrow = c(1, 3))
for (i in 1:n) {
  
  X_i <- X[,,i]
  
  X_mean <- MCMC_samp$init$mean
  
  mean_samp <- MCMC_samp$output$mean  # [iter, k, 2]
  
  mean_low <- apply(mean_samp, c(2,3), quantile, probs = 0.025)
  mean_high <- apply(mean_samp, c(2,3), quantile, probs = 0.975)
  mean_samp_mean <- apply(mean_samp, c(2,3), mean)
  
  plot(mean_samp_mean, asp = 1, pch = 21, bg = "darkblue",
       main = "Estimated Mean + CI")
  
  lines(c(mean_samp_mean[,1], mean_samp_mean[1,1]),
        c(mean_samp_mean[,2], mean_samp_mean[1,2]),
        col = "blue", lwd = 2)
  
  for (k in 1:nrow(X_mean)) {
    
    polygon(
      x = c(mean_low[k,1], mean_high[k,1],
            mean_high[k,1], mean_low[k,1]),
      y = c(mean_low[k,2], mean_low[k,2],
            mean_high[k,2], mean_high[k,2]),
      border = NA,
      col = rgb(0, 0, 1, 0.15)
    )
  }
  
  plot(X_i, asp = 1, pch = 21, bg = "darkred",
       main = "i-th unit")
  
  lines(c(X_i[,1], X_i[1,1]),
        c(X_i[,2], X_i[1,2]),
        col = "red", lwd = 2)
  
  samp_i <- array(NA, dim = c(n_iter, k_l, 2))
  for (j in 1:n_iter) {
    samp_i[j,,] <- MCMC_samp$output$mean_i[[j]][,,i]
    
  }
  
  
  #samp_i <- MCMC_samp$output$mean_i[[i]]
  mean_samp_i <- apply(samp_i, c(2,3), mean)
  low_i <- apply(samp_i, c(2,3), quantile, 0.025)
  high_i <- apply(samp_i, c(2,3), quantile, 0.975)
  plot(mean_samp_i, asp = 1, pch = 21, bg = "black",
       main = "Estimated config + CI")
  
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
  
  # post_samp_i <- apply(posterior_pred, c(1,2), mean)
  # low_post_i <- apply(posterior_pred, c(1,2), quantile, 0.025)
  # high_post_i <- apply(posterior_pred, c(1,2), quantile, 0.975)
  # plot(post_samp_i, asp = 1, pch = 21, bg = "black",
  #      main = "Posterior Predictive")
  # lines(c(post_samp_i[,1], post_samp_i[1,1]),
  #       c(post_samp_i[,2], post_samp_i[1,2]),
  #       col = "brown", lwd = 2)
}
nsamp <- dim(output$lambdas)[1]
mean_orig <- array(NA, c(k, 2, nsamp))
mean_noloc <- array(NA, c(k, 2, nsamp))
mean_norot <- array(NA, c(k, 2, nsamp))
mean_nosize <- array(NA, c(k, 2, nsamp))

sigma_orig <- array(NA, c(2, 2, nsamp))
sigma_noloc <- array(NA, c(2, 2, nsamp))
sigma_norot <- array(NA, c(2, 2, nsamp))
sigma_nosize <- array(NA, c(2, 2, nsamp))

for (iobs in 1:120) {
  
  # 1. Feedback reale sulla console R
  cat(sprintf("\rElaborazione Unità: %3d / 100 ... ", iobs))
  flush.console()
  
  krot <- pair[2]
  for (isamp in 1:nsamp){
    ## Creo oggetti originali
    rrr <- Rmat(MCMC_samp$output$lambdas[isamp, iobs])
    sigma_orig[, , isamp] <- rrr %*% MCMC_samp$output$Sigma[isamp, , ] %*% t(rrr) * MCMC_samp$output$alphas[isamp, iobs]^2
    mean_orig[, , isamp] <- MCMC_samp$output$mean_i[[isamp]][, , iobs]
    
    # Metto una coordinata a zero (Rimuovo traslazione)
    mean_noloc[, , isamp] <- mean_orig[, , isamp] - matrix(mean_orig[pair[1], , isamp], nrow = k, ncol = 2, byrow = T)
    sigma_noloc[, , isamp] <- sigma_orig[, , isamp]
    
    # Tolgo la rotazione
    angle <- atan2(mean_noloc[krot, 2, isamp], mean_noloc[krot, 1, isamp])
    mean_norot[, , isamp] <- mean_noloc[, , isamp] %*% t(Rmat(angle))
    sigma_norot[, , isamp] <- Rmat(angle) %*% sigma_noloc[, , isamp] %*% t(Rmat(angle))
    
    # Stabilizzo elevando la media campionaria
    mean_norot[, , isamp] <- mean_norot[, , isamp] - matrix(colMeans(mean_norot[, , isamp]), nrow = k, ncol = 2, byrow = T)
    normx <- norm(mean_norot[, , isamp])
    
    mean_nosize[, , isamp] <- mean_norot[, , isamp] / normx
    sigma_nosize[, , isamp] <- sigma_norot[, , isamp] / normx^2
  }
  
  # Calcoliamo la forma media a posteriori per disegnare una linea guida pulita
  mean_shape_posterior <- apply(mean_nosize, c(1, 2), mean)
  
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
      
      plot(mean_nosize[iii, ip, ], type = "l", col = "steelblue",
           xlab = "Iterazioni MCMC", ylab = coord_label,
           main = paste("Media", iii, "-", coord_label))
      
      # Linea rossa tratteggiata sul valore medio stimato per controllare la convergenza
      abline(h = mean_shape_posterior[iii, ip], col = "red", lty = 2, lwd = 1.5)
    }
  }
  
  # ====================================================================
  # PAGINA FINALE: I Traceplot della matrice di covarianza Sigma_nosize
  # ====================================================================
  # Creiamo una griglia 2x2 perfetta appositamente per la matrice 2x2 di Sigma
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
gg_mcmc_diagnostics(as.vector(MCMC_samp$output$phi), param_name = "phi", real_values = NA)
dev.off()
library(LaplacesDemon)
par(mfrow = c(1, 1))
ESS(output_identificato$eta[,3,1])
plot(output_identificato$eta[,3,1], type = "l")
ESS(output_identificato$Sigma[,2,2])
plot(output_identificato$Sigma[,2,2])
par(mfrow = c(1,1))
plot(X[,,2])
points(X[1,1,2],X[1,2,2],  col = "red", pch = 16)
points(X[12,1,2],X[12,2,2],  col = "red", pch = 16)
dim(output_identificato$mean)
plot(output_identificato$mean[1,,])
plot(output$mean[1,,])


par(mfrow = c(1, 2))
for (i in 1:100) {
  X_new_i <- output_identificato$mean[i,,]
  X_mean <- output$mean[i,,]
  
  # plot(rot, asp = 1, pch = 21, bg = "darkgreen", main = "Latent mu")
  # lines(c(rot[,1], rot[1,1]), c(rot[,2], rot[1,2]), col = "forestgreen", lwd = 2)
  
  plot(X_mean, asp = 1, pch = 21, bg = "darkblue", main = "Not id")
  lines(c(X_mean[,1], X_mean[1,1]), c(X_mean[,2], X_mean[1,2]), col = "blue", lwd = 2)
  plot(X_new_i, asp = 1, pch = 21, bg = "black", main = "Ident")
  lines(c(X_new_i[,1], X_new_i[1,1]), c(X_new_i[,2], X_new_i[1,2]), col = "brown", lwd = 2)
  
  
}

