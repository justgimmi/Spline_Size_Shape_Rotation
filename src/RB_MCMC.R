source("src/Packages.R", keep.source = TRUE)
source("src/Utils.R", keep.source = TRUE)
source("src/sample_tau.R", keep.source = TRUE)
source("src/sample_thetas.R", keep.source = TRUE)
source("src/sample_lambdas.R", keep.source = TRUE)
source("src/sample_Sigma.R", keep.source = TRUE)
source("src/sample_alphas.R", keep.source = TRUE)
source("src/sample_betas.R", keep.source = TRUE)
source("src/sample_eta.R", keep.source = TRUE)
source("src/sample_phi.R", keep.source = TRUE)
source("src/sample_psi.R", keep.source = TRUE)
source("src/sample_sigma2.R", keep.source = TRUE)
source("src/sample_xi.R", keep.source = TRUE)
source("src/sample_z.R", keep.source = TRUE)
total_iter = 100000
burnin = 30000
thinning  = 10
par(mfrow = c(1, 3))
# MCMC Sampler for all the parameters -----
RB_MCMC <- function(total_iter, burnin, thinning, 
                    X, init, hyper, plot_file = "mcmc.pdf"){
  # total_iter : --> total number of iterations for the sampler
  # burnin : --> burnin for the chain
  # thinning: --> how much thinning do we apply
  # X: --> array of size K x p x n (where K is the number of lanmdarks, p the dimension and n the number of observations)
  # init_param --> starting values for the parameters (it depends on how we decide to build the order)
  # hyper --> list of all the hyperparameters 
  
  init_env <- as.environment(init)
  hyper_env <- as.environment(hyper)
  hyper_env$Sigma_gamma_inv <-  (1/init_env$tau[1])*hyper_env$K_gamma1 + 
    (1/init_env$tau[2])*hyper_env$K_gamma2
  hyper_env$U_lambda <- t(chol(hyper_env$Sigma_gamma_inv))
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
  mean_constructor(hyper_env$n_basis, hyper_env$degree, init_env, hyper_env, X) # construct the mean configuration

  burn  = TRUE
  hyper_env$B_sim <- Basis_Construction(init_env$thetas, hyper_env$n_basis, hyper_env$degree)
  hyper_env$log_lik <- log_density(init_env)
  dev.off()
  par(mfrow = c(1, 3))
  #init_env$z <- 1
  # z_temp <- numeric(thin)
  for (i in 1:mcmc_iter) {
    
    for (j in 1:thin) {
      #print(j)
      hyper_env$t <- hyper_env$t + 1
      cli_progress_update(id = pbar, inc = 1)
      #print(init_env$mean)
      sample_tau(init_env, hyper_env, burn) 
      sample_beta(init_env, hyper_env, X, k_l)
      sample_theta_logit_normal(init_env, hyper_env, X)
      # sample_z(init_env, hyper_env, X)
      sample_sigma2_ig(init_env, hyper_env, X)
      sample_psi(init_env, hyper_env, X)
      # check_state(init_env, hyper_env, X, tag = "NA")
      # #sample_sigma2(init_env, hyper_env, X, burn = burn)
      # check_state(init_env, hyper_env, X, tag = "NA")
      hyper_env$log_lik <- log_density(init_env)
      
      #print(init_env$mean)
      for (indice in 1:n) {
        #sample_eta(init_env, hyper_env, indice, X) 
        sample_lambda_eta(init_env, hyper_env, indice, X, burn)
      }
      hyper_env$log_lik <- log_density(init_env)
      
      for (y in 1:n) {
        hyper_env$log_lik_a[y] <- eval_log_lik_alpha(log(init_env$alphas[y]), hyper_env, init_env, k_l, X, y)$log_lik

      }
      sample_alpha(init_env, hyper_env, X, k_l, burn)
      # hyper_env$a_sigma_post <- hyper_env$a_sigma + (k_l - 3)
      # hyper_env$b_sigma_post <- hyper_env$b_sigma + tr(t(init_env$W) %*% hyper_env$S_mat %*% init_env$W)/2
      # init_env$sigma_2 <- rinvgamma(n = 1, shape = hyper_env$a_sigma_post, scale = hyper_env$b_sigma_post)
      #sample_phi(init_env, hyper_env, X, burn = TRUE)
      sample_Sigma(init_env, hyper_env, k_l)
      # check_state(init_env, hyper_env, X, tag = "NA")
      #print(init_env$z)
      #print(init_env$sigma_2)
      # print(init_env$psi)
      # check <- eval_log_lik_theta_joint(init_env$thetas, init_env, hyper_env, X)
      # cat("hyper$log_lik =", hyper_env$log_lik,
      #     " recomputed at current theta =", check$log_lik - check$val_3,
      #     " diff =", check$log_lik - check$val_3 - hyper_env$log_lik, "\n")
      
      # X_mean <- init_env$mean
      # plot(X_mean, asp = 1, pch = 21, bg = "darkblue", main = "Latent mu")
      # lines(c(X_mean[,1], X_mean[1,1]), c(X_mean[,2], X_mean[1,2]), col = "blue", lwd = 2)
      # X_new_i <- init_env$mean_i[,,1]
      # plot(X_new_i, asp = 1, pch = 21, bg = "black", main = "Latent mu unit i")
      # lines(c(X_new_i[,1], X_new_i[1,1]), c(X_new_i[,2], X_new_i[1,2]), col = "brown", lwd = 2)
      # 
      # plot(X[,,1], asp = 1, pch = 21, bg = "black", main = "Observed Unit i")
      # lines(c(X[,1,1], X[1,1,1]), c(X[,2,1], X[1,2,1]), col = "brown", lwd = 2)
      #print(init_env$z)
      #z_temp[j] <- init_env$z
      if (i == 1 & j %% 100 == 0) {
        # cli_inform(paste0("\n[Burn-in - Iterazione globale: ", hyper$t, "]"))
        #print(hyper_env$a_aver_a/hyper_env$t)
        # print(init$Sigma)
        # print(init_env$tau)
        #print(init_env$betas)
        # print(hyper_env$Eigen_vector_null%*%init_env$betas)
        # print(init_env$Sigma)
        # print(init_env$alphas)
        # print(init_env$thetas)
        # print(init_env$C)
        # #print(init_env$thetas_sorted)
        # #print(init_env$lambdas)
        # #print(init_env$alphas)
        # print(hyper_env$tau_acc/j)
        # print(init_env$tau)
        # print(init_env$Sigma)
        # #print(init_env$eta)
        # print(hyper_env$lambda_acc/j)
        # print(hyper_env$a_aver_a/j)
        # print(hyper_env$phi_acc/j)
        # print(init_env$phi)
        # print(paste("sigma2", init_env$sigma_2))
        # print(hyper_env$sigma_acc_cp/j)
        # print(hyper_env$sigma_acc_ncp/j)
        # print(init_env$W[,,1])
  
        
        tryCatch(plot_status(init_env, hyper_env, X, plot_file, hyper_env$t, i_unit = 1),
                 error = function(e) message("plot_status fallito: ", conditionMessage(e)))
        tryCatch(plot_status(init_env, hyper_env, X, plot_file, hyper_env$t, i_unit = 10),
                 error = function(e) message("plot_status fallito: ", conditionMessage(e)))
        # print(hyper_env$S_mat %*% hyper_env$D_p %*% solve(t(hyper_env$D_p)%*%hyper_env$S_mat %*% hyper_env$D_p) %*% 
        #                  t(hyper_env$D_p) %*% t(hyper_env$S_mat ))
      }
    
      
      
    }
    if(i == 1){ # for the first iteration we are doing burnin, then we are doing thinning
    thin <- thinning
    burn <- FALSE
    }
    
    if (i %% 50 == 0) {
      tryCatch(plot_status(init_env, hyper_env, X, plot_file, hyper_env$t),
               error = function(e) message("plot_status fallito: ", conditionMessage(e)))
    }
    output$tau[i,] <- init_env$tau
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
    # output$psi[i,,] <- init_env$psi
    output$psi[[i]] <- init_env$psi
    output$W[[i]] <- init_env$W
    output$sigma_2[i, ] <- init_env$sigma_2
    # output$z[i,] <- init_env$z
    #setTxtProgressBar(pb, value = i, title = "Ziocan")
    #
    
  }
  cli_progress_done(id = pbar)
  out <- list()
  out$output <- output
  out$init <- init_env
  return(out)
  
}
# check_state <- function(init, hyper, X, tag){
#   chk <- eval_log_lik_beta(init$gammas, init, hyper, X, init$k_l)
#   d <- c(lik    = abs((chk$log_lik - chk$val_3) - log_density(init)),
#          prior  = abs(chk$val_3 - hyper$traces),
#          mean_i = max(abs(chk$mean_i - init$mean_i)),
#          W      = max(abs(chk$W - init$W)))
#   if (any(d > 1e-6)) message(sprintf("[%s] ", tag), paste(names(d), signif(d, 3), collapse = " | "))
#   invisible(d)
# }

# hyper <- hyperparameters(2, 0.5, 8, 3, width_theta = 1, m = 6, nu = 10,
#                          psi = psi,n = dim(X)[3], a = 0.01, b = 0.01, a_phi = a_phi, b_phi = b_phi, X = X)
X <- sam$X/100
X <- mat_not_demented
X <- mat_demented/100
X <- final_array/100
X <- final_array

prova <- init_mu_star (X,n_basis = 14, degree = 3, start = "gpa", ridge   = 1e-1, plot = TRUE,
                       Sigma_eta = diag(10, 2))
init <- prova$init
hyper <- prova$hyper

X <- prova$X
init$sigma_2
hyper$a_sigma
init$C
hyper$nu
init$Sigma
init$sigma_2
init$C
hyper$Sigma_a
X <- prova$X
hyper$b_phi
hyper$Sigma_eta_new
hyper$a_sigma
hyper$degree
hyper$a_phi
# save(sam, X, hyper, init, file = "No_warp.RData")
# save(sam, X, hyper, init, file = "warp_indep.RData")
load("warp.RData")
save(sam, X, hyper, init, file = "warp.RData")
# save(X, hyper, init, file = "No_Demente.RData")
save(X, hyper, init, file = "Demente.RData")
save(X, hyper, init, file = "foglie.RData")
save(X, hyper, init, file = "pesci.RData")
# load("No_Demente.RData")
# load("No_warp.RData")
load("foglie.RData")
hyper$log_lik
hyper$K
hyper$L
init$psi
init$W
init$sigma_2
init$chi
# init$alphas <- init$alphas/init$alphas[1]
# init$lambdas <- ((init$lambdas - init$lambdas[1])%% (2*pi)) 
# init$eta <- init$eta - matrix(init$eta[1,], nrow = n, ncol = 2, byrow = TRUE)
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
hyper$a_sigma
init$Z
tic()
par(mfrow = c(1, 3))
library(profvis)
profvis(expr = {
  MCMC_samp <- RB_MCMC(total_iter = 2e5, burnin = 50000, thinning  = 10, X = X, 
                     init = init, hyper = hyper)
})
  toc()
hyper$b_phi
init$mean
mean(MCMC_samp$output$z)
plot(MCMC_samp$output$z)
gcinfo(FALSE)
total_iter = 100000
burnin = 30000
thinning  = 20
X <- sam$X
thetas <- MCMC_samp$output$theta[1,]
k_l <- length(thetas)
mat <- matrix(NA, nrow = k_l, ncol = k_l)
d <- covariance_mat(mat, thetas)
sam$theta

sum(d >= 1)/2
18*18/2
sample_posterior_predictive_identita <- function(output_mcmc, subject_idx) {
  
  n_saved_iter <- dim(output_mcmc$mean)[1]
  k_l <- dim(output_mcmc$mean)[2] # Numero di landmark (K)
  p <- 2                           # Dimensione spaziale (x, y)
  
  X_pred <- array(0, dim = c(k_l, p, n_saved_iter))
  mat_dist <- matrix(0, nrow = k_l, ncol = k_l)
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
    mat_dist <- covariance_mat(mat_dist, output_mcmc$theta[i,])
    C <- exp(-mat_dist/output_mcmc$phi[i])
    # 5. Campiona la matrice di rumore bianco standard Z (K x 2)
    Z <- matrix(rnorm(k_l * p), nrow = k_l, ncol = p)
    
    # 6. Trasforma il rumore solo a destra (poiché L_C è l'identità)
    # E_i = Z %*% t(L_Sigma_spatial)
    
    E_i <- chol(C) %*% Z %*% t(L_Sigma_spatial)
    
    # Configurazione finale predetta per l'iterazione i
    X_pred[, , i] <- M_subject + E_i
  }
  
  return(X_pred)
}

MCMC_samp$output$t
posterior_pred <- sample_posterior_predictive_identita(MCMC_samp$output, 1)

dim(posterior_pred)
length(MCMC_samp$output$mean_i[[1]])
k_l = 18
mat_dist_true <- matrix(0, nrow = k_l, ncol = k_l)
mat_dist_true <- covariance_mat(mat_dist_true, sam$theta)
n <- 100
B_sim <- Basis_Construction(sam$theta, 10, degree) # build the basis with the current value of theta
r <- exp(B_sim %*% as.vector(sam$betas))# compute the current value of the radius
mu_mean_x <- r*cos(sam$theta)
mu_mean_y <- r*sin(sam$theta)
mean <- cbind(mu_mean_x, mu_mean_y) # average configuration 
R <- array(0, dim = c(2, 2, n))
mu_i_true <- array(0, dim = c(k_l, 2, n))
Q_R <- R
X <- sam$X
X <- X_new
for (i in 1:n) { # compute rotation matrix 
  R[,,i] <- matrix(c(cos(sam$lambda[i]), -sin(sam$lambda[i]), 
                                sin(sam$lambda[i]), cos(sam$lambda[i])), 
                              byrow = T, nrow = 2, ncol = 2)
  eta_matrix <- matrix(sam$eta[i,], nrow = length(sam$theta), ncol = 2, byrow = T)
  mu_i_true[,,i] <- sam$alphas[i]*mean%*%R[,,i] + eta_matrix # compute the mean configuration for every unit
  Q_R[,,i] <- (1/(sam$alphas[i]^2)) *t(R[,,i])%*%solve(sam$Sigma_e)%*%R[,,i] 

}

par(mfrow = c(1,1))
phi_grid <- seq(0.01, 2, length.out = 1000)
loglik_grid <- sapply(phi_grid, function(p) {
  C_p <- exp(-mat_dist_true / p)
  chol_p <- chol(C_p)
  val <- 0
  for (i in 1:n) {
    res <- backsolve(chol_p, X[,,i] - mu_i_true[,,i], transpose = TRUE)
    val <- val + sum(Q_R[,,i] * (t(res) %*% res))
  }
  -0.5*val - n*log(det(C_p)) -2*k_l*sum(log(sam$alphas)) -((n*k_l)/2)*log(det(sam$Sigma_e))
})
plot(phi_grid, loglik_grid, type = "l")
abline(v = sam$phi, col = "red")




save(MCMC_samp,X, file = "Fish_31_August.RData")
load("Torino.RData")

par(mfrow = c(1, 1))
plot((MCMC_samp$output$lambdas[,5] - MCMC_samp$output$lambdas[,1])%%(2*pi), type ="l")
library(LaplacesDemon)
ESS(MCMC_samp$output$lambdas[,3] - MCMC_samp$output$lambdas[,1])
plot((MCMC_samp$output$theta[,1] - MCMC_samp$output$lambdas[,10]) %%(2*pi), type = "l")
ESS((MCMC_samp$output$theta[,1] - MCMC_samp$output$theta[,12]) %%(2*pi))
plot((MCMC_samp$output$theta[,1] - MCMC_samp$output$theta[,3]) %%(2*pi), type ="l")
plot(MCMC_samp$output$phi, type = "l")
ESS(MCMC_samp$output$phi)
plot(density(MCMC_samp$output$phi))

ESS(MCMC_samp$output$theta[,5] - MCMC_samp$output$theta[,2])
sam$theta[3] - sam$theta[1]
ESS(MCMC_samp$output$alphas[,2]/MCMC_samp$output$alphas[,1])
sam$lambda[2] - sam$lambda[1]
sam$alphas[2]/sam$alphas[1]

plot(MCMC_samp$output$betas[,13] - exp(MCMC_samp$output$alphas[,1]))
plot(MCMC_samp$output$gammas[,10], type = "l")
plot(MCMC_samp$output$theta[,3] - MCMC_samp$output$theta[,2], type = "l")
ESS(MCMC_samp$output$theta[,3] - MCMC_samp$output$theta[,2])
gg_mcmc_diagnostics(MCMC_samp$output$tau[,1], param_name = "tau", real_values = NA)
gg_mcmc_diagnostics(MCMC_samp$output$tau[,2], param_name = "tau", real_values = NA)
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

pdf("24Settembre.pdf")
par(mfrow = c(1, 3))
dim(X)
# n = 120
# n <- 38
n_iter <- dim(MCMC_samp$output$mean)[1]
#i = 1
k_l <- dim(MCMC_samp$output$mean)[2]
par(mfrow = c(1, 2))
n = dim(sam$X)[3]
#i = 1
for (i in 1:n) {
  
  X_i <- sam$mu
  
  X_mean <- MCMC_samp$init$mean

  mean_samp <- MCMC_samp$output$mean  # [iter, k, 2]

  mean_low <- apply(mean_samp, c(2,3), quantile, probs = 0.025)
  mean_high <- apply(mean_samp, c(2,3), quantile, probs = 0.975)
  mean_samp_mean <- apply(mean_samp, c(2,3), mean)
  # 
  plot(mean_samp_mean, asp = 1, pch = 21, bg = "darkblue",
       main = "Estimated Mean + CI", xlim = range(mean_samp_mean[,1]) + c(-0.1, 0.1),
       ylim = range(mean_samp_mean[,2]) + c(-0.1, 0.1))

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
  # 
  plot(X_i, asp = 1, pch = 21, bg = "darkred",
       main = "i-th unit", xlim = range(X_i[,1]) + c(-0.1, 0.1),
       ylim = range(X_i[,2]) + c(-0.1, 0.1) )

  lines(c(X_i[,1], X_i[1,1]),
        c(X_i[,2], X_i[1,2]),
        col = "red", lwd = 2)
  # 
  # samp_i <- array(NA, dim = c(n_iter, k_l, 2))
  # for (j in 1:n_iter) {
  #   samp_i[j,,] <- MCMC_samp$output$mean_i[[j]][,,i]
  # 
  # }
  # 
  # 
  # #samp_i <- MCMC_samp$output$mean_i[[i]]
  # mean_samp_i <- apply(samp_i, c(2,3), mean)
  # low_i <- apply(samp_i, c(2,3), quantile, 0.025)
  # high_i <- apply(samp_i, c(2,3), quantile, 0.975)
  # plot(mean_samp_i, asp = 1, pch = 21, bg = "black",
  #      main = "Estimated config + CI",  xlim = range(mean_samp_i[,1]) + c(-0.1, 0.1),
  #      ylim = range(mean_samp_i[,2]) + c(-0.1, 0.1))
  # 
  # lines(c(mean_samp_i[,1], mean_samp_i[1,1]),
  #       c(mean_samp_i[,2], mean_samp_i[1,2]),
  #       col = "brown", lwd = 2)
  # for (k in 1:k_l) {
  # 
  #   polygon(
  #     x = c(low_i[k,1], high_i[k,1],
  #           high_i[k,1], low_i[k,1]),
  #     y = c(low_i[k,2], low_i[k,2],
  #           high_i[k,2], high_i[k,2]),
  #     col = rgb(0,0,1,0.40),
  #     border = NA
  #   )
  # 
  # }
  posterior_pred <- sample_posterior_predictive_identita(MCMC_samp$output, i)
  post_samp_i <- apply(posterior_pred, c(1,2), mean)
  low_post_i <- apply(posterior_pred, c(1,2), quantile, 0.025)
  high_post_i <- apply(posterior_pred, c(1,2), quantile, 0.975)
  plot(post_samp_i, asp = 1, pch = 21, bg = "black",
        main = "Posterior Predictive", xlim = range(post_samp_i[,1]) + c(-0.1, 0.1), ylim = range(post_samp_i[,2]) + c(-0.1, 0.1))
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
gg_mcmc_diagnostics(MCMC_samp$output$phi[,1], param_name = "phi", real_values = NA)
gg_mcmc_diagnostics(MCMC_samp$output$phi[,1], param_name = "phi", real_values = NA)
gg_mcmc_diagnostics(MCMC_samp$output$tau[,2], param_name = "tau", real_values = NA)
gg_mcmc_diagnostics(MCMC_samp$output$theta, param_name = "theta", real_values = NA, TRUE)
gg_mcmc_diagnostics(MCMC_samp$output$lambda, param_name = "lambda", real_values = NA, TRUE)
Sigma_samp <- cbind(MCMC_samp$output$Sigma[,1,1], MCMC_samp$output$Sigma[,1,2], MCMC_samp$output$Sigma[,1,2], MCMC_samp$output$Sigma[,2,2])
gg_mcmc_diagnostics(Sigma_samp, param_name = "Sigma", real_values = c(NA), TRUE)
gg_mcmc_diagnostics(MCMC_samp$output$alphas, param_name = "alpha", real_values = NA, TRUE)
gg_mcmc_diagnostics(MCMC_samp$output$betas, param_name = "betas", real_values = NA, TRUE)

save(MCMC_samp,sam, X,  file = "Torino.RData")
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

X <- sam$X

output <- MCMC_samp$output

X0 <- X[, , 2]
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
par(mfrow = c(1, 1))
n <- 100
for (i in 1:n) {
  
  X_i <- X[,,i]
  
  # X_mean <- MCMC_samp$init$mean
  # 
  # mean_samp <- MCMC_samp$output$mean  # [iter, k, 2]
  # 
  # mean_low <- apply(mean_samp, c(2,3), quantile, probs = 0.025)
  # mean_high <- apply(mean_samp, c(2,3), quantile, probs = 0.975)
  # mean_samp_mean <- apply(mean_samp, c(2,3), mean)
  # 
  # plot(mean_samp_mean, asp = 1, pch = 21, bg = "darkblue",
  #      main = "Estimated Mean + CI")
  # 
  # lines(c(mean_samp_mean[,1], mean_samp_mean[1,1]),
  #       c(mean_samp_mean[,2], mean_samp_mean[1,2]),
  #       col = "blue", lwd = 2)
  # 
  # for (k in 1:nrow(X_mean)) {
  #   
  #   polygon(
  #     x = c(mean_low[k,1], mean_high[k,1],
  #           mean_high[k,1], mean_low[k,1]),
  #     y = c(mean_low[k,2], mean_low[k,2],
  #           mean_high[k,2], mean_high[k,2]),
  #     border = NA,
  #     col = rgb(0, 0, 1, 0.15)
  #   )
  # }
  # 
  # plot(X_i, asp = 1, pch = 21, bg = "darkred",
  #      main = "i-th unit")
  # 
  # lines(c(X_i[,1], X_i[1,1]),
  #       c(X_i[,2], X_i[1,2]),
  #       col = "red", lwd = 2)
  # 
  # samp_i <- array(NA, dim = c(n_iter, k_l, 2))
  # for (j in 1:n_iter) {
  #   samp_i[j,,] <- MCMC_samp$output$mean_i[[j]][,,i]
  #   
  # }
  # 
  # 
  # #samp_i <- MCMC_samp$output$mean_i[[i]]
  # mean_samp_i <- apply(samp_i, c(2,3), mean)
  # low_i <- apply(samp_i, c(2,3), quantile, 0.025)
  # high_i <- apply(samp_i, c(2,3), quantile, 0.975)
  # plot(mean_samp_i, asp = 1, pch = 21, bg = "black",
  #      main = "Estimated config + CI")
  # 
  # lines(c(mean_samp_i[,1], mean_samp_i[1,1]),
  #       c(mean_samp_i[,2], mean_samp_i[1,2]),
  #       col = "brown", lwd = 2)
  # for (k in 1:k_l) {
  #   
  #   polygon(
  #     x = c(low_i[k,1], high_i[k,1],
  #           high_i[k,1], low_i[k,1]),
  #     y = c(low_i[k,2], low_i[k,2],
  #           high_i[k,2], high_i[k,2]),
  #     col = rgb(0,0,1,0.40),
  #     border = NA
  #   )
  #   
  # }
  # 
  post_samp_i <- apply(posterior_pred, c(1,2), mean)
  low_post_i <- apply(posterior_pred, c(1,2), quantile, 0.025)
  high_post_i <- apply(posterior_pred, c(1,2), quantile, 0.975)
  plot(post_samp_i, asp = 1, pch = 21, bg = "black",
       main = "Posterior Predictive")
  lines(c(post_samp_i[,1], post_samp_i[1,1]),
        c(post_samp_i[,2], post_samp_i[1,2]),
        col = "brown", lwd = 2)
}
pdf("Risulati24.pdf")
dev.off()
output <- MCMC_samp$output
gg_mcmc_diagnostics(as.vector(MCMC_samp$output$tau[,1]), param_name = "tau1", real_values = sam$tau[1])
gg_mcmc_diagnostics(as.vector(MCMC_samp$output$tau[,2]), param_name = "tau2", real_values = sam$tau[2])
gg_mcmc_diagnostics(as.vector(MCMC_samp$output$sigma_2), param_name = "phi", real_values = NA)
for (i in 1:14) {
  
  for (j in 1:2) {
    print(gg_mcmc_diagnostics(as.vector(MCMC_samp$output$psi[,i,j]), param_name = paste(paste("psi", sep = "_", i), sep = "", j), real_values = NA))
    
  }
  
}
dev.off()
gg_mcmc_diagnostics(as.vector(MCMC_samp$output$sigma_2), param_name = "sigma2", real_values = 0)

nsamp <- dim(output$lambdas)[1]
mean_orig <- array(NA, c(k, 2, nsamp))
mean_noloc <- array(NA, c(k, 2, nsamp))
mean_norot <- array(NA, c(k, 2, nsamp))
mean_nosize <- array(NA, c(k, 2, nsamp))

sigma_orig <- array(NA, c(2, 2, nsamp))
sigma_noloc <- array(NA, c(2, 2, nsamp))
sigma_norot <- array(NA, c(2, 2, nsamp))
sigma_nosize <- array(NA, c(2, 2, nsamp))
# pair <- c(1, 10)
for (iobs in 1:n) {
  
  # 1. Feedback reale sulla console R
  cat(sprintf("\rElaborazione Unità: %3d / 120 ... ", iobs))
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
      print(ESS(mean_nosize[iii, ip, ]))
      
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


apply(MCMC_samp$output$psi, c(2,3), mean)
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

