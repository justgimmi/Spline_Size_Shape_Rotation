source("src/Packages.R")
source("src/Utils.R")
source("src/sample_tau.R")
source("src/sample_thetas.R")
source("src/sample_lambdas.R")
source("src/sample_Sigma.R")
source("src/sample_alphas.R")
source("src/sample_betas.R")

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
  init_env$gammas <- hyper_env$A_bar %*% init_env$betas
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
  #hyper$log_theta <- log_density_theta(init) # compute the starting value for the log-density of theta
  burn  = TRUE
  hyper_env$B_sim <- Basis_Construction(init_env$thetas, hyper_env$n_basis, hyper_env$degree)
  hyper_env$log_lik <- log_density(init_env)
  for (i in 1:mcmc_iter) {
    
    for (j in 1:thin) {
      #print(j)
      cli_progress_update(id = pbar, inc = 1)
      
      sample_lambda(init_env, hyper_env, X, k_l)
      
      sample_theta(init_env, hyper_env, X)
      
      sample_beta(init_env, hyper_env, X, k_l)
      
      hyper_env$a_new <- hyper_env$L/2 + hyper_env$a_tau # shape parameter
      #hyper_env$b_new <- (t(init_env$betas)%*%hyper_env$K%*%init_env$betas)/2+ hyper_env$b_tau # scale parameter
      hyper_env$b_new <- (t(init_env$gammas)%*%init_env$gammas)/2+ hyper_env$b_tau
      sample_tau(init_env, a_new = hyper_env$a_new, b_new = hyper_env$b_new) # sample new tau^2

      for (y in 1:n) {
        hyper_env$log_lik_a[y] <- eval_log_lik_alpha(log(init_env$alphas[y]), hyper_env, init_env, k_l, X, y)$log_lik

      }
      sample_alpha(init_env, hyper_env, X, k_l, burn)
      #mean_constructor(hyper_env$n_basis, hyper_env$degree, init_env, X)
      #hyper_env$log_lik <- log_density(init_env)

      # Tau update ----
      # alpha update ------
      
      sample_Sigma(init_env, hyper_env, k_l)
      
      #sample_theta(init_env, hyper_env, X)
      # init <- app_lambda$init
      # hyper <- app_lambda$hyper
      # init <- init_env
      # hyper <- hyper_env
      # output$theta[j, ] <- init_env$thetas
      # output$betas[j, ] <- init_env$betas
      
      if (i == 1 & j %% 100 == 0) {
        # cli_inform(paste0("\n[Burn-in - Iterazione globale: ", hyper$t, "]"))
        #print(hyper_env$a_aver_a/hyper_env$t)
        # print(init$Sigma)
        # print(init_env$tau)
        print(init_env$betas)
        # print(hyper_env$Eigen_vector_null%*%init_env$betas)
        # print(init_env$Sigma)
        # print(init_env$alphas)
        print(init_env$thetas)
        print(init_env$lambdas)
        print(init_env$alphas)
        print(init_env$tau)
        print(init_env$Sigma)
      }
      
    }
    if(i == 1){ # for the first iteration we are doing burnin, then we are doing thinning
    thin <- thinning
    burn <- FALSE
    }
    output$tau[i] <- init_env$tau
    output$theta[i,] <- init_env$thetas
    output$lambdas[i,] <- init_env$lambdas
    output$Sigma[i,,] <- init_env$Sigma
    output$alphas[i, ] <- init_env$alphas
    output$betas[i, ] <- init_env$betas
    #setTxtProgressBar(pb, value = i, title = "Ziocan")
    #
    
  }
  cli_progress_done(id = pbar)
  out <- list()
  out$output <- output
  out$init <- init_env
  return(out)
  
}

library(shapes)
X <- sam$X
X_centered <- array(NA, dim = dim(X))

for (i in 1:100) {
  X_centered[,,i] <- X[,,i] - matrix(colMeans(X[,,i]), nrow = 25, ncol = 2, byrow = TRUE)
}

gpa <- procGPA(X_centered, scale = FALSE)
mean_shape <- gpa$mshape
theta_init <- atan2(mean_shape[,2], mean_shape[,1])
theta_init <- ifelse(theta_init < 0, theta_init + 2*pi, theta_init)
theta_init <- sort(theta_init)

r_emp <- sqrt(mean_shape[,1]^2 + mean_shape[,2]^2)
log_r_emp <- log(r_emp)
B_init <- Basis_Construction(theta_init, 10, 3)
beta_init <- solve(t(B_init) %*% B_init + 0.01 * diag(ncol(B_init)), t(B_init) %*% log_r_emp)

centroid_sizes <- apply(sam$X, 3, function(mat) {
  centroid <- colMeans(mat)
  sqrt(sum(sweep(mat, 2, centroid, "-")^2))
})

alpha_init <- centroid_sizes / mean(centroid_sizes)
sam$alphas

#alpha_init <- apply(raggi_individuali, 2, mean) / raggio_medio_globale


hyper <- hyperparameters(0.01, 0.01, 10, 3, width_theta = 0.5, m = 12, nu = 4,
                         psi = diag(0.01,2),n <- dim(sam$X)[3], a = 0.01, b = 0.01)

lambda_init <- rep(0, 100)
init <- init_param(0.01,lambda_init, theta_init,beta_init,
                   diag(0.0001, nrow = 2),
                   sam$eta,alpha_init, n <- dim(sam$X)[3])
alpha_init - sam$alphas
alpha_init
init <- init_param(0.01, runif(n = 100)*2*pi, sort(runif(n = 25)*2*pi),rmvnorm_rd(n = 1, mu = 0, Precision =  (1/0.01)*hyper$K, tol = 1e-6),
                   diag(0.01, nrow = 2),
                   sam$eta,rgamma(n = 100, 1, 1), n <- dim(sam$X)[3])
init <- init_param(0.01, runif(n = 100)*2*pi, sam$theta,sam$betas, diag(0.01, nrow = 2),
                   sam$eta,rgamma(n = 100, 1, 1), n <- dim(sam$X)[3])

init$betas

hyper$Eigen_vector_null%*%sam$betas
sam$betas
sam$tau^2
sam$Sigma_e
tic()
MCMC_samp <- RB_MCMC(total_iter = 100000, burnin = 10000, thinning  = 20, X <- sam$X, 
                     init = init, hyper = hyper)
toc()
gcinfo(FALSE)
total_iter = 100000
burnin = 10000
thinning  = 10
X <- sam$X
init$Q_R
sam$theta
sam$lambda

20*0.31
2*pi
gg_mcmc_diagnostics(MCMC_samp$output$tau, param_name = "tau", real_values = sam$tau^2)
gg_mcmc_diagnostics(MCMC_samp$output$theta, param_name = "theta", real_values = sam$theta, TRUE)
gg_mcmc_diagnostics(MCMC_samp$output$lambda, param_name = "lambda", real_values = sam$lambda, TRUE)
Sigma_samp <- cbind(MCMC_samp$output$Sigma[,1,1], MCMC_samp$output$Sigma[,1,2], MCMC_samp$output$Sigma[,1,2], MCMC_samp$output$Sigma[,2,2])
gg_mcmc_diagnostics(Sigma_samp, param_name = "Sigma", real_values = c(sam$Sigma_e), TRUE)
gg_mcmc_diagnostics(MCMC_samp$output$alphas, param_name = "alpha", real_values = sam$alphas, TRUE)
gg_mcmc_diagnostics(MCMC_samp$output$betas, param_name = "betas", real_values = sam$betas, TRUE)

par(mfrow = c(1,1))
acf(MCMC_samp$output$theta[seq(1, 4000, 50),1])

apply(Sigma_samp, MARGIN = 2, quantile, probs = c(0.025, 0.975))
sam$Sigma_e
Sigma_samp <- cbind(MCMC_samp$output$Sigma[,1,1], MCMC_samp$output$Sigma[,1,2], MCMC_samp$output$Sigma[,1,2], MCMC_samp$output$Sigma[,2,2])
dim(MCMC_samp$output$Sigma)

par(mfrow = c(1, 4))
for (i in 1:n) {
  
  rot <- sam$mu
  mu_i <- sam$mu_i[,,i]
  X_i <- sam$X[,,i]
  X_new_i <- init_env$mean_i[,,i]
  X_mean <- init_env$mean

  plot(rot, asp = 1, pch = 21, bg = "darkgreen", main = "Latent mu")
  lines(c(rot[,1], rot[1,1]), c(rot[,2], rot[1,2]), col = "forestgreen", lwd = 2)
  
  plot(X_mean, asp = 1, pch = 21, bg = "darkblue", main = "Estimated Mean")
  lines(c(X_mean[,1], X_mean[1,1]), c(X_mean[,2], X_mean[1,2]), col = "blue", lwd = 2)
  
  plot(mu_i, asp = 1, pch = 21, bg = "darkred", main = "mu for the i-th unit")
  lines(c(mu_i[,1], mu_i[1,1]), c(mu_i[,2], mu_i[1,2]), col = "red", lwd = 2)
  
  plot(X_new_i, asp = 1, pch = 21, bg = "black", main = "estimated Conf for unit i")
  lines(c(X_new_i[,1], X_new_i[1,1]), c(X_new_i[,2], X_new_i[1,2]), col = "brown", lwd = 2)
  
  
}

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

pdf("prova.pdf")

gg_mcmc_diagnostics(MCMC_samp$output$tau, param_name = "tau", real_values = sam$tau^2)
gg_mcmc_diagnostics(MCMC_samp$output$theta, param_name = "theta", real_values = sam$theta, TRUE)
gg_mcmc_diagnostics(MCMC_samp$output$lambda, param_name = "lambda", real_values = sam$lambda, TRUE)
Sigma_samp <- cbind(MCMC_samp$output$Sigma[,1,1], MCMC_samp$output$Sigma[,1,2], MCMC_samp$output$Sigma[,1,2], MCMC_samp$output$Sigma[,2,2])
gg_mcmc_diagnostics(Sigma_samp, param_name = "Sigma", real_values = c(sam$Sigma_e), TRUE)
gg_mcmc_diagnostics(MCMC_samp$output$alphas, param_name = "alpha", real_values = sam$alphas, TRUE)
gg_mcmc_diagnostics(MCMC_samp$output$betas, param_name = "betas", real_values = sam$betas, TRUE)

save(MCMC_samp,sam , file = "prova.RData")
par(mfrow = c(1, 4))
for (i in 1:n) {
  
  rot <- sam$mu
  mu_i <- sam$mu_i[,,i]
  X_i <- sam$X[,,i]
  X_new_i <- init_env$mean_i[,,i]
  X_mean <- init_env$mean
  
  plot(rot, asp = 1, pch = 21, bg = "darkgreen", main = "Latent mu")
  lines(c(rot[,1], rot[1,1]), c(rot[,2], rot[1,2]), col = "forestgreen", lwd = 2)
  
  plot(X_mean, asp = 1, pch = 21, bg = "darkblue", main = "Estimated Mean")
  lines(c(X_mean[,1], X_mean[1,1]), c(X_mean[,2], X_mean[1,2]), col = "blue", lwd = 2)
  
  plot(mu_i, asp = 1, pch = 21, bg = "darkred", main = "mu for the i-th unit")
  lines(c(mu_i[,1], mu_i[1,1]), c(mu_i[,2], mu_i[1,2]), col = "red", lwd = 2)
  
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




init_env$lambdas
sam$lambda
X_new_i <- MCMC_samp$init$mean_i[,,3]
mu_i <- init_env$mean
par(mfrow = c(1, 4))
plot(rot, asp = 1, pch = 21, bg = "darkgreen", main = "Latent mu")
lines(c(rot[,1], rot[1,1]), c(rot[,2], rot[1,2]), col = "forestgreen", lwd = 2)

plot(mu_i, asp = 1, pch = 21, bg = "darkred", main = "mu for the i-th unit")
lines(c(mu_i[,1], mu_i[1,1]), c(mu_i[,2], mu_i[1,2]), col = "red", lwd = 2)

plot(X_i, asp = 1, pch = 21, bg = "darkblue", main = "Observed Conf for unit i")
lines(c(X_i[,1], X_i[1,1]), c(X_i[,2], X_i[1,2]), col = "blue", lwd = 2)

plot(X_new_i, asp = 1, pch = 21, bg = "black", main = "estimated Conf for unit i")
lines(c(X_new_i[,1], X_new_i[1,1]), c(X_new_i[,2], X_new_i[1,2]), col = "brown", lwd = 2)
