source("src/Packages.R")
source("src/Utils.R")
source("src/sample_tau.R")
source("src/sample_thetas.R")
source("src/sample_lambdas.R")
source("src/sample_Sigma.R")
source("src/sample_alphas.R")

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
  mcmc_iter <- floor((total_iter - burnin)/thinning) # how many iterations do we need to save?
  total_steps <- burnin + (mcmc_iter - 1) * thinning
  #pbar <- cli_progress_bar(total = total_steps)
  k_l <- dim(X)[1] # number of landmarks
  n <- init$n
  output <- create_output(mcmc_iter, k_l, init$n) # create the output list
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
  for (i in 1:mcmc_iter) {
    
    for (j in 1:thin) {
      cli_progress_update(id = pbar, inc = 1)
      
      # Tau update ----
      hyper_env$a_new <- hyper_env$L/2 + hyper_env$a_tau # shape parameter
      hyper_env$b_new <- (t(init_env$betas)%*%hyper_env$K%*%init_env$betas)/2+ hyper_env$b_tau # scale parameter 
      sample_tau(init_env, a_new = hyper_env$a_new, b_new = hyper_env$b_new) # sample new tau^2
      # alpha update ------
      for (y in 1:n) {
        hyper_env$log_lik_a[y] <- eval_log_lik_alpha(log(init_env$alphas[y]), hyper_env, init_env, k_l, X, y)$log_lik
        
      }
      sample_alpha(init_env, hyper_env, X, k_l, burn)
      # init <- app_alpha$init
      # hyper <- app_alpha$hyper
      
      
      # Sigma update ---- 
      sample_Sigma(init_env, hyper_env, k_l)
      # init <- app_sigma$init
      # hyper <- app_sigma$hyper
      # print(init$thetas)
      # theta update ----
      hyper_env$log_lik <- log_density(init_env)
      sample_theta(init_env, hyper_env, X)
      # init <- app_theta$init
      # hyper <- app_theta$hyper
      # lambda update -----
      sample_lambda(init_env, hyper_env, X, k_l)
      # init <- app_lambda$init
      # hyper <- app_lambda$hyper
      
      
      if (i == 1 & j %% 500 == 0) {
        # cli_inform(paste0("\n[Burn-in - Iterazione globale: ", hyper$t, "]"))
        # print(hyper$a_aver_a/hyper$t)
        # print(init$Sigma)
        print(init_env$thetas)
        
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
    #setTxtProgressBar(pb, value = i, title = "Ziocan")
    #
    
  }
  cli_progress_done(id = pbar)
  out <- list()
  out$output <- output
  out$init <- init_env
  return(out)
  
}




init <- init_param(0.01, runif(n = 100)*2*pi, sort(runif(n = 25)*2*pi), sam$betas, diag(2, nrow = 2),
                   sam$eta,rgamma(n = 100, 1, 1), n <- dim(sam$X)[3])
# init <- init_param(0.01, sam$lambda, sam$theta, sam$betas, diag(1, nrow = 2),
#                    sam$eta,sam$alphas, n <- dim(sam$X)[3])
hyper <- hyperparameters(0.01, 0.01, 10, 3, width_theta = 1, m = 6, nu = 4,
                         psi = diag(0.01,2),n <- dim(sam$X)[3], a = 1, b = 1)


prova <- init
hyper$mu_a
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

apply(Sigma_samp, MARGIN = 2, quantile, probs = c(0.025, 0.975))
sam$Sigma_e
Sigma_samp <- cbind(MCMC_samp$output$Sigma[,1,1], MCMC_samp$output$Sigma[,1,2], MCMC_samp$output$Sigma[,1,2], MCMC_samp$output$Sigma[,2,2])
dim(MCMC_samp$output$Sigma)

X_new_i <- MCMC_samp$init$mean_i[,,3]
par(mfrow = c(1, 4))
plot(rot, asp = 1, pch = 21, bg = "darkgreen", main = "Latent mu")
lines(c(rot[,1], rot[1,1]), c(rot[,2], rot[1,2]), col = "forestgreen", lwd = 2)

plot(mu_i, asp = 1, pch = 21, bg = "darkred", main = "mu for the i-th unit")
lines(c(mu_i[,1], mu_i[1,1]), c(mu_i[,2], mu_i[1,2]), col = "red", lwd = 2)

plot(X_i, asp = 1, pch = 21, bg = "darkblue", main = "Observed Conf for unit i")
lines(c(X_i[,1], X_i[1,1]), c(X_i[,2], X_i[1,2]), col = "blue", lwd = 2)

plot(X_new_i, asp = 1, pch = 21, bg = "black", main = "estimated Conf for unit i")
lines(c(X_new_i[,1], X_new_i[1,1]), c(X_new_i[,2], X_new_i[1,2]), col = "brown", lwd = 2)
