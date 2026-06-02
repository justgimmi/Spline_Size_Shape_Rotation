source("src/Packages.R")
source("src/Utils.R")
source("src/sample_tau.R")
source("src/sample_thetas.R")
source("src/sample_lambdas.R")
# MCMC struct useful for the sample ------


# MCMC Sampler for all the parameters -----

RB_MCMC <- function(total_iter, burnin, thinning, 
                    X, init, hyper){
  # total_iter : --> total number of iterations for the sampler
  # burnin : --> burnin for the chain
  # thinning: --> how much thinning do we apply
  # X: --> array of size K x p x n (where K is the number of lanmdarks, p the dimension and n the number of observations)
  # init_param --> starting values for the parameters (it depends on how we decide to build the order)
  # hyper --> list of all the hyperparameters 
  
  mcmc_iter <- floor((total_iter - burnin)/thinning) # how many iterations do we need to save?
  k_l <- dim(X)[1] # number of landmarks
  output <- create_output(mcmc_iter, k_l, init$n)
  pb <- txtProgressBar(min = 1, max = mcmc_iter, style = 3, char = "*")
  thin <- burnin
  init <- mean_constructor(hyper$n_basis, hyper$degree, init, X)
  hyper$log_theta <- log_density_theta(init)
  for (i in 1:mcmc_iter) {
    
    for (j in 1:thin) {
      
      # Tau update ----
      hyper$a_new <- hyper$L/2 + hyper$a_tau
      hyper$b_new <- (t(init$betas)%*%hyper$K%*%init$betas)/2+ hyper$b_tau
      init$tau <- sample_tau(hyper$a_new, hyper$b_new)
      
      # theta update ----
      app_theta <- sample_theta(init, hyper, X)
      init <- app_theta$init
      hyper <- app_theta$hyper
      
      # lambda update -----
      app_lambda <- sample_lambda(init, hyper, X)
      init <- app_lambda$init
      hyper <- app_lambda$hyper
      # print(init$thetas)
      
    }
    if(i == 1){ # for the first iteration we are doing burnin, then we are doing thinning
    thin <- thinning
    }
    output$tau[i] <- init$tau
    output$theta[i,] <- init$thetas
    output$lambdas[i,] <- init$lambdas
    setTxtProgressBar(pb, value = i, title = "Ziocan")
    
  }
  out <- list()
  out$output <- output
  out$init <- init
  return(out)
  
}
init$Q_R
init <- init_param(0.01, runif(n = 100)*2*pi, runif(n = 25)*2*pi, sam$betas, sam$Sigma_e, sam$eta,sam$alphas, n <- dim(sam$X)[3])
hyper <- hyperparameters(0.01, 0.01, 10, 3, width_theta = 2*pi)
tic()
MCMC_samp <- RB_MCMC(total_iter = 15000, burnin = 100, thinning  = 30, X <- sam$X, 
                     init = init, hyper = hyper)
toc()
gcinfo(FALSE)
total_iter = 10000
burnin = 10
thinning  = 10
X <- sam$X
init$Q_R
sam$theta
sam$lambda

gg_mcmc_diagnostics(MCMC_samp$output$tau, param_name = "tau", real_values = sam$tau^2)
gg_mcmc_diagnostics(MCMC_samp$output$theta, param_name = "theta", real_values = sam$theta)
gg_mcmc_diagnostics(MCMC_samp$output$lambda, param_name = "lambda", real_values = sam$lambda)
