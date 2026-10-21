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
total_iter = 1000
burnin = 1
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
  if (!is.null(dev.list())) {
    dev.off() 
  }
  par(mfrow = c(1, 3))
  #init_env$z <- 1
  # z_temp <- numeric(thin)
  candidate <- list()
  candidate$B_cand <- hyper_env$B_sim
  candidate$r_cand <- init_env$r
  candidate$mean_cand <- init_env$mean
  candidate$S_mat_cand <- hyper_env$S_mat
  candidate$D_p_cand <- hyper_env$D_p
  candidate$D_P_S_cand <- candidate$S_mat_cand %*% candidate$D_p_cand
  candidate$W_cand <- init_env$W
  candidate$SW_cand <- hyper_env$SW
  candidate$app_cand <- sweep(candidate$SW_cand, MARGIN = c(1, 2), STATS = candidate$mean_cand, FUN = "+")
  candidate$val <- 0
  candidate <- as.environment(candidate)
  X_pred <- X
  for (i in 1:mcmc_iter) {
    
    for (j in 1:thin) {
      #print(j)
      hyper_env$t <- hyper_env$t + 1
      cli_progress_update(id = pbar, inc = 1)
      #print(init_env$mean)
      sample_tau(init_env, hyper_env, burn) 
      sample_beta(init_env, hyper_env,candidate,X)
      sample_theta_logit_normal(init_env, hyper_env, candidate, X)
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
      
      # for (y in 1:n) {
      #   hyper_env$log_lik_a[y] <- eval_log_lik_alpha(log(init_env$alphas[y]), hyper_env, init_env, k_l, X, y)$log_lik
      # 
      # }
      sample_alpha(init_env, hyper_env, X, k_l, burn)
      # hyper_env$a_sigma_post <- hyper_env$a_sigma + (k_l - 3)
      # hyper_env$b_sigma_post <- hyper_env$b_sigma + tr(t(init_env$W) %*% hyper_env$S_mat %*% init_env$W)/2
      # init_env$sigma_2 <- rinvgamma(n = 1, shape = hyper_env$a_sigma_post, scale = hyper_env$b_sigma_post)
      #sample_phi(init_env, hyper_env, X, burn = TRUE)
      sample_Sigma(init_env, hyper_env, k_l)

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
        print(init_env$sigma_2/hyper_env$V)
        print(init_env$sigma_2*hyper_env$V/(norm(init_env$mean, type = "F")^2))

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
    # save all the objects 
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
    output$psi[[i]] <- init_env$psi
    output$W[[i]] <- init_env$W
    output$S[[i]] <- hyper_env$S_mat
    output$SW[[i]] <- hyper_env$SW
    output$sigma_2[i, ] <- init_env$sigma_2
    output$sigma_2_epurated[i,] <- init_env$sigma_2/hyper_env$V
    app <- compute_proc_dist(init_env, X, X_pred)
    output$dis[i, ] <- app$dis
    output$pred_i[[i]] <- app$X_pred
  }
  cli_progress_done(id = pbar)
  out <- list()
  out$output <- output
  out$init <- init_env
  out$hyper <- hyper_env
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
X <- X/100
X
prova <- init_mu_star (X,n_basis = 12, degree = 3, start = "gpa", ridge   = 1e-1, plot = TRUE,
                       Sigma_eta = diag(10, 2))

prova <- init_mu_star_ordered(X,n_basis = 12, degree = 3,
                                 nu = 4, a_tau = 2, b_tau = 0.5,
                                 pen = 1e-3,            # weak smoothness penalty for the fit of log r
                                 tau0 = c(0.1, 0.1), phi0 = NULL,
                                 Sigma_eta = diag(1, 2),
                                 d0 = 0.1,             # initial individual warp sd (fraction of shape size)
                                 plot = TRUE) 

init <- prova$init
hyper <- prova$hyper

X <- prova$X
hyper$b_s2 <- 100
hyper$a_s2 <- 1
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
load("pesci.RData")
save(sam, X, hyper, init, file = "warp.RData")
# save(X, hyper, init, file = "No_Demente.RData")
save(X, hyper, init, file = "Demente.RData")
save(X, hyper, init, file = "foglie.RData")
save(X, hyper, init, file = "pesci.RData")
load("pesci_MCMC_I.RData")
# load("No_Demente.RData")
# load("No_warp.RData")
load("Warp_MCMC_I.RData")
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
  MCMC_samp <- RB_MCMC(total_iter = 2e5, burnin = 1e5, thinning  = 100, X = X, 
                     init = init, hyper = hyper)
})
  toc()
output <- MCMC_samp$output
output
dev.off()
dev.list()
options(error = NULL)
try(visualize_results(output, X, "diobestia", sam))
load("risultati_new.RData")
output <- MCMC_samp$output
#pair <- c(1, 12)
data(gorf.dat)
dati <- schizophrenia.dat
dati
plot(dati[,,1], type = "l")
plotshapes(dati)
