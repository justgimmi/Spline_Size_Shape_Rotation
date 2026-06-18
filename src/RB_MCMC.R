source("src/Packages.R")
source("src/Utils.R")
source("src/sample_tau.R")
source("src/sample_thetas.R")
source("src/sample_lambdas.R")
source("src/sample_Sigma.R")
source("src/sample_alphas.R")
source("src/sample_betas.R")
source("src/sample_eta.R")


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

  burn  = TRUE
  hyper_env$B_sim <- Basis_Construction(init_env$thetas, hyper_env$n_basis, hyper_env$degree)
  hyper_env$log_lik <- log_density(init_env)
  for (i in 1:mcmc_iter) {
    
    for (j in 1:thin) {
      #print(i)
      #print(j)
      cli_progress_update(id = pbar, inc = 1)
      
      sample_lambda(init_env, hyper_env, X, k_l)
      
      sample_theta(init_env, hyper_env, X)
      
      # init = init_env
      # hyper =  hyper_env
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
      
      for (indice in 1:n) {
        sample_eta(init_env, hyper_env, indice, X) 
      }
      hyper_env$log_lik <- log_density(init_env)
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
        print(init_env$eta)
        # X_mean <- init_env$mean
        # plot(X_mean, asp = 1, pch = 21, bg = "darkblue", main = "Estimated Mean")
        # lines(c(X_mean[,1], X_mean[1,1]), c(X_mean[,2], X_mean[1,2]), col = "blue", lwd = 2)
        # X_new_i <- init_env$mean_i[,,1]
        # plot(X_new_i, asp = 1, pch = 21, bg = "black", main = "estimated Conf for unit i")
        # lines(c(X_new_i[,1], X_new_i[1,1]), c(X_new_i[,2], X_new_i[1,2]), col = "brown", lwd = 2)
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
    output$eta[i,,] <- init_env$eta
    output$mean_i[[i]] <- init_env$mean_i
    output$mean[i, ,] <- init_env$mean
    if (i %% 100 == 0) {
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
      X_mean <- init_env$mean
      plot(X_mean, asp = 1, pch = 21, bg = "darkblue", main = "Estimated Mean")
      lines(c(X_mean[,1], X_mean[1,1]), c(X_mean[,2], X_mean[1,2]), col = "blue", lwd = 2)
    }
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

X <- sam/100
hyper <- hyperparameters(0.01, 0.01, 10, 3, width_theta = 0.1, m = 60, nu = 4,
                         psi = diag(0.01,2),n <- dim(X)[3], a = 0.01, b = 0.01)
X_centered <- array(NA, dim = dim(X))

for (i in 1:dim(X)[3]) {
  X_centered[,,i] <- X[,,i] - matrix(colMeans(X[,,i]), nrow = dim(X)[1], ncol = 2, byrow = TRUE)
}

gpa <- procGPA(X_centered, scale = FALSE)
#?procGPA
mean_shape <- gpa$mshape
theta_init <- atan2(mean_shape[,2], mean_shape[,1])
theta_init <- ifelse(theta_init < 0, theta_init + 2*pi, theta_init)
ord <- order(theta_init)
theta_init <- theta_init[ord]
mean_shape_ordered <- mean_shape[ord, ]


X_new <- array(NA, dim = dim(X))
for(i in 1:dim(X)[3]) {
  X_new[,,i] <- X[ord, , i] # Ordiniamo i landmark per il sarago i
}

X <- X_new
#theta_init <- sort(theta_init)

r_emp <- sqrt(mean_shape_ordered[,1]^2 + mean_shape_ordered[,2]^2)
log_r_emp <- log(r_emp)
B_init <- Basis_Construction(theta_init, 10, 3)
beta_init <- solve(t(B_init) %*% B_init + 0.01 * diag(ncol(B_init)), t(B_init) %*% log_r_emp)

gamma_init_ridotto <- as.vector(hyper$A_bar %*% beta_init)
beta_init <- as.vector(hyper$C_bar %*% gamma_init_ridotto)

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
  R_i <- svd_decomp$u %*% t(svd_decomp$v)
  
  lambda_sync[i] <- atan2(R_i[2,1], R_i[1,1])
  alpha_sync[i]  <- sum((X_i_centered %*% R_i) * mu_pure_centered) / sum(mu_pure_centered^2)
}
lambda_sync <- ifelse(lambda_sync < 0, lambda_sync + 2*pi, lambda_sync)


#lambda_init <- rep(0, 120)
init <- init_param(0.01,lambda_sync, theta_init,beta_init,
                   diag(0.001, 2),
                   eta_sync,alpha_sync, n = dim(X)[3])
hyper$L
init$lambdas


# init <- init_param(0.01, runif(n = 100)*2*pi, sort(runif(n = 25)*2*pi),rmvnorm_rd(n = 1, mu = 0, Precision =  (1/0.01)*hyper$K, tol = 1e-6),
#                    diag(0.01, nrow = 2),
#                    sam$eta,rgamma(n = 100, 1, 1), n <- dim(sam$X)[3])
# init <- init_param(0.01, runif(n = 100)*2*pi, sam$theta,sam$betas, diag(0.01, nrow = 2),
#                    sam$eta,rgamma(n = 100, 1, 1), n <- dim(sam$X)[3])

init$betas
hyper$n_basis

hyper$Eigen_vector_null%*%sam$betas
sam$betas
sam$tau^2
sam$Sigma_e
tic()
MCMC_samp <- RB_MCMC(total_iter = 100000, burnin = 20000, thinning  = 20, X = X, 
                     init = init, hyper = hyper)
toc()
gcinfo(FALSE)
total_iter = 100000
burnin = 30000
thinning  = 20
X <- sam$X
init$Q_R
sam$theta
sam$lambda
sam$betas
sam$theta
20*0.31
2*pi
save(MCMC_samp, file = "risultati_new.RData")
load("risultati_new.RData")
gg_mcmc_diagnostics(MCMC_samp$output$tau, param_name = "tau", real_values = NA)
gg_mcmc_diagnostics(MCMC_samp$output$theta, param_name = "theta", real_values = NA, TRUE)
gg_mcmc_diagnostics(MCMC_samp$output$lambda, param_name = "lambda", real_values = NA, TRUE)
Sigma_samp <- cbind(MCMC_samp$output$Sigma[,1,1], MCMC_samp$output$Sigma[,1,2], MCMC_samp$output$Sigma[,1,2], MCMC_samp$output$Sigma[,2,2])
gg_mcmc_diagnostics(Sigma_samp, param_name = "Sigma", real_values = c(NA), TRUE)
gg_mcmc_diagnostics(MCMC_samp$output$alphas, param_name = "alpha", real_values = NA, TRUE)
gg_mcmc_diagnostics(MCMC_samp$output$betas, param_name = "betas", real_values = NA, TRUE)
gg_mcmc_diagnostics(MCMC_samp$output$eta[,,1], param_name = "etas", real_values = NA, TRUE)
par(mfrow = c(1,1))
acf(MCMC_samp$output$theta[seq(1, 4000, 50),1])

apply(Sigma_samp, MARGIN = 2, quantile, probs = c(0.025, 0.975))
sam$Sigma_e
Sigma_samp <- cbind(MCMC_samp$output$Sigma[,1,1], MCMC_samp$output$Sigma[,1,2], MCMC_samp$output$Sigma[,1,2], MCMC_samp$output$Sigma[,2,2])
dim(MCMC_samp$output$Sigma)

par(mfrow = c(1, 4))
MCMC_samp$init$mean_i[]
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

pdf("Risultati_Not_Ident.pdf")
par(mfrow = c(1, 3))
n = 120
n <- 120
n_iter <- dim(MCMC_samp$output$mean)[1]
#i = 1
k_l <- dim(MCMC_samp$output$mean)[2]
for (i in 1:n) {
  
  X_i <- X[,,i]
  
  X_mean <- MCMC_samp$init$mean
  
  # ==========================================================
  # CI per la mean shape (se hai posterior samples)
  # ==========================================================
  
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
}


gg_mcmc_diagnostics(MCMC_samp$output$tau, param_name = "tau", real_values = NA)
gg_mcmc_diagnostics(MCMC_samp$output$theta, param_name = "theta", real_values = NA, TRUE)
gg_mcmc_diagnostics(MCMC_samp$output$lambda, param_name = "lambda", real_values = NA, TRUE)
Sigma_samp <- cbind(MCMC_samp$output$Sigma[,1,1], MCMC_samp$output$Sigma[,1,2], MCMC_samp$output$Sigma[,1,2], MCMC_samp$output$Sigma[,2,2])
gg_mcmc_diagnostics(Sigma_samp, param_name = "Sigma", real_values = c(NA), TRUE)
gg_mcmc_diagnostics(MCMC_samp$output$alphas, param_name = "alpha", real_values = NA, TRUE)
gg_mcmc_diagnostics(MCMC_samp$output$betas, param_name = "betas", real_values = NA, TRUE)

#save(MCMC_samp,sam , file = "prova.RData")

# for (i in 1:n) {
#   
#   #rot <- sam$mu
#   #mu_i <- sam$mu_i[,,i]
#   X_i <- X[,,i]
#   X_new_i <- MCMC_samp$init$mean_i[,,i]
#   X_mean <- MCMC_samp$init$mean
#   
#   # plot(rot, asp = 1, pch = 21, bg = "darkgreen", main = "Latent mu")
#   # lines(c(rot[,1], rot[1,1]), c(rot[,2], rot[1,2]), col = "forestgreen", lwd = 2)
# 
#   plot(X_mean, asp = 1, pch = 21, bg = "darkblue", main = "Estimated Mean")
#   lines(c(X_mean[,1], X_mean[1,1]), c(X_mean[,2], X_mean[1,2]), col = "blue", lwd = 2)
#   
#   plot(X_i, asp = 1, pch = 21, bg = "darkred", main = "i-th unit")
#   lines(c(X_i[,1], X_i[1,1]), c(X_i[,2], X_i[1,2]), col = "red", lwd = 2)
#   
#   # plot(mu_i, asp = 1, pch = 21, bg = "darkred", main = "mu for the i-th unit")
#   # lines(c(mu_i[,1], mu_i[1,1]), c(mu_i[,2], mu_i[1,2]), col = "red", lwd = 2)
#   
#   plot(X_new_i, asp = 1, pch = 21, bg = "black", main = "estimated Conf for unit i")
#   lines(c(X_new_i[,1], X_new_i[1,1]), c(X_new_i[,2], X_new_i[1,2]), col = "brown", lwd = 2)
#   
#   
# }

dev.off()


cor(thetas[!is.na(thetas)], betas[!is.na(betas)])
plot(thetas, betas)

for (i in 1:init_env$k_l) {
  
  thetas <- output$theta[,i]
  betas <- output$betas[,i]
  plot(thetas, betas)
  
  
}



MCMC_samp$output$mean[1, ,]


# Identifiability ------ 
# Funzione per allineare una forma X alla target
# ==============================================================================
# FUNZIONE DI ALLINEAMENTO DI PROCRUSTE (SOLO POST-PROCESSING)
# ==============================================================================
output <- MCMC_samp$output
# ============================================================
# 1. scegli pair (massima distanza sulla forma iniziale media)
# ============================================================

X0 <- X[, , 1]

pair <- c(1, 2)
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
align_prof_method <- function(output, pair, third_landmark = NULL) {
  
  n_iter <- dim(output$mean)[1]
  k <- dim(output$mean)[2]
  
  mean_aligned <- array(0, dim = c(n_iter, k, 2))
  eta_aligned <- array(0, dim = dim(output$eta))
  
  alphas_aligned  <- matrix(0, n_iter, ncol(output$alphas))
  betas_aligned   <- matrix(0, n_iter, ncol(output$betas))
  lambdas_aligned <- matrix(0, n_iter, ncol(output$lambdas))
  theta_aligned   <- matrix(0, n_iter, ncol(output$theta))
  
  Sigma_aligned <- array(0, dim(output$Sigma))
  
  for (i in 1:n_iter) {
    
    M_raw <- output$mean[i,,]
    
    # ======================================================
    # 1. TRASLAZIONE
    # ======================================================
    
    shift <- M_raw[pair[1],]
    M <- sweep(M_raw, 2, shift, "-")
    
    eta_aligned[i,,1] <- output$eta[i,,1] + shift[1]
    eta_aligned[i,,2] <- output$eta[i,,2] + shift[2]
    
    # ======================================================
    # 2. ROTAZIONE
    # ======================================================
    
    v <- M[pair[2],]
    
    phi <- atan2(v[2], v[1])
    
    R <- matrix(c(cos(phi), -sin(phi),
                  sin(phi),  cos(phi)),
                2,2,byrow=TRUE)
    
    M <- M %*% t(R)
    
    # update angular params (CHECK SIGN CONSISTENCY MODEL-WISE)
    lambdas_aligned[i,] <- output$lambdas[i,] - phi
    theta_aligned[i,]   <- (output$theta[i,] + phi) %% (2*pi)
    
    # ======================================================
    # 3. SCALA (NORMA = 1)
    # ======================================================
    
    scale <- sqrt(sum(M^2))
    M <- M / scale
    
    alphas_aligned[i,] <- output$alphas[i,] * scale
    betas_aligned[i,]  <- output$betas[i,] - log(scale)
    
    # ======================================================
    # 4. RIFLESSIONE (CRUCIALE)
    # ======================================================
    
    # default: scegli landmark più informativo se non dato
    if (is.null(third_landmark)) {
      third_landmark <- which.max(abs(M[,2]))
      if (third_landmark %in% pair) third_landmark <- setdiff(1:k, pair)[1]
    }
    
    flip <- FALSE
    
    if (M[third_landmark,2] < 0) {
      M[,2] <- -M[,2]
      flip <- TRUE
    }
    
    mean_aligned[i,,] <- M
    
    # ======================================================
    # 5. SIGMA
    # ======================================================
    
    Sigma_i <- (R %*% output$Sigma[i,,] %*% t(R)) / (scale^2)
    
    if (flip) {
      Fmat <- diag(c(1,-1))
      Sigma_i <- Fmat %*% Sigma_i %*% Fmat
    }
    
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

# ==============================================================================
# UTILIZZO
# ==============================================================================
# Finito l'MCMC, ti basta lanciare:
output_identificato <- align_prof_method(output, pair = pair)
gg_mcmc_diagnostics(MCMC_samp$output$tau, param_name = "tau", real_values = NA)
gg_mcmc_diagnostics(output_identificato$theta, param_name = "theta", real_values = NA, TRUE)
gg_mcmc_diagnostics(output_identificato$lambda, param_name = "lambda", real_values = NA, TRUE)
Sigma_samp <- cbind(output_identificato$Sigma[,1,1], output_identificato$Sigma[,1,2],output_identificato$Sigma[,1,2], output_identificato$Sigma[,2,2])
gg_mcmc_diagnostics(Sigma_samp, param_name = "Sigma", real_values = c(NA), TRUE)
gg_mcmc_diagnostics(output_identificato$alphas, param_name = "alpha", real_values = NA, TRUE)
gg_mcmc_diagnostics(output_identificato$betas, param_name = "betas", real_values =NA, TRUE)

plot(X[,,1])
points(X[19,1,1],X[19,2,1],  col = "red", pch = 16)
points(X[10,1,1],X[10,2,1],  col = "red", pch = 16)
