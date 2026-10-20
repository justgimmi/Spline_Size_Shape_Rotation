###### Theta Sampler #####
# In this Rcode, you can find the sampler for thetas:
# The basic idea is to not directly model thetas, but to impose a logit normal prior on the consecutive angular differences
# basically, we have derived an elliptical slice sampler
# the only problem of this type of sampler is that we do not know how many calls we will make to the function
# an interesting idea that can be used is to define first all the objects that do not depend on theta on a separate function
# so we just need to compute the objects that actually depend on theta 

eval_log_lik_theta_joint <- function(theta_cand, init, hyper, candidate, X) {
  
  k_l <- init$k_l # landmarks
  n <- init$n # sample size
  k_free <- hyper$k_free # free dimensions
  
  candidate$B_cand <- Basis_Construction(theta_cand, hyper$n_basis, hyper$degree) # construct the candidate basis
  candidate$B_sim <- candidate$B_cand %*% hyper$C_bar # apply the subspace representation
  candidate$r_cand <- exp(candidate$B_sim %*% as.vector(init$gammas)) # propose a new radius
  # Construct the global latent mean shape and updaye the relevant objects depending on \mu^*
  candidate$mean_cand <- cbind(candidate$r_cand*cos(theta_cand), candidate$r_cand*sin(theta_cand)) # proposed mean
  candidate$S_mat_cand <- S_mat_constructor(candidate$mean_cand) # proposed S(mu^*) 
  candidate$D_p_cand <- D_p_constructor(candidate$mean_cand, k_l) # proposed D_p
  candidate$D_P_S_cand <- candidate$S_mat_cand %*% candidate$D_p_cand # proposed S x D_p (k x k-3)
  #candidate$W_cand <- tensor(candidate$D_p_cand, init$psi, 2, 1) # create the candidate value W
  #candidate$SW_cand <- tensor(candidate$S_mat_cand, candidate$W_cand, 2, 1)# candidate value of SW
  #mean_i_array <- init$mean_i # candidate new mean i
  #residual_array <- init$residual # candidate new residual 
  # this function is used to do mean_cand + SW_i applied over every i
  # candidate$app_cand <- sweep(candidate$SW_cand, MARGIN = c(1, 2), STATS = candidate$mean_cand, FUN = "+")
  # candidate$val <- 0
  # for (i in 1:n) {
  #   hyper$mean_i_array[,,i] <- init$alphas[i] * candidate$app_cand[,,i] %*% init$R[,,i] + init$eta_array[,,i]
  #   hyper$residual_array[,,i] <- crossprod(X[,,i] - hyper$mean_i_array[,,i])
  #   candidate$val <- candidate$val + sum(init$Q_R[,,i] * hyper$residual_array[,,i] )
  # }
  candidate$val <-  units_loglik(candidate$mean_cand, candidate$D_p_cand, candidate$D_P_S_cand,
                                 init$psi, init$R, init$eta_array, X, init$Q_R,
                                 as.numeric(init$alphas))$val
  candidate$val_2 <- -2 * k_l * sum(log(init$alphas)) - ((n * k_l) / 2) * log(det(init$Sigma))
  K <- t(candidate$D_p_cand) %*% candidate$D_P_S_cand
  Lk <- tryCatch(chol(K), error = function(e) NULL)
  if (is.null(Lk) || any(!is.finite(Lk)) || max(diag(Lk)) / min(diag(Lk)) > 1e7) {
    candidate$log_lik <- -Inf
    candidate$val_3   <- 0
    return(candidate)                                      # reject: the slice shrinks
  }
  # Lk <- tryCatch(chol(K), error = function(e) NULL)
  # if (is.null(Lk)) return(list(log_lik = -Inf, val_3 = 0))
  Kinv    <- chol2inv(Lk)
  candidate$V <-  median(diag(candidate$D_P_S_cand %*% Kinv %*% t(candidate$D_P_S_cand)))
  #lp_sig <- -log(V) - log(1 + (init$sigma_2/(V*hyper$a_sigma))^2)
  lp_sig  <- hyper$a_s2 * log(candidate$V) - hyper$b_s2 * candidate$V / init$sigma_2 # elements of the prior distributions of sigma_2 depending on \mu^*
  logdetK <- n*(chol_det(Lk) - k_free*log(init$sigma_2)) # log determinant 
  
  candidate$val_3 <- -0.5 * (1/init$sigma_2) * trace_function(init$psi,  K) +
    logdetK  + lp_sig
  candidate$log_lik <- -0.5 * candidate$val +  candidate$val_2 + candidate$val_3

}


sample_theta_logit_normal <- function(init, hyper, candidate, X){
  # this function implements the logit-normal reparametrization of the gaps into the angles

  k_l <- init$k_l
  omega_samp <- c(rnorm(k_l), 0)
  u <- runif(n = 1)
  thresh <- hyper$log_lik + hyper$traces + log(u)
  #thresh <- hyper$log_lik + hyper$traces + log(u)
  angle <- runif(n = 1)*2*pi
  angle_min <- angle - 2*pi
  angle_max <- angle
  omega_cand <- init$omega_theta * cos(angle) + omega_samp*sin(angle)
  gaps_cand <- exp(omega_cand)/(sum(exp(omega_cand))) 
  theta_cand <- cumsum(gaps_cand)[-c(k_l + 1)]*(2*pi)
  eval_log_lik_theta_joint(theta_cand, init, hyper, candidate, X) 
  candidate$val <- 0
  while(candidate$log_lik <=   thresh){ # elliptical slice sampling step

    if (angle >=  0) { # bracketing 
      angle_max <- angle
    }
    else{
      angle_min <- angle
    }
    
    angle <- runif(n = 1, min = angle_min, max = angle_max)
    omega_cand <- init$omega_theta * cos(angle) + omega_samp*sin(angle)
    gaps_cand <- exp(omega_cand)/(sum(exp(omega_cand))) 
    theta_cand <- cumsum(gaps_cand)[-c(k_l + 1)]*(2*pi)
    eval_log_lik_theta_joint(theta_cand, init, hyper, candidate, X) 
    candidate$val <- 0
  }
  hyper$log_lik <- candidate$log_lik - candidate$val_3 # update likelihood based on the saved values
  hyper$traces <- candidate$val_3 # update the component depending on sigma^2
  # update all the relevant parameters of interest modified in the function
  res <-  units_loglik(candidate$mean_cand, candidate$D_p_cand, candidate$D_P_S_cand,
                       init$psi, init$R, init$eta_array, X, init$Q_R,
                       as.numeric(init$alphas), keep = TRUE)
  init$mean_i   <- res$mean_i
  init$residual <- res$residual
  hyper$SW      <- res$SW
  init$W        <- res$W
  init$omega_theta <- omega_cand # omegas
  init$thetas <- theta_cand # transformed thetas
  # init$residual <- hyper$residual_array # new residual
  # init$mean_i <- hyper$mean_i_array # new mean for each unit
  init$mean <- candidate$mean_cand # new star shaped mean
  hyper$B_sim <- candidate$B_cand # new Basis matrix 
  hyper$S_mat <- candidate$S_mat_cand # new S matrix
  hyper$D_p <- candidate$D_p_cand # new basis matrix for the warping
  hyper$V <- candidate$V # new version of the scalar dividing sigma_2
  # hyper$SW <- candidate$SW_cand
  # init$W <- candidate$W_cand # new warp 
  init$r <- candidate$r_cand # new radius

}
