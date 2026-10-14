###### Theta Sampler #####
# In this Rcode, you can find the sampler for thetas:
# The basic idea is to not directly model thetas, but to impose a logit normal prior on the consecutive angular differences
# basically, we have derived an elliptical slice sampler
# the only problem of this type of sampler is that we do not know how many calls we will make to the function
# an interesting idea that can be used is to define first all the objects that do not depend on theta on a separate function
# so we just need to compute the objects that actually depend on theta 

eval_log_lik_theta_joint <- function(theta_cand, init, hyper, X) {
  
  k_l <- init$k_l # landmarks
  n <- init$n # sample size
  k_free <- hyper$k_free # free dimensions
  
  B_cand <- Basis_Construction(theta_cand, hyper$n_basis, hyper$degree) # construct the candidate basis
  B_sim <- B_cand %*% hyper$C_bar # apply the subspace representation
  r_cand <- exp(B_sim %*% as.vector(init$gammas)) # propose a new radius
  
  # Construct the global latent mean shape and updaye the relevant objects depending on \mu^*
  mean_cand <- cbind(r_cand*cos(theta_cand), r_cand*sin(theta_cand)) # proposed mean
  S_mat_cand <- S_mat_constructor(mean_cand) # proposed S(mu^*) 
  D_p_cand <- D_p_constructor(mean_cand, k_l) # proposed D_p
  D_P_S_cand <- S_mat_cand %*% D_p_cand # proposed S x D_p (k x k-3)
  W_cand <- tensor(D_p_cand, init$psi, 2, 1) # create the candidate value W
  SW_cand <- tensor(S_mat_cand, W_cand, 2, 1)# candidate value of SW
  mean_i_array <- init$mean_i # candidate new mean i
  residual_array <- init$residual # candidate new residual 
  # this function is used to do mean_cand + SW_i applied over every i
  app_cand <- sweep(SW_cand, MARGIN = c(1, 2), STATS = mean_cand, FUN = "+")
  val <- 0
  for (i in 1:n) {
    mean_i_cand <- init$alphas[i] * app_cand[,,i] %*% init$R[,,i] + init$eta_array[,,i]
    res_cand <- X[,,i] - mean_i_cand
    res_mat <- crossprod(res_cand) # t(res_cand) %*% res_cand
    mean_i_array[,,i] <- mean_i_cand
    residual_array[,,i] <- res_mat
    val <- val + sum(init$Q_R[,,i] * res_mat)
  }
  
  val_2 <- -2 * k_l * sum(log(init$alphas)) - ((n * k_l) / 2) * log(det(init$Sigma))
  K <- t(D_p_cand) %*% D_P_S_cand
  Lk <- chol(K)
  # Lk <- tryCatch(chol(K), error = function(e) NULL)
  # if (is.null(Lk)) return(list(log_lik = -Inf, val_3 = 0))
  Kinv    <- chol2inv(Lk)
  V <-  median(diag(D_P_S_cand %*% Kinv %*% t(D_P_S_cand)))
  #lp_sig <- -log(V) - log(1 + (init$sigma_2/(V*hyper$a_sigma))^2)
  lp_sig  <- hyper$a_s2 * log(V) - hyper$b_s2 * V / init$sigma_2 # elements of the prior distributions of sigma_2 depending on \mu^*
  logdetK <- n*(chol_det(Lk) - k_free*log(init$sigma_2)) # log determinant 
  
  val_3 <- -0.5 * (1/init$sigma_2) * trace_function(W_cand, S_mat_cand) +
    logdetK  + lp_sig
  log_lik_total <- -0.5 * val + val_2 + val_3
  return(list(log_lik = log_lik_total, 
              mean = mean_cand, 
              mean_i = mean_i_array, 
              residual = residual_array,
              B_sim = B_cand,
              S_mat = S_mat_cand,
              W = W_cand,
              SW = SW_cand,
              D_p = D_p_cand,
              val_3 = val_3,
              r = r_cand, 
              V = V))
}


sample_theta_logit_normal <- function(init, hyper, X){
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
  res_cand <-eval_log_lik_theta_joint(theta_cand, init, hyper, X) 
  while(res_cand$log_lik <=   thresh){ # elliptical slice sampling step

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
    res_cand <- eval_log_lik_theta_joint(theta_cand, init, hyper, X) 
  }
  hyper$log_lik <- res_cand$log_lik - res_cand$val_3 # update likelihood based on the saved values
  hyper$traces <- res_cand$val_3 # update the component depending on sigma^2
  # update all the relevant parameters of interest modified in the function
  
  init$omega_theta <- omega_cand # omegas
  init$thetas <- theta_cand # transformed thetas
  init$residual <- res_cand$residual # new residual
  init$mean_i <- res_cand$mean_i # new mean for each unit
  init$mean <- res_cand$mean # new star shaped mean
  hyper$B_sim <- res_cand$B_sim # new Basis matrix 
  hyper$S_mat <- res_cand$S_mat # new S matrix
  hyper$D_p <- res_cand$D_p # new basis matrix for the warping
  hyper$V <- res_cand$V # new version of the scalar dividing sigma_2
  hyper$SW <- res_cand$SW
  init$W <- res_cand$W # new warp 
  init$r <- res_cand$r # new radius

}
