# Function to sample betas -----
# we have decided to implement elliptical slice sampling for the betas

eval_log_lik_beta <- function(gamma_cand, init, hyper, candidate, X){
  n      <- init$n
  k_l    <- init$k_l
  k_free <- hyper$k_free
  
  candidate$r_cand      <- exp(candidate$B_sim %*% as.vector(gamma_cand))
  candidate$mean_cand   <- cbind(candidate$r_cand * cos(init$thetas),
                                 candidate$r_cand * sin(init$thetas))
  candidate$S_mat_cand  <- S_mat_constructor(candidate$mean_cand)
  candidate$D_p_cand    <- D_p_constructor(candidate$mean_cand, k_l)
  candidate$D_P_S_cand  <- candidate$S_mat_cand %*% candidate$D_p_cand
  
  ul <- units_loglik(candidate$mean_cand, candidate$D_p_cand, candidate$D_P_S_cand,
                     init$psi, init$R, init$eta_array, X, init$Q_R,
                     as.numeric(init$alphas))
  
  candidate$val   <- ul$val
  candidate$val_2 <- -2 * k_l * sum(log(init$alphas)) -
    ((n * k_l) / 2) * log(det(init$Sigma))
  
  SD   <- candidate$D_P_S_cand
  K    <- crossprod(candidate$D_p_cand, SD)
  Lk <- tryCatch(chol(K), error = function(e) NULL)
  if (is.null(Lk) || any(!is.finite(Lk)) || max(diag(Lk)) / min(diag(Lk)) > 1e7) {
    candidate$log_lik <- -Inf
    candidate$val_3   <- 0
    return(candidate)                                      
  }
  Kinv <- chol2inv(Lk)
  candidate$V <- median(rowSums((SD %*% Kinv) * SD))
  lp_sig  <- hyper$a_s2 * log(candidate$V) - hyper$b_s2 * candidate$V / init$sigma_2
  logdetK <- n * (chol_det(Lk) - k_free * log(init$sigma_2))
  
  candidate$val_3   <- -0.5 / init$sigma_2 * ul$trace + logdetK + lp_sig
  candidate$log_lik <- -0.5 * candidate$val + candidate$val_2 + candidate$val_3
}


sample_beta <- function(init, hyper, candidate, X) {
  # B_sim depends only on theta, constant within this step: build it once
  candidate$B_sim <- hyper$B_sim %*% hyper$C_bar
  
  gamma_samp <- solve(t(hyper$U_lambda), rnorm(hyper$L))
  thresh     <- hyper$log_lik + hyper$traces + log(runif(1))
  theta      <- runif(1) * 2 * pi
  theta_min  <- theta - 2 * pi
  theta_max  <- theta
  
  gamma_cand <- init$gammas * cos(theta) + gamma_samp * sin(theta)
  eval_log_lik_beta(gamma_cand, init, hyper, candidate, X)
  
  while (candidate$log_lik <= thresh) {
    if (theta >= 0) theta_max <- theta else theta_min <- theta
    theta      <- runif(1, theta_min, theta_max)
    gamma_cand <- init$gammas * cos(theta) + gamma_samp * sin(theta)
    eval_log_lik_beta(gamma_cand, init, hyper, candidate, X)
  }
  
  # accepted: build the n-unit arrays once
  full <- units_loglik(candidate$mean_cand, candidate$D_p_cand, candidate$D_P_S_cand,
                       init$psi, init$R, init$eta_array, X, init$Q_R,
                       as.numeric(init$alphas), keep = TRUE)
  
  hyper$log_lik <- candidate$log_lik - candidate$val_3
  hyper$traces  <- candidate$val_3
  init$gammas   <- gamma_cand
  init$betas    <- hyper$C_bar %*% gamma_cand
  init$mean     <- candidate$mean_cand
  init$r        <- candidate$r_cand
  hyper$S_mat   <- candidate$S_mat_cand
  hyper$D_p     <- candidate$D_p_cand
  hyper$V       <- candidate$V
  init$W        <- full$W
  hyper$SW      <- full$SW
  init$mean_i   <- full$mean_i
  init$residual <- full$residual
}
