# Sample the smoothness/symmetry parameter ------

# sample_tau <- function(init, a_new, b_new){
#   
#   init$tau <- rinvgamma(n = 1, shape = a_new, scale = b_new)
# }



sample_tau <- function(init, hyper, burn){
  log_tau_proposal <- log(init$tau) + t(chol(hyper$lambda_tau* hyper$Sigma_tau))%*%rnorm(n = 2)
  tau_proposal <- exp(log_tau_proposal)
  log_tau_curr <- log(init$tau)
  tau_curr <- init$tau
  
  Sigma_gamma_proposal <- (1/tau_proposal[1])*hyper$K_gamma1 + 
    (1/tau_proposal[2])*hyper$K_gamma2
  U_lambda_proposal <-t(chol(Sigma_gamma_proposal))
  log_cand <- sum(log_tau_proposal) + sum(log(diag(U_lambda_proposal))) - 
    0.5*t(init$gammas)%*%Sigma_gamma_proposal%*%init$gammas - 
    log_tau_proposal[1]*(hyper$a_tau + 1) - hyper$b_tau/tau_proposal[1] - 
    log_tau_proposal[2]*(hyper$a_tau + 1) - hyper$b_tau/tau_proposal[2]
  
  log_curr <- sum(log_tau_curr) + sum(log(diag(hyper$U_lambda))) - 
    0.5*t(init$gammas)%*%hyper$Sigma_gamma_inv%*%init$gammas - 
    log_tau_curr[1]*(hyper$a_tau + 1) - hyper$b_tau/tau_curr[1] - 
    log_tau_curr[2]*(hyper$a_tau + 1) - hyper$b_tau/tau_curr[2]
  
  log_alpha <- min(0, log_cand - log_curr)
  thresh <- log(runif(n = 1))
  
  if (thresh <= log_alpha) {
    init$tau <- tau_proposal
    hyper$U_lambda <- U_lambda_proposal
    hyper$Sigma_gamma_inv <- Sigma_gamma_proposal
    hyper$tau_acc <-  hyper$tau_acc + 1
  }
  
  if (burn == TRUE) {
    hyper$gamma_a <- min(0.01, hyper$t^(-0.5))
    #hyper$lambda_a[i] <- hyper$lambda_a[i]*exp(hyper$gamma_a *(hyper$a_aver_a[i]/hyper$t - hyper$a_opt))
    hyper$lambda_tau <- hyper$lambda_tau*exp(hyper$gamma_a *(exp(log_alpha) - hyper$a_opts))
    hyper$Sigma_tau <- hyper$Sigma_tau + hyper$gamma_a *((log(init$tau) - hyper$mu_tau)%*% t(log(init$tau) - hyper$mu_tau)-
                                                           hyper$Sigma_tau)
    hyper$mu_tau <- hyper$mu_tau + hyper$gamma_a *(log(init$tau)- hyper$mu_tau)
    
  }
  
}
