# Sigma sampler -----

sample_Sigma <- function(init, hyper, k_l){
  
  n <- init$n
  
  hyper$nu_post <- hyper$nu + k_l*n
  R_11 <- init$R[1, 1, ]
  R_12 <- init$R[1, 2, ]
  R_21 <- init$R[2, 1, ]
  R_22 <- init$R[2, 2, ]
  
  res_11 <- init$residual[1, 1, ]
  res_12 <- init$residual[1, 2, ]
  res_21 <- init$residual[2, 1, ]
  res_22 <- init$residual[2, 2, ]
  
  Int_11 <- R_11 * res_11 + R_12 * res_21
  Int_12 <- R_11 * res_12 + R_12 * res_22
  Int_21 <- R_21 * res_11 + R_22 * res_21
  Int_22 <- R_21 * res_12 + R_22 * res_22
  
  M_11 <- Int_11 * R_11 + Int_12 * R_12
  M_12 <- Int_11 * R_21 + Int_12 * R_22
  M_21 <- Int_21 * R_11 + Int_22 * R_12
  M_22 <- Int_21 * R_21 + Int_22 * R_22
  
  w_alpha <- 1 / (init$alphas^2)
  
  psi_update <- matrix(c(
    sum(w_alpha * M_11), sum(w_alpha * M_12),
    sum(w_alpha * M_21), sum(w_alpha * M_22)
  ), nrow = 2, ncol = 2, byrow = TRUE)
  
  
  hyper$psi_post <- hyper$psi + psi_update
  init$Sigma <- riwish(hyper$nu_post,  hyper$psi_post)
  init$Sigma_inv <- solve(init$Sigma)
  
  for (i in 1:n) {
    init$Q_R[,,i] <- (1/(init$alphas[i]^2)) *t(init$R[,,i])%*%init$Sigma_inv%*%init$R[,,i] 
  }
  
  # out <- list()
  # out$init <- init
  # out$hyper <- hyper
  # return(out)
  
}
