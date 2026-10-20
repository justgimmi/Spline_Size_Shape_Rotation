# Sigma sampler -----

# sample_Sigma <- function(init, hyper, k_l){# this function has been checked
#   # this is the gibbs step for the Sigma parameter
#   
#   n <- init$n
#   
#   hyper$nu_post <- hyper$nu + k_l*n
#   R_11 <- init$R[1, 1, ]
#   R_12 <- init$R[1, 2, ]
#   R_21 <- init$R[2, 1, ]
#   R_22 <- init$R[2, 2, ]
#   
#   res_11 <- init$residual[1, 1, ]
#   res_12 <- init$residual[1, 2, ]
#   res_21 <- init$residual[2, 1, ]
#   res_22 <- init$residual[2, 2, ]
#   
#   Int_11 <- R_11 * res_11 + R_12 * res_21
#   Int_12 <- R_11 * res_12 + R_12 * res_22
#   Int_21 <- R_21 * res_11 + R_22 * res_21
#   Int_22 <- R_21 * res_12 + R_22 * res_22
#   
#   M_11 <- Int_11 * R_11 + Int_12 * R_12
#   M_12 <- Int_11 * R_21 + Int_12 * R_22
#   M_21 <- Int_21 * R_11 + Int_22 * R_12
#   M_22 <- Int_21 * R_21 + Int_22 * R_22
#   
#   w_alpha <- 1 / (init$alphas^2)
#   
#   psi_update <- matrix(c(
#     sum(w_alpha * M_11), sum(w_alpha * M_12),
#     sum(w_alpha * M_21), sum(w_alpha * M_22)
#   ), nrow = 2, ncol = 2, byrow = TRUE)
#   
#   
#   hyper$psi_post <- hyper$psi + psi_update
#   init$Sigma <- riwish(hyper$nu_post,  hyper$psi_post)
#   init$Sigma_inv <- solve(init$Sigma)
#   
#   for (i in 1:n) {
#     init$Q_R[,,i] <- (1/(init$alphas[i]^2)) *t(init$R[,,i])%*%init$Sigma_inv%*%init$R[,,i] 
#   }
#   hyper$log_lik <- log_density(init)
#   
# }

sample_Sigma <- function(init, hyper, k_l) {
  n <- init$n
  hyper$nu_post <- hyper$nu + k_l * n
  
  cs <- init$R[1, 1, ]; sn <- init$R[2, 1, ]          # cos(lambda_i), sin(lambda_i)
  r11 <- init$residual[1, 1, ]
  r12 <- init$residual[1, 2, ]
  r22 <- init$residual[2, 2, ]
  
  # M_i = R_i res_i R_i'  (res_i symmetric)
  M11 <- cs^2 * r11 - 2 * cs * sn * r12 + sn^2 * r22
  M22 <- sn^2 * r11 + 2 * cs * sn * r12 + cs^2 * r22
  M12 <- cs * sn * (r11 - r22) + (cs^2 - sn^2) * r12
  
  a2 <- init$alphas^2
  w  <- 1 / a2
  psi_update <- matrix(c(sum(w * M11), sum(w * M12),
                         sum(w * M12), sum(w * M22)), 2, 2)
  hyper$psi_post <- hyper$psi + psi_update
  
  # Sigma ~ IW(nu, psi)  <=>  Sigma^-1 ~ Wishart(nu, psi^-1)
  init$Sigma_inv <- rWishart(1, hyper$nu_post, solve(hyper$psi_post))[, , 1]
  init$Sigma     <- solve(init$Sigma_inv)
  
  # Q_R_i = R_i' Sigma^-1 R_i / alpha_i^2, vectorised over i
  Si  <- init$Sigma_inv
  t11 <- cs^2 * Si[1, 1] + 2 * cs * sn * Si[1, 2] + sn^2 * Si[2, 2]
  t22 <- sn^2 * Si[1, 1] - 2 * cs * sn * Si[1, 2] + cs^2 * Si[2, 2]
  t12 <- cs * sn * (Si[2, 2] - Si[1, 1]) + (cs^2 - sn^2) * Si[1, 2]
  init$Q_R[1, 1, ] <- t11 * w
  init$Q_R[2, 2, ] <- t22 * w
  init$Q_R[1, 2, ] <- t12 * w
  init$Q_R[2, 1, ] <- t12 * w
  
  hyper$log_lik <- log_density(init)
}