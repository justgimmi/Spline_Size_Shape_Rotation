# Sample the angle parameters -----


eval_log_lik_l <- function(theta_candidate, hyper, init, l, Q_R_11, Q_R_22, Q_R_12, X) {
  # this function is just a wrap up to compute the contribution of the new element inside the log-lik
  theta_mod <- theta_candidate %% (2*pi)
  #theta_mod <- theta_candidate
  B_raw <- Basis_Construction(theta_mod, hyper$n_basis, hyper$degree)
  r_raw <- as.numeric(exp(B_raw %*% as.vector(init$betas)))
  mean_raw <- c(r_raw * cos(theta_mod), r_raw * sin(theta_mod))
  
  out_x <- mean_raw[1] * init$R[1, 1, ] + mean_raw[2] * init$R[2, 1, ]
  out_y <- mean_raw[1] * init$R[1, 2, ] + mean_raw[2] * init$R[2, 2, ]
  
  m_i_1 <- init$alphas * out_x + init$eta[, 1]
  m_i_2 <- init$alphas * out_y + init$eta[, 2]
  
  v_x <- X[l, 1, ] - m_i_1
  v_y <- X[l, 2, ] - m_i_2
  quad <- v_x^2 * Q_R_11 + 2 * v_x * v_y * Q_R_12 + v_y^2 * Q_R_22
  
  return(list(log_lik = -0.5 * sum(quad), mean_i_1 = m_i_1, mean_i_2 = m_i_2, 
              r = r_raw, mean = mean_raw, v_x = v_x, v_y = v_y))
}
# # # 
sample_theta <- function(init, hyper, X){
  # This function build the slice sampler for the angle parameters
  n <- init$n # number of samples
  k <- dim(X)[1] # number of landmarks
  w <- hyper$width_theta # estimated size of the interval
  m <- hyper$m
  # In this way, for each unit, we are decomposing the R^TSigma^-1R in three different elements
  Q_R_11 <- init$Q_R[1,1,]
  Q_R_12 <- init$Q_R[1,2,]
  Q_R_22 <- init$Q_R[2,2,]

  # We want to implement a stepping out procedure as in Neal 2003

  for (l in 1:k) {
    e <- rexp(1)
    #log_thresh <- hyper$log_theta - e
    v <- runif(1)
    left <- init$thetas[l] - v * w
    right <- left + w
    bool <- TRUE
    theta_l <- init$thetas[l]
    v_old_x <- X[l, 1, ] - init$mean_i[l, 1, ] # compute the old value for the residuals x
    v_old_y <- X[l, 2, ] - init$mean_i[l, 2, ] # compute the old value for the residuals x
    quad_old <- v_old_x^2 * Q_R_11 +
      2 * v_old_x * v_old_y * Q_R_12 +
      v_old_y^2 * Q_R_22 # this is the contribution of the l-th landmark on the log-likelihood
    # it is very easy to evaluate and it is the only part changing in the likelihood to be honest
    log_lik_l_old <- -0.5 * sum(quad_old)
    log_thresh <- log_lik_l_old - e
    u <- runif(1)
    J <- floor(m*u)
    up <- (m-1) - J

    # STEPPING OUT
    while ((right - left) < 2*pi && eval_log_lik_l(left, hyper, init, l, Q_R_11, Q_R_22, Q_R_12, X)$log_lik >= log_thresh && J > 0) {
      left <- left - w
       J <- J - 1
    }
    while ((right - left) < 2*pi && eval_log_lik_l(right, hyper, init,l, Q_R_11, Q_R_22, Q_R_12, X)$log_lik >= log_thresh && up > 0) {

      right <- right + w
      up <- up - 1
    }

    if ((right - left) >= 2*pi) {
      right <- left + 2*pi
    }
    #print(c(left, right))
    while (bool == TRUE) {
      theta_raw <- left + runif(n = 1, min = 0, max = 1)*(right - left)
      res_cand <-  eval_log_lik_l(theta_raw, hyper, init, l, Q_R_11, Q_R_22, Q_R_12, X)
      log_theta_prop <- res_cand$log_lik

      if (log_theta_prop >= log_thresh) {
        theta_cand <- theta_raw %% (2*pi)
        init$thetas[l] <- theta_cand
        init$r[l] <- res_cand$r
        init$mean[l, ] <- res_cand$mean

        init$residual[1, 1, ] <- init$residual[1, 1, ] - v_old_x^2 + res_cand$v_x^2
        init$residual[1, 2, ] <- init$residual[1, 2, ] - v_old_x * v_old_y + res_cand$v_x * res_cand$v_y
        init$residual[2, 1, ] <- init$residual[1, 2, ]
        init$residual[2, 2, ] <- init$residual[2, 2, ] - v_old_y^2 + res_cand$v_y^2

        init$mean_i[l, 1, ] <- res_cand$mean_i_1
        init$mean_i[l, 2, ] <- res_cand$mean_i_2
        hyper$log_lik <- hyper$log_lik - log_lik_l_old + res_cand$log_lik

        bool <- FALSE
      } else {
        if (theta_raw < theta_l) {
          left <- theta_raw
        } else {
          right <- theta_raw
        }
      }
    }
  }
  hyper$B_sim <- Basis_Construction(init$thetas, hyper$n_basis, hyper$degree)

  # return(list(init = init, hyper = hyper))
}


# # 
# sample_theta <- function(init, hyper, X){
#   n <- init$n
#   k <- dim(X)[1]
#   w <- hyper$width_theta
#   m <- hyper$m
# 
#   Q_R_11 <- init$Q_R[1,1,]
#   Q_R_12 <- init$Q_R[1,2,]
#   Q_R_22 <- init$Q_R[2,2,]
# 
#   for (l in 1:k) {
# 
#     if (l == 1) {
# 
#       theta_min <- 0
#       theta_max <- init$thetas[2]
#     } else if (l == k) {
# 
#       theta_min <- init$thetas[k - 1]
#       theta_max <-  2*pi
#     } else {
# 
#       theta_min <- init$thetas[l - 1]
#       theta_max <- init$thetas[l + 1]
#     }
# 
#     # -----------------------------------------------------------------
# 
#     e <- rexp(1)
#     v <- runif(1)
# 
#     # Inizializzazione della finestra dello slice sampler
#     left <- init$thetas[l] - v * w
#     right <- left + w
# 
#     # --- NUOVO: Forziamo la finestra iniziale dentro i confini dei vicini ---
#     left <- max(left, theta_min)
#     right <- min(right, theta_max)
#     # -----------------------------------------------------------------------
# 
#     bool <- TRUE
#     theta_l <- init$thetas[l]
# 
#     v_old_x <- X[l, 1, ] - init$mean_i[l, 1, ]
#     v_old_y <- X[l, 2, ] - init$mean_i[l, 2, ]
#     quad_old <- v_old_x^2 * Q_R_11 + 2 * v_old_x * v_old_y * Q_R_12 + v_old_y^2 * Q_R_22
#     log_lik_l_old <- -0.5 * sum(quad_old)
#     log_thresh <- log_lik_l_old - e
# 
#     u <- runif(1)
#     J <- floor(m * u)
#     up <- (m - 1) - J
# 
#     # --- MODIFICATO: STEPPING OUT VINCOLATO ---
#     # Espandiamo a sinistra solo SE non tocchiamo il "muro" del landmark precedente
#     while (left > theta_min &&
#            eval_log_lik_l(left, hyper, init, l, Q_R_11, Q_R_22, Q_R_12, X)$log_lik >= log_thresh &&
#            J > 0) {
#       left <- max(left - w, theta_min) # non scavalcare mai theta_min
#       J <- J - 1
#     }
#     # Espandiamo a destra solo SE non tocchiamo il "muro" del landmark successivo
#     while (right < theta_max &&
#            eval_log_lik_l(right, hyper, init, l, Q_R_11, Q_R_22, Q_R_12, X)$log_lik >= log_thresh &&
#            up > 0) {
#       right <- min(right + w, theta_max) # non scavalcare mai theta_max
#       up <- up - 1
#     }
#     # ------------------------------------------
# 
#     # Rimosso il controllo "if ((right - left) >= 2*pi)" poiché ora l'intervallo
#     # è strutturalmente limitato dai vicini e non supererà mai 2pi.
# 
#     while (bool == TRUE) {
#       theta_raw <- left + runif(n = 1, min = 0, max = 1) * (right - left)
#       res_cand <- eval_log_lik_l(theta_raw, hyper, init, l, Q_R_11, Q_R_22, Q_R_12, X)
#       log_theta_prop <- res_cand$log_lik
# 
#       if (log_theta_prop >= log_thresh) {
#         # theta_cand <- theta_raw %% (2 * pi)
#         theta_cand <- theta_raw
#         init$thetas[l] <- theta_cand
#         init$r[l] <- res_cand$r
#         init$mean[l, ] <- res_cand$mean
# 
#         init$residual[1, 1, ] <- init$residual[1, 1, ] - v_old_x^2 + res_cand$v_x^2
#         init$residual[1, 2, ] <- init$residual[1, 2, ] - v_old_x * v_old_y + res_cand$v_x * res_cand$v_y
#         init$residual[2, 1, ] <- init$residual[1, 2, ]
#         init$residual[2, 2, ] <- init$residual[2, 2, ] - v_old_y^2 + res_cand$v_y^2
# 
#         init$mean_i[l, 1, ] <- res_cand$mean_i_1
#         init$mean_i[l, 2, ] <- res_cand$mean_i_2
#         hyper$log_lik <- hyper$log_lik - log_lik_l_old + res_cand$log_lik
# 
#         bool <- FALSE
#       } else {
#         if (theta_raw < theta_l) {
#           left <- theta_raw
#         } else {
#           right <- theta_raw
#         }
#       }
#     }
#   }
# 
#   #
#   hyper$B_sim <- Basis_Construction(init$thetas, hyper$n_basis, hyper$degree)
# }
