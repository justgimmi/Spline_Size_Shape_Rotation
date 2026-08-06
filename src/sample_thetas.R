# # # Sample the angle parameters -----
# # 
# # 
# eval_log_lik_l <- function(theta_candidate, hyper, init, l, Q_R_11, Q_R_22, Q_R_12, X) {
#   # this function is just a wrap up to compute the contribution of the new element inside the log-lik
#   theta_mod <- theta_candidate %% (2*pi)
#   #theta_mod <- theta_candidate
#   B_raw <- Basis_Construction(theta_mod, hyper$n_basis, hyper$degree)
#   r_raw <- as.numeric(exp(B_raw %*% as.vector(init$betas)))
#   mean_raw <- c(r_raw * cos(theta_mod), r_raw * sin(theta_mod))
# 
#   out_x <- mean_raw[1] * init$R[1, 1, ] + mean_raw[2] * init$R[2, 1, ]
#   out_y <- mean_raw[1] * init$R[1, 2, ] + mean_raw[2] * init$R[2, 2, ]
# 
#   m_i_1 <- init$alphas * out_x + init$eta[, 1]
#   m_i_2 <- init$alphas * out_y + init$eta[, 2]
# 
#   v_x <- X[l, 1, ] - m_i_1
#   v_y <- X[l, 2, ] - m_i_2
#   quad <- v_x^2 * Q_R_11 + 2 * v_x * v_y * Q_R_12 + v_y^2 * Q_R_22
# 
#   return(list(log_lik = -0.5 * sum(quad), mean_i_1 = m_i_1, mean_i_2 = m_i_2,
#               r = r_raw, mean = mean_raw, v_x = v_x, v_y = v_y))
# }
# # # # #
# # sample_theta <- function(init, hyper, X){
# #   # This function build the slice sampler for the angle parameters
# #   n <- init$n # number of samples
# #   k <- dim(X)[1] # number of landmarks
# #   w <- hyper$width_theta # estimated size of the interval
# #   m <- hyper$m
# #   # In this way, for each unit, we are decomposing the R^TSigma^-1R in three different elements
# #   Q_R_11 <- init$Q_R[1,1,]
# #   Q_R_12 <- init$Q_R[1,2,]
# #   Q_R_22 <- init$Q_R[2,2,]
# # 
# #   # We want to implement a stepping out procedure as in Neal 2003
# # 
# #   for (l in 1:k) {
# #     e <- rexp(1)
# #     #log_thresh <- hyper$log_theta - e
# #     v <- runif(1)
# #     left <- init$thetas[l] - v * w
# #     right <- left + w
# #     bool <- TRUE
# #     theta_l <- init$thetas[l]
# #     v_old_x <- X[l, 1, ] - init$mean_i[l, 1, ] # compute the old value for the residuals x
# #     v_old_y <- X[l, 2, ] - init$mean_i[l, 2, ] # compute the old value for the residuals x
# #     quad_old <- v_old_x^2 * Q_R_11 +
# #       2 * v_old_x * v_old_y * Q_R_12 +
# #       v_old_y^2 * Q_R_22 # this is the contribution of the l-th landmark on the log-likelihood
# #     # it is very easy to evaluate and it is the only part changing in the likelihood to be honest
# #     log_lik_l_old <- -0.5 * sum(quad_old)
# #     log_thresh <- log_lik_l_old - e
# #     u <- runif(1)
# #     J <- floor(m*u)
# #     up <- (m-1) - J
# # 
# #     # STEPPING OUT
# #     while ((right - left) < 2*pi && eval_log_lik_l(left, hyper, init, l, Q_R_11, Q_R_22, Q_R_12, X)$log_lik >= log_thresh && J > 0) {
# #       left <- left - w
# #        J <- J - 1
# #     }
# #     while ((right - left) < 2*pi && eval_log_lik_l(right, hyper, init,l, Q_R_11, Q_R_22, Q_R_12, X)$log_lik >= log_thresh && up > 0) {
# # 
# #       right <- right + w
# #       up <- up - 1
# #     }
# # 
# #     if ((right - left) >= 2*pi) {
# #       right <- left + 2*pi
# #     }
# #     #print(c(left, right))
# #     while (bool == TRUE) {
# #       theta_raw <- left + runif(n = 1, min = 0, max = 1)*(right - left)
# #       res_cand <-  eval_log_lik_l(theta_raw, hyper, init, l, Q_R_11, Q_R_22, Q_R_12, X)
# #       log_theta_prop <- res_cand$log_lik
# # 
# #       if (log_theta_prop >= log_thresh) {
# #         theta_cand <- theta_raw %% (2*pi)
# #         init$thetas[l] <- theta_cand
# #         init$r[l] <- res_cand$r
# #         init$mean[l, ] <- res_cand$mean
# # 
# #         init$residual[1, 1, ] <- init$residual[1, 1, ] - v_old_x^2 + res_cand$v_x^2
# #         init$residual[1, 2, ] <- init$residual[1, 2, ] - v_old_x * v_old_y + res_cand$v_x * res_cand$v_y
# #         init$residual[2, 1, ] <- init$residual[1, 2, ]
# #         init$residual[2, 2, ] <- init$residual[2, 2, ] - v_old_y^2 + res_cand$v_y^2
# # 
# #         init$mean_i[l, 1, ] <- res_cand$mean_i_1
# #         init$mean_i[l, 2, ] <- res_cand$mean_i_2
# #         hyper$log_lik <- hyper$log_lik - log_lik_l_old + res_cand$log_lik
# # 
# #         bool <- FALSE
# #       } else {
# #         if (theta_raw < theta_l) {
# #           left <- theta_raw
# #         } else {
# #           right <- theta_raw
# #         }
# #       }
# #     }
# #   }
# #   hyper$log_lik <- log_density(init)
# #   hyper$B_sim <- Basis_Construction(init$thetas, hyper$n_basis, hyper$degree)
# # 
# #   # return(list(init = init, hyper = hyper))
# # }
# 
# #
# 
# 
# sample_theta_gaps <- function(init, hyper, X){
#   n <- init$n
#   k <- dim(X)[1]
#   w <- hyper$width_theta
#   m <- hyper$m
#   Q_R_11 <- init$Q_R[1,1,]
#   Q_R_12 <- init$Q_R[1,2,]
#   Q_R_22 <- init$Q_R[2,2,]
#   
#   # Converti a gap
#   delta <- numeric(k)
#   delta[1] <- init$thetas[1]  # delta[1] = theta[1] - 0
#   for (l in 2:k) {
#     delta[l] <- init$thetas[l] - init$thetas[l-1]
#   }
#   
#   # Campiona ogni gap
#   for (l in 1:k) {
#     
#     # Likelihood corrente
#     v_old_x <- X[l, 1, ] - init$mean_i[l, 1, ]
#     v_old_y <- X[l, 2, ] - init$mean_i[l, 2, ]
#     quad_old <- v_old_x^2 * Q_R_11 + 2 * v_old_x * v_old_y * Q_R_12 + v_old_y^2 * Q_R_22
#     log_lik_l_old <- -0.5 * sum(quad_old)
#     
#     # Slice sampler SUL GAP (senza vincoli espliciti)
#     e <- rexp(1)
#     log_thresh <- log_lik_l_old - e
#     
#     v <- runif(1)
#     w_actual <- w
#     left <- delta[l] - v * w_actual
#     right <- left + w_actual
#     
#     bool <- TRUE
#     delta_l <- delta[l]
#     max_iter <- 1000
#     iter <- 0
#     
#     left <- max(left, 1e-6)   # delta deve essere positivo
#     right <- max(right, 1e-6)
#     
#     # Stepping out
#     J <- floor(m * runif(1))
#     up <- (m - 1) - J
#     
#     while (J > 0) {
#       # Calcola theta proposto da questo gap
#       theta_prop <- if(l == 1) left else sum(delta[1:(l-1)]) + left
#       res_left <- eval_log_lik_l(theta_prop, hyper, init, l, Q_R_11, Q_R_22, Q_R_12, X)
#       if (res_left$log_lik < log_thresh) break
#       left <- left - w_actual
#       left <- max(left, 1e-6)
#       J <- J - 1
#     }
#     
#     while (up > 0) {
#       theta_prop <- if(l == 1) right else sum(delta[1:(l-1)]) + right
#       res_right <- eval_log_lik_l(theta_prop, hyper, init, l, Q_R_11, Q_R_22, Q_R_12, X)
#       if (res_right$log_lik < log_thresh) break
#       right <- right + w_actual
#       up <- up - 1
#     }
#     
#     # Shrinkage
#     while (bool && iter < max_iter) {
#       iter <- iter + 1
#       delta_raw <- left + runif(1) * (right - left)
#       delta_raw <- max(delta_raw, 1e-6)
#       
#       theta_prop <- if(l == 1) delta_raw else sum(delta[1:(l-1)]) + delta_raw
#       
#       res_cand <- eval_log_lik_l(theta_prop, hyper, init, l, Q_R_11, Q_R_22, Q_R_12, X)
#       log_theta_prop <- res_cand$log_lik
#       
#       if (log_theta_prop >= log_thresh) {
#         # Accettato
#         delta[l] <- delta_raw
#         init$thetas[l] <- if(l == 1) delta[l] else sum(delta[1:l])
#         init$r[l] <- res_cand$r
#         init$mean[l, ] <- res_cand$mean
#         
#         init$residual[1, 1, ] <- init$residual[1, 1, ] - v_old_x^2 + res_cand$v_x^2
#         init$residual[1, 2, ] <- init$residual[1, 2, ] - v_old_x * v_old_y + res_cand$v_x * res_cand$v_y
#         init$residual[2, 1, ] <- init$residual[1, 2, ]
#         init$residual[2, 2, ] <- init$residual[2, 2, ] - v_old_y^2 + res_cand$v_y^2
#         init$mean_i[l, 1, ] <- res_cand$mean_i_1
#         init$mean_i[l, 2, ] <- res_cand$mean_i_2
#         
#         hyper$log_lik <- hyper$log_lik - log_lik_l_old + res_cand$log_lik
#         bool <- FALSE
#       } else {
#         if (delta_raw < delta_l) {
#           left <- delta_raw
#         } else {
#           right <- delta_raw
#         }
#       }
#     }
#   }
#   
#   hyper$log_lik <- log_density(init)
#   hyper$B_sim <- Basis_Construction(init$thetas, hyper$n_basis, hyper$degree)
# }


# # #
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
#     # 1. Definizione rigorosa dei muri di contenimento (Prior Ordinata)
#     if (l == 1) {
#       theta_min <- 0
#       theta_max <- init$thetas[2]
#     } else if (l == k) {
#       theta_min <- init$thetas[k - 1]
#       theta_max <- 2*pi
#     } else {
#       theta_min <- init$thetas[l - 1]
#       theta_max <- init$thetas[l + 1]
#     }
# 
#     e <- rexp(1)
#     v <- runif(1)
# 
#     # Inizializzazione STANDARD della finestra (SENZA tagliare con max/min)
#     left <- init$thetas[l] - v * w
#     right <- left + w
# 
#     bool <- TRUE
#     theta_l <- init$thetas[l]
# 
#     # Calcolo locale della likelihood corrente per il landmark l
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
#     # --- STEPPING OUT VINCOLATO (Teoricamente Corretto) ---
#     # Se left supera il confine sinistro, assegniamo -Inf senza sprecare calcoli
#     # nella funzione di likelihood. Il loop si interrompe immediatamente.
#     while (J > 0) {
#       log_lik_left <- if (left < theta_min) -Inf else eval_log_lik_l(left, hyper, init, l, Q_R_11, Q_R_22, Q_R_12, X)$log_lik
#       if (log_lik_left < log_thresh) break
#       left <- left - w
#       J <- J - 1
#     }
# 
#     while (up > 0) {
#       log_lik_right <- if (right > theta_max) -Inf else eval_log_lik_l(right, hyper, init, l, Q_R_11, Q_R_22, Q_R_12, X)$log_lik
#       if (log_lik_right < log_thresh) break
#       right <- right + w
#       up <- up - 1
#     }
# 
#     # --- SHRINKAGE SYSTEM ---
#     while (bool == TRUE) {
#       theta_raw <- left + runif(n = 1, min = 0, max = 1) * (right - left)
# 
#       # Se il candidato estratto cade fuori dai vincoli di ordinamento,
#       # la sua likelihood totale è -Inf (rifiutato e usato per restringere)
#       if (theta_raw < theta_min || theta_raw > theta_max) {
#         log_theta_prop <- -Inf
#       } else {
#         res_cand <- eval_log_lik_l(theta_raw, hyper, init, l, Q_R_11, Q_R_22, Q_R_12, X)
#         log_theta_prop <- res_cand$log_lik
#       }
# 
#       if (log_theta_prop >= log_thresh) {
#         # Accettato! Aggiorna lo stato locale
#         theta_cand <- theta_raw
#         init$thetas[l] <- theta_cand
#         init$r[l] <- res_cand$r
#         init$mean[l, ] <- res_cand$mean
# 
#         # Sottrazione del vecchio e addizione del nuovo (Vettorizzato su n soggetti)
#         init$residual[1, 1, ] <- init$residual[1, 1, ] - v_old_x^2 + res_cand$v_x^2
#         init$residual[1, 2, ] <- init$residual[1, 2, ] - v_old_x * v_old_y + res_cand$v_x * res_cand$v_y
#         init$residual[2, 1, ] <- init$residual[1, 2, ]
#         init$residual[2, 2, ] <- init$residual[2, 2, ] - v_old_y^2 + res_cand$v_y^2
# 
#         init$mean_i[l, 1, ] <- res_cand$mean_i_1
#         init$mean_i[l, 2, ] = res_cand$mean_i_2
#         hyper$log_lik <- hyper$log_lik - log_lik_l_old + res_cand$log_lik
# 
#         bool <- FALSE
#       } else {
#         # Restringimento della finestra dello slice sampler
#         if (theta_raw < theta_l) {
#           left <- theta_raw
#         } else {
#           right <- theta_raw
#         }
#       }
#     }
#   }
# 
#   # A fine ciclo di tutti i landmark, rinfreschiamo le matrici globali per sample_beta
#   hyper$log_lik <- log_density(init)
#   hyper$B_sim <- Basis_Construction(init$thetas, hyper$n_basis, hyper$degree)
# }

sample_theta <- function(init, hyper, X){
  n <- init$n
  k <- dim(X)[1]
  
  # Precision matrix components for local likelihood calculations
  Q_R_11 <- init$Q_R[1,1,]
  Q_R_12 <- init$Q_R[1,2,]
  Q_R_22 <- init$Q_R[2,2,]
  
  # 1. INITIALIZE GAMMAS (If they don't exist yet)
  # Recovers latent gammas from the initial thetas to ensure a seamless chain continuation
  if (is.null(init$gam)) {
    theta_extended <- c(0, init$thetas)
    w <- diff(theta_extended) / (2 * pi)
    w <- pmax(w, 1e-8)  # Prevent log(0)
    w <- w / sum(w)
    init$gam <- log(w)
  }
  
  gamma_current <- init$gam
  
  # 2. COMPUTE CURRENT TOTAL LOG-LIKELIHOOD
  # Under the independence assumption across landmarks, the total log-likelihood 
  # is simply the sum of individual landmark likelihoods
  log_lik_current <- 0
  for (l in 1:k) {
    res <- eval_log_lik_l(init$thetas[l], hyper, init, l, Q_R_11, Q_R_22, Q_R_12, X)
    log_lik_current <- log_lik_current + res$log_lik
  }
  log_thresh <- log_lik_current + log(runif(1))
  
  # 3. ELLIPTICAL SLICE SAMPLING SETUP
  # Define prior variance for your independent gammas (defaults to 1 if not specified)
  sigma_gamma <- 1.0
  nu <- rnorm(k, mean = 0, sd = sigma_gamma) # Draw an ideal sample from the Gaussian prior
  
  # Set up the elliptical search bracket
  phi <- runif(1, 0, 2 * pi)
  phi_min <- phi - 2 * pi
  phi_max <- phi
  
  # 4. ESS SHRINKAGE LOOP
  while (TRUE) {
    # Generate proposal on the ellipse
    gamma_prop <- gamma_current * cos(phi) + nu * sin(phi)
    
    # Softmax transformation to obtain ordered increments mapping to [0, 2*pi]
    w_prop <- exp(gamma_prop) / sum(exp(gamma_prop))
    theta_prop <- 2 * pi * cumsum(w_prop)
    
    # Evaluate global log-likelihood for the proposed vector arrangement
    log_lik_prop <- 0
    res_prop_list <- list()
    for (l in 1:k) {
      res <- eval_log_lik_l(theta_prop[l], hyper, init, l, Q_R_11, Q_R_22, Q_R_12, X)
      log_lik_prop <- log_lik_prop + res$log_lik
      res_prop_list[[l]] <- res
    }
    
    # Check Acceptance Condition
    if (log_lik_prop >= log_thresh) {
      # --- ACCEPT STEP ---
      init$gam <- gamma_prop
      init$thetas <- theta_prop
      
      # Unpack and update landmark-specific fields inside init
      for (l in 1:k) {
        res_cand <- res_prop_list[[l]]
        init$r[l] <- res_cand$r
        init$mean[l, ] <- res_cand$mean
        init$mean_i[l, 1, ] <- res_cand$mean_i_1
        init$mean_i[l, 2, ] <- res_cand$mean_i_2
      }
      
      # 2. Vectorized Residual Reconstruction (No zeros, no loops, no redundant self-adding)
      # This stacks the length-n vectors from all k landmarks into k x n matrices
      v_x_mat <- do.call(rbind, lapply(res_prop_list, function(x) x$v_x))
      v_y_mat <- do.call(rbind, lapply(res_prop_list, function(x) x$v_y))
      
      # colSums collapses the k landmarks instantly into a clean vector of length n
      init$residual[1, 1, ] <- colSums(v_x_mat^2)
      init$residual[1, 2, ] <- colSums(v_x_mat * v_y_mat)
      init$residual[2, 2, ] <- colSums(v_y_mat^2)
      init$residual[2, 1, ] <- init$residual[1, 2, ] # Mirror the covariance
      
      break # Exit the slice sampling loop
    } else {
      # --- SHRINKAGE STEP ---
      if (phi < 0) {
        phi_min <- phi
      } else {
        phi_max <- phi
      }
      phi <- runif(1, phi_min, phi_max)
    }
  }
  
  # 5. REFRESH GLOBAL MATRICES (Preserved from original script for downstream steps)
  hyper$log_lik <- log_density(init)
  hyper$B_sim <- Basis_Construction(init$thetas, hyper$n_basis, hyper$degree)
  
  # Return statements (if your MCMC loop does not pass init/hyper as environments)
  # return(list(init = init, hyper = hyper))
}

# ##### new code ##### 
# eval_log_lik_global <- function(thetas_prop, hyper, init, X) {
#   thetas_sorted <- sort(thetas_prop %% (2*pi))
#   
#   B_raw <- Basis_Construction(thetas_sorted, hyper$n_basis, hyper$degree)
#   r_raw <- as.numeric(exp(B_raw %*% as.vector(init$betas)))
#   mean_raw_x <- r_raw * cos(thetas_sorted)
#   mean_raw_y <- r_raw * sin(thetas_sorted)
#   
#   n <- init$n
#   k <- length(thetas_prop)
#   Q_R_11 <- init$Q_R[1,1,]
#   Q_R_12 <- init$Q_R[1,2,]
#   Q_R_22 <- init$Q_R[2,2,]
#   
#   quad_total <- 0
#   mean_i_1_all <- array(0, dim = c(k, n))
#   mean_i_2_all <- array(0, dim = c(k, n))
#   v_x_all <- array(0, dim = c(k, n))
#   v_y_all <- array(0, dim = c(k, n))
#   residual_mat <- array(0, dim = c(2, 2, n))
#   
#   for (i in 1:n) {
#     out_x <- mean_raw_x * init$R[1,1,i] + mean_raw_y * init$R[2,1,i]
#     out_y <- mean_raw_x * init$R[1,2,i] + mean_raw_y * init$R[2,2,i]
#     m_i_1 <- init$alphas[i] * out_x + init$eta[i,1]
#     m_i_2 <- init$alphas[i] * out_y + init$eta[i,2]
#     
#     v_x <- X[, 1, i] - m_i_1
#     v_y <- X[, 2, i] - m_i_2
#     
#     quad <- v_x^2*Q_R_11[i] + 2*v_x*v_y*Q_R_12[i] + v_y^2*Q_R_22[i]
#     quad_total <- quad_total + sum(quad)
#     
#     mean_i_1_all[,i] <- m_i_1
#     mean_i_2_all[,i] <- m_i_2
#     v_x_all[,i] <- v_x
#     v_y_all[,i] <- v_y
#     
#     # Pre-compute aggregated residuals here
#     residual_mat[1,1,i] <- sum(v_x^2)
#     residual_mat[1,2,i] <- sum(v_x * v_y)
#     residual_mat[2,1,i] <- residual_mat[1,2,i]
#     residual_mat[2,2,i] <- sum(v_y^2)
#   }
#   
#   val_2 <- -2*init$k_l*sum(log(init$alphas)) - ((n*init$k_l)/2)*log(det(init$Sigma))
#   list(log_lik = -0.5*quad_total + val_2, 
#        mean_i_1 = mean_i_1_all, 
#        mean_i_2 = mean_i_2_all,
#        r = r_raw, 
#        mean = cbind(mean_raw_x, mean_raw_y),
#        thetas_sorted = thetas_sorted, 
#        v_x = v_x_all, 
#        v_y = v_y_all,
#        residual = residual_mat)  # ← Return pre-computed
# }
# 
# sample_theta <- function(init, hyper, X){
#   n <- init$n
#   k <- dim(X)[1]
#   w <- hyper$width_theta
#   m <- hyper$m
#   
#   for (l in 1:k) {
#     e <- rexp(1)
#     v <- runif(1)
#     left <- init$thetas[l] - v * w
#     right <- left + w
#     bool <- TRUE
#     theta_l <- init$thetas[l]
#     
#     log_lik_old <- hyper$log_lik
#     log_thresh <- log_lik_old - e
#     
#     u <- runif(1)
#     J <- floor(m * u)
#     up <- (m - 1) - J
#     
#     # STEPPING OUT (optimized)
#     while ((right - left) < 2*pi && J > 0) {
#       thetas_left <- init$thetas
#       thetas_left[l] <- left
#       log_lik_left <- eval_log_lik_global(thetas_left, hyper, init, X)$log_lik
#       
#       if (log_lik_left < log_thresh) break
#       left <- left - w
#       J <- J - 1
#     }
#     
#     while ((right - left) < 2*pi && up > 0) {
#       thetas_right <- init$thetas
#       thetas_right[l] <- right
#       log_lik_right <- eval_log_lik_global(thetas_right, hyper, init, X)$log_lik
#       
#       if (log_lik_right < log_thresh) break
#       right <- right + w
#       up <- up - 1
#     }
#     
#     if ((right - left) >= 2*pi) {
#       right <- left + 2*pi
#     }
#     
#     # SHRINKAGE SYSTEM
#     while (bool == TRUE) {
#       theta_raw <- left + runif(n = 1, min = 0, max = 1) * (right - left)
#       
#       thetas_cand_vec <- init$thetas
#       thetas_cand_vec[l] <- theta_raw
#       
#       res_cand <- eval_log_lik_global(thetas_cand_vec, hyper, init, X)
#       log_theta_prop <- res_cand$log_lik
#       
#       if (log_theta_prop >= log_thresh) {
#         # Accept: update all state at once
#         init$thetas[l] <- theta_raw %% (2*pi)
#         init$thetas_sorted <- res_cand$thetas_sorted
#         init$r <- res_cand$r
#         init$mean <- res_cand$mean
#         init$residual <- res_cand$residual  # ← Use pre-computed
#         
#         # Update mean_i efficiently
#         for(i in 1:n) {
#           init$mean_i[,,i] <- cbind(res_cand$mean_i_1[,i], res_cand$mean_i_2[,i])
#         }
#         
#         hyper$log_lik <- log_theta_prop
#         bool <- FALSE
#       } else {
#         # Shrinkage
#         if (theta_raw < theta_l) {
#           left <- theta_raw
#         } else {
#           right <- theta_raw
#         }
#       }
#     }
#   }
#   
#   # Update B_sim with sorted thetas
#   hyper$B_sim <- Basis_Construction(init$thetas_sorted, hyper$n_basis, hyper$degree)
# }



###### new idea ###### 



eval_log_lik_theta_joint <- function(theta_cand, init, hyper, X) {
  k_l <- init$k_l
  n <- init$n
  
  B_cand <- Basis_Construction(theta_cand, hyper$n_basis, hyper$degree)
  r_cand <- as.numeric(exp(B_cand %*% as.vector(init$betas)))
  
  # 3. Construct the global latent mean shape
  mean_cand <- cbind(r_cand * cos(theta_cand), r_cand * sin(theta_cand))
  
  # 4. Evaluate the aggregate fit across all subjects
  val <- 0
  mean_i_array <- array(0, dim = c(k_l, 2, n))
  residual_array <- array(0, dim = c(2, 2, n))
  mat_cand <- covariance_mat(init$mat_dist, theta_cand)
  C_cand <- exp(-mat_cand/init$phi)
  chol_c_cand <- chol(C_cand)
  for (i in 1:n) {
    eta_matrix <- matrix(init$eta[i, ], nrow = k_l, ncol = 2, byrow = TRUE)
    mean_i_cand <- init$alphas[i] * (mean_cand %*% init$R[,,i]) + eta_matrix
    res_cand <- X[,,i] - mean_i_cand
    res_cand <- backsolve(chol_c_cand,res_cand,  transpose = TRUE)
    res_mat <- t(res_cand) %*% res_cand
    
    mean_i_array[,,i] <- mean_i_cand
    residual_array[,,i] <- res_mat
    val <- val + sum(init$Q_R[,,i] * res_mat)
  }
  
  val_2 <- -2 * k_l * sum(log(init$alphas)) - ((n * k_l) / 2) * log(det(init$Sigma)) - n*log(det(C_cand))
  log_lik_total <- -0.5 * val + val_2
  
  return(list(log_lik = log_lik_total, 
              mean = mean_cand, 
              mean_i = mean_i_array, 
              residual = residual_array,
              B_sim = B_cand,
              C = C_cand, 
              chol_c = chol_c_cand,
              mat_cand = mat_cand))
}


sample_theta_logit_normal <- function(init, hyper, X){
  # this function tries to implement the logit-normal reparametrization of the gaps  
  k_l <- init$k_l
  #beta_samp <- sqrt(init$tau) * hyper$U_lambda%*%beta_samp
  omega_samp <- c(rnorm(k_l), 0)
  u <- runif(n = 1)
  thresh <- hyper$log_lik + log(u)
  angle <- runif(n = 1)*2*pi
  angle_min <- angle - 2*pi
  angle_max <- angle
  omega_cand <- init$omega_theta * cos(angle) + omega_samp*sin(angle)
  gaps_cand <- exp(omega_cand)/(sum(exp(omega_cand))) 
  theta_cand <- cumsum(gaps_cand)[-c(k_l + 1)]*(2*pi)
  res_cand <-eval_log_lik_theta_joint(theta_cand, init, hyper, X) 
  while(res_cand$log_lik <=   thresh){
    #print(theta)
    #print(res_cand$log_lik)
    #print(thresh)
    if (angle >=  0) {
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
  
  hyper$log_lik <- res_cand$log_lik
  init$omega_theta <- omega_cand 
  init$thetas <- theta_cand
  init$residual <- res_cand$residual
  init$mean_i <- res_cand$mean_i
  init$mean <- res_cand$mean
  init$mat_dist <- res_cand$mat_cand
  init$C <- res_cand$C
  init$chol_c <- res_cand$chol_c
  hyper$B_sim <- res_cand$B_sim
  for (i in 1:k_l) {
    for (j in 1:k_l) {
      init$mat_dist[i, j] <- min(2*pi - abs(init$thetas[i] - init$thetas[j]),
                                  abs(init$thetas[i] - init$thetas[j]))
      init$mat_dist[j, i] <- init$mat_dist[i, j]
    }
  }
  init$C <- exp(-init$mat_dist/init$phi)

}
