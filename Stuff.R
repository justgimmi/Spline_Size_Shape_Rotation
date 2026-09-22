
# Source file -----
setwd("C:/Users/gmsan/Documenti/GitHub/Spline_Size_Shape_Rotation")
source("src/Packages.R")
source("src/Utils.R")

x <- seq(0.0001, 0.2, by = 0.0001)
library(MASS)
?manipulate
manipulate(curve(dinvgamma(x, a, b), from = 0.0001, to = 1), a = slider(0.00001, 10), b = slider(0.0001, 10))
plot(x, dinvgamma(x, 1, 0.01))
plot(curve(dinvgamma(x, 1, 0.01), from = 0.0001, to = 1))
# install.packages("R.matlab")
# 
# # Caricamento della libreria
library(R.matlab)
#
# # Lettura del file .mat
dati <- readMat("LANDMARKS.mat")
class <- readMat("class.mat")
class$class
dati$imm.name
dati$ll

plot(dati$Landmarks[,,20], type = "l")
plot(dati$Landmarks[,,20], asp = 1, pch = 21, bg = "darkgreen", main = "Unità 20")
lines(c(dati$Landmarks[,1,20], dati$Landmarks[1,1,20]), c(dati$Landmarks[,2,20], dati$Landmarks[1,2,20]), col = "forestgreen", lwd = 2)
text(ind[ord, ])
table(class$class)/3
ind <- dati$Landmarks[,,1]
ind <- sam$X[,,5]
r <- sqrt(ind[,1]^2 + ind[, 2]^2)
angles <- atan2(ind[,2], ind[,1])
ord <- order(angles)
angles
plot(angles, r)
angles_unwrapped <- angles
d <- diff(angles_unwrapped)
d[d >  pi] <- d[d >  pi] - 2*pi
d[d < -pi] <- d[d < -pi] + 2*pi
angles_unwrapped <- c(angles_unwrapped[1], angles_unwrapped[1] + cumsum(d))

plot(angles_unwrapped, type = "b", xlab = "landmark index (boundary order)", ylab = "unwrapped angle")
is_monotone <- all(diff(angles_unwrapped) > 0) || all(diff(angles_unwrapped) < 0)
is_monotone

plot(angles[ord], r[ord], type = "l")

# Try stuff ----
x <- runif(50) * 2* pi
x <- c(x)
par(mfrow = c(1,1))
Basis <- Basis_Construction(x, 17)
Basis |> matplot(type = "l") # L = 10



mat <- runif(n = 10, min = 0, max = 2*pi)
r <- runif(n = 10, min = 0, max = 10) 

mu = cbind(r*cos(mat), r*sin(mat))

plot(mu)
x <- -1000:1000
values <- x%%(2*pi)
values_bis <- x - 2*pi*floor(x/(2*pi))
m <- cbind(2*cos(values), 2*sin(values))
plot(x, values)
plot(m)
?floor

# data <- read.csv("src/saraghi_final_dataset.csv")
# data|>
#   filter(species == "D.sargus") -> Sargus_data
# final_array <- array(NA, dim = c(19, 2, nrow(Sargus_data)))
# for (i in 1:nrow(Sargus_data)) {
#   row_data <- unlist(Sargus_data[i, 8:45])
#   
#   final_array[,,i] <- matrix(row_data, nrow = 19, ncol = 2, byrow = TRUE)
# }
# library(geomorph)
# raw_fish <- final_array[,,1]
# y.gpa<-gpagen(final_array[,,1:3],ProcD = F)
# raw_fish <- y.gpa$consensus
# thetas_raw <- atan2(raw_fish[,2], raw_fish[,1])
# 
# thetas_raw <- ifelse(thetas_raw < 0, thetas_raw + 2*pi, thetas_raw)
# radii_raw  <- sqrt(raw_fish[,1]^2 + raw_fish[,2]^2)
# sort_idx         <- order(thetas_raw)
# empirical_thetas <- thetas_raw[sort_idx]
# empirical_radius <- radii_raw[sort_idx]        # Fixed: Now aligns with sorted thetas
# 
# fish_sorted    <- raw_fish[sort_idx, ]



# Simulation ------
# A <- hyper$Eigen_vector_null
# base <- t(A) %*% A
# eig <- eigen(base)
# A_bar <- t(eig$vectors[, eig$values < 1e-10])
# matrice <- rbind(A, A_bar)
# det(matrice)
# inv_matrice <- solve(matrice)
# A_bar %*% t(A)
n_points <- 21
# Create evenly spaced angles from 0 to almost 2*pi
thetas <- seq(0 + 0.0001, 2 * pi - 0.0001, length.out = n_points + 1)[-(n_points + 1)]
# 
# B_sim <- Basis_Construction(empirical_thetas, L = L_intervals, degree = degree)
# tau = 0.2# smoothness
# K1 = K1_construction(n_basis)
# K2 = K2_construction(n_basis)
# #P = (1/(tau^2))*(K1 + 2*K2)
# P <- (1/(tau^2))*(K1 + 2*K2)

#### almost zero warping ##### 
n_points <- 21
thetas <- seq(0 + 0.0001, 2*pi - 0.0001, length.out = n_points + 1)[-(n_points+1)]

sam <- in_model_sample(n = 100, K_l = n_points, thetas = thetas,
                        n_int_knots = 6, degree = 3,
                        tau = c(0.05, 0.05),   # tight-ish, smooth radial profile
                        phi = 0.5,
                        sigma_2 = 1e-8,        # effectively kills the warp
                        Sigma_e = diag(c(9e-4, 9e-4)))  # fixed, not riwish-drawn

sam <- in_model_sample(n = 100, K_l = n_points, thetas = thetas,
                        n_int_knots = 8, degree = 3,
                        tau = c(0.15, 0.35),   # loosen the spline
                        phi = 0.5,
                        sigma_2 = 0,        # same warp level as rung 2
                        Sigma_e = diag(c(9e-4, 9e-4)))


# beta_values <- c(0.73331465, 0.27866933, -0.04889285, -0.28155420, -0.53854766, -0.46705397, 
#                     -0.32019879, -0.11568639 ,-0.10122263, 0.29319813, 0.56797438)
# n_points <- 19
# # Create evenly spaced angles from 0 to almost 2*pi
# degree = 3
# n_points <- 21
# # Create evenly spaced angles from 0 to almost 2*pi
# thetas <- seq(0 + 0.0001, 2 * pi - 0.0001, length.out = n_points + 1)[-(n_points + 1)]
# sam <- in_model_sample(n = 100, K_l = length(thetas), thetas = thetas, n_int_knots = 14, degree = 3, tau = c(0.25, 0.3), phi = 0.5, sigma_2 = 1e-3,
#                        beta_values = NA)


rot <- sam$mu
mu_i <- sam$mu_i[,,1]
X_i <- sam$X[,,1]
sam$Sigma_e
sam$betas
# rot <- mu_mean%*%R
par(mfrow = c(1, 3))
plot(rot, asp = 1, pch = 21, bg = "darkgreen", main = "Latent mu")
lines(c(rot[,1], rot[1,1]), c(rot[,2], rot[1,2]), col = "forestgreen", lwd = 2)

plot(mu_i, asp = 1, pch = 21, bg = "darkred", main = "mu for the i-th unit")
lines(c(mu_i[,1], mu_i[1,1]), c(mu_i[,2], mu_i[1,2]), col = "red", lwd = 2)


plot(X_i, asp = 1, pch = 21, bg = "darkred", main = "mu for the i-th unit")
lines(c(X_i[,1], X_i[1,1]), c(X_i[,2], X_i[1,2]), col = "red", lwd = 2)
k = 18
k <- dim(sam$X)[1] # recuperiamo il numero di landmark
prova <- array(0, dim  = c(100, k, 2))
for (i in 1:100) {
  eta_matrix <- matrix(sam$eta[i, ], nrow = k, ncol = 2, byrow = TRUE)
  mu_i <- sam$X[,,i] - eta_matrix
  
  mu_i <- mu_i / sam$alphas[i]
  
  R_inv <- Rmat(sam$lambda[i])
  mu_i <- mu_i %*% R_inv
  
  # 4. PLOT
  prova[i,,] <- mu_i
  plot(mu_i, asp = 1, pch = 21, bg = "darkred", 
       main = paste("Mu ricostruito per l'unità", i),
       xlim = c(-2, 2), ylim = c(-1.5, 1.5)) 
  lines(c(mu_i[,1], mu_i[1,1]), c(mu_i[,2], mu_i[1,2]), col = "red", lwd = 2)
}
par(mfrow = c(1, 1))
mean_behav <- apply(prova, MARGIN = c(2, 3), FUN = mean)
plot(mean_behav, xlim = c(-1, 1.5), pch = 16, ylim = c(-0.9, 0.9))
#plot()
for (i in 1:100) {
  points(prova[i,,], asp = 1, pch = 16, col = "red")
  
}
points(apply(prova, MARGIN = c(2, 3), FUN = mean), asp = 1, pch = 16,
     xlim = c(-0.5, 0.5), ylim = c(-1, 1))
plot(apply(prova, MARGIN = c(2, 3), FUN = mean))
i = 1
sam$X - sam$mu_i
(X_i )
sum(sam$betas)
sam$betas
length(sam$betas)
R<- matrix(c(cos(sam$lambda[i]), -sin(sam$lambda[i]), sin(sam$lambda[i]), cos(sam$lambda[i])), 
                 byrow = T, nrow = 2, ncol = 2)
sam$alphas[i]^2 *(t(R) %*%sam$Sigma_e %*%R)

data <- read.csv("src/saraghi_final_dataset.csv")
data|>
  filter(species == "D.sargus", age == "Adult") -> Sargus_data

data|>
  filter(species == "D.sargus") -> Sargus_data
final_array <- array(NA, dim = c(19, 2, nrow(Sargus_data)))
for (i in 1:nrow(Sargus_data)) {
  row_data <- unlist(Sargus_data[i, 8:45])

  final_array[,,i] <- matrix(row_data, nrow = 19, ncol = 2, byrow = TRUE)
}
sam <- final_array



k <- 19
prova <- array(0, dim  = c(120, k, 2))
for (i in 1:120) {
  eta_matrix <- matrix(init_env$eta[i, ], nrow = k, ncol = 2, byrow = TRUE)
  mu_i <- X[,,i] - eta_matrix
  
  mu_i <- mu_i / init_env$alphas[i]
  
  R_inv <- Rmat(init_env$lambdas[i])
  mu_i <- mu_i %*% R_inv
  
  # 4. PLOT
  prova[i,,] <- mu_i
  plot(mu_i, asp = 1, pch = 21, bg = "darkred", 
       main = paste("Mu ricostruito per l'unità", i),
       xlim = c(-2, 2), ylim = c(-1.5, 1.5)) 
  lines(c(mu_i[,1], mu_i[1,1]), c(mu_i[,2], mu_i[1,2]), col = "red", lwd = 2)
}
par(mfrow = c(1, 1))
mean_behav <- apply(prova, MARGIN = c(2, 3), FUN = mean)
plot(mean_behav, xlim = c(-2, 1.5), pch = 16, ylim = c(-1, 1), asp = 1)
#plot()
for (i in 1:120) {
  points(prova[i,,], asp = 1, pch = 16, col = "red")
  
}
points(apply(prova, MARGIN = c(2, 3), FUN = mean), asp = 1, pch = 16,
       xlim = c(-0.5, 0.5), ylim = c(-1, 1))
plot(apply(prova, MARGIN = c(2, 3), FUN = mean))

load("Sim2.RData")

X_i <- final_array[,,1]
ind <- final_array[,,1]
r <- sqrt(ind[,1]^2 + ind[, 2]^2)
angles <- atan2(ind[,2], ind[,1])
angles
par(mfrow = c(1,1))
plot(angles, r)
points(angles, r, type = "l")



cc <- c(1, 1)
S = matrix(NA, nrow = nrow(X_i), ncol = nrow(X_i))
X_i <-  final_array[,,1]
for (i in 1:nrow(X_i)) {
  for (j in 1:nrow(X_i)) {
    normm <- sum((X_i[i, ] - X_i[j, ])^2)
    S[i, j] <- ifelse(normm  > 0,normm*log(sqrt(normm)) ,0)
    
  }
  
}

plot(S)
##### new ####
k <- nrow(X_i)

# TPS kernel matrix (this part was already correct)
S <- matrix(0, nrow = k, ncol = k)
for (i in 1:k) {
  for (j in 1:k) {
    normm <- sum((X_i[i, ] - X_i[j, ])^2)   # = r^2
    S[i, j] <- ifelse(normm > 0, normm * log(sqrt(normm)), 0)  # r^2 log(r)
  }
}

# Constraint matrix: enforce 1_k' W = 0 and X_i' W = 0
A_constr <- rbind(matrix(1, nrow = 1, ncol = k), t(X_i))   # 3 x k
mat <- t(A_constr)%*%A_constr
eigs <- eigen(mat)
eg_vec <- eigs$vectors[,eigs$values< 1e-10, drop = FALSE] 

matts <- rbind(A_constr, t(eg_vec))
inv <- solve(matts)
inv

C_bar_w <- inv[, (nrow(A_constr) + 1):k]   # k x (k-3), orthonormal basis of ker(A_constr)
C_bar_w <- eg_vec
# Free unconstrained parameter
0.03^2
L <- t(chol(t(C_bar_w)%*% C_bar_w))
W1 <- solve(L, rnorm((k - 3), sd = 0.001))
W2 <- solve(L, rnorm((k - 3), sd = 0.001))
gamma_W <- cbind(W1, W2)

W <- C_bar_w %*% gamma_W   # k x 2, guaranteed to satisfy both constraints
apply(W, MARGIN = 2, FUN = sum)
cc <- c(1, 1)
A  <- matrix(rnorm(4, sd = 1), nrow = 2)
ones <- matrix(1, nrow = k)

Y <- S %*% W + X_i

par(mfrow = c(1, 2))
r <- sqrt(Y[,1]^2 + Y[,2]^2)
angles <- atan2(Y[,2], Y[,1])
ord <- order(angles)
plot(angles[ord], r[ord], type = "b")
plot(Y[ord, ], type = "l")
text(Y[ord, ])
Y_ord <- Y[ord, ]
plot(Y_ord,asp = 1, pch = 21, bg = "darkgreen", main = "Non star-shaped")
lines(c(Y_ord[,1], Y_ord[1,1]), 
      c(Y_ord[,2], Y_ord[1,2]), col = "forestgreen", lwd = 2)
text(Y[ord, ])
X_ord = X_i[ord, ]
plot(X_ord,asp = 1, pch = 21, bg = "darkgreen", main = "Mean star-shaped")
lines(c(X_i[,1], X_i[1,1]), 
      c(X_i[,2], X_i[1,2]), col = "forestgreen", lwd = 2)
text(X_ord)



angles_unwrapped <- angles
d <- diff(angles_unwrapped)
d[d >  pi] <- d[d >  pi] - 2*pi
d[d < -pi] <- d[d < -pi] + 2*pi
angles_unwrapped <- c(angles_unwrapped[1], angles_unwrapped[1] + cumsum(d))

plot(angles_unwrapped, type = "b", xlab = "landmark index (boundary order)", ylab = "unwrapped angle")
is_monotone <- all(diff(angles_unwrapped) > 0) || all(diff(angles_unwrapped) < 0)
is_monotone


###### effect of S ######

# ---- fix: compute ord/r/angles right after Y is built, before any indexing uses them ----
r <- sqrt(Y[,1]^2 + Y[,2]^2)
angles <- atan2(Y[,2], Y[,1])
ord <- order(angles)

# ---- TPS evaluation at arbitrary query points (not just the k landmarks) ----
tps_warp <- function(query, mu_star, W, cc = c(0,0), Aff = diag(2)) {
  # query: m x 2 matrix of points to warp
  # mu_star: k x 2 control points (the star-shaped template landmarks)
  # W: k x 2 non-affine warp coefficients (already constraint-satisfying)
  m <- nrow(query)
  k <- nrow(mu_star)
  S_query <- matrix(0, nrow = m, ncol = k)
  for (a in 1:m) {
    for (b in 1:k) {
      d2 <- sum((query[a, ] - mu_star[b, ])^2)
      S_query[a, b] <- ifelse(d2 > 0, d2 * log(sqrt(d2)), 0)
    }
  }
  S_query %*% W + query
}

# ---- build a regular grid covering the shape's bounding box (with margin) ----
rng_x <- range(X_i[,1]); rng_y <- range(X_i[,2])
pad <- 0.25 * max(diff(rng_x), diff(rng_y))
gx <- seq(rng_x[1] - pad, rng_x[2] + pad, length.out = 25)
gy <- seq(rng_y[1] - pad, rng_y[2] + pad, length.out = 25)
grid_pts <- as.matrix(expand.grid(x = gx, y = gy))

grid_warped <- tps_warp(grid_pts, X_i, W)   # X_i is your mu_star here

# ---- plot: original grid+shape vs warped grid+shape, side by side ----
par(mfrow = c(1, 2))
X_ord <- X_i[ord, ]
plot(grid_pts, pch = ".", col = "grey70", asp = 1, main = "Star Shaped Object",
     xlab = "", ylab = "")
for (xi in gx) lines(rep(xi, length(gy)), gy, col = "grey70")
for (yi in gy) lines(gx, rep(yi, length(gx)), col = "grey70")
points(X_i, pch = 21, bg = "darkgreen")
lines(c(X_i[,1], X_i[1,1]), 
      c(X_i[,2], X_i[1,2]), col = "forestgreen", lwd = 2)
text(X_ord, pos = 2)

#lines(X_i[c(ord, ord[1]), ], col = "forestgreen", lwd = 2)
Y_ord <- Y[ord, ]
plot(grid_warped, pch = ".", col = "grey70", asp = 1, main = "Warped grid of (S %*% W)",
     xlab = "", ylab = "")
grid_warped_mat <- matrix(grid_warped, ncol = 2)
gw_x <- matrix(grid_warped_mat[,1], nrow = length(gx), ncol = length(gy))
gw_y <- matrix(grid_warped_mat[,2], nrow = length(gx), ncol = length(gy))
for (j in 1:length(gy)) lines(gw_x[, j], gw_y[, j], col = "grey70")
for (i in 1:length(gx)) lines(gw_x[i, ], gw_y[i, ], col = "grey70")
points(Y, pch = 21, bg = "darkgreen")
lines(Y[c(ord, ord[1]), ], col = "forestgreen", lwd = 2)
text(Y_ord, pos = 2)
