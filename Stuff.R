
# Source file -----
source("src/Packages.R")
source("src/Utils.R")

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



n_points <- 25
# Create evenly spaced angles from 0 to almost 2*pi
thetas <- seq(0, 2 * pi, length.out = n_points + 1)[-(n_points + 1)]
# 
# B_sim <- Basis_Construction(empirical_thetas, L = L_intervals, degree = degree)
# tau = 0.2# smoothness
# K1 = K1_construction(n_basis)
# K2 = K2_construction(n_basis)
# #P = (1/(tau^2))*(K1 + 2*K2)
# P <- (1/(tau^2))*(K1 + 2*K2)

beta_values <- c(0.73331465, 0.27866933, -0.04889285, -0.28155420, -0.53854766, -0.46705397, 
                    -0.32019879, -0.11568639 ,-0.10122263, 0.29319813, 0.56797438)
n_points <- 25
# Create evenly spaced angles from 0 to almost 2*pi
degree = 3
sam <- in_model_sample(n = 10, K_l = length(thetas), thetas = thetas, n_int_knots = 8, degree = 3, tau = 0.2, 
                       beta_values = beta_values)

rot <- sam$mu
mu_i <- sam$mu_i[,,3]
X_i <- sam$X[,,3]
sam$Sigma_e
# rot <- mu_mean%*%R
par(mfrow = c(1, 3))
plot(rot, asp = 1, pch = 21, bg = "darkgreen", main = "Latent mu")
lines(c(rot[,1], rot[1,1]), c(rot[,2], rot[1,2]), col = "forestgreen", lwd = 2)

plot(mu_i, asp = 1, pch = 21, bg = "darkred", main = "mu for the i-th unit")
lines(c(mu_i[,1], mu_i[1,1]), c(mu_i[,2], mu_i[1,2]), col = "red", lwd = 2)

plot(X_i, asp = 1, pch = 21, bg = "darkblue", main = "Observed Conf for unit i")
lines(c(X_i[,1], X_i[1,1]), c(X_i[,2], X_i[1,2]), col = "blue", lwd = 2)
