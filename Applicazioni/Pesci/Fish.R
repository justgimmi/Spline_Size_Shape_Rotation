### ---- Import data ----- #### 
require(here)
source(here("src/Packages.R"))
source(here("src/Utils.R"))
source(here("src/RB_MCMC.R"))
data <- read.csv(here("Applicazioni/Pesci/saraghi_final_dataset.csv"))
data|>
  filter(species == "D.sargus") -> Sargus_data
unique(data$species)
final_array <- array(NA, dim = c(19, 2, nrow(Sargus_data)))
for (i in 1:nrow(Sargus_data)) {
  row_data <- unlist(Sargus_data[i, 8:45])
  
  final_array[,,i] <- matrix(row_data, nrow = 19, ncol = 2, byrow = TRUE)
}
sam <- final_array
# save(sam, file = here("Applicazioni/Pesci/Sargus/Fishes.RData"))
prova <- init_mu_star(sam,n_basis = 16, degree = 3, start = "circle", ridge   = 1e-2, plot = TRUE,
                      Sigma_eta = diag(10, 2)) # we start from the gpa mean

### ---- Fit the model ----- ###### 
init <- prova$init # define init object 
hyper <- prova$hyper # define hyper object 
X <- prova$X # define X 
MCMC_fish <- RB_MCMC(total_iter = 10000, burnin = 1000, thinning  = 10, X = X, 
                       init = init, hyper = hyper) 

load(file = here("Applicazioni/Pesci/Sargus/Fishes_MCMC.RData"))
output <- MCMC_samp$output
hypers <- MCMC_samp$hyper
try(visualize_results(output, hypers, X, "Applicazioni/Pesci/Sargus/Fishes"))

# how about the different sizes of the fishes
identifier <- Sargus_data$age # Adult or Juvenile
alphas_identified <- t(sapply(1:nrow(output$alphas), function(i){
  log(output$alphas[i,]) - mean(log(output$alphas[i,]))}
))
dim(output$alphas)
dim(alphas_identified)
Adult_alphas <- alphas_identified[,identifier == "Adult"]
Juvenile_alphas <- alphas_identified[,identifier != "Adult"]
# If you want to visualize that those chains are identified, you can look up at this 
# show(gg_mcmc_diagnostics(Adult_alphas, param_name = "alphas", real_values = NA, TRUE))
# show(gg_mcmc_diagnostics(Juvenile_alphas, param_name = "alphas", real_values = NA, TRUE))
df_all <- data.frame(
  alpha = c(as.vector(Adult_alphas), as.vector(Juvenile_alphas)),
  group = rep(c("Adult", "Juvenile"),
              times = c(length(Adult_alphas), length(Juvenile_alphas)))
)

p <- ggplot(df_all, aes(x = alpha, fill = group, colour = group)) +
  geom_density(alpha = 0.35, adjust = 1.2) +
  labs(x = expression(log(alpha[i]) - mean(log(alpha))),
       y = "Density", fill = NULL, colour = NULL) +
  theme_minimal()
# ggsave(filename = here("Applicazioni/Pesci/Sargus/Fishes_Sizes.png"), dpi = 600, width = 18, height = 10, unit = "cm")

# it is quite clear that there is an effective separation between adults and juveniles in the Sargus data 


###### On the other hand, we would like to reproduce the analysis of those papers ######
# Does morphology predict trophic niche differentiation? Relationship 
# between feeding habits and body shape in four co-occurring juvenile species (Pisces: Perciformes, Sparidae)

# Are developmental shifts the main driver of phenotypic evolution in Diplodus spp. (Perciformes: Sparidae)?
# we will fit the model for all the 4 species in the juvenile regime
rm(list = ls()) # clean all the environment
source(here("src/Utils.R"))
source(here("src/RB_MCMC.R"))
data <- read.csv(here("Applicazioni/Pesci/saraghi_final_dataset.csv"))
unique(data$species)
data|>
  filter(species == "D.sargus" & age == "Juvenile") -> Sargus_data
final_array_sargus <- array(NA, dim = c(19, 2, nrow(Sargus_data)))
for (i in 1:nrow(Sargus_data)) {
  row_data <- unlist(Sargus_data[i, 8:45])
  
  final_array_sargus[,,i] <- matrix(row_data, nrow = 19, ncol = 2, byrow = TRUE)
}
sam <- final_array_sargus
# save(sam, file = here("Applicazioni/Pesci/Juveniles/Sargus.RData"))

prova <- init_mu_star(sam,n_basis = 16, degree = 3, start = "circle", ridge   = 1e-2, plot = TRUE,
                      Sigma_eta = diag(10, 2)) # we start from the gpa mean
data|>
  filter(species == "D.puntazzo" & age == "Juvenile") -> Puntazzo_data

final_array_puntazzo <- array(NA, dim = c(19, 2, nrow(Puntazzo_data)))
for (i in 1:nrow(Puntazzo_data)) {
  row_data <- unlist(Puntazzo_data[i, 8:45])
  
  final_array_puntazzo[,,i] <- matrix(row_data, nrow = 19, ncol = 2, byrow = TRUE)
}
sam <- final_array_puntazzo
# save(sam, file = here("Applicazioni/Pesci/Juveniles/Puntazzo.RData"))
prova <- init_mu_star(sam,n_basis = 16, degree = 3, start = "circle", ridge   = 1e-2, plot = TRUE,
                      Sigma_eta = diag(10, 2)) # we start from the gpa mean

data|>
  filter(species == "D.annularis" & age == "Juvenile") -> Annularis_data

final_array_annularis <- array(NA, dim = c(19, 2, nrow(Annularis_data)))
for (i in 1:nrow(Annularis_data)) {
  row_data <- unlist(Annularis_data[i, 8:45])
  
  final_array_annularis[,,i] <- matrix(row_data, nrow = 19, ncol = 2, byrow = TRUE)
}
sam <- final_array_annularis
# save(sam, file = here("Applicazioni/Pesci/Juveniles/Annularis.RData"))
prova <- init_mu_star(sam,n_basis = 16, degree = 3, start = "circle", ridge   = 1e-2, plot = TRUE,
                      Sigma_eta = diag(10, 2)) # we start from the gpa mean
data|>
  filter(species == "D.vulgaris" & age == "Juvenile") -> Vulgaris_data

final_array_vulfaris <- array(NA, dim = c(19, 2, nrow(Vulgaris_data)))
for (i in 1:nrow(Vulgaris_data)) {
  row_data <- unlist(Vulgaris_data[i, 8:45])
  
  final_array_vulfaris[,,i] <- matrix(row_data, nrow = 19, ncol = 2, byrow = TRUE)
}
sam <- final_array_vulfaris
# save(sam, file = here("Applicazioni/Pesci/Juveniles/Vulgaris.RData"))
prova <- init_mu_star(sam,n_basis = 16, degree = 3, start = "circle", ridge   = 1e-2, plot = TRUE,
                      Sigma_eta = diag(10, 2)) # we start from the gpa mean



# init <- prova$init # define init object 
# hyper <- prova$hyper # define hyper object 
# X <- prova$X # define X 
# MCMC_fish <- RB_MCMC(total_iter = 10000, burnin = 1000, thinning  = 10, X = X, 
#                      init = init, hyper = hyper) 


##### I need also the adult to be included ######## 
rm(list = ls()) # clean all the environment
source(here("src/Utils.R"))
source(here("src/RB_MCMC.R"))
data <- read.csv(here("Applicazioni/Pesci/saraghi_final_dataset.csv"))
unique(data$species)
data|>
  filter(species == "D.sargus" & age == "Adult") -> Sargus_data
final_array_sargus <- array(NA, dim = c(19, 2, nrow(Sargus_data)))
for (i in 1:nrow(Sargus_data)) {
  row_data <- unlist(Sargus_data[i, 8:45])
  
  final_array_sargus[,,i] <- matrix(row_data, nrow = 19, ncol = 2, byrow = TRUE)
}
sam <- final_array_sargus
save(sam, file = here("Applicazioni/Pesci/Adults/Adults_Sargus.RData"))

prova <- init_mu_star(sam,n_basis = 16, degree = 3, start = "circle", ridge   = 1e-2, plot = TRUE,
                      Sigma_eta = diag(10, 2)) # we start from the gpa mean
data|>
  filter(species == "D.puntazzo" & age == "Adult") -> Puntazzo_data

final_array_puntazzo <- array(NA, dim = c(19, 2, nrow(Puntazzo_data)))
for (i in 1:nrow(Puntazzo_data)) {
  row_data <- unlist(Puntazzo_data[i, 8:45])
  
  final_array_puntazzo[,,i] <- matrix(row_data, nrow = 19, ncol = 2, byrow = TRUE)
}
sam <- final_array_puntazzo
save(sam, file = here("Applicazioni/Pesci/Adults/Adults_Puntazzo.RData"))
prova <- init_mu_star(sam,n_basis = 16, degree = 3, start = "circle", ridge   = 1e-2, plot = TRUE,
                      Sigma_eta = diag(10, 2)) # we start from the gpa mean

data|>
  filter(species == "D.annularis" & age == "Adult") -> Annularis_data

final_array_annularis <- array(NA, dim = c(19, 2, nrow(Annularis_data)))
for (i in 1:nrow(Annularis_data)) {
  row_data <- unlist(Annularis_data[i, 8:45])
  
  final_array_annularis[,,i] <- matrix(row_data, nrow = 19, ncol = 2, byrow = TRUE)
}
sam <- final_array_annularis
save(sam, file = here("Applicazioni/Pesci/Adults/Adults_Annularis.RData"))
prova <- init_mu_star(sam,n_basis = 16, degree = 3, start = "circle", ridge   = 1e-2, plot = TRUE,
                      Sigma_eta = diag(10, 2)) # we start from the gpa mean
data|>
  filter(species == "D.vulgaris" & age == "Adult") -> Vulgaris_data

final_array_vulfaris <- array(NA, dim = c(19, 2, nrow(Vulgaris_data)))
for (i in 1:nrow(Vulgaris_data)) {
  row_data <- unlist(Vulgaris_data[i, 8:45])
  
  final_array_vulfaris[,,i] <- matrix(row_data, nrow = 19, ncol = 2, byrow = TRUE)
}
sam <- final_array_vulfaris
save(sam, file = here("Applicazioni/Pesci/Adults/Adults_Vulgaris.RData"))
prova <- init_mu_star(sam,n_basis = 16, degree = 3, start = "circle", ridge   = 1e-2, plot = TRUE,
                      Sigma_eta = diag(10, 2)) # we start from the gpa mean



# init <- prova$init # define init object 
# hyper <- prova$hyper # define hyper object 
# X <- prova$X # define X 
# MCMC_fish <- RB_MCMC(total_iter = 10000, burnin = 1000, thinning  = 10, X = X, 
#                      init = init, hyper = hyper) 

### ---- Let's analyze the results among fishes ---- ####
# Juveniles 
load(here("Applicazioni/Pesci/Juveniles/Annularis/Annularis_MCMC.RData"))
MCMC_annularis_juvenile <- MCMC_samp
output <- MCMC_samp$output
hypers <- MCMC_samp$hyper
try(visualize_results(output, hypers, X, "Applicazioni/Pesci/Juveniles/Annularis/Annularis_Juvenile"))

load(here("Applicazioni/Pesci/Juveniles/Vulgaris_MCMC.RData"))
MCMC_vulgaris_juvenile <- MCMC_samp
output <- MCMC_samp$output
hypers <- MCMC_samp$hyper
try(visualize_results(output, hypers, X, "Applicazioni/Pesci/Juveniles/Vulgaris_Juvenile"))


load(here("Applicazioni/Pesci/Juveniles/Sargus/Sargus_MCMC.RData"))
MCMC_sargus_juvenile <- MCMC_samp
output <- MCMC_samp$output
hypers <- MCMC_samp$hyper
try(visualize_results(output, hypers, X, "Applicazioni/Pesci/Juveniles/Sargus/Sargus_Juvenile"))

load(here("Applicazioni/Pesci/Juveniles/Puntazzo/Puntazzo_MCMC.RData"))
MCMC_puntazzo_juvenile <- MCMC_samp
output <- MCMC_samp$output
hypers <- MCMC_samp$hyper
try(visualize_results(output, hypers, X, "Applicazioni/Pesci/Juveniles/Puntazzo/Puntazzo_Juvenile"))


load(here("Applicazioni/Pesci/Adults/Puntazzo/Puntazzo_Adults_MCMC.RData"))
MCMC_puntazzo_Adults <- MCMC_samp
output <- MCMC_samp$output
hypers <- MCMC_samp$hyper
try(visualize_results(output, hypers, X, "Applicazioni/Pesci/Adults/Puntazzo/Puntazzo_Adults"))


load(here("Applicazioni/Pesci/Adults/Sargus/Sargus_Adults_MCMC.RData"))
MCMC_sargus_Adults <- MCMC_samp
output <- MCMC_samp$output
hypers <- MCMC_samp$hyper
try(visualize_results(output, hypers, X, "Applicazioni/Pesci/Adults/Sargus/Sargus_Adults"))

load(here("Applicazioni/Pesci/Adults/Annularis/Annularis_Adults_MCMC.RData"))
MCMC_Annularis_Adults <- MCMC_samp
output <- MCMC_samp$output
hypers <- MCMC_samp$hyper
try(visualize_results(output, hypers, X, "Applicazioni/Pesci/Adults/Annularis/Annularis_Adults"))


load(here("Applicazioni/Pesci/Adults/Vulgaris/Vulgaris_Adults_MCMC.RData"))
MCMC_Vulgaris_Adults <- MCMC_samp
output <- MCMC_samp$output
hypers <- MCMC_samp$hyper
try(visualize_results(output, hypers, X, "Applicazioni/Pesci/Adults/Vulgaris/Vulgaris_Adults"))
# To understand if there are any differences, we should first of all identify the mean w.r.t the same
# pairs of landmarks 



compute_identified_mean <- function(output, pair){
  nsamp <- dim(output$lambdas)[1]
  k <- dim(output$theta)[2]
  n <- dim(output$lambdas)[2]
  mean_orig <- array(NA, c(k, 2, nsamp))
  mean_noloc <- array(NA, c(k, 2, nsamp))
  mean_noloc_i <- array(NA, c(k, 2, nsamp))
  mean_norot <- array(NA, c(k, 2, nsamp))
  mean_norot_i <- array(NA, c(k, 2, nsamp))
  mean_nosize <- array(NA, c(k, 2, nsamp))
  mean_nosize_i <- array(NA, c(k, 2, nsamp))
  # warp_mean <- array(NA, c(k, 2, n))
  warping_i <- array(NA, dim = c(k,2,nsamp))
  
  warp_mean <- list()
  for (iobs in 1:n) {
    
    # 1. Feedback reale sulla console R
    cat(sprintf("\rElaborazione Unità: %3d", iobs))
    flush.console()
    
    krot <- pair[2]
    for (isamp in 1:nsamp){
      mean_orig[, , isamp] <- output$mean[isamp, ,] 
      S <- S_mat_constructor(mean_orig[, , isamp])
      Wbar <- apply(output$W[[isamp]], c(1, 2), mean) 
      SWbar <- matrix(apply(S %*% Wbar, MARGIN = 2, FUN = "mean"), ncol = 2, nrow = k, byrow = T)
      # Remove the mean of the Warping for every unit
      warping <- S%*%output$W[[isamp]][,,iobs] - SWbar
      # Remove Location effect 
      mean_noloc[, , isamp] <- mean_orig[, , isamp] - matrix(mean_orig[pair[1], , isamp], nrow = k, ncol = 2, byrow = T) + SWbar
      warping <- warping
      mean_noloc_i[, , isamp] <- mean_noloc[, , isamp] + warping
      
      # Remove Rotation effect 
      angle <- atan2(mean_noloc[krot, 2, isamp], mean_noloc[krot, 1, isamp])
      mean_norot[, , isamp] <- mean_noloc[, , isamp] %*% Rmat(angle)
      warping <- warping %*% Rmat(angle)
      mean_norot_i[, , isamp] <- mean_norot[, , isamp] +  warping
      
      # Remove scale effect
      normx <- norm(sweep(mean_orig[, , isamp], 2, colMeans(mean_orig[, , isamp])), "F")
      mean_nosize[, , isamp] <- mean_norot[, , isamp] / normx
      warping_i[,,isamp] <- warping/normx
      mean_nosize_i[, , isamp] <- mean_nosize[, , isamp] + warping_i[,,isamp]
    }
    
    warp_mean[[iobs]] <- warping_i
  }
  
  return(list(mean = mean_nosize , warping = warp_mean ))
}
Annularis_adult <- compute_identified_mean(MCMC_Annularis_Adults$output, c(9, 19))
Annularis_juvenile <- compute_identified_mean(MCMC_annularis_juvenile$output, c(9, 19))
Sargus_adult <- compute_identified_mean(MCMC_sargus_Adults$output, c(9, 19))
Sargus_juvenile <- compute_identified_mean(MCMC_sargus_juvenile$output, c(9, 19))
Vulgaris_adult <- compute_identified_mean(MCMC_Vulgaris_Adults$output, c(9, 19))
Vulgaris_juvenile <- compute_identified_mean(MCMC_vulgaris_juvenile$output, c(9, 19))
Puntazzo_adult <- compute_identified_mean(MCMC_puntazzo_Adults$output, c(9, 19))
Puntazzo_juvenile <- compute_identified_mean(MCMC_puntazzo_juvenile$output, c(9, 19))
save(Annularis_adult, Annularis_juvenile, Sargus_adult, Sargus_juvenile, Vulgaris_adult, 
     Vulgaris_juvenile, Puntazzo_adult, Puntazzo_juvenile, file = here("Applicazioni/Pesci/Mean_Identified.RData"))

# Let's compare the mean 

# long data frame of the closed polygons; amp = 1 gives the plain means,
# amp > 1 amplifies the deviation from the grand mean
build_df <- function(amp = 1) {
  do.call(rbind, lapply(1:8, function(g) {
    M <- grand + amp * (Mrot[, , g] - grand)
    M <- rbind(M, M[1, ])                                  # close the contour
    data.frame(x = M[, 1], y = M[, 2], lm = c(1:k, 1),
               group   = names(groups)[g],
               species = factor(sp[g],  levels = names(cols)),
               age     = factor(age[g], levels = c("adult", "juvenile")))
  }))
}

plot_means <- function(df, title) {
  ggplot(df, aes(x, y, colour = species, linetype = age, group = group)) +
    geom_path(linewidth = 0.9) +
    geom_point(size = 0.8) +
    scale_colour_manual(values = cols, drop = FALSE) +
    scale_linetype_manual(values = c(adult = "solid", juvenile = "dashed"), drop = FALSE) +
    coord_equal(xlim = rg, ylim = rg) +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))
}

make_three <- function(amp = 1, label = "Posterior mean shapes") {
  df <- build_df(amp)
  p_juv <- plot_means(subset(df, age == "juvenile"), "Juveniles")
  p_ad  <- plot_means(subset(df, age == "adult"),    "Adults")
  p_all <- plot_means(df,                            "Combined")
  (p_juv + p_ad + p_all) + plot_layout(guides = "collect") +
    plot_annotation(title = label)
}

ggsave("Applicazioni/Pesci/means_overlay.png", make_three(1), width = 16, height = 5.5, dpi = 600)





idx_ad <- c(1, 3, 5, 7) # juvenile index
idx_juv <- c(2, 4, 6, 8) # adult index
nd <- min(sapply(groups, function(g) dim(g$mean)[3])) # number of mcmc iterations
Dd <- array(0, c(8, 8, nd))
for (t in 1:nd) for (a in 1:7) for (b in (a + 1):8)
  Dd[a, b, t] <- Dd[b, a, t] <- sin(riemdist(groups[[a]]$mean[, , t], groups[[b]]$mean[, , t])) # compute
# the full Procrustes distance between each pair of means across MCMC_sample
stage  <- sapply(1:4, function(s) Dd[idx_ad[s], idx_juv[s], ]) # here we just consider pairwise couple, namely inside the same species
# stage represents the within mean 
ub     <- function(M) mean(M[upper.tri(M)])
sp_ad  <- apply(Dd[idx_ad,  idx_ad,  , drop = FALSE], 3, ub) # average distance across species within same juvenile stage
sp_juv <- apply(Dd[idx_juv, idx_juv, , drop = FALSE], 3, ub) # average distance across species within same adult stage
ratio  <- rowMeans(stage) / ((sp_ad + sp_juv) / 2)          # stage effect >> species effect
quantile(ratio, c(.025, .5, .975))
# on average the differences within the same species are double the differences across different species
colMeans(stage) # this tell us which fishes have the greater change across stages
colnames(stage) <- c("Annularis", "Sargus", "Vulgaris", "Puntazzo")
rbind(mean = colMeans(stage), apply(stage, 2, quantile, c(.025, .975)))

# which species has the larger change? posterior probability for each pair
pp <- sapply(combn(4, 2, simplify = FALSE), function(p) mean(stage[, p[1]] > stage[, p[2]]))
names(pp) <- combn(colnames(stage), 2, paste, collapse = " > ")
round(pp, 2)
