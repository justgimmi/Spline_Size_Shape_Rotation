require(here)
source(here("src/Packages.R"))
source(here("src/Utils.R"))
source(here("src/RB_MCMC.R"))
data <- read.csv(here("Applicazioni/Pesci/saraghi_final_dataset.csv"))
data|>
  filter(species == "D.sargus") -> Sargus_data

final_array <- array(NA, dim = c(19, 2, nrow(Sargus_data)))
for (i in 1:nrow(Sargus_data)) {
  row_data <- unlist(Sargus_data[i, 8:45])
  
  final_array[,,i] <- matrix(row_data, nrow = 19, ncol = 2, byrow = TRUE)
}
sam <- final_array
save(sam, file = here("Applicazioni/Pesci/Fishes.RData"))
prova <- init_mu_star(sam,n_basis = 16, degree = 3, start = "circle", ridge   = 1e-2, plot = TRUE,
                      Sigma_eta = diag(10, 2)) # we start from the gpa mean

init <- prova$init # define init object 
hyper <- prova$hyper # define hyper object 
X <- prova$X # define X 
MCMC_fish <- RB_MCMC(total_iter = 10000, burnin = 1000, thinning  = 10, X = X, 
                       init = init, hyper = hyper) 
