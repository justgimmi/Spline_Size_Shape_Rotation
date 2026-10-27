##### Import data ###### 
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


###### Fit the model ##### 
init <- prova$init # define init object 
hyper <- prova$hyper # define hyper object 
X <- prova$X # define X 
MCMC_fish <- RB_MCMC(total_iter = 10000, burnin = 1000, thinning  = 10, X = X, 
                       init = init, hyper = hyper) 

load(file = here("Applicazioni/Pesci/Fishes_MCMC.RData"))
output <- MCMC_samp$output
hypers <- MCMC_samp$hyper
try(visualize_results(output, hypers, X, "Applicazioni/Pesci/Fishes"))

# how about the different sizes of the fishes
identifier <- Sargus_data$age # Adult or Juvenile
alphas_identified <- t(sapply(1:nrow(output$alphas), function(i){
  log(output$alphas[i,]) - mean(log(output$alphas[i,]))}
))
dim(output$alphas)
dim(alphas_identified)
Adult_alphas <- alphas_identified[,identifier == "Adult"]
Juvenile_alphas <- alphas_identified[,identifier != "Adult"]

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
ggsave(filename = here("Applicazioni/Pesci/Fishes_Sizes.png"), dpi = 600, width = 18, height = 10, unit = "cm")
