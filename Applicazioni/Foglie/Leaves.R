###### Source file ###### 
require(here)
source(here("src/Packages.R"))
source(here("src/Utils.R"))
source(here("src/RB_MCMC.R"))
set.seed(250601)
# import data
data <- read.csv(file = "df_landmarks.csv")
n <- length(unique(data$treeID)) # number of unique trees
total_n <- 94820/20

total_leaves <- array(NA, dim = c(20, 2,94820/20 ))
type_of_tree <- numeric(94820/20)
for (i in 1:dim(total_leaves)[3]){
  total_leaves[,,i] <- as.matrix(data[((i - 1)*20 + 1):(i*20), c("x_bar", "y_bar")])
  type_of_tree[i] <- as.vector(data[(i - 1)*20 + 1, c("treeID")])
}

table(data$treeID)/20 # we have several leaves for each of the trees
sum(table(data$treeID))/20
leaf_df <- data.frame(idx = seq_along(type_of_tree), treeID = type_of_tree)
s <- leaf_df %>% group_by(treeID) %>% slice_sample(prop = 500 / nrow(leaf_df))

final_array <- total_leaves[,,s$idx]
type_of_tree_sampled <- type_of_tree[s$idx]
prop.table(table(type_of_tree_sampled))
prop.table(table(data$treeID)/20)

# save(total_leaves, final_array, type_of_tree_sampled,file = here("Leaves.RData"))
# load(here("Leaves.RData"))
prova <- init_mu_star(final_array,n_basis = 18, degree = 3, start = "circle", ridge   = 1e-2, plot = TRUE,
                      Sigma_eta = diag(10, 2)) # we start from a circle

init <- prova$init # define init object 
hyper <- prova$hyper # define hyper object 
X <- prova$X # define X 
MCMC_leaves <- RB_MCMC(total_iter = 2e5, burnin = 1e5, thinning  = 100, X = X, 
                      init = init, hyper = hyper) # it takes a while, if you trust me you can find everything on the 
# RData, otherwise you can easily wait.
