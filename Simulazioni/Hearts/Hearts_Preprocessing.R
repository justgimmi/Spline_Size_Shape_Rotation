require(here)
source(here("src/Packages.R"))
source(here("src/Utils.R"))
set.seed(1984)
### ---- If you want to visualize ---- ####

hearts %T>%                    # A toy dataset
  stack() %>%                  # Take a family picture of raw outlines
  fgProcrustes() %>%           # Full generalized Procrustes alignment
  coo_slide(ldk = 2) %T>%      # Redefine a robust 1st point between the cheeks
  stack() %>%                  # Another picture of aligned outlines
  efourier(6, norm=FALSE) %>%  # Elliptical Fourier Transforms
  PCA() %T>%                   # Principal Component Analysis
  plot_PCA(~aut) %>%           # A PC1:2 plot
  LDA(~aut) %>%                # Linear Discriminant Analysis
  plot_CV()                    # And the confusion matrix after leave one out cross validation

### ---- Find the mean shape configuration ----- ##### 

dati <- hearts # hearts dataset inside Momocs library
hearts_list <- unlist(dati, recursive = FALSE) # unlist data
hearts_list <- hearts_list[sapply(hearts_list, is.matrix)] # transform each of them into matrix

X <- simplify2array(hearts_list) # stack into a k x 2 x n array

gpa <- procGPA(X, scale = TRUE, reflect = FALSE)
mean_shape_full <- gpa$mshape
mean_shape <- gpa$mshape          # k x 2 Procrustes mean shape --> we consider it as an average configuration
par(mfrow = c(1,1))
plotshapes(mean_shape, joinline = 1:nrow(mean_shape))
M <- mean_shape           
k <- nrow(M) # number of landmakrs (80 too many =( )

# Look at the landmark numbers and let's pick some anchors
plot(M, type = "b", asp = 1) 
text(M, labels = 1:k, pos = 3, cex = 0.6)
# we argue that two important landmarks are 59 and 18

top_idx <- 59
bot_idx <- 18


# Select equispaced landmarks between the two anchors
sel <- pick_landmarks(M, top_idx, bot_idx, n_side = 11)   # 24 landmarks

plot(M, type = "l", asp = 1)
points(M[sel, ], pch = 19, col = "red")
text(M[sel, ], labels = seq_along(sel), pos = 3, cex = 0.7)

X_simm   <- X[sel, , ] # this can be regarded as the symmetric case --> 240 cases are more than enough
gpa <- procGPA(X_simm, scale = TRUE, reflect = FALSE)
mean_shape <- gpa$mshape # we will consider this as the star shaped object 
plotshapes(mean_shape)
### ----- Sample for the symmetric case  --- #####
# The basic idea is that small warping and high warping as an absolue concept has no proper meaning
# for this reason, it has to be define relative to the distance between landmarks
d <- dist(mean_shape)
consecutive_dist <- numeric(dim(mean_shape)[1] - 1)
for (i in 1:(dim(mean_shape)[1] - 1)) {
  consecutive_dist[i] <- d[i, i+1]
  
}

sam_no_warp <- sample_from_mean(mean_shape = mean_shape, n = 100,sigma_2 = 0)
sam_poco_warping <- sample_from_mean(mean_shape = mean_shape, n = 100,sigma_2 = min(consecutive_dist)/5)
sam_tanto_warping <- sample_from_mean(mean_shape = mean_shape, n = 100,sigma_2 = min(consecutive_dist)/3)

# save(sam_no_warp, file = here("Simulazioni/Hearts/no_warp_sym.RData"))
# save(sam_poco_warping, file = here("Simulazioni/Hearts/little_warp_sym.RData"))
# save(sam_tanto_warping, file = here("Simulazioni/Hearts/more_warp_sym.RData"))


load(here("Simulazioni/Hearts/more warp/asym/more_warp_gpa_asym_MCMC.RData"))
dev.off()
dev.list()
options(error = NULL)
output <- MCMC_samp$output
hypers <- MCMC_samp$hyper
try(visualize_results(output, hypers, X, "Simulazioni/Hearts/more warp/asym/more_warp_gpa_asym", sam)) 
### ------ Sample for the asymmetric case ------ 
delta <- 0.5                                  # asymmetry strength relative to unit size
mean_asym <- mean_shape                          # the symmetric 24-landmark mean
mean_asym[, 1] <- mean_asym[, 1] + delta * mean_asym[, 2] # try to have the slope of the second part stronger
plotshapes(mean_asym, joinline = 1:nrow(mean_asym)) # it works!
d <- dist(mean_asym)
consecutive_dist <- numeric(dim(mean_asym)[1] - 1)
for (i in 1:(dim(mean_asym)[1] - 1)) {
  consecutive_dist[i] <- d[i, i+1]
  
}

sam_no_warp <- sample_from_mean(mean_shape = mean_asym, n = 100,sigma_2 = 0)
sam_poco_warping <- sample_from_mean(mean_shape = mean_asym, n = 100,sigma_2 = min(consecutive_dist)/5)
sam_tanto_warping <- sample_from_mean(mean_shape = mean_asym, n = 100,sigma_2 = min(consecutive_dist)/3)


# save(sam_no_warp, file = here("Simulazioni/Hearts/no_warp_asym.RData"))
# save(sam_poco_warping, file = here("Simulazioni/Hearts/little_warp_asym.RData"))
# save(sam_tanto_warping, file = here("Simulazioni/Hearts/more_warp_asym.RData"))



load(here("more warp/sym/more_warp_gpa_sym_MCMC.RData"))