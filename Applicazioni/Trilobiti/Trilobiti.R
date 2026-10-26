# Source file -----
require(here)
source(here("src/Packages.R"))
source(here("src/Utils.R"))
source(here("src/RB_MCMC.R"))
# Import data 
# We consider the trilo dataset available on the momocs package
# Arranged from: https://folk.universitetetioslo.no/ 
# The original data included 51 outlines and 5 ontogenetic stages,
# but one of them has just a single outline thas has been removed.


# Outline coordinates of 50 cephalic outlines of trilobite 

dati <- trilo # trilo dataset inside Momocs library
trilo_list <- unlist(dati, recursive = FALSE) # unlist data
trilo_list <- trilo_list[sapply(trilo_list, is.matrix)] # transform each of them into matrix

X <- simplify2array(trilo_list) 
X
gpa <- procGPA(X, scale = TRUE, reflect = FALSE)
mean_shape_full <- gpa$mshape
mean_shape <- gpa$mshape          # k x 2 Procrustes mean shape --> we consider it as an average configuration
par(mfrow = c(1,1))
plotshapes(mean_shape, joinline = 1:nrow(mean_shape))
M <- mean_shape           
k <- nrow(M) # number of landmakrs (80 too many =( )
plot(M, type = "b", asp = 1) 
text(M, labels = 1:k, pos = 3, cex = 0.6)



# Look at the landmark numbers and let's pick some anchors
plot(M, type = "b", asp = 1) 
text(M, labels = 1:k, pos = 3, cex = 0.6)
# we argue that two important landmarks are 1 and 34

top_idx <- 1
bot_idx <- 34

sel <- pick_landmarks(M, top_idx, bot_idx, n_side = 7)   # 16 landmarks
plot(M, type = "l", asp = 1)
points(M[sel, ], pch = 19, col = "red")
text(M[sel, ], labels = seq_along(sel), pos = 3, cex = 0.7)  # it seems a reasonable amount of landmarks

X_final   <- X[sel, , ]
aligned <- lapply(1:n, function(i)
((X_final[, , i])))
par(mfrow = c(1,1))
rng <- range(unlist(aligned))
plot(NA, xlim = rng, ylim = rng, asp = 1, xlab = "", ylab = "")
for (A in aligned){
  points(A, col = adjustcolor("grey40", 0.3), pch = 16, cex = 0.6)
  lines(A, type = "l")
}

#### ----  Fit the model ---- ####

prova <- init_mu_star(X_final/100,n_basis = 18, degree = 3, start = "gpa", ridge   = 1e-2, plot = TRUE,
                                               Sigma_eta = diag(10, 2)) # we start from the gpa mean


init <- prova$init # define init object 
hyper <- prova$hyper # define hyper object 
X <- prova$X # define X 
MCMC_trilo <- RB_MCMC(total_iter = 1e5, burnin = 30000, thinning  = 30, X = X, 
                                           init = init, hyper = hyper)
save(X,MCMC_trilo, file = here("Applicazioni/Trilobiti/Trilo_MCMC.RData")) # here you can find the MCMC object if u do not want to wait =)
dev.off()
dev.list()
options(error = NULL)
output <- MCMC_trilo$output
hypers <- MCMC_trilo$hyper
try(visualize_results(output, hypers, X, "Trilobite")) # this create several plots in the same folder
