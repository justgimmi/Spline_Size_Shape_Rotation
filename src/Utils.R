source(file = "Packages.R")

# Basis Construction ------
# Here we define a function to define the cubic splines 

# 
# library(mgcv)
# x <- runif(50) * 2* pi
# 
# sx <- s(x, bs = "bs", k = 16, m = 3)
# smooth <- smoothCon(sx, data = data.frame(x = x), absorb.cons = F)
# smooth[[1]]$knots # prints the knots
# par(mfrow = c(1, 1))
# smooth[[1]]$X |> matplot(type = "l") # shows the basis functions
# 


Basis_Construction <- function(thetas, L, degree = 3){
  # I do not like to use already implemented function because they are based
  # on the range of the vector to define the knots so every iteration of the MCMC
  # we could have different basis following this approach
  
  # thetas: --> vector of theta values where to evaluate the basis 
  # L: --> how many internal knots do we want?
  # degree: --> spline degree
  
  # The idea of this function is the following. Given the fact that I want the starting value in 0 
  # and the end point in 2 \pi, I first define equispaced knots on this interval. Then, we habe to
  # add padding on the left and on the right. A conservative choice would be to add degree times 0 before 
  # and degree times 2\pi after. Talking with Johannes, I understood that the best approach to have 
  # meaningfull penalties is to have equispaced knots also on the left and right part of the interval
  delta <- (2 * pi) / L # equispaced shift 
  
  knots_base <- seq(0, 2 * pi, by = delta)[-c(1, L+1)] # base knots
  
  total_knots <- c(seq(0 - degree*delta, 0, by = delta), knots_base, seq(2*pi , 2*pi + degree*delta, 
                                                                         by = delta))
  Basis <- splineDesign(total_knots, thetas, ord = degree + 1, outer.ok = FALSE) # Basis construction
  # At the End we obtain a matrix of size length(theta) x (L + degree)
  return(Basis)
}



# Useful function for K1 and K2 -----
K1_construction <- function(J){ # smoothness constraint
  # J: --> number of basis 
  I <- diag(1, J) # define identity matrix
  D2  <-  matrix(0, J-2, J) # define empty matrix
  for (i in 1:nrow(D2)) {
    D2[i, i:(i+2)] <- c(-1, 2, -1)
  }
  return(t(D2) %*% D2)
  
}

K2_construction <- function(J){ # symmetry constraint
  # J: --> number of basis 
  I <- diag(1, J) # define identity matrix
  R_J <- matrix(0, J, J) # reverse identity matrix 
  for (i in 1:nrow(R_J)) {
    R_J[i, J - (i-1)] <- 1
  } 
  return(I - R_J)
}

# Simulation Function ----
