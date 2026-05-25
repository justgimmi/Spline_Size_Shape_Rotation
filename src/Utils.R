source(file = "Packages.R")

# Simulation Function ----- 
# Here we define a function to sample from the underlying process 




# Basis Construction ------
# Here we define a function to define the cubic splines 
# x <- runif(100) * 2* pi
# equi <- seq(0, 2*pi, by = 2*pi/20)
# basis <- bs(x, knots = equi[2:20], Boundary.knots = c(0, 2*pi))



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

K1_construction(21)
K2_construction(21)