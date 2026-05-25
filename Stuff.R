J <- 10
I <- diag(1, J)
D2  <-  matrix(0, J-2, J)
for (i in 1:nrow(D2)) {
  D2[i, i:(i+2)] <- c(-1, 2, -1)
  
}
R_J <- matrix(0, J, J)
for (i in 1:nrow(R_J)) {
  R_J[i, J - (i-1)] <- 1
  
}

K <- t(D2)%*%(D2) + (I- R_J)
det(K)
prod(eigen(K)$values[eigen(K)$values > 1e-2])
install.packages("maotai")
library(maotai)
pdeterminant(K)



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
