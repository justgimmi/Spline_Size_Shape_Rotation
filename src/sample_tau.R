# Sample the smoothness/symmetry parameter ------

sample_tau <- function(a_new, b_new){
  
  return(rinvgamma(n = 1, shape = a_new, scale = b_new))
}