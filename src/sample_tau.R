# Sample the smoothness/symmetry parameter ------

sample_tau <- function(init, a_new, b_new){
  
  init$tau <- rinvgamma(n = 1, shape = a_new, scale = b_new)
}
