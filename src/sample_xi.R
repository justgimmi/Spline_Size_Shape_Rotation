# sample_xi <- function(init, hyper){
#   ##### this function is useful for sampling the inclusion probility of the warping function #### 
#   
#   hyper$a_xi_post <- hyper$a_chi + init$z
#   hyper$b_xi_post <- hyper$b_chi + 1 - init$z
#   init$chi <- rbeta(n = 1, shape1 =hyper$a_xi_post, shape2 =  hyper$b_xi_post)
# }
