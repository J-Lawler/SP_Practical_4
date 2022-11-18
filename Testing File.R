

source("proj4.R")


###############################################################################




# Test Functions from Task Sheet

# Function to minimise
test_func <- function(param,j){ 
  j*(param[2]-param[1]^2)^2 + (1-param[1])^2
}

# Gradient
test_grad <- function(param,j){ 
  c(-2*(1-param[1])-j*4*param[1]*(param[2]-param[1]^2),j*2*(param[2]-param[1]^2))
}

# test_hess<- function(theta, j){
#   Inf
# }

# Hessian
test_hess <- function(param,j){ 
  h <- matrix(0,2,2) 
  h[1,1] <- 2-j*2*(2*(param[2]-param[1]^2) - 4*param[1]^2) 
  h[2,2] <- 2*j 
  h[1,2] <- h[2,1] <- -4*j*param[1] 
  h
}


test_inputs <- list(rb = test_func,
                    gb = test_grad,
                    hb = test_hess)

# Initial Parameters


init_param <- c(2,2)

########## Use In-built Function ##########


estimate_nlm <- nlm(f= test_inputs$rb, p = init_param,j=2)$estimate


########## Use Hand-built Function ##########


# Estimate

estimate_newt <- newt(theta = init_param,func = test_inputs$rb,
                      grad = test_inputs$gb,hess = test_inputs$hb, j = 2)

# Estimate when Hessian is not supplied

estimate_newt <- newt(theta = init_param,func = test_inputs$rb,
                      grad = test_inputs$gb, j = 40, max.half = 10)


####################

# Compare Estimates

print("Estimate from in-built R function:")

estimate_nlm


print("Estimate from newt function:")

estimate_newt$theta




################################################################################

# testing the approximate_hess function

testing_parameters <- c(4,2)

hess_approx <- approximate_hess(test_grad(testing_parameters), testing_parameters, test_grad,1e-6)

hess_true <- test_hess(testing_parameters)

print(hess_approx)

print(hess_true)
