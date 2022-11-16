

source("Initial Workings.R")


###############################################################################




# Test Functions from Task Sheet

# Function to minimise
test_func <- function(th,k=2){ 
  k*(th[2]-th[1]^2)^2 + (1-th[1])^2
}

# Gradient
test_grad <- function(th,k=2){ 
  c(-2*(1-th[1])-k*4*th[1]*(th[2]-th[1]^2),k*2*(th[2]-th[1]^2))
}

# Hessian
test_hess <- function(th,k=2){ 
  h <- matrix(0,2,2) 
  h[1,1] <- 2-k*2*(2*(th[2]-th[1]^2) - 4*th[1]^2) 
  h[2,2] <- 2*k 
  h[1,2] <- h[2,1] <- -4*k*th[1] 
  h
}


test_inputs <- list(rb = test_func,
                    gb = test_grad,
                    hb = test_hess)

# Initial Parameters


init_param <- c(2,2)

########## Use In-built Function ##########


estimate_nlm <- nlm(f= test_inputs$rb, p = init_param)$estimate


########## Use Hand-built Function ##########


# Estimate

estimate_newt <- newt(theta = init_param,func = test_inputs$rb,
                      grad = test_inputs$gb,hess = test_inputs$hb)


####################


print("Estimate from in-built R function:")

estimate_nlm


print("Estimate from newt function:")

estimate_newt[,1]


