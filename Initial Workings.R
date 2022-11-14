

###############################################################################


#The code below implements a crude version of the newt function

#### Problems still to solve:

# At each step need to check that the hessian is +ve definite and perturb it
# to be if it is not

# What if Hessian is not supplied - need to create using finite differences

# What if gradient is not supplied - do we need to consider this case?

# Need to improve on checking convergence - now only uses gradient < tolerance

# Efficiency - Cholesky decomposition to find inverse instead of solve function

# What does fscale do? How do we use it?

# Need to make sure the ... argument is working properly

# Implement appropriate errors and warnings




###############################################################################




# Test Functions from Task Sheet

# Function to minimise
rb <- function(th,k=2){ 
  k*(th[2]-th[1]^2)^2 + (1-th[1])^2
  }

# Gradient
gb <- function(th,k=2){ 
  c(-2*(1-th[1])-k*4*th[1]*(th[2]-th[1]^2),k*2*(th[2]-th[1]^2))
  }

# Hessian
hb <- function(th,k=2){ 
  h <- matrix(0,2,2) 
  h[1,1] <- 2-k*2*(2*(th[2]-th[1]^2) - 4*th[1]^2) 
  h[2,2] <- 2*k 
  h[1,2] <- h[2,1] <- -4*k*th[1] 
  h
  }




########## Use In-built Function ##########


init_nlm <- c(2,2)

estimate_nlm <- nlm(f= rb, p = init_nlm)$estimate




########## Use Hand-built Function ##########



# Parameters
init <- c(2,2)



newt <- function(theta,func,grad,hess=NULL,...,tol=1e-8,fscale=1,maxit=100,max.half=20,eps=1e-6){
  
  # Starting point
  step_prev <- init
  
  # Gradient at starting point
  gradient <- grad(init)
  
  # Count iterations 
  iter <-  1
  
  # Loop through steps
  # Stop when gradient is (very near) 0 or we exceed the maximum number of iterations
  while(abs(sum(gradient))>=tol & iter <= maxit ){ # hessian should be +ve def also
    
    # hessian at current step
    hessian <- hess(step_prev)
    
    # gradient at current step
    gradient <- grad(step_prev)
    
    # the change to find the next step - see lecture notes page 66
    delta <- - solve(hessian)%*%gradient
    
    # next step
    step <-  step_prev+delta
    
    # Check step decreases function, if not half the delta and try again, repeat a maximum of max.half times
    
    for(j in 1:max.half+1){
      
      # if we half the delta too many times
      if(j == max.half+1){stop("Oops")}
      
      # Check function has decreased
      if(func(step_prev)-func(step)>0){
        break
      }
      else{
        delta <- delta/2
        step <- step_prev+delta
      }
    }
    
    step_prev <- step
    iter <- iter + 1
  }
  
  
  final <- step
  
  final
  
}




estimate_newt <- newt(theta = init,func = rb,grad = gb,hess = hb,
     tol=1e-8,fscale=1,maxit=100,max.half=20,eps=1e-6)


####################


print("Estimate from in-built R function:")

estimate_nlm


print("Estimate from newt function:")

estimate_newt[,1]

