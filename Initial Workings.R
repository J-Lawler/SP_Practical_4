

###############################################################################


#The code below implements a crude version of the newt function

#### Problems still to solve:

# At each step need to check that the hessian is +ve definite and perturb it
# to be if it is not

# What if Hessian is not supplied - need to create using finite differences

# Need to improve on checking convergence - now only uses gradient < tolerance

# Need to make sure the ... argument is working properly

# Implement appropriate errors and warnings

# What if gradient is not supplied - do we need to consider this case?







######## Function Code ########



newt <- function(theta,func,grad,hess=NULL,...,tol=1e-8,fscale=1,maxit=100,max.half=20,eps=1e-6){
  
  # Starting point
  step_prev <- theta
  
  # Gradient at starting point
  gradient <- grad(theta)
  
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
    R <- chol(hessian)
    # Solve Hessian^-1 * gradient
    delta <- - backsolve(R,forwardsolve(t(R),gradient))
    
    # next step
    step <-  step_prev+delta
    
    # Check step decreases function, if not half the delta and try again, repeat a maximum of max.half times
    
    for(j in 1:max.half+1){
      
      # if we half the delta too many times
      if(j == max.half+1){stop("Oops")}
      
      # Check function has decreased (and the objective function is not infinte)
      if(is.finite(func(step_prev)-func(step)) & func(step_prev)-func(step)>0){
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










