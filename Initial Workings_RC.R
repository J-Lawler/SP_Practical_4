

###############################################################################


#The code below implements a crude version of the newt function

#### Problems still to solve:

# At each step need to check that the hessian is +ve definite and perturb it
# to be if it is not

# What if Hessian is not supplied - need to create using finite differences

# What if gradient is not supplied - do we need to consider this case?

# Need to improve on checking convergence - now only uses gradient < tolerance**

# Efficiency - Cholesky decomposition to find inverse instead of solve function

# What does fscale do? How do we use it?

# Need to make sure the ... argument is working properly

# Implement appropriate errors and warnings**







######## Function Code ########



newt <- function(theta,func,grad,hess=NULL,...,tol=1e-8,fscale=1,maxit=100,max.half=20,eps=1e-6){
  
  # stop if the objective or derivatives are not finite at the initial theta
  
  if( is.infinite(grad(theta)) || is.infinite(fun(theta))){ stop("infinite at theta")}
  
  # Starting point
  step_prev <- theta
  
  # Gradient at starting point
  gradient <- grad(theta)
  
  # Count iterations 
  iter <-  1
  
  # Loop through steps
  # Stop when gradient is (very near) 0 or we exceed the maximum number of iterations
  while(judge(gradient, tol, fscale, func(step)) == FALSE & iter <= maxit){ # hessian should be +ve def also
    
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
      
      # is this for error 3
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
      
    # stop if maxit is reached without convergence, adding try code here since 
    #inside while loop means not convergence
    if(iter == maxit){
      warning("reached maxit without convergence")
    }
      
    }
    
    step_prev <- step
    
    # if iter not reach maxit, add 1 on iter
    if ( iter < maxit ){
    iter <- iter + 1
    }
    
    # warning if Hessian is not +ve def
    e <- try(chol(hessian), silent=TRUE)
    if(inherits(e, "try-error")) warning("not +ve def hessian matrix at convergence ")
  }
  
  
  final <- step
  
  final
  
}

judge <- function(grad, tol, fscale, f){
  con <- tol * (abs(f) + fscale)
  
  #get gradient
  abs_value <- abs(grad)
  
  #get judge value, if it is 0, then condition satisfied
  result <- sum(abs_value >= con)
  
  if( result != 0 ){
    FALSE
  }
  else {
    TRUE
  }
  
  
}




