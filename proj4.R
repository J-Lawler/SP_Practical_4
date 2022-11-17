







newt <- function(theta,func,grad,hess=NULL,...,tol=1e-8,fscale=1,maxit=100,max.half=20,eps=1e-6){
  
  # stop if the objective or derivatives are not finite at the initial theta
  #if( is.infinite(grad(theta)) || is.infinite(func(theta))){ stop("infinite at theta")}
  
  # Starting point
  step_prev <- theta
  
  # Gradient at starting point
  gradient <- grad(theta)
  
  # Count iterations 
  iter <-  1
  stop <-  FALSE
  
  # Loop through steps
  # Stop when gradient is (very near) 0 or we exceed the maximum number of iterations
  while(stop == FALSE & iter <= maxit){ # hessian should be +ve def also
    
    # gradient at current step
    gradient <- grad(step_prev)
    
    # hessian at current step
    if(is.null(hess)){
      hessian <- approximate_hess(gradient, step_prev, grad,eps)
    }
    else{
      hessian <- hess(step_prev)
    }
    
    # Propose a descent direction
    
    delta <- propose_descent(hessian,gradient)
    
    # next step
    step <-  step_prev+delta
    
    # Check step decreases function, if not half the delta and try again, repeat a maximum of max.half times
    
    for(j in 1:max.half+1){
      
      # if we half the delta too many times
      if(j == max.half+1){warning("Fails to reduce by harlfing")}
      
      # Check function has decreased
      if(func(step_prev)-func(step)>0){
        break
      }
      else{
        delta <- delta/2
        step <- step_prev+delta
      }
      
      # stop if maxit is reached without convergence, adding try code here since 
      # inside while loop means not convergence
      if(iter == maxit){
        warning("reached maxit without convergence")
      }
      
    }
    
    step_prev <- step
    
    stop <- judge(func,grad, step, tol, fscale)
    
    # if iter not reach maxit, add 1 on iter
    if ( iter < maxit ){
      iter <- iter + 1
    }
    
  }
  
  
  # warning if Hessian is not +ve def
  e <- try(chol(hessian), silent=TRUE)
  if(inherits(e, "try-error")) warning("not +ve def hessian matrix at convergence ")
  
  final <- step
  
  final
  
}


#########################

# Stopping Criteria

judge <- function(func,grad, step, tol, fscale){
  
  con <- tol * (abs(func(step)) + fscale)
  
  #get gradient
  abs_value <- abs(grad(step))
  
  #get judge value, if it is 0, then condition satisfied
  result <- sum(abs_value >= con)
  
  if( result != 0 ){
    FALSE
  }
  else {
    TRUE
  }
  
  
}


#########################

# Propose a descent direction

propose_descent = function(hessian_matrix, gradient_vector){
  
  R <- NULL
  
  # Check Hessian Positive Definite
  try(R <- chol(hessian_matrix), silent = TRUE)
  
  # If Not Perturb it to Be
  if(is.null(R)){
    min_eigen <- min(eigen(hessian_matrix)$values)
    
    perturb <- diag(nrow=nrow(hessian_matrix))*(1e-6-min_eigen)
    
    R <- chol(hessian_matrix+perturb)
  }
  
  # Solve Hessian^-1 * gradient
  
  descent <- - backsolve(R,forwardsolve(t(R),gradient_vector))
  
}


#########################


# Approximate Hessian


approximate_hess <- function(gradient_value, step_prev, grad_fun,eps){
  
  grad_len <- length(gradient_value)
  
  # variable to allocate space for hessian
  hess_fd <- matrix(0,nrow = grad_len, ncol = grad_len)
  
  # loop through each parameter in range 1 to length theta nought 
  for (i in 1:grad_len) {
    # Set theta 1 to theta 0
    th1 <- step_prev
    # for each element in theta 1, update theta 1 by the value of eps
    th1[i] <- th1[i] + eps
    # gradient vector at theta 1
    grad1 <- grad_fun(th=th1, 2)
    # hessian approximated using the difference between grad1 and 0 over eps
    hess_fd[i,] <- (grad1 - gradient_value)/eps
  }
  # Ensure symmetry
  hess_approx <- (t(hess_fd)+hess_fd)/2
}
