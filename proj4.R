



newt <- function(theta,func,grad,hess=NULL,...,tol=1e-8,fscale=1,maxit=100,max.half=20,eps=1e-6){
  

  #if( is.infinite(grad(theta)) || is.infinite(func(theta))){ stop("infinite at theta")}
  
  # Starting point
  step_prev <- theta
  
  # Function at starting point
  func_start <- func(theta,...)
  
  # Gradient at starting point
  gradient <- grad(theta,...)
  
  # hessian at starting point
  hessian_start <- hessian_func(theta,gradient,hess,grad,eps,...)
  
  # stop if the objective or derivatives are not finite at the initial theta
  if(sum(is.infinite(c(func_start,gradient,hessian_start)))!=0){
    stop(stop_warn_text(1,stop = TRUE))
  }
  
  # Count iterations 
  iter <-  1
  stop <-  FALSE
  stop_warn_message = NULL
  
  # Loop through steps
  # Stop when gradient is (very near) 0 or we exceed the maximum number of iterations
  while(stop == FALSE){ # hessian should be +ve def also
    
    # gradient at current step
    gradient <- grad(step_prev,...)
    
    # hessian at current step
    hessian <- hessian_func(step_prev,gradient,hess,grad,eps,...)
    
    # Propose a descent direction
    
    delta <- propose_descent(hessian,gradient)
    
    # next step
    step <-  step_prev+delta
    
    # Check step decreases function, if not half the delta and try again, repeat a maximum of max.half times
    
    for(j in 1:max.half){
      
      # Check function has decreased
      if(func(step_prev,...)-func(step,...)>0){
        break
      }
      else{
        delta <- delta/2
        step <- step_prev+delta
      }
      
      # if we half the delta too many times
      if(j == max.half){stop(stop_warn_text(2,stop = TRUE))}
      
      }
    
    step_prev <- step
    
    stop <- judge(func,grad, step, tol, fscale,...)
    
    # if iter not reach maxit, add 1 on iter
    if ( iter <= maxit ){
      iter <- iter + 1
    } else{
      stop_warn_message <- stop_warn_text(3,stop_warn_message)
      break
    }
  }
  
  # Collate elements to return
  
  fun_at_stop <- func(step,...)
  
  grad_at_stop <- grad(step,...)

  hess_at_stop <- hessian_func(step,grad_at_stop,hess,grad,eps,...)
  
  inv_hess_at_stop <- NULL
  
  try(inv_hess_at_stop <- chol2inv(chol(hess_at_stop)), silent = TRUE)
  
  if(is.null(inv_hess_at_stop)){stop_warn_message <- stop_warn_text(4,stop_warn_message)}
  
  iterations <- iter-1
  
  result <- list(f = fun_at_stop,
                 theta = step,
                 iter = iterations,
                 g = grad_at_stop,
                 Hi = inv_hess_at_stop)
  
  if(!is.null(stop_warn_message)){
    warning(stop_warn_message)
  }
  return(result)
}


#########################

# Stopping Criteria

judge <- function(func,grad, step, tol, fscale,...){
  
  con <- tol * (abs(func(step,...)) + fscale)
  
  #get gradient
  abs_value <- abs(grad(step,...))
  
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


approximate_hess <- function(gradient_value, step_prev, grad_fun,eps,...){
  
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
    grad1 <- grad_fun(th1,...)
    # hessian approximated using the difference between grad1 and 0 over eps
    hess_fd[i,] <- (grad1 - gradient_value)/eps
  }
  # Ensure symmetry
  hess_approx <- (t(hess_fd)+hess_fd)/2
}



##############################


stop_warn_text = function(i,stop_warn_message=NULL,stop = FALSE){
  
  if(stop == TRUE){
    stop_warn_message <- "Minimum not found."
  }
  
  if(is.null(stop_warn_message)& stop == FALSE){
    stop_warn_message <- "Minimum likely not found."
  }
  
  if(i == 1){
    message = paste(stop_warn_message,"The given function or its derivatives are not finite at the initial theta.")
  }
  else if(i == 2){
    message = paste(stop_warn_message,"\nDescent direction failed to reduce the function after max.half step halvings.")
  }
  else if(i == 3){
    message = paste(stop_warn_message,"\nFailed to find convergence within maxit iterations.")
  }
  else{
    message = paste(stop_warn_message,"\nHessian at convergence is not positive definite.")
  }
  message
}


###############################


hessian_func <- function(step_point,grad_at_step,hess,grad,eps,...){
  if(is.null(hess)){
    hessian <- approximate_hess(grad_at_step, step_point, grad,eps,...)
  }
  else{
    hessian <- hess(step_point,...)
  }
  
}
