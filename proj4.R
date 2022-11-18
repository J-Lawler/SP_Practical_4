#################### Code for Statistical Programming Practical 4 ####################




#################### Names ####################

#
# Runchen Geng (s2332611)
#
# Sulaiman Zaman (s2405379)
#
# Jake Lawler (s2451377)


#################### Github Repository ####################

# https://github.com/J-Lawler/SP_Practical_4


#################### Contributions ####################

#
# Runchen Geng - Stopping criterion, errors and warnings, compiling, comments & copy-editing
#
# Sulaiman Zaman - Approximating hessian, compiling, comments & copy-editing
#
# Jake Lawler - Perturbing hessian, compiling, comments & copy-editing
#
# Else - All together.
#
# Proportion of work performed by each team member roughly equal.



#################### Task Description ####################

# In the code below, we define a function newt that takes as input a target 
# function and (minimally) its gradient and a starting point, and 
# returns the minimum of that target function.

# It implements Newton's method for minimisation.
# From the initial starting point, we fit a quadratic at that point, locally 
# approximating the target function.
# We minimise the quadratic to propose a descent direction.
# If the proposed step in the descent direction decreases the function,
# we accept the step and repeat the process.
# Otherwise we half the step and check if that decreases the function, repeating
# until we accept a step that decreases the function.
# This loop is repeated until convergence.
# Convergence occurs when the gradient at the proposed step is tolerably near 0.

# The function takes optional inputs:
  # hess - a function that returns the Hessian matrix of the target function
  # tol - the tolerance for checking if newt has converged to a minimum
  # fscale - an estimate of the size of the target function at its minimum
  # maxit - the maximum number of iterations attempted to find convergence
  # max.half - the maximum number of step halvings to atttempt.
  # eps - if the Hessian is not provided, the finite difference interval to use
          # in approximating it.

# The function returns a list containing:
  # f - the value of the target function at its minimum
  # theta - the value of the parameters at the minimum
  # iter - the number of steps taken to reach the minimum
  # g - the gradient vector at the minimum
  # Hi - the inverse of the Hessian at the minimum




#################### Over-arching Optimisation Function ####################


# This section defines the overarching function.
# It calls several functions defined below, but 
# should be independently understandable

##### Function Definition #####


# Set function input defaults
newt <- function(theta,func,grad,hess=NULL,...,tol=1e-8,fscale=1,maxit=100,max.half=20,eps=1e-6){
  
  # Starting point
  step_prev <- theta
  
  ### Check function and derivatives are well-defined at start
  
  # Function at starting point
  func_start <- func(theta,...)
  
  # Gradient at starting point
  gradient <- grad(theta,...)
  
  # Hessian at starting point
  hessian <- hessian_func(theta,gradient,hess,grad,eps,...)
  
  # Error if the objective or derivatives are not finite at the initial theta
  if(sum(is.infinite(c(func_start,gradient,hessian)))!=0){
    stop(stop_warn_text(1,stop = TRUE))
  }
  
  # Count iterations 
  iter <-  1
  # Initialise convergence criteria
  stop <-  FALSE
  # Initialise warning message object
  stop_warn_message = NULL
  
  # Loop through steps
  # The stop object is true when convergence criteria are met
  while(stop == FALSE){
    
    # Propose a descent direction
    delta <- propose_descent(hessian,gradient)
    
    # Next step proposal
    step <-  step_prev+delta
    
    # Loop to find proposal step that decreases function
    for(j in 1:max.half+1){
      
      # Error if function has not decreased after max.half halving attempts
      if(j > max.half){stop(stop_warn_text(2,stop = TRUE))}
      
      # Has function decreased?
      if(func(step_prev,...)-func(step,...)>0){
        break
      }
      # If not, half the delta and try again, repeat a maximum of max.half times
      else{
        delta <- delta/2
        step <- step_prev+delta
      }
      
      }
    
    # Current step becomes previous step
    step_prev <- step
    
    # Gradient at current step
    gradient <- grad(step,...)
    
    # Hessian at current step
    hessian <- hessian_func(step,gradient,hess,grad,eps,...)
    
    # Check whether convergence criteria are met
    stop <- judge(func,grad, step, tol, fscale,...)
    
    # Break loop with warning if we have reached maximum iterations
    if ( iter >= maxit ){
      stop_warn_message <- stop_warn_text(3,stop_warn_message)
      break
    }
    
    iter <- iter + 1
  }
  
  #### Collate elements to return
  
  # Function at convergence (or maxit)
  fun_at_stop <- func(step,...)
  
  # Calculate inverse Hessian
  inv_hess_at_stop <- NULL
  try(inv_hess_at_stop <- chol2inv(chol(hessian)), silent = TRUE)
  # If Cholesky decomposition fails, warn that Hessian is not positive definite
  if(is.null(inv_hess_at_stop)){stop_warn_message <- stop_warn_text(4,stop_warn_message)}
  
  # Iterations tried
  iterations <- iter
  
  # Results list
  result <- list(f = fun_at_stop,
                 theta = step,
                 iter = iterations,
                 g = gradient,
                 Hi = inv_hess_at_stop)
  
  # Return warnings if there are any
  if(!is.null(stop_warn_message)){
    warning(stop_warn_message)
  }
  # Return result
  return(result)
}




#################### Helper Functions ####################


########## Test for Convergence Criteria  ##########

# Function returns TRUE if convergence criteria met, otherwise FALSE

judge <- function(func,grad, step, tol, fscale,...){
  
  # Convergence tolerance (as independent of function scale as possible)
  con <- tol * (abs(func(step,...)) + fscale)
  
  # Absolute value of gradient
  abs_value <- abs(grad(step,...))
  
  # Is each element of absolute gradient within tolerance?
  result <- sum(abs_value >= con)
  
  if( result != 0 ){
    FALSE
  }
  else {
    TRUE
  }
  
  
}




########## Propose a Descent Direction  ##########

# This function returns a delta. The function decreases at some distance along
# the direction of this delta

propose_descent = function(hessian_matrix, gradient_vector){
  
  R <- NULL
  
  # Check Hessian positive definite
  try(R <- chol(hessian_matrix), silent = TRUE)
  
  # If Not perturb it to be positive definite
  if(is.null(R)){
    
    # One or more eigenvalues will be 0 or negative
    min_eigen <- min(eigen(hessian_matrix)$values)
    
    # Create a diagonal matrix with entries the inverse of the smallest eigenvalue
    # (plus a little bit, in case the smallest eigenvalue is 0)
    perturb <- diag(nrow=nrow(hessian_matrix))*(1e-6-min_eigen)
    
    # Adding this perturbation to the hessian guarantees positive definiteness
    R <- chol(hessian_matrix+perturb)
  }
  
  # Solve Hessian^-1 * gradient - this is our proposed descent direction
  # It minimisies the local quadratic approximating the target function
  descent <- - backsolve(R,forwardsolve(t(R),gradient_vector))
  
}




########## Approximate Hessian if None Supplied  ##########

# Approximate Hessian

# Takes as input the value of the gradient at the current step,
# the value of the current step, the supplied gradient function,
# and the finite difference size input eps

# Returns approximate hessian

approximate_hess <- function(gradient_value, step_prev, grad_fun,eps,...){
  
  # Dimension we are working in
  grad_len <- length(gradient_value)
  
  # Allocate space for hessian
  hess_fd <- matrix(0,nrow = grad_len, ncol = grad_len)
  
  # Loop through parameter for each dimension 
  for (i in 1:grad_len) {
    # Set theta 1 to current step
    th1 <- step_prev
    # For each element in parameter vector, update theta 1 by the value of eps
    th1[i] <- th1[i] + eps
    # Gradient vector at theta 1
    grad1 <- grad_fun(th1,...)
    # Hessian approximated using the difference between gradients over eps
    hess_fd[i,] <- (grad1 - gradient_value)/eps
  }
  # Ensure symmetry
  hess_approx <- (t(hess_fd)+hess_fd)/2
}


########## Function Wrapper for (Approximate) Hessian  ##########

# If Hessian supplied, returns Hessian at a given step
# If Hessian not supplied, returns approximate Hessian at a given step
# Does this by calling approximate_hess function above

hessian_func <- function(step_point,grad_at_step,hess,grad,eps,...){
  if(is.null(hess)){
    hessian <- approximate_hess(grad_at_step, step_point, grad,eps,...)
  }
  else{
    hessian <- hess(step_point,...)
  }
  
}




########## Warning and Error Messages ##########

# Provided text for required error and warning messages.
# If error message, input stop here will be TRUE

stop_warn_text = function(i,stop_warn_message=NULL,stop = FALSE){
  
  # If we need an error - state that minimum was not found
  if(stop == TRUE){
    stop_warn_message <- "Minimum not found."
  }
  
  # If we only need a warning, check that there is not already a warning
  # If not, state that minimum was likely not found
  if(is.null(stop_warn_message)& stop == FALSE){
    stop_warn_message <- "Minimum likely not found."
  }
  
  # Appropriate error or warning message appended to existing warning message
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



