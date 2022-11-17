
test_grad <- function(th,k=2){ 
  c(-2*(1-th[1])-k*4*th[1]*(th[2]-th[1]^2),k*2*(th[2]-th[1]^2))
}

if hess=NULL {
  
  # test parameter values
  th0 <- c(2, 2) 
  #remove
  eps <- 1e-6
  # gradient values at theta0
  grad0 <- test_grad(th=th0, 2) 
  
  # variable to allocate space for hessian
  hess_fd <- matrix(0,2,2) 
  
  # loop through each parameter in range 1 to length theta nought 
  for (i in 1:length(th0)) {
    # Set theta 1 to theta 0
    th1 <- th0
    # for each element in theta 1, update theta 1 by the value of eps
    th1[i] <- th1[i] + eps
    # gradient vector at theta 1
    grad1 <- test_grad(th=th1, 2)
    # hessian approximated using the difference between grad1 and 0 over eps
    hess_fd[i,] <- (grad1 - grad0)/eps
    # Assigns our approximated hessian matrix to our hess variable
    hess <- hess_fd
    # Returns our hessian matrix
    hess
  }
}


test_hess <- function(th,k=2){ 
  h <- matrix(0,2,2) 
  h[1,1] <- 2-k*2*(2*(th[2]-th[1]^2) - 4*th[1]^2) 
  h[2,2] <- 2*k 
  h[1,2] <- h[2,1] <- -4*k*th[1] 
  h
}

hess; test_hess(th0, 2)
