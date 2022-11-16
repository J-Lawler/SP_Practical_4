test_gradient <- c(2,3,21,-23,-4,-5,-7,9)
test_gradient_success <- c(1,1,1,1,1,1)

#set test tol value
test_tol <- 2
fscale <- 1
test_f <- 2

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


judge(test_gradient, test_tol, fscale, test_f)
judge(test_gradient_success, test_tol, fscale, test_f)
