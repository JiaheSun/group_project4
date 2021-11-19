# Group number: 7  (Jiahe Sun)  (Lingda Wang)  (Qi Yu)
# Github: 


######################################################################################
###### 1. Compute Hessian matrix #####################################################
######################################################################################
compute_numerical_Hessian <- function(theta, f, ...) {
  # description: Calculate approximate Hessian matrix using finite differencing
  
  # arguments:
  # theta: the initial values for the optimization parameters
  # f: the objective function to minimize
  
  # return: Hessian matrix
  
  n <- length(theta)
  # finite difference interval
  eps <- 1e-7
  # Create empty Hessian matrix (filled by NA)
  H <- matrix(rep(NA, n^2), ncol = n)
  # loop over parameters
  for (i in 1:n) {
    for (j in 1:n) {    
      ei <- rep(0, n)
      ej <- rep(0, n)
      ei[i] <- 1
      ej[j] <- 1
      f1 <-  f(theta + eps * ei + eps * ej, T, ...)
      f2 <-  f(theta + eps * ei - eps * ej, T, ...)
      f3 <-  f(theta - eps * ei + eps * ej, T, ...)
      f4 <-  f(theta - eps * ei - eps * ej, T, ...)
      numdiff = (f1-f2-f3+f4)/(4*eps*eps)
      H[i,j] <- numdiff
    }
  }
  # fix the Hessian to be symmetric
  (H + t(H))/2
}


######################################################################################
###### 2. Search for step length #####################################################
######################################################################################
line_search <- function(theta,f,...,d=p,c1 = 0,c2 = 0.9) {
  # DESCRIPTION: 
  #   Use wolfe conditions and bisection method to find suitable step "alpha".
  #   1) Find an interval containing desirable step lengths, 
  #   2) Use bisection computes a good step length within this interval.
  
  # ARGUMENT:
  #   theta: is a vector of initial values for the optimization parameters
  #   f: is the objective function to minimize. Its first argument is the vector of optimization parameters. 
  #      Its second argument is a logical indicating whether or not gradients of the objective 
  #      w.r.t. the parameters should be computed. Remaining arguments will be passed from bfgs using '...'. 
  #      The scalar value returned by f will have a gradient attribute if the second argument to f is TRUE.
  #   ...: any arguments of f after the first two (the parameter vector, and gradient logical indicator) are passed using this.
  #   d: Search direction given theta
  #   c1: Wolfe condition(a) parameter
  #   c2: Wolfe condition(b) parameter ( 0<c1<c2<1 )
  
  # RETURN (containing): 
  #   alpha: The value of step length given theta
  #   theta_alpha: The new parameter value
  #   f_alpha: The new value of the objective function
  #   g_alpha: The new gradient vector 
  #   fail: Whether the search process failed (ie: 0 means success)
  
  # Initialization -------------------------------------
  n <- length(theta)
  # original parameter
  theta0 <- theta
  # original result of objective function
  out <- f(theta0, T, ...)
  # original scale value of objective function
  f0 <- as.numeric(out)
  # original gradient vector
  grd0 <- matrix(attr(out, "gradient"), ncol = 1)
  # original (gradient vector)%*%(search direction)
  gtd0 <- t(grd0) %*% d
  
  # Create the objects to be output --------------------
  # (The specific meaning is given at the beginning of this function)
  alpha <- 0
  theta_alpha <- 0
  f_alpha <- f0
  g_alpha <- grd0
  fail <- 0
  
  # Suppose an interval containing desirable step lengths
  # (i.e. the max of step lengths is beta)
  beta <- Inf
  
  # When the scale value of objective function is negative infinity,
  # stopping searching value
  fvalquit <- -Inf
  
  # first try for the step length 
  t <- 1
  
  # The maximum number of times to keep taking the median
  maxit_bisection <- 1e2
  n_bisect <- 0
  # The maximum number of times to increase step length
  maxit_increas <- 1e2
  n_increas <- 0
  
  # Loop begin to find suitable step length
  done <- F
  while (!done) {
    # The (next point) is equal to 
    # the (current point) plus (the changed distance multiplied by the direction)
    theta <- theta0 + t * d
    out <- f(theta, T, ...)
    fun <- as.numeric(out)
    grd <- matrix(attr(out, "gradient"), ncol = 1)
    # In order to reduce the amount of calculation, create (g^T)*(d)
    gtd <- t(grd) %*% d 
    
    # When the value of the objective function is small enough, 
    # it can be output directly
    if (fun < fvalquit) {
      return(list( alpha = t, theta_alpha = theta, f_alpha = fun,
                   g_alpha = grd, fail = fail ) ) 
    }
    
    
    # Check Wolfe conditions
    if (fun > f0 + c1 * t * gtd0) {
      # Wolfe condition (a) is violated
      beta <- t # decrease step length
    } else if (gtd < c2 * gtd0 || is.na(gtd)) {
      # Wolfe condition (b) is violated
      alpha <- t # increase step length
    } else {
      # Both conditions are satisfied
      return(list( alpha = t, theta_alpha = theta, f_alpha = fun,
                   g_alpha = grd, fail = fail ) )
    }
    
    
    # For next checking
    if (beta < Inf) {
      if (n_bisect < maxit_bisection) {
        # There is an interval for the step size at this time, 
        # take the middle value
        n_bisect <- n_bisect + 1
        t <- (alpha + beta) / 2
      } else {
        # Stop after reaching the maximum number of loops
        done <- T
      }
    } 
    else {
      # beta is Inf
      if (n_increas < maxit_increas) {
        # There is no interval at this time, 
        # the original step length is directly multiplied by two
        n_increas <- n_increas + 1
        t <- 2 * alpha
      } else {
        # Stop after reaching the maximum number of loops
        done <- T
      }
    }
    
  } # end while
  
  stop("Linesearch: Failed to satisfy Wolfe conditions.")
}


######################################################################################
###### 3. Combined all ###############################################################
######################################################################################
bfgs <- function(theta,f,...,tol = 1e-5,fscale = 1,maxit = 100) {
  # DESCRIPTION: use BFGS method to calculate the optimal value of function f, use the initial values theta
  
  # ARGUMENT:
  #  theta: is a vector of initial values for the optimization parameters
  #  f: is the objective function to minimize. Its first argument is the vector of optimization parameters. 
  #     Its second argument is a logical indicating whether or not gradients of the objective 
  #     w.r.t. the parameters should be computed. Remaining arguments will be passed from bfgs using '...'. 
  #     The scalar value returned by f will have a gradient attribute if the second argument to f is TRUE.
  #  ...: any arguments of f after the first two (the parameter vector, and gradient logical indicator) are passed using this.
  #  tol: the convergence tolerance.
  #  fscale: a rough estimate of the magnitude of f at the optimum - used in convergence testing.
  #  maxit: the maximum number of BFGS iterations to try before giving up.
  
  # RETURN (a list containing): 
  #  f the scalar value of the objective function at the minimum. theta the vector of values of the parameters at the minimum.
  #  iter the number of iterations taken to reach the minimum.
  #  g the gradient vector at the minimum (so the user can judge closeness to numerical zero).
  #  H the approximate Hessian matrix (obtained by finite differencing) at the minimum.
  
  # theta_new is the value of the parameters at the minimum
  n <- length(theta)
  theta <- as.matrix(theta)
  
  # initial value for approx. inverse Hessian, i.e., B0 is assumed to be I  (in the note)
  B0 <- diag(1, n) 
  
  # The result of objective funtion
  out <- f(theta, T, ...)
  
  # When the supplied f does not supply a gradient attribute
  # Compute the gradient by finite differencing
  if (is.null(attr(out, "gradient"))) {
    # finite difference interval
    eps <- 1e-7
    fold <- f
    f <- function(theta, T, ...) {
      n <- length(theta)
      out <- fold(theta, T, ...)
      grad <- rep(NA, n)
      # loop over parameters
      for (i in 1:n) {
        theta_new <- theta
        theta_new[i] <- theta_new[i] + eps # increase theta_new[i] by eps
        grad[i] <- (fold(theta_new, T, ...) - fold(theta, T, ...)) / eps }
      attr(out, "gradient") <- grad
      out
    }
  }
  
  
  # Initializations ---------------------------------------------------------
  # approx. inverse Hessian
  B <- B0
  # objective function result
  out <- f(theta, T, ...)
  # scalar value of objective function
  fval <- as.numeric(out)
  # the original scalar value of objective function
  fval_init <- fval 
  # gradient vactor 
  g <- matrix(attr(out, "gradient"),ncol=1)
  
  
  # Create a function "is.nainf()"  for checking na or infinite --------------
  is.nainf <- function(x)
    any(is.na(x) | is.infinite(x))
  # Checking If the objective or derivatives are not finite at the initial theta
  if (is.nainf(fval) || is.nainf(g)) {
    mes <-"Function or its gradient is not well defined at the initual value"
    stop(mes) }
  
  # Begin main loop ----------------------------------------------------------
  # The maximum number of loops is "maxit"
  for (iter in 1:maxit) {
    # Quasi-Newton step from theta_old
    p <- -B %*% g 
    # for BFGS update
    gprev <- g 
    
    # Find step length by quoting function(created above) line_search
    out_line_search <- line_search(theta, f, ..., d = p, c1 = 0, c2 = 0.9)
    alpha <- out_line_search$alpha  # actual step size
    theta <- as.matrix(out_line_search$theta_alpha)
    fval <- out_line_search$f_alpha
    g <- out_line_search$g_alpha
    
    #  To judge whether the gradient vector is close enough to zero
    if (max(abs(g)) < (abs(fval) + fscale) * tol) {
      return( list( theta = theta, f = fval, iter = iter, g = g,
                    H = compute_numerical_Hessian(theta, f, ...) ) ) }
    
    # BFGS update (the equation used accroding to the notes)
    s <- alpha * p ; y <- g - gprev ; sty <- c(t(s) %*% y)
    rho <- 1 / sty
    rhoByst <- rho * (B %*% y) %*% t(s)
    B <- B - t(rhoByst) - rhoByst + rho * s %*% (t(y) %*% rhoByst) + rho * s %*% t(s)
  } # end for
  warning("The maxit is reached without convergence.")
  
  # If convergence has not occurred (objective function value cannnot decrease)
  if(fval > fval_init) {
    warnings("We fail to reduce the objective but convergence has not occurred.")}
  
  # return a list required
  return(list( theta = theta, f = fval, iter = iter, g = as.vector(g),
               H = compute_numerical_Hessian(theta, f, ...) ) )
}


######################################################################################
###### 4. Use Rosenbrock function test ###############################################
######################################################################################
rb <-function(theta, getg = F, k = 10) {
  ## Rosenbrock objective function, suitable for use by 'bfgs'
  z <- theta[1] ; x <- theta[2]
  f <- k * (z - x ^ 2) ^ 2 + (1 - x) ^ 2 + 1
  if (getg) {
    attr(f,"gradient") <- c(2*k*(z-x^2), -4*k*x*(z-x^2) -2*(1-x))
    }
  f
} 
out <- bfgs(c(-1, 2),rb,k = 10,tol = 1e-5,fscale = 1,maxit = 100)
print(out)
print(optim(c(-1, 2), rb, method = "BFGS", hessian = T))

