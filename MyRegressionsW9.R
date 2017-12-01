# MyRegressionsW9


##################################
## Function 1: QR decomposition ##
##################################

myQR <- function(A){
  
  ## Perform QR decomposition on the matrix A
  ## Input: 
  ## A, an n x m matrix
  
  ########################
  ## FILL IN CODE BELOW ##
  ########################  
  
  n <- nrow(A)
  m <- ncol(A)
  Q <- diag(n)
  R <- A
  
  for(k in 1:(m - 1)){
    x      <- rep(0, n)
    x[k:n] <- R[k:n, k]
    s      <- -1 * sign(x[k])
    v      <- x
    v[k]   <- x[k] - s * norm(x, type = "2")
    u      <- v / norm(v, type = "2")
    
    R <- R - 2 * u %*% t(u) %*% R
    Q <- Q - 2 * u %*% t(u) %*% Q
    
  }
  
  ## Function should output a list with Q.transpose and R
  ## Q is an orthogonal n x n matrix
  ## R is an upper triangular n x m matrix
  ## Q and R satisfy the equation: A = Q %*% R
  return(list("Q" = t(Q), "R" = R))
  
}



#################################
## Function 2: Sweep operation ##
#################################

mySweep <- function(A, m){
  
  # Perform a SWEEP operation on A with the pivot element A[m,m].
  # 
  # A: a square matrix.
  # m: the pivot element is A[m, m].
  # Returns a swept matrix.
  
  ## Leave this function as is unless you want to make it 
  ## more efficient!
  
  n <- nrow(A)
  
  for(k in 1:m){ 
    for(i in 1:n)     
      for(j in 1:n)   
        if(i != k  & j != k)     
          A[i,j] <- A[i,j] - A[i,k]*A[k,j]/A[k,k]    
        
        for(i in 1:n) 
          if(i != k) 
            A[i,k] <- A[i,k]/A[k,k]  
          
          for(j in 1:n) 
            if(j != k) 
              A[k,j] <- A[k,j]/A[k,k]
            
            A[k,k] <- - 1/A[k,k]
  }
  
  return(A)
  
}




###############################################
## Function 3: Linear regression based on QR ##
###############################################

myLM <- function(X, Y){
  
  ## Perform the linear regression of Y on X
  ## Input: 
  ## X is an n x p matrix of explanatory variables
  ## Y is an n dimensional vector of responses
  ## Use myQR inside of this function
  
  ########################
  ## FILL IN CODE BELOW ##
  ########################  
  
  n <- nrow(X)
  p <- ncol(X)
  
  ## Stack (X, Y) and solve it by QR decomposition
  X1 <- cbind(rep(1,n),X)
  Z <- cbind(X1, Y)
  R <- myQR(Z)$R
  
  R1 <- R[1:(p + 1), 1:(p + 1)]
  Y1 <- R[1:(p + 1), p + 2]
  
  
  beta_ls <- solve(R1) %*% Y1
  err <- (sum((Y - X1%*%beta_ls)^2)/n)
  
  
  ## Function returns the 1 x (p + 1) vector beta_ls, 
  ## the least squares solution vector
  return(list("beta" = beta_ls, "error" = err))
  
}


##################################
## Function 4: PCA based on QR  ##
##################################

myEigen_QR <- function(A, numIter = 1000){
  
  ## Perform PCA on matrix A using your QR function, myQRC.
  ## Input:
  ## A: Square matrix
  ## numIter: Number of iterations
  
  ########################
  ## FILL IN CODE BELOW ##
  ########################  
  
  r <- nrow(A)
  c <- ncol(A)
  V <- matrix(runif(r * r), r, r)
  
  for(i in 1:numIter){
    Q <- myQR(V)$Q
    V <- A %*% Q
  }
  
  eigen_out <- myQR(V)
  Q <- eigen_out$Q
  R <- eigen_out$R
  
  ## Function should output a list with D and V
  ## D is a vector of eigenvalues of A
  ## V is the matrix of eigenvectors of A (in the 
  ## same order as the eigenvalues in D.)
  return(list("D" = diag(R), "V" = Q))
  
}








######################################
## Function 5: Logistic regression  ##
######################################

## Expit/sigmoid function
expit <- function(x){
  1 / (1 + exp(-x))
}

myLogistic <- function(X, Y, epsilon = 1e-6){
  
  
  n <- nrow(X)
  X = cbind(rep(1,n),X)
  p <- ncol(X)    
  
  beta <- matrix(rep(0, p), nrow = p)
  epsilon <- 1e-6
  repeat
  {
    eta <- X%*%beta
    pr <- expit(eta)
    w <- pr*(1-pr)
    Z <- eta + (Y-pr)/w
    sw <- sqrt(w)
    mw <- matrix(sw, n, p)
    Xwork <- mw*X
    Ywork <- sw*Z
    Z1 = cbind( Xwork, Ywork)
    R = myQR(Z1)$R
    R1 = R[1:(p), 1:(p)]
    Y1 = R[1:(p), p+1]
    beta_new = solve(R1, Y1)
    err <- sum(abs(beta_new-beta))
    beta <- beta_new
    if (err<epsilon)
      break
  }
  beta_logistic = beta
  
  ## Function returns beta_logistic, the solution to 
  ## the logistic regression problem
  return(beta_logistic)
  
  
}


##################################
## Function 6: Ridge regression ##
##################################

myRidge <- function(X, Y, lambda){
  
  # Perform ridge regression of Y on X.
  # 
  # X: an n x p matrix of explanatory variables.
  # Y: an n vector of dependent variables. Y can also be a 
  # matrix, as long as the function works.
  # lambda: regularization parameter (lambda >= 0)
  # Returns beta, the ridge regression solution.

  n <- dim(X)[1]
  p <- dim(X)[2]
  
  ## Define cross product matrix
  Z <- cbind(rep(1, n), X, Y)
  A <- t(Z) %*% Z
  
  ## Add ridge to cross product matrix
  D <- diag(rep(lambda, p + 2))
  D[p + 2, p + 2] <- 0
  
  ## Don't regularize intercept
  D[1, 1] <- 0
  A <- A + D
  S <- mySweep(A, p + 1)
  beta_ridge <- S[1:(p + 1), p + 2]
  
  ## Function should output the vector beta_ridge, the 
  ## solution to the ridge regression problem. beta_ridge
  ## should have p + 1 elements.
  return(beta_ridge)
  
}


####################################################
## Function 7: Piecewise linear spline regression ##
####################################################


mySpline <- function(x, Y, lambda, p = 100){
  
  # Perform spline regression of Y on X.
  # 
  # x: An n x 1 vector or matrix of explanatory variables.
  # Y: An n x 1 vector of dependent variables. Y can also be a 
  # matrix, as long as the function works.
  # lambda: regularization parameter (lambda >= 0)
  # p: Number of cuts to make to the x-axis.
  
  ##################################
  ## FILL IN THIS SECTION OF CODE ##
  ##################################
  
  n <- length(x)
  X <- matrix(x, nrow = n)
  for(k in (1:(p-1))/p){
    X <- cbind(X, (x > k) * (x - k))
  }
  beta_spline <- myRidge(X, Y, lambda)
  Yhat <- cbind(rep(1,n), X) %*% beta_spline
  # plot(x, Y, ylim = c(-.2, 1.2), col = "red")
  # par(new = TRUE)
  # plot(x, Yhat, ylim = c(-.2, 1.2), type = 'l', col = "green")
  # 
  
  ## Function should a list containing two elements:
  ## The first element of the list is the spline regression
  ## beta vector, which should be p + 1 dimensional (here, 
  ## p is the number of cuts we made to the x-axis).
  ## The second element is y.hat, the predicted Y values
  ## using the spline regression beta vector. This 
  ## can be a numeric vector or matrix.
  output <- list(beta_spline = beta_spline, predicted_y = Yhat)
  return(output)
  
}



#####################################
## Function 8: Lasso solution path ##
#####################################

myLasso <- function(X, Y, lambda_all){
  
  # Find the lasso solution path for various values of 
  # the regularization parameter lambda.
  # 
  # X: n x p matrix of explanatory variables.
  # Y: n dimensional response vector
  # lambda_all: Vector of regularization parameters. Make sure 
  # to sort lambda_all in decreasing order for efficiency.
  #
  # Returns a matrix containing the lasso solution vector 
  # beta for each regularization parameter.
  
  n = length(Y)
  X1 = cbind(rep(1,n),X)
  p = ncol(X1)
  
  
  beta = matrix(rep(0, p), nrow = p)
  L = length(lambda_all)
  beta_all = matrix(rep(0, p*L), nrow = p)
  lambda_all = sort(lambda_all, decreasing = TRUE)
  
  
  
  R = Y
  ss = rep(0, p)
  for (j in 1:p)
    ss[j] = sum(X1[, j]^2)
  for (l in 1:L)
  {
    lambda = lambda_all[l]
    for (t in 1:100)
    {
      for (j in 1:p)
      {
        db = sum(R*X1[, j])/ss[j]
        b = beta[j]+db
        b = sign(b)*max(0, abs(b)-lambda/ss[j])
        db = b - beta[j]
        R = R - X1[, j]*db   # + - *
        beta[j] = b
      }
    }
    beta_all[, l] = beta
    colnames(beta_all) = lambda_all
  }
  
  
  err_all = matrix(rep(0, L),ncol = 1)
  for(l in 1:L){
    err_all[l,1] <- (sum((Y - X1%*%beta_all[,l])^2)/n)#^0.5
  }
  
  
  ## Function should output the matrix beta_all, the 
  ## solution to the lasso regression problem for all
  ## the regularization parameters. 
  ## beta_all is (p+1) x length(lambda_all)
  return(list("beta_all" = beta_all, "error_all" = err_all))
  
  
}







package.skeleton(name = "MyRegressionsW9")



