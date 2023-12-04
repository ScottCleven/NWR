#' NWR -  A function to run network weighted regression models by performing optimization of the posterior distribution.
#' @param formula A formula line usually denoted y ~ X or outcome ~ predictors.
#' @param data The dataset the variables from the formula line come from. 
#' @param y If you would rather input an outcome vector directly. Must be used in conjunction with X and without a formula.
#' @param X If you would rather input a predictor matrix directly. DO NOT include an intercept column, the intercept is modeled differently from the other predictors. Must be used in conjunction with y and without a formula.
#' @param A An adjacency matrix describing your network. It will be row-normalized in the function if it is not already.
#' @param k The number of important eigenvectors of A. The true value of this number is probably unknown to you but you can test for optimal k using the model selection method of your choice and the NWR_kstep function. The larger your k, the longer the function will take to run (the more parameters there will need to be estimated). A good place to start might be the number of predictors.
#' @param gamma_prior The prior for each gamma vector of the model. Requires list or function inputs. Currently requires all vectors to have the same prior. Default is a multivariate normal with mean=0 and sigma=I where I is the identity matrix. If you supply a list, you must have one value named 'mean' which is the mean of the multivariate normal and one named 'sigma' which is the variance-covariance matrix. If mean is a single integer, that mean will be used for every element in the mean vector and if sigma is a single integer, it will multiply that integer by the identity matrix. Otherwise, both mean and sigma will be what you specify them to be. Keep in mind, if you're adjacency matrix is not symmetrical (when it isn't row-normalized), you will need the vectors to be 2(k+1)p in dimension. And if it is symmetrical, they need to be (k+1)p in dimension. You can also supply your own prior by supplying a function name of the prior function you would like to use.  
#' @param beta0_prior The prior for the intercept of the model. Requires list or function inputs. Default is a N(0,1). If you supply a list, you must have one value named 'mean' which is the mean of the normal and one named 'sd' which is the standard deviation with both values of length=1. You can also supply your own prior by supplying a function name of the prior function you would like to use.
#' @param s2_prior The prior for the variance of the model. Requires list or function inputs. Default is an inverse-gamma(0.5,0.5). If you supply a list, you must have one value named 'shape' which is the shape of the inverse-gamma and one named 'scale' which is the scale with both values of length=1. You can also supply your own prior by supplying a function name of the prior function you would like to use.
#' @param rho_prior The prior for the correlation of the adjacency matrix of the model. Requires list or function inputs. Default is a N(0.36,0.49) (Dittrich et al. 2017). If you supply a list, you must have one value named 'mean' which is the mean of the normal and one named 'sd' which is the standard deviation with both values of length=1. You can also supply your own prior by supplying a function name of the prior function you would like to use.
#' @param ... Values to be passed into the optim() function. Currently defaults to: optim(par=c(rep(0, length(gamma_mat)),1,0.2), fn=posterior_function, gr=NULL, method="L-BFGS-B", lower=c(rep(-15, length(gamma_mat)), 1e-3, -.85), upper=c(rep(15, length(gamma_mat)), 20, .85), control=list(fnscale=-1), hessian=TRUE). I don't allow you to change fn or hessian.  
#' @return A list containing values:
#' @return \code{beta0} The estimated intercept of the model.
#' @return \code{beta} The estimated coefficient matrix of the model.
#' @return \code{s2} The estimated variance of the model.
#' @return \code{rho} The estimated correlation of the adjacency matrix of the model.
#' @return \code{gamma} The estimated gamma parameters of the model.
#' @return \code{U} The Orthonormal matrix of k important eigenvectors of the spectral or QR decomposition of the model (Plus intercept column).
#' @return \code{A} The inputted adjacency matrix (row-normalized).
#' @return \code{H} The hessian matrix of the optimization.
#' @return \code{Hinv} The inverse of the hessian matrix.
#' @return \code{logpost} The log of the posterior density.
#' @return \code{y} The inputted outcome vector.
#' @return \code{X} The inputted prediction matrix. 
#' @export
NWR <- function(formula, data, y=FALSE, X=FALSE, A,
                k = 5,
                gamma_prior = list(mean=0, sigma = 1),
                beta0_prior = list(mean=0, sd=1),
                s2_prior = list(shape=1/2, scale=1/2),
                rho_prior = list(mean=0.36, sd=0.7^2),
                ...){
  ######### Getting the Data Prepared for optimization #############
  if(y==0 & X==0){
    X = model.matrix(formula, data=data)[,-1]
    y = unlist(data[as.character(formula[[2]])])
  }else if(y==0 & X==1){
    stop("You need an outcome y")
  }else if(X==0){
    stop("You need covariates X. If you want an intercept only model use lm.")
  }
  
  rownorm <- function(networkmat){
    networkmat/ifelse(rowSums(networkmat) == 0,1,rowSums(networkmat))
  }
  n = NROW(y)
  
  if(isSymmetric(A)){
    U = Re(as.matrix(cbind(rep(1,n), eigen(A)$vectors[,1:k])))
    A = rownorm(A)
    k <- ncol(U)
  }else if(sum(rowSums(A) %in% c(0,1)) == nrow(A)){
    A_backtrans = diag(rowSums(A > 0)) %*% A
    rownames(A_backtrans) <- colnames(A_backtrans)
    if(isSymmetric(A_backtrans)){  
      U <- Re(as.matrix(cbind(rep(1,n), eigen(A_backtrans)$vectors[,1:k])))
      k <- ncol(U)
    }else{
      svdU <- svd(A_backtrans)
      U <- cbind(rep(1,n), svdU$u[,1:k], svdU$v[,1:k])
      k <- ncol(U)
    }
  }else{
    svdU <- svd(A)
    U <- cbind(rep(1,n), svdU$u[,1:k], svdU$v[,1:k])
    k <- ncol(U)
  }
  
  
  NWR_dat = list(X=X, y=y, A=A, U=U)
  
  ########## LOG priors checking and cleaning ########################
  
  if(typeof(gamma_prior) == "list"){
    if(length(gamma_prior$mean) == 1){
      gmean <- rep(gamma_prior$mean, k)
    }else{
      gmean <- gamma_prior$mean
    }
    
    if(length(gamma_prior$sigma) == 1){
      gsigma <- diag(k) * gamma_prior$sigma
    }else{
      gsigma <- gamma_prior$sigma
    }
    
    gprior <- function(gamma){
      -1/2 * log(det(gsigma)) - 1/2 * t(gamma - gmean) %*% qr.solve(gsigma) %*% gamma
    }
  }else{
    gprior <- gamma_prior
  }
  
  if(typeof(beta0_prior) == "list"){
    b0prior <- function(beta0){
      -1 * log(beta0_prior$sd) - 1/2 * (beta0 - beta0_prior$mean)^2 / beta0_prior$sd^2
    }
  }else{
    b0prior <- beta0_prior
  }
  
  if(typeof(s2_prior) == "list"){
    s2prior <- function(s2){
      (-s2_prior$shape - 1)*s2 - s2_prior$scale / s2
    }
  }else{
    s2prior <- s2_prior
  }
  
  if(typeof(rho_prior) == "list"){
    rhoprior <- function(rho){
      -log(rho_prior$sd) - 1/2 * (rho - rho_prior$mean)^2 / rho_prior$sd^2 - log(pnorm(1, mean = 0, sd = rho_prior$sd) - pnorm(-1, mean = 0, sd = rho_prior$sd))
    }
  }else{
    rhoprior <- rho_prior
  }
  
  ######## Function that creates the log posterior and optimizes it ###########
  
  NWR_opt <- function(NWR_dat){
    y = NWR_dat$y
    X = NWR_dat$X
    U = NWR_dat$U
    A = NWR_dat$A
    
    n <- length(y)
    I <- diag(n)
    k <- ncol(U)
    p <- ncol(X)
    Ut <- t(U)
    b1 <- matrix(rep(1, n))
    At <- t(A)
    
    
    #theta will be parameterized as the first p*(k+1) rows are the gamma_j vectors: 
    # i.e. 1:k will be gamma_1 and (k+1):2k will be gamma_2 and so on.
    # then beta0, sigma^2, and rho. So I will have p*k + 3 rows.
    # For the record, k is now the number of columns of U, not the k from the input. 
    # I did this for simplicity and because it's a function of the k from the input. 
    logpost <- function(theta){
      #Likelihood*(1:p gamma_j priors)*(beta0 prior)*(sigma^2 prior)*(rho prior)
      
      # The prior for gamma is based on each column of gamma so I go through them individually
      priorg = 0
      for(i in 1:p){
        priorg = priorg + gprior(theta[((i-1)*k+1):(i*k)])
      }
      nbeta = U %*% matrix(theta[1:(k*p)], nrow = k, ncol = p)
      m = theta[p*k + 1] + rowSums(X * nbeta)
      
      #loglik
      mvtnorm::dmvnorm(y,
                       m,
                       (theta[p*k + 2])*tcrossprod(qr.solve(I - theta[p*k + 3]*A)),
                       log=T) +
        #gamma
        priorg +
        #beta0
        b0prior(theta[p*k + 1]) + 
        #s2
        s2prior(theta[p*k + 2]) +
        #rho
        rhoprior(theta[p*k + 3])
    }
    
    
    optim_pars <- function(...){
      
      mylist <- list(...)
      
      if(is.null(mylist$par)){
        mylist$par <- c(rep(0, p*k+1),1,0.2)
      }
      if(is.null(mylist$lower)){
        mylist$lower <- c(rep(-15, p*k+1), 1e-3, -.85)
      }
      if(is.null(mylist$upper)){
        mylist$upper <- c(rep(15, p*k+1), 20, .85)
      }
      if(is.null(mylist$control)){
        mylist$control <- list(fnscale=-1)
      }
      if(is.null(mylist$method)){
        mylist$method <- "L-BFGS-B"
      }
      mylist
    }
    
    opars <- optim_pars(...)
    
    
    myopt <- optim(par = opars$par,
                   fn = logpost, gr=opars$gr, method = opars$method, control = opars$control, lower = opars$lower, upper = opars$upper, hessian = TRUE)
    
    if(opars$method == "L-BFGS-B"){
      # Just checking if you hit a boundary
      if(1 %in% (myopt$par == opars$lower)){
        index <- which(1 == (myopt$par == opars$lower))
        if(length(index) > 1){
          if(sum(index > k*p) == length(index)){
            warning("You hit the lower boundary on multiple of beta0, sigma^2, or rho. Use/change the 'lower=' argument.")
          }else if(sum(index <= k*p) == length(index)){
            warning("You hit the lower boundary on multiple gamma parameters. Use/change the 'lower=' argument.")
          }else{
            warning("You hit the lower boundary on many different parameters. Use/change the 'lower=' argument.")
          }
        } else if(index == k*p + 1){
          warning("You hit the lower boundary on the intercept parameter. Use/change the 'lower=' argument.")
        }else if(index == k*p + 2){
          warning("You hit the lower boundary on the sigma^2 parameter. Use/change the 'lower=' argument.")
        }else if(index == k*p + 3){
          warning("You hit the lower boundary on the rho parameter. Use/change the 'lower=' argument.")
        }else{
          warning("You hit the lower boundary on one of the gamma parameters. Use/change the 'lower=' argument.")
        }
      }
      
      
      if(1 %in% (myopt$par == opars$upper)){
        index <- which(1 == (myopt$par == opars$upper))
        if(length(index) > 1){
          if(sum(index > k*p) == length(index)){
            warning("You hit the upper boundary on multiple of beta0, sigma^2, or rho. Use/change the 'upper=' argument.")
          }else if(sum(index <= k*p) == length(index)){
            warning("You hit the upper boundary on multiple gamma parameters. Use/change the 'upper=' argument.")
          }else{
            warning("You hit the upper boundary on many different parameters. Use/change the 'upper=' argument.")
          }
        } else if(index == k*p + 1){
          warning("You hit the upper boundary on the intercept parameter. Use/change the 'upper=' argument.")
        }else if(index == k*p + 2){
          warning("You hit the upper boundary on the sigma^2 parameter. Use/change the 'upper=' argument.")
        }else if(index == k*p + 3){
          warning("You hit the upper boundary on the rho parameter. Use/change the 'upper=' argument.")
        }else{
          warning("You hit the upper boundary on one of the gamma parameters. Use/change the 'upper=' argument.")
        }
      }
    }  
    
    beta = U %*% matrix(myopt$par[1:(k*p)], ncol=p)
    colnames(beta) <- colnames(X)
    
    list(beta0 = myopt$par[k*p + 1], beta = beta, s2 = myopt$par[k*p + 2],
         rho = myopt$par[k*p + 3], gamma = matrix(myopt$par[1:(k*p)], ncol=p), U = U, A = A, H = myopt$hessian, Hinv = qr.solve(myopt$hessian), logpost = myopt$value,
         y = y, X = X)
  }
  
  out <- NWR_opt(NWR_dat)
  class(out) <- "NWR"
  out
}





