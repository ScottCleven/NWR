
NWR <- function(formula, data, y=FALSE, X=FALSE, A,
                k = 5,
                c = rep(0, k),
                D = diag(k+1),
                s = 1,
                m0 = 1,
                a = 1,
                b = 1,
                mr = 0,
                sr = 1){
  if(y==0 & X==0){
    X = model.matrix(formula, data=data)[,-1]
    y = unlist(data[as.character(formula[[2]])])
  }
  else if(y==0 & X==1){
    stop("You need an outcome y")
  }
  else if(X==0){
    stop("You need covariates X. If you want an intercept only model use lm.")
  }
  
  rownorm <- function(networkmat){
    networkmat/ifelse(rowSums(networkmat) == 0,1,rowSums(networkmat))
  }
  n = NROW(y)
  
  if(isSymmetric(A)){
    U = Re(as.matrix(cbind(rep(1,n), eigen(A)$vectors[,1:k])))
    A = rownorm(A)
  }
  else if(sum(rowSums(A) %in% c(0,1)) == nrow(A)){
    A_backtrans = diag(rowSums(A > 0)) %*% A
    rownames(A_backtrans) <- colnames(A_backtrans)
    if(isSymmetric(A_backtrans)){  
      U <- Re(as.matrix(cbind(rep(1,n), eigen(A_backtrans)$vectors[,1:k])))
    }
    else{
      svdU <- svd(A_backtrans)
      U <- cbind(rep(1,n), svdU$u[,1:k], svdU$v[,1:k])
    }
  }
  else{
    svdU <- svd(A)
    U <- cbind(rep(1,n), svdU$u[,1:k], svdU$v[,1:k])
  }
  
  NWR_dat = list(X=X, y=y, A=A, U=U, c=c, D=D, s=s, m0=m0, a=a, b=b, mr=mr, sr=sr)
  
  NWR_CD <- function(NWR_dat){
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
    
    #priors
    c <- rep(0, k)
    D <- diag(k)
    s <- 1
    m0 <- 1
    a <- 1
    b <- 1
    mr <- 0
    sr <- 1
    
    
    
    #theta will be parameterized as the first p*k rows are the gamma_j vectors: 
    # i.e. 1:k will be gamma_1 and (k+1):2K will be gamma_2 and so on.
    # then beta0, sigma^2, and rho. So I will have p*k + 3 rows. 
    logpost <- function(theta){
      #Likelihood*(1:p gamma_j priors)*(beta0 prior)*(sigma^2 prior)*(rho prior)
      
      nbeta = matrix(NA, nrow = n, ncol = p)
      priorg = 0
      for(i in 1:p){
        nbeta[,i] = (U %*% theta[((i-1)*k+1):(i*k)])
        priorg = priorg + -1/2 * log(det(D)) - 1/2 * t(theta[((i-1)*k+1):(i*k)] - c) %*% qr.solve(D) %*% theta[((i-1)*k+1):(i*k)]
      }
      m = theta[p*k + 1] + rowSums(X * nbeta)
      #loglik
      mvtnorm::dmvnorm(y,
                       m,
                       (theta[p*k + 2])*tcrossprod(qr.solve(I - theta[p*k + 3]*A)),
                       log=T) +
        #log(det((theta[p*k + 2])*tcrossprod(qr.solve(I - theta[p*k + 3]*A)))^(-1/2)) -
        #(1/(2*theta[p*k+2]))* (crossprod(matrix(y - m), crossprod(I - theta[p*k + 3]*A)) %*% (y-m)) + 
        #gamma
        priorg +
        #beta0
        -1 * log(s) - 1/2 * (theta[p*k + 1] - m0)^2 / s^2 +
        #s2
        (-a/2 - 1)*theta[p*k + 2] - b/2 / theta[p*k + 2] +
        #rho
        -log(sr) - 1/2 * (theta[p*k + 3])^2 / sr^2 - log(pnorm(1, mean = 0, sd = sr) - pnorm(-1, mean = 0, sd = sr))
    }
    # loggrad <- function(theta){
    #   out = c()
    #   OR = (I - theta[p*k+3]*t(A)) %*% (I - theta[p*k+3]*A)
    #   OE = (1/theta[p*k+2])*OR
    #   nbeta = matrix(NA, nrow = n, ncol = p)
    #   for(i in 1:p){
    #     nbeta[,i] = (U %*% theta[((i-1)*k+1):(i*k)])
    #   }
    #   elemp = X * nbeta
    #   m = theta[p*k + 1] + rowSums(elemp)
    #   av = y-m
    #   avt = t(av)
    #   
    #   
    #   for(i in 1:p){
    #     DX = diag(X[,i])
    #     out[((i-1)*k+1):(i*k)] = - Ut %*% DX %*% OE %*% DX %*% U %*% theta[((i-1)*k+1):(i*k)] +
    #       Ut %*% DX %*% OE %*% (y - m + elemp[,i])  +
    #       -1/2 * (qr.solve(t(D)) + qr.solve(D)) %*% theta[((i-1)*k+1):(i*k)] + qr.solve(t(D)) %*% c
    #   }
    #   out[p*k + 1] = -theta[p*k+1] * t(b1) %*% OR %*% b1 - t(b1) %*% OR %*% (m - y - theta[p*k + 1]) -
    #     theta[p*k + 1]/s^2 - m0/2
    #   out[p*k + 2] = -n/(2*theta[p*k+2]) + 1/(2*theta[p*k+2]^2) * avt %*% OR %*% av +
    #     b/(2*theta[p*k+2]^2) - a/2 - 1
    #   out[p*k + 3] = -psych::tr(A %*% qr.solve(I- theta[p*k+3]*A)) + 
    #     1/(2*out[p*k + 2]) * (avt %*% At %*% av + avt %*% A %*% av - 2*theta[p*k + 3] %*% avt %*% At %*% A %*% av) +
    #     theta[p*k + 3]/sr^2
    #   out
    # }
    myopt <- optim(par = c(rep(0, p*k+1),
                           1,
                           0.2),
                   fn = logpost, method = "L-BFGS-B", control = list(fnscale=-1), lower = c(rep(-15, p*k+1), 1e-3, -.85), upper = c(rep(15, p*k+1), 20, .85), hessian = TRUE)
    
    beta = U %*% matrix(myopt$par[1:(k*p)], ncol=p)
    colnames(beta) <- colnames(X)
    
    list(beta0 = myopt$par[k*p + 1], beta = beta, s2 = myopt$par[k*p + 2],
         rho = myopt$par[k*p + 3], gamma = matrix(myopt$par[1:(k*p)], ncol=p), U = U, A = A, H = myopt$hessian, Hinv = qr.solve(myopt$hessian), loglik = myopt$value,
         y = y, X = X)
  }
  
  out <- NWR_CD(NWR_dat)
  class(out) <- "NWR"
  out
}





