
# TO DO: 
  # Figure out good prior parameters for b_j, sigmasq_t. 
  


NWR_sampler <- function(y, X, A, n_draws=1e4, mu_t, sigmasq_t, a_j, b_j, a_e, b_e, mu_rhoj, var_rhoj, mu_rhoe, var_rhoe, tpj, tpe){
  library(Matrix) 
  #if(class(A == "this")) A <- that
  #"network"
  #"igraph"
  #"matrix"
  p <- ncol(X)
  n <- nrow(X)
  I <- diag(x=1,nrow=n, ncol = n)
  
  ## Set up posterior draws objects:
  if(missing(a_e)) a_e <- 1
  if(missing(b_e)) b_e <- var(y) * (a - 2) # E(\sigma^2) = a/2 / (b/2 - 1)
  if(missing(a_j)) a_j <- rep(1, p)
  #if(missing(b_j)) b_j <- var(y)a_j 
  

  
  beta_draws = array(0.0, dim=c(n, p, n_draws+1))
  beta_tilde_draws = matrix(0.0, n_draws+1, p)
  sigma_ep_draws = rep(1, n_draws+1)  
  sigma_j_draws = matrix(1, n_draws+1, p) 
  rho_ep_draws = numeric(n_draws+1)
  rho_j_draws = matrix(0.0, n_draws+1, p)
  Omega_j_draws = array(0.0, dim=c(n,n,p,n_draws+1))
  Omega_e_draws = array(0.0, dim=c(n,n,n_draws+1))
  
  
  # Start with initial values of Beta and Beta_tilde
    # Beta tilde is the OLS estimates and each subject's Beta will be the leave one out estimator
  lmfit = lm(y~X)
  beta_tilde_draws[1,] = coef(lmfit)
  e = resid(lmfit)
  
  for (i in 1:n){
    XtXinv <- qr.solve(crossprod(X))
    hi = X[i,] %*% tcrossprod(XtXinv, X[i,])
    leave_out = beta_tilde_draws[1,] - t(tcrossprod(XtXinv, X[i,])*(e[i]/(1-hi))) 
    beta_draws[i,,1] = n*beta_tilde_draws[1,] - (n-1)*leave_out
  }
  
  
  for (i in 2:(n_draws+1)){  
    
    
    # I need to simulate each rho_e -> sigma_e ->  and then each rho_j -> sigma_j -> beta_tilde_j -> beta_j 
    
    beta_draws[,,i] = beta_draws[,,i-1]
    beta_tilde_draws[i,] = beta_tilde_draws[i-1,]
    sigma_ep_draws[i] = sigma_ep_draws[i-1]
    sigma_j_draws[i,] = sigma_j_draws[i-1,]
    rho_ep_draws[i] = rho_ep_draws[i-1]
    rho_j_draws[i,] = rho_j_draws[i-1,] 
    
    M <- rowsum(X * beta_draws[,,i])
    
    # Draw all rhos using MH
    rho = c(rho_ep_draws[i], rho_j_draws[i,])
    
    
    ## Draw rho_ep:
    
    sqrtOmege = qr.solve(I - rho_ep_draws[i] * A)
    Sigma_e = tcrossprod(sqrtOmege)
    Omega_e_draws[,,i] = qr.solve(Sigma_e)
    
    rhoestar = truncdist::rtrunc(1, spec="norm", a=-1, b=1, mean = rho_ep_draws[i], sd=tpe) 
    uj = runif(1)
    Sigmarhoestar = tcrossprod(qr.solve(I - rhoestar*A))
    
    checkfuncnum = mvtnorm::dmvnorm(y, mean = M, 
                                    sigma = sigma_ep_draws[i]*Sigmarhoestar)*
      truncdist::dtrunc(rhoestar, spec = "norm", a=-1, b=1, sd = sqrt(var_rhoe))*
      truncdist::dtrunc(rho_ep_draws[i], spec = "norm", a=-1, b=1, mean = rhoestar, sd = tpe)
    
    checkfuncden = mvtnorm::dmvnorm(y, mean = M, 
                                    sigma = sigma_ep_draws[i]*Sigma_e)*
      truncdist::dtrunc(rho_ep_draws[i], spec = "norm", a=-1, b=1, sd = sqrt(var_rhoe))*
      truncdist::dtrunc(rhoestar, spec = "norm", a=-1, b=1, mean = rho_ep_draws[i], sd = tpe)
    
    if(uj < checkfuncnum/checkfuncden){
      rho_ep_draws[i] = rhoestar
      Omega_e_draws[,,i] = qr.solve(Sigmarhoestar)
    }
    
    
    
    ## Draw sigma_ep:
    
    sigmamu = as.matrix(y - M)
    
    sigma_ep_draws[i] = invgamma::rinvgamma(1, 
                                            shape = (a_e+n)/2,
                                            scale = .5*(b_e + crossprod(sigmamu, Omega_e_draws[,,i]) %*% sigmamu)
    )
    
    
    
    for (j in 1:p){
      
      ## Draw rho_j:
      
      sqrtOmegj = qr.solve(I - rho_j_draws[i,j] * A)
      Sigma_j = tcrossprod(sqrtOmegj)
      Omega_j[,,j,i] = qr.solve(Sigma_j)
      
      rhojstar = truncdist::rtrunc(1, spec="norm", a=-1, b=1, mean = rho_j_draws[i,j], sd = tpj) 
      uj = runif(1)
      Sigmarhojstar = tcrossprod(qr.solve(I - rhojstar*A))
      
      checkfuncnum = mvtnorm::dmvnorm(beta_draws[,j,i], mean = beta_tilde_draws[i,j]*rep(1,n), 
                                    sigma = sigma_j_draws[i,j]*Sigmarhojstar)*
        truncdist::dtrunc(rhojstar, spec = "norm", a=-1, b=1, sd = sqrt(var_rhoj[j]))*
        truncdist::dtrunc(rho_j_draws[i,j], spec = "norm", a=-1, b=1, mean = rhojstar, sd = tpj)
      
      checkfuncden = mvtnorm::dmvnorm(beta_draws[,j,i], mean = beta_tilde_draws[i,j]*rep(1,n), 
                                      sigma = sigma_j_draws[i,j]*Sigma_j)*
        truncdist::dtrunc(rho_j_draws[i,j], spec = "norm", a=-1, b=1, sd = sqrt(var_rhoj[j]))*
        truncdist::dtrunc(rhojstar, spec = "norm", a=-1, b=1, mean = rho_j_draws[i,j], sd = tpj)
      
      if(uj < checkfuncnum/checkfuncden){
        rho_j_draws[i,j] = rhojstar
        Omega_j[,,j,i] = qr.solve(Sigmarhojstar)
      }
         
    
    
    
    ## Draw sigma_j:
   
      sigma_j_draws[i,j] = invgamma::rinvgamma(1, 
                                              shape = (a_j[j]+n)/2,
                                              scale = .5*(b_j[j] + crossprod(as.matrix(beta_draws[,j,i] - beta_tilde_draws[i,j]), Omega_j[,,j,i]) %*% as.matrix(beta_draws[,j,i] - beta_tilde_draws[i,j]))
      )
    
    
    ## Draw Beta_tilde_j:
    
      sigmma_tilde = 1/((1/sigma_j_draws[i,j])*
                             crossprod(as.matrix(rep(1,n)), Omega_j_draws[,,j,i]) %*% as.matrix(rep(1,n)) +
                          1/sigmasq_t)
      mu_tilde = sigmma_tilde * ((1/sigma_j_draws[i,j])*crossprod(beta_draws[,j,i], Omega_j_draws[,,j,i]) %*% 
                                   rep(1,n) + (1/sigmasq_t)*mu_t)
      beta_tilde_draws[i,j] = rnorm(mu_tilde, sigmma_tilde)
    
    ## Draw Beta_j:
    
      y_check = y 
      for (k in 1:p){
        if (k != j){
          y_check = y_check - X[,k]*beta_draws[,k,i]
        }
      }
      oneoversige = (1/sigma_ep_draws[i])
      oneoversigj = (1/sigma_j_draws[i,j])
      Omega_sig_ei = oneoversige*Omega_e_draws[,,i]
      Omega_sig_ji = oneoversigj*Omega_j_draws[,,j,i]
      diagxj = diag(x=X[,j])
      diagxj_Omega = diagxj %*% Omega_sig_ei
      
      sigma_j = qr.solve(diagxj_Omega %*% diagxj - Omega_sig_ji)
      mu_j = sigma_j %*% (diagxj_Omega %*% as.matrix(y_check) - 
                            Omega_sig_ji %*% as.matrix(beta_tilde_draws[i,j] * rep(1, n)))
      beta_draws[,j,i] = mvtnorm::rmvnorm(n=n, mean = mu_j, sigma = sigma_j)
    }
  }
  list(beta = beta_draws[,,2:(n_draws + 1)], betatilde = beta_tilde_draws[2:(n_draws + 1),], sigmaj = sigma_j_draws[2:(n_draws + 1),], sigmae = sigma_ep_draws[2:(n_draws + 1)], rhoj = rho_j_draws[2:(n_draws + 1),], rhoe = rho_ep_draws[2:(n_draws + 1)])
}





