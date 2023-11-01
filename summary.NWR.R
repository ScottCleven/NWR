summary.NWR <- function(NWR_dat, ci = 0.95){
  if(class(NWR_dat) != "NWR"){
    stop("Only takes NWR Objects")
  }
  n <- NROW(NWR_dat$beta)
  p <- NCOL(NWR_dat$beta)
  k <- NROW(NWR_dat$gamma)
  
  pecent <- paste0(as.character(100*ci), '%')
  
  z <- -qnorm((1-ci)/2)
  
  hinvd <- sqrt(-diag(NWR_dat$Hinv))
  
  c1 <- matrix(rep(1, n), ncol = 1)
  
  betase <- c()
  for (i in 1:p){
    dims <- ((i-1)*k+1):(i*k)
    betase[i] <- sqrt(diag(-(1/(n^2)) * crossprod(c1, NWR_dat$U) %*% tcrossprod(NWR_dat$Hinv[dims,dims], NWR_dat$U) %*% c1))
  }
  
  
  partab <- cbind(c(NWR_dat$beta0, NWR_dat$rho, NWR_dat$s2),
                  c(NWR_dat$beta0 - z*hinvd[k*p + 1],
                    NWR_dat$rho - z*hinvd[k*p + 3],
                    NWR_dat$s2 - z*hinvd[k*p + 2]),
                  c(NWR_dat$beta0 + z*hinvd[k*p + 1],
                    NWR_dat$rho + z*hinvd[k*p + 3],
                    NWR_dat$s2 + z*hinvd[k*p + 2]))
  colnames(partab) <- c("Estimate", paste(pecent, "LB"), paste(pecent, "UB"))
  rownames(partab) <- c("Intercept", "Rho", "Sigma^2")
  
  
  gtab = matrix(NA, nrow = NROW(NWR_dat$gamma), ncol = 1)
  for(i in 1:p){
    gtab = cbind(gtab, NWR_dat$gamma[,i], 
                 NWR_dat$gamma[,i] - z*hinvd[((i-1)*k+1):((i-1)*k+k)],
                 NWR_dat$gamma[,i] + z*hinvd[((i-1)*k+1):((i-1)*k+k)])
  }
  gtab <- gtab[,-1]
  
  #fullbetatab <- NWR_dat$U %*% gtab
  #betatab <- matrix(colMeans(fullbetatab), ncol = 3, byrow = TRUE)
  
  betatab <- matrix(c(colMeans(NWR_dat$beta), colMeans(NWR_dat$beta) - 1.96*betase, colMeans(NWR_dat$beta) + 1.96*betase), nrow = p)
  rownames(betatab) <- colnames(NWR_dat$X)
  
  
  partab <- round(rbind(partab[1,], betatab, partab[2:3,]), 5)
  rownames(partab)[1] <- "(Intercept)"
  
  
  numpars <- (k*p + 3)
  
  BIC <- -2*NWR_dat$loglik + numpars*log(n) 
  AICc <- -2*NWR_dat$loglik + 2*numpars + (2*numpars*(numpars+1))/(n-numpars+1)
  
  list(Parameters = partab, loglik = NWR_dat$loglik, AICc = AICc, BIC = BIC)
}