#' summary.NWR -  A function that summarizes the results of an NWR object. 
#' @param NWR An NWR object.
#' @param ci The level of the credible interval.
#' @return A list with values:
#' @return \code{Parameters} The estimated parameters (with mean of each coefficient vector) as well as their 100*ci% credible intervals.
#' @return \code{logpost} The log of the posterior of the estimated parameters. 
#' @return \code{BIC} The BIC of the model. 
#' @return \code{AICc} The corrected AIC of the model. Corrected AIC is suggested to be used when comparing the NWR model to other models in order to account for the bias of AIC which doesn't penalize models with too many parameters as much as it should. Doing so is important for a model with a matrix of parameters of unspecified size.  
#' @export
summary.NWR <- function(NWR, ci = 0.95){
  if(class(NWR) != "NWR"){
    stop("Only takes NWR Objects")
  }
  n <- NROW(NWR$beta)
  p <- NCOL(NWR$beta)
  k <- NROW(NWR$gamma)
  
  pecent <- paste0(as.character(100*ci), '%')
  
  z <- -qnorm((1-ci)/2)
  
  hinvd <- sqrt(-diag(NWR$Hinv))
  
  c1 <- matrix(rep(1, n), ncol = 1)
  
  betase <- c()
  for (i in 1:p){
    dims <- ((i-1)*k+1):(i*k)
    betase[i] <- sqrt(diag(-(1/(n^2)) * crossprod(c1, NWR$U) %*% tcrossprod(NWR$Hinv[dims,dims], NWR$U) %*% c1))
  }
  
  
  partab <- cbind(c(NWR$beta0, NWR$rho, NWR$s2),
                  c(NWR$beta0 - z*hinvd[k*p + 1],
                    NWR$rho - z*hinvd[k*p + 3],
                    NWR$s2 - z*hinvd[k*p + 2]),
                  c(NWR$beta0 + z*hinvd[k*p + 1],
                    NWR$rho + z*hinvd[k*p + 3],
                    NWR$s2 + z*hinvd[k*p + 2]))
  colnames(partab) <- c("Estimate", paste(pecent, "LB"), paste(pecent, "UB"))
  rownames(partab) <- c("Intercept", "Rho", "Sigma^2")
  
  
  gtab = matrix(NA, nrow = NROW(NWR$gamma), ncol = 1)
  for(i in 1:p){
    gtab = cbind(gtab, NWR$gamma[,i], 
                 NWR$gamma[,i] - z*hinvd[((i-1)*k+1):((i-1)*k+k)],
                 NWR$gamma[,i] + z*hinvd[((i-1)*k+1):((i-1)*k+k)])
  }
  gtab <- gtab[,-1]
  
  #fullbetatab <- NWR$U %*% gtab
  #betatab <- matrix(colMeans(fullbetatab), ncol = 3, byrow = TRUE)
  
  betatab <- matrix(c(colMeans(NWR$beta), colMeans(NWR$beta) - 1.96*betase, colMeans(NWR$beta) + 1.96*betase), nrow = p)
  rownames(betatab) <- colnames(NWR$X)
  
  
  partab <- round(rbind(partab[1,], betatab, partab[2:3,]), 5)
  rownames(partab)[1] <- "(Intercept)"
  
  
  numpars <- (k*p + 3)
  
  BIC <- -2*NWR$logpost + numpars*log(n) 
  AICc <- -2*NWR$logpost + 2*numpars + (2*numpars*(numpars+1))/(n-numpars+1)
  
  list(Parameters = partab, logpost = NWR$logpost, AICc = AICc, BIC = BIC)
}