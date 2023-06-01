#######
## MOST OF THIS CODE IS TAKEN FROM:
##    https://github.com/cran/ordinalNet/tree/master/R
## with modifications to the functions that end with ".mine" postfix.
##
## The most critical modification is captured by the function "ordinalNet.mine()", 
## where we added three new parameters to allow for more tailored optimization schemes.
#######


## Coordinate descent inner loop function
# Note: should not need to check for improvement because each coordinate step necessarily improves the approximate objective
cdIn <- function(wtsum, betaHat, score, info, alpha, lambdaMod, positiveID, threshIn, maxiterIn)
{
  # Update Active Set
  activeSet <- which((betaHat!=0 | lambdaMod==0) & diag(info)!=0)
  # check lambdaMod==0 because if data are balanced, an intercept could be initialized exactly to zero
  betaHat[-activeSet] <- 0
  betaHatActive <- betaHat[activeSet]
  scoreActive <- score[activeSet]
  infoActive <- info[activeSet, activeSet, drop=FALSE]
  infoInactive <- info[-activeSet, activeSet, drop=FALSE]
  lambdaModActive <- lambdaMod[activeSet]
  lambdaModInactive <- lambdaMod[-activeSet]
  positiveIDActive <- positiveID[activeSet]
  positiveIDInactive <- positiveID[-activeSet]
  
  # softThreshTerms vector does not change during the inner loop, even if active set changes
  softThreshTerms <- c(info[, activeSet, drop=FALSE] %*% betaHatActive) + score
  # softThreshTerms = I(beta^(r)) %*% beta^(r) + U(beta^(r)) = I(beta^(r)) %*% beta^(r+1)
  softThreshTermsActive <- softThreshTerms[activeSet]
  softThreshTermsInactive <- softThreshTerms[-activeSet]
  
  # Initialize quadratic approximation to the log-likelihood and objective
  loglikApprox <- getLoglikApprox(betaHatActive, scoreActive, infoActive)
  penalty <- getPenalty(betaHat, lambdaMod, alpha)
  objApprox <- -loglikApprox/wtsum + penalty
  
  iterIn <- 0
  kktAll <- FALSE
  while (!kktAll && iterIn<maxiterIn)
  {
    conv <- FALSE
    while (!conv && iterIn<maxiterIn)
    {
      iterIn <- iterIn + 1
      for (i in seq_along(activeSet))
      {
        numTerm <- softThreshTermsActive[i] - sum(infoActive[i, -i, drop=FALSE] * betaHatActive[-i])
        denTerm <- infoActive[i, i]
        penTerm <- wtsum * lambdaModActive[i]
        betaHatActive[i] <- softThresh(numTerm, penTerm*alpha) / (denTerm + penTerm*(1-alpha))
        if (positiveIDActive[i]) betaHatActive[i] <- max(0, betaHatActive[i])
      }
      
      betaHatOld <- betaHat
      betaHat[activeSet] <- betaHatActive
      loglikApprox <- getLoglikApprox(betaHatActive, scoreActive, infoActive)
      penalty <- getPenalty(betaHat, lambdaMod, alpha)
      objApproxOld <- objApprox
      objApprox <- -loglikApprox / wtsum + penalty
      dif <- abs((objApproxOld - objApprox) / (abs(objApproxOld) + 1e-100))
      conv <- dif < threshIn
      
    }  # end while (!conv && iterIn<maxiterIn)
    
    kkt <- rep(TRUE, length(betaHat))
    kktInactiveTerms <- softThreshTermsInactive - c(infoInactive %*% betaHatActive)
    kktInactiveTerms[!positiveIDInactive] <- abs(kktInactiveTerms[!positiveIDInactive])
    kktInactive <- kktInactiveTerms <= wtsum * lambdaModInactive * alpha
    kkt[-activeSet] <- kktInactive
    kktAll <- all(kkt)
    if (!kktAll)
    {
      iterIn <- iterIn - 1  # repeat the iteration if kkt conditions are not satisfied
      activeSet <- union(activeSet, which(!kkt))
      betaHatActive <- betaHat[activeSet]
      scoreActive <- score[activeSet]
      infoActive <- info[activeSet, activeSet, drop=FALSE]
      infoInactive <- info[-activeSet, activeSet, drop=FALSE]
      softThreshTermsActive <- softThreshTerms[activeSet]
      softThreshTermsInactive <- softThreshTerms[-activeSet]
      lambdaModActive <- lambdaMod[activeSet]
      lambdaModInactive <- lambdaMod[-activeSet]
      positiveIDActive <- positiveID[activeSet]
      positiveIDInactive <- positiveID[-activeSet]
    }
    
  }  # end while (!kktAll && iterIn<maxiterIn)
  
  list(betaHat=betaHat, iterIn=iterIn)
}  # end cdIn




# Coordinate descent outer loop function
cdOut <- function(betaHat, lambdaIndex, lambdaNum, lambdaMod,
                  xList, xMat, yMat, alpha, positiveID, linkfun,
                  pMin, threshOut, threshIn, maxiterOut, maxiterIn,
                  printIter, printBeta)
{
  nObs <- nrow(yMat)
  wts <- if (is.null(attr(yMat, "wts"))) rowSums(yMat) else attr(yMat, "wts")
  wtsum <- if (is.null(attr(yMat, "wtsum"))) sum(wts) else attr(yMat, "wtsum")
  nLev <- ncol(yMat)
  nVar <- ncol(xMat)
  
  # Could carry over epsMat, pMat, and loglik from previous cdOut iteration
  betaNonzeroIndex <- which(betaHat != 0)
  etaMat <- matrix(xMat[, betaNonzeroIndex, drop=FALSE] %*% betaHat[betaNonzeroIndex], nrow=nObs, byrow=TRUE)
  pMat <- do.call(rbind, lapply(1:nrow(etaMat), function(i) linkfun$h(etaMat[i, ])))
  loglik <- getLoglik(pMat, yMat)
  penalty <- getPenalty(betaHat, lambdaMod, alpha)
  obj <- -loglik/wtsum + penalty
  
  conv <- FALSE
  iterOut <- 0
  while (!conv && iterOut<maxiterOut)
  {
    iterOut <- iterOut + 1
    
    # Update score and info
    si <- getScoreInfo(xList, yMat, pMat, pMin, linkfun)
    score <- si$score
    info <- si$info
    
    # Run coordinate descent inner loop
    betaHatOld <- betaHat
    cdInResult <- cdIn(wtsum, betaHat, score, info, alpha, lambdaMod, positiveID, threshIn, maxiterIn)
    betaHat <- cdInResult$betaHat
    iterIn <- cdInResult$iterIn
    betaNonzeroIndex <- which(betaHat != 0)
    etaMat <- matrix(xMat[, betaNonzeroIndex, drop=FALSE] %*% betaHat[betaNonzeroIndex], nrow=nObs, byrow=TRUE)
    pMat <- do.call(rbind, lapply(1:nrow(etaMat), function(i) linkfun$h(etaMat[i, ])))
    
    # Update log-likelihood and objective
    loglikOld <- loglik
    objOld <- obj
    loglik <- getLoglik(pMat, yMat)
    penalty <- getPenalty(betaHat, lambdaMod, alpha)
    obj <- -loglik / wtsum + penalty
    
    # Take half steps if obj does not improve. Loglik is set to -Inf
    # if any fitted probabilities are negative, which can happen for
    # the nonparallel or semiparallel cumulative probability model.
    nhalf <- 0
    while (obj > objOld && nhalf < 10) {
      nhalf <- nhalf + 1
      betaHat <- (betaHat + betaHatOld) / 2
      betaNonzeroIndex <- which(betaHat != 0)
      etaMat <- matrix(xMat[, betaNonzeroIndex, drop=FALSE] %*% betaHat[betaNonzeroIndex], nrow=nObs, byrow=TRUE)
      pMat <- do.call(rbind, lapply(1:nrow(etaMat), function(i) linkfun$h(etaMat[i, ])))
      loglik <- getLoglik(pMat, yMat)
      penalty <- getPenalty(betaHat, lambdaMod, alpha)
      obj <- -loglik / wtsum + penalty
    }
    dif <- (objOld - obj) / (abs(objOld) + 1e-100)
    conv <- dif < threshOut
    
    # Convergence is declared if objective worsens. In this case, set betaHat
    # to previous value. (Typically means model is saturated.)
    objImproved <- obj <= objOld
    if (!objImproved)
    {
      betaHat <- betaHatOld
      loglik <- loglikOld
    }
    
    # Print iteration info if printIter=TRUE
    if (printIter)
    {
      if (iterOut == 1) cat("\nLambda", lambdaIndex, " of ", lambdaNum, '\n')
      cat("outer iteration ", iterOut, ":  ", iterIn,
          " inner iterations, relative change in objective: ",
          signif(dif, 2), '\n', sep='')
    }
    
    # Print betaHat if printBeta=TRUE
    if (printBeta)
    {
      if (!printIter)
      {
        if (iterOut==1) cat("\nLambda", lambdaIndex, " of ", lambdaNum, '\n', sep='')
        cat("outer iteration ", iterOut, " ", '\n', sep='')
      }
      cat("betaHat: ", signif(betaHat, 2), '\n\n')
    }
    
  }  # end while (!conv && iterOut<maxiterOut)
  
  # Opting not to return penalty or obj because they depend on covariate scaling.
  list(betaHat=betaHat, loglik=loglik, iterOut=iterOut, iterIn=iterIn, dif=dif)
}







# Each of these functions creates a list of functions consisting of
# g (link), h (inverse link), and getQ (Jacobian of inverse link, dh/ddeta^T).
# delta is a transformation of p to which the element-wise link function is applied (e.g. logit)
# tt(p) = delta, ttinv(delta) = p, and ttinvprime(delta) = dttinv/ddelta^T

# Wrapper for link function construction functions
makeLinkfun <- function(family, link)
{
  if (family == "cumulative") {
    linkfun <- makeLinkCumulative(link)
  } else if (family == "sratio") {
    linkfun <- makeLinkRatio(link, stopping=TRUE)
  } else if (family == "cratio") {
    linkfun <- makeLinkRatio(link, stopping=FALSE)
  } else if (family == "acat") {
    linkfun <- makeLinkACAT(link)
  }
  
  linkfun
}

# Cumulative probability family
makeLinkCumulative <- function(link)
{
  lf <- stats::make.link(link)
  
  tt <- function(p) cumsum(p)
  
  ttinv <- function(delta) delta - c(0, delta[-length(delta)])
  
  ttinvprime <- function(delta)
  {
    k <- length(delta)
    ttipDiag <- diag(1, nrow=k)
    ttipOffDiag <- cbind(-ttipDiag[, -1], 0)
    ttip <- ttipDiag + ttipOffDiag
    ttip
  }
  
  g <- function(p) lf$linkfun(tt(p))
  h <- function(eta) ttinv(lf$linkinv(eta))
  getQ <- function(eta) rep(lf$mu.eta(eta), each=length(eta)) * ttinvprime(lf$linkinv(eta))
  list(g=g, h=h, getQ=getQ)
}

# (Forward) stopping and continuation ratio families
# stopping is logical argument, if FALSE, then make continuation ratio link
makeLinkRatio <- function(link, stopping)
{
  lf <- stats::make.link(link)
  
  tt <- function(p)
  {
    k <- length(p)
    cp <- c(1, 1-cumsum(p[-k]))  # 1-cumulative probabilites
    delta <- p / cp
    delta
  }
  
  ttinv <- function(delta)
  {
    k <- length(delta)
    p <- rep(NA, k)
    p[1] <- delta[1]
    cp <- 1 - p[1]  # 1-cumulative probability dummy variable
    if (k >= 2) {
      for (i in 2:k) {
        p[i] <- delta[i] * cp
        cp <- cp - p[i]
      }
    }
    p
  }
  
  ttinvprime <- function(delta)
  {
    k <- length(delta)
    p <- ttinv(delta)
    cp <- c(1, 1-cumsum(p))  # 1-cumulative probabilities
    ttip <- matrix(nrow=k, ncol=k)
    ttip[1, ] <- c(1, rep(0, k-1))
    cs <- ttip[1, ]  # cumulative row sum dummy variable
    if (k >= 2) {
      for (i in 2:k) {
        ttip[i, ] <- -delta[i] * cs
        ttip[i, i] <- cp[i]
        cs <- cs + ttip[i, ]
      }
    }
    ttip
  }
  
  if (stopping)
  {
    g <- function(p) lf$linkfun(tt(p))
    h <- function(eta) ttinv(lf$linkinv(eta))
    getQ <- function(eta) rep(lf$mu.eta(eta), each=length(eta)) * ttinvprime(lf$linkinv(eta))
  } else
  {
    g <- function(p) lf$linkfun(1 - tt(p))
    h <- function(eta) ttinv(1 - lf$linkinv(eta))
    getQ <- function(eta) rep(lf$mu.eta(eta), each=length(eta)) * -ttinvprime(1-lf$linkinv(eta))
  }
  list(g=g, h=h, getQ=getQ)
}

## Adjacent category
makeLinkACAT <- function(link)
{
  lf <- stats::make.link(link)
  
  tt <- function(p)
  {
    pp <- c(p[-1], 1-sum(p))
    delta <- pp / (p + pp)
    delta
  }
  
  ttinv <- function(delta)
  {
    k <- length(delta)
    ar <- delta / (1-delta)  # adjacent ratios p2/p1, p3/p2, ... , pk+1/pk
    r1 <- cumprod(ar)  # p relative to p1:  p2/p1, p3/p1, ..., pk+1/p1
    p1 <- 1 / (1 + sum(r1))
    p <- c(p1, p1*r1[-k])
    p
  }
  
  ttinvprime <- function(delta)
  {
    k <- length(delta)
    p <- ttinv(delta)
    p1 <- p[1]
    cp <- cumsum(p)  # cumulative probabilities
    ar <- delta / (1-delta)  # adjacent ratios p2/p1, p3/p2, ... , pk+1/pk
    dttdar <- matrix(nrow=k, ncol=k)  # jacobian of ttinv with respect to ar
    dttdar[1, ] <- -p1 * (1-cp) / ar
    if (k >= 2) {
      for (i in 2:k) {
        dttdar[i, ] <- ar[i-1] * dttdar[i-1, ]
        dttdar[i, i-1] <- dttdar[i, i-1] + p[i-1]
      }
    }
    darddelta <- 1 / (1-delta)^2  # jacobian of ar with respect to delta (diagonal matrix)
    ttip <- rep(darddelta, each=k) * dttdar  # multiply each row of dttdar by darddelta
    ttip
  }
  
  g <- function(p) lf$linkfun(tt(p))
  h <- function(eta) ttinv(lf$linkinv(eta))
  getQ <- function(eta) rep(lf$mu.eta(eta), each=length(eta)) * ttinvprime(lf$linkinv(eta))
  list(g=g, h=h, getQ=getQ)
}







## General optimization algorithm for multinomial regression models via coordinate descent
mirlsNet <- function(xList, yMat, alpha, penaltyFactors, positiveID, linkfun, betaStart,
                     lambdaVals, nLambda, lambdaMinRatio, includeLambda0, alphaMin,
                     pMin, stopThresh, threshOut, threshIn, maxiterOut, maxiterIn,
                     printIter, printBeta)
{
  wtsum <- if (is.null(attr(yMat, "wtsum"))) sum(yMat) else attr(yMat, "wtsum")
  nObs <- nrow(yMat)
  xMat <- do.call(rbind, xList)
  if (!is.null(lambdaVals)) lambdaVals <- sort(lambdaVals, decreasing=TRUE)
  lambdaNum <- if (is.null(lambdaVals)) nLambda + includeLambda0 else length(lambdaVals)
  fits <- vector("list", length=lambdaNum)
  
  # If lambdaVals=NULL, need to determine the minimum lambda value that sets all penalized coefficients to zero
  if (is.null(lambdaVals))
  {
    lambdaMod <- ifelse(penaltyFactors==0, 0, Inf)  # to find solution with only unpenalized terms
    fits[[1]] <- cdOut(betaHat=betaStart, lambdaIndex=1, lambdaNum, lambdaMod,
                       xList, xMat, yMat, max(alpha, alphaMin), positiveID, linkfun,
                       pMin, threshOut, threshIn, maxiterOut, maxiterIn,
                       printIter, printBeta)
    betaStart <- fits[[1]]$betaHat
    
    # Calculate starting lambda value
    etaMat <- matrix(xMat %*% betaStart, nrow=nObs, byrow=TRUE)
    pMat <- do.call(rbind, lapply(1:nrow(etaMat), function(i) linkfun$h(etaMat[i, ])))
    si <- getScoreInfo(xList, yMat, pMat, pMin, linkfun)
    # betaHat is zero for all penalized terms, so the soft threshold argument is just the score function
    penID <- penaltyFactors != 0
    lambdaMaxVals <- si$score[penID] / (wtsum * max(alpha, alphaMin) * penaltyFactors[penID])
    lambdaMaxVals[positiveID[penID]] <- pmax(0, lambdaMaxVals[penID & positiveID])
    lambdaMaxVals <- abs(lambdaMaxVals)
    lambdaMax <- max(lambdaMaxVals)
    lambdaMin <- lambdaMax * lambdaMinRatio
    lambdaVals <- exp(seq(log(lambdaMax), log(lambdaMin), length.out=nLambda))
    if (includeLambda0) lambdaVals <- c(lambdaVals, 0)
  }
  
  # If alpha < alphaMin, the model needs to be re-fit the first lambda value
  # using alpha instead of alphaMin
  if (alpha < alphaMin)
    fits[1] <- list(NULL)
  
  # fits[[1]] is NULL if alpha < alphaMin or if lambdaVals is specified by user.
  llik <- if (is.null(fits[[1]])) -Inf else fits[[1]]$loglik
  for (i in (1+!is.null(fits[[1]])):lambdaNum)
  {
    # If relative change in loglik is < stopThresh, then use the current fit
    # for all remaining lambda. Do not stop if loglik stays the same, because
    # this can happen if the first several lambda values produce null models,
    # e.g. in cross validation.
    if ((i > 2) && llikOld != llik && (abs((llikOld - llik) / llikOld) < stopThresh))
    {
      fits[[i]] <- fits[[i-1]]
    } else
    {
      lambdaMod <- lambdaVals[i] * penaltyFactors
      lambdaMod <- ifelse(penaltyFactors==0, 0, lambdaVals[i] * penaltyFactors)
      fits[[i]] <- cdOut(betaHat=betaStart, lambdaIndex=i, lambdaNum, lambdaMod,
                         xList, xMat, yMat, alpha, positiveID, linkfun,
                         pMin, threshOut, threshIn, maxiterOut, maxiterIn,
                         printIter, printBeta)
      betaStart <- fits[[i]]$betaHat
      llikOld <- llik
      llik <- fits[[i]]$loglik
    }
  }
  
  iterOut <- sapply(fits, function(x) x$iterOut)
  iterIn <- sapply(fits, function(x) x$iterIn)
  dif <- sapply(fits, function(x) x$dif)
  
  if (any(iterOut==maxiterOut))
    warning(paste0("Reached outer loop iteration limit before convergence ",
                   "for at least one lambda value. Consider increasing maxiterOut."))
  # if (any(iterIn==maxiterIn & dif>0))
  #     warning(paste0("Outer loop converged upon inner loop reaching iteration ",
  #                    "limit for at least one lambda value. Consider increasing maxiterIn."))
  
  betaHat <- t(sapply(fits, function(f) f$betaHat))
  loglik <- sapply(fits, function(f) f$loglik)
  list(lambdaVals=lambdaVals, betaHat=betaHat, loglik=loglik, iterOut=iterOut, iterIn=iterIn, dif=dif)
}






# returns wtsum * lambda * (alpha*|betaHat|_1 + (1-alpha)/2*|betaHat|_2^2)
getPenalty <- function(betaHat, lambdaMod, alpha)
{
  lambdaMod[betaHat == 0] <- 0  # 0*Inf=0
  l1 <- sum(abs(betaHat) * lambdaMod)
  l22 <- sum(betaHat^2 * lambdaMod)
  pen1 <- alpha * l1
  pen2 <- .5 * (1-alpha) * l22
  pen <- pen1 + pen2
  pen
}

getLoglik <- function(pMat, yMat)
{
  pkplusone <- 1 - rowSums(pMat)
  pMatFull <- cbind(pMat, pkplusone)
  
  if (any(pMatFull < 0)) return(-Inf)

  llMat <- yMat * log(pMatFull)
  llMat[yMat==0] <- 0  # -Inf*0 = 0
  llik <- sum(llMat)
  llik
}

getMisclass <- function(pMat, yMat)
{
  pkplusone <- 1 - rowSums(pMat)
  pMatFull <- cbind(pMat, pkplusone)
  predClass <- apply(pMatFull, 1, which.max)
  nMisclass <- sapply(1:nrow(yMat), function(i) sum(yMat[i, -predClass[i]]))
  misclass <- sum(nMisclass) / sum(yMat)
  misclass
}

getBrier <- function(pMat, yMat)
{
  pkplusone <- 1 - rowSums(pMat)
  pMatFull <- cbind(pMat, pkplusone)
  n <- rowSums(yMat)
  brier <- sum(yMat * (1 - pMatFull)^2 + (n - yMat) * pMatFull^2) / sum(n)
  brier
}

# Returns approximate log-likelihood (as a function of beta, up to a constant)
getLoglikApprox <- function(betaHatActive, scoreActive, infoActive)
{
  -sum(betaHatActive * (infoActive %*% betaHatActive)) + sum(scoreActive * betaHatActive)
}

getLoglikNull <- function(yMat)
{
  pHatNull <- colSums(yMat) / sum(yMat)
  llmat0 <- yMat * rep(log(pHatNull), each=nrow(yMat))
  llmat0[yMat == 0] <- 0  # 0*-Inf=0
  loglik0 <- sum(llmat0)
  loglik0
}

invLogit <- function(x) 1/(1+exp(-x))

softThresh <- function(z, g) sign(z)*max(0, abs(z)-g)

yFactorToMatrix <- function(y)
{
  nObs <- length(y)
  nLev <- length(levels(y))
  yMat <- matrix(0, nrow=nObs, ncol=nLev, dimnames=list(NULL, levels(y)))
  yInt <- as.integer(y)
  yMat[cbind(1:nObs, yInt)] <- 1
  yMat
}


# Enhanced subroutine for makexMat function,
# to allow enforcing a particular non-parallel odds fit 
# (with specific cumulative odds categories chose)
makeNonparallelBlock.mine <- function(x, nLev, nonparallel.categ)
{
  nVar <- length(x)
  zeros <- rep(0, nVar)
  block <- do.call(rbind, lapply(c(1:(nLev-1)), function(i)
  {
    c(rep(zeros, i-1), x, rep(zeros, nLev-1-i))
  }))
  
  needed.col.ind <- sapply(c(1:(nLev-1))[nonparallel.categ], function(x) (x-1)*nVar + c(1:(nVar)))
  block <- block[,needed.col.ind]
}

# at least one of parallelTerms and nonparellelTerms should be TRUE
makexList.mine <- function(x, nLev, parallelTerms, nonparallelTerms, nonparallel.factor, nonparallel.categ, parallel.factor)
{
  x1 <- diag(nLev-1)  # intercept columns
  xList <- lapply(1:nrow(x), function(i)
  {
    xi <- x[i, ]
    x2 <- if (!parallelTerms) NULL else rbind(xi[parallel.factor])[rep(1, nLev-1), , drop=FALSE]
    x3 <- if (!nonparallelTerms) NULL else  makeNonparallelBlock.mine(xi[nonparallel.factor], nLev, nonparallel.categ)
    xListi <- cbind(x1, x2, x3)
    rownames(xListi) <- NULL
    xListi
  })
  xList
}

getScoreInfo <- function(xList, yMat, pMat, pMin, linkfun)
{
  nObs <- length(xList)
  nCoef <- ncol(xList[[1]])
  nLev <- ncol(yMat)
  
  # Fitted probabilities less than pMin are set to pMin; everything is rescaled to sum to 1
  pMatFull <- cbind(pMat, 1-rowSums(pMat))
  pMatFull[pMatFull < pMin] <- pMin
  pMatFull <- pMatFull / rowSums(pMatFull)
  pkplusone <- pMatFull[, nLev]  # vector of fitted probabilities for class K+1
  pMat <- pMatFull[, -nLev, drop=FALSE]
  
  d <- yMat[, nLev] / pkplusone
  d[yMat[, nLev] == 0] <- 0  # 0/0 = 0
  uMat <- yMat[, -nLev, drop=FALSE] / pMat
  uMat[yMat[, -nLev] == 0] <- 0  # 0/0 = 0
  uminusdMat <- uMat - d
  
  wts <- if (is.null(attr(yMat, "wts"))) rowSums(yMat) else attr(yMat, "wts")
  wpMat <- wts / pMat
  wpkplusone <- wts / pkplusone
  
  score <- rep(0, nCoef)
  info <- matrix(0, nrow=nCoef, ncol=nCoef)
  for (i in 1:nObs)
  {
    # Compute score term
    x <- xList[[i]]
    uminusd <- uminusdMat[i, ]
    eta <- linkfun$g(pMat[i, ])  # compute eta from p, capped at pMin
    q <- linkfun$getQ(eta)
    score <- score + crossprod(x, crossprod(q, uminusd))
    
    # Compute info term
    sigInv <- diag(wpMat[i, ], nrow=nLev-1) + wpkplusone[i]
    w <- crossprod(q, sigInv) %*% q
    info <- info + crossprod(x, w %*% x)
  }
  
  score <- c(score)
  list(score=score, info=info)
}





# Function to get column names for coefficient matrices and predictions.
getDeltaNames <- function(family, reverse, nLev)
{
  index <- if (reverse) nLev:2 else 1:(nLev-1)
  deltaNames <- sapply(index, function(i)
  {
    if (family=="cumulative") {
      if (reverse) {
        paste0("P[Y>=", i, "]")  # P(Y>=i)
      } else {
        paste0("P[Y<=", i, "]")  # P(Y<=i)
      }
    } else if (family=="sratio") {
      if (reverse) {
        paste0("P[Y=", i, "|Y<=", i, "]")  # P(Y=i|Y<=i)
      } else {
        paste0("P[Y=", i, "|Y>=", i, "]")  # P(Y=i|Y>=i)
      }
    } else if (family=="cratio") {
      if (reverse) {
        paste0("P[Y<", i, "|Y<=", i, "]")  # P(Y<i|Y<=i)
      } else {
        paste0("P[Y>", i, "|Y>=", i, "]")  # P(Y>i|Y>=i)
      }
    } else if (family=="acat") {
      if (reverse) {
        paste0("P[Y=", i, "|", i, "<=Y<=", i+1, "]")  # P(Y=i|i<=Y<=i+1)
      } else {
        paste0("P[Y=", i+1, "|", i, "<=Y<=", i+1, "]")  # P(Y=i+1|i<=Y<=i+1)
      }
    }
  })
  
  deltaNames
}

#' Summary method for an "ordinalNet" object.
#'
#' Provides a data frame which summarizes the model fit at each lambda value in
#' the solution path.model fit summary as a data frame.
#'
#' @param object An "ordinalNet" S3 object
#' @param ... Not used. Additional summary arguments.
#'
#' @return A data frame containing a record for each lambda value in the solution
#' path. Each record contains the following fields: lambda value, degrees of freedom
#' (number of nonzero parameters), log-likelihood, AIC, BIC, and percent deviance explained.
#'
#' @seealso
#' \code{\link{ordinalNet}}
#'
#' @examples
#' # See ordinalNet.mine() documentation for examples.
#'
#' @export
summary.ordinalNet.mine <- function(object, ...)
{
  with(object, data.frame(lambdaVals, nNonzero, loglik, devPct, aic, bic))
}

#' Print method for an "ordinalNet" object.
#'
#' Prints the data frame returned by the \code{summary.ordinalNet.mine()} method.
#'
#' @param x An "ordinalNet" S3 object
#' @param ... Not used. Additional plot arguments.
#'
#' @seealso
#' \code{\link{ordinalNet}}
#'
#' @examples
#' # See ordinalNet.mine() documentation for examples.
#'
#' @export
print.ordinalNet.mine <- function(x, ...)
{
  cat("\nSummary of fit:\n\n")
  print(summary.ordinalNet.mine(x))
  cat("\n")
  invisible(x)
}

#' Method to extract fitted coefficients from an "ordinalNet" object.
#'
#' @param object An "ordinalNet" S3 object.
#' @param matrix Logical. If \code{TRUE}, coefficient estimates are returned in
#' matrix form. Otherwise a vector is returned.
#' @param whichLambda Optional index number of the desired \code{lambda} within
#' the sequence of \code{lambdaVals}. By default, the solution with the best AIC
#' is returned.
#' @param criteria Selects the best \code{lambda} value by AIC or BIC. Only used
#' if \code{whichLambda=NULL}.
#' @param ... Not used. Additional coef arguments.
#'
#' @return The object returned depends on \code{matrix}.
#'
#' @seealso
#' \code{\link{ordinalNet}}
#'
#' @examples
#' # See ordinalNet.mine() documentation for examples.
#'
#' @export
coef.ordinalNet.mine <- function(object, matrix=FALSE, whichLambda=NULL, criteria=c("aic", "bic"), ...)
{
  if (!is.null(whichLambda) && length(whichLambda)!=1)
    stop("whichLambda should be a single value, or NULL.")
  criteria <- match.arg(criteria)
  if (is.null(whichLambda)) whichLambda <- which.min(object[[criteria]])
  betaHat <- object$coefs[whichLambda, ]
  
  if (!matrix) return(betaHat)
  # Extract variables from ordinalNet object
  nLev <- object$nLev
  nVar <- object$nVar
  
  parallelTerms <- object$args$parallelTerms
  nonparallelTerms <- object$args$nonparallelTerms
  parallel.factor <- object$args$parallel.factor
  nonparallel.factor <- object$args$nonparallel.factor
  nonparallel.categ <- object$args$nonparallel.categ
  
  
  family <- object$args$family
  reverse <- object$args$reverse
  xNames <- object$xNames
  link <- if (!is.null(object$args$customLink)) "g" else object$args$link
  # Create coefficient matrix
  intercepts <- betaHat[1:(nLev-1)]
  
  nonintercepts <- matrix(0, nrow=nVar, ncol=nLev-1)
  
  # Making sure to properly set up the coefficient matrices accounting for the enhanced functionality
  if (parallelTerms) nonintercepts[parallel.factor, ] <- nonintercepts[parallel.factor, ] + betaHat[nLev:(nLev-1+sum(parallel.factor))]
  if (nonparallelTerms) nonintercepts[nonparallel.factor, nonparallel.categ] <- nonintercepts[nonparallel.factor, nonparallel.categ] + betaHat[-(1:(nLev-1+sum(parallel.factor)*parallelTerms))]
  
  
  betaMat <- rbind(intercepts, nonintercepts)
  
  rownames(betaMat) <- c("(Intercept)", xNames)
  deltaNames <- getDeltaNames(family, reverse, nLev)
  colnames(betaMat) <- paste0(link, "(", deltaNames, ")")
  betaMat
}

#' Predict method for an "ordinalNet" object
#'
#' Obtains predicted probabilities, predicted class, or linear predictors.
#'
#' @param object An "ordinalNet" S3 object.
#' @param newx Optional covariate matrix. If NULL, fitted values will be obtained
#' for the training data, as long as the model was fit with the argument
#' \code{keepTrainingData=TRUE}.
#' @param whichLambda Optional index number of the desired lambda value within
#' the solution path sequence.
#' @param criteria Selects the best lambda value by AIC or BIC. Only used
#' if \code{whichLambda=NULL}.
#' @param type The type of prediction required.  Type "response" returns a
#' matrix of fitted probabilities. Type "class" returns a vector containing the
#' class number with the highest fitted probability. Type "link" returns a
#' matrix of linear predictors.
#' @param ... Not used. Additional predict arguments.
#'
#' @return The object returned depends on \code{type}.
#'
#' @seealso
#' \code{\link{ordinalNet}}
#'
#' @examples
#' # See ordinalNet.mine() documentation for examples.
#'
#' @export
predict.ordinalNet.mine <- function(object, newx=NULL, whichLambda=NULL, criteria=c("aic", "bic"),
                                    type=c("response", "class", "link"), ...)
{
  criteria <- match.arg(criteria)
  type <- match.arg(type)
  if (is.null(newx) && !object$args$keepTrainingData)
    stop(paste0("Model was fit with keepTrainingData=FALSE, so training data was not saved. ",
                "A newx argument is required."))
  # Extract variables from ordinalNet object
  nLev <- object$nLev
  xNames <-  object$xNames
  parallelTerms <- object$args$parallelTerms
  nonparallelTerms <- object$args$nonparallelTerms
  
  nonparallel.factor <- object$args$nonparallel.factor
  parallel.factor <- object$args$parallel.factor
  nonparallel.categ <- object$args$nonparallel.categ
  
  
  
  family <- object$args$family
  link <- object$args$link
  reverse <- object$args$reverse
  linkfun <- if (is.null(object$args$customLink)) makeLinkfun(family, link) else object$args$customLink
  x <- if (is.null(newx)) object$args$x else newx
  
  betaMat <- coef.ordinalNet.mine(object, matrix=TRUE, whichLambda=whichLambda, criteria=criteria)
  # Compute prediction values
  etaMat <- cbind(1, x) %*% betaMat
  
  deltaNames <- getDeltaNames(family, reverse, nLev)
  colnames(etaMat) <- colnames(betaMat)
  if (type == "link") return(etaMat)
  probMat <- do.call(rbind, lapply(1:nrow(etaMat), function(i) linkfun$h(etaMat[i, ])))
  
  probMat <- cbind(probMat, 1-rowSums(probMat))
  if (reverse) probMat <- probMat[, nLev:1]
  colnames(probMat) <- paste0("P[Y=", 1:nLev, "]")
  if (type == "response") return(probMat)
  class <- c(apply(probMat, 1, which.max))
  class
}








#' Ordinal regression models with elastic net penalty
#'
#' Fits ordinal regression models with elastic net penalty by coordinate descent.
#' Supported model families include cumulative probability, stopping ratio, continuation ratio,
#' and adjacent category. These families are a subset of vector glm's which belong to a model
#' class we call the elementwise link multinomial-ordinal (ELMO) class. Each family
#' in this class links a vector of covariates to a vector of class probabilities.
#' Each of these families has a parallel form, which is appropriate for ordinal response
#' data, as well as a nonparallel form that is appropriate for an unordered categorical
#' response, or as a more flexible model for ordinal data. The parallel model
#' has a single set of coefficients, whereas the nonparallel model has a set of coefficients
#' for each response category except the baseline category. It is also possible
#' to fit a model with both parallel and nonparallel terms, which we call the semi-parallel
#' model. The semi-parallel model has the flexibility of the nonparallel model,
#' but the elastic net penalty shrinks it toward the parallel model.
#'
#' @param x Covariate matrix. It is recommended that categorical covariates are
#' converted to a set of indicator variables with a variable for each category
#' (i.e. no baseline category); otherwise the choice of baseline category will
#' affect the model fit.
#' @param y Response variable. Can be a factor, ordered factor, or a matrix
#' where each row is a multinomial vector of counts. A weighted fit can be obtained
#' using the matrix option, since the row sums are essentially observation weights.
#' Non-integer matrix entries are allowed.
#' @param alpha The elastic net mixing parameter, with \code{0 <= alpha <= 1}.
#' \code{alpha=1} corresponds to the lasso penalty, and \code{alpha=0} corresponds
#' to the ridge penalty.
#' @param standardize If \code{standardize=TRUE}, the predictor variables are
#' scaled to have unit variance. Coefficient estimates are returned on the
#' original scale.
#' @param penaltyFactors Optional nonnegative vector of penalty factors with
#' length equal to the number of columns in \code{x}. If this argument is used,
#' then the penalty for each variable is scaled by its corresponding factor.
#' If \code{NULL}, the penalty factor is one for each coefficient.
#' @param positiveID Logical vector indicating whether each coefficient should
#' be constrained to be non-negative. If \code{NULL}, the default value is \code{FALSE}
#' for all coefficients.
#' @param family Specifies the type of model family. Options are "cumulative"
#' for cumulative probability, "sratio" for stopping ratio, "cratio" for continuation ratio,
#' and "acat" for adjacent category.
#' @param reverse Logical. If \code{TRUE}, then the "backward" form of the model
#' is fit, i.e. the model is defined with response categories in reverse order.
#' For example, the reverse cumulative model with \eqn{K+1} response categories
#' applies the link function to the cumulative probabilities \eqn{P(Y \ge 2),
#' \ldots, P(Y \ge K+1)}, rather then \eqn{P(Y \le 1), \ldots, P(Y \le K)}.
#' @param link Specifies the link function. The options supported are logit,
#' probit, complementary log-log, and cauchit. Only used if \code{customLink=NULL}.
#' @param customLink Optional list containing a vectorized link function \code{g},
#' a vectorized inverse link \code{h}, and the Jacobian function of the inverse link
#' \code{getQ}. The Jacobian should be defined as \eqn{\partial h(\eta) / \partial \eta^T}
#' (as opposed to the transpose of this matrix).
#' @param parallelTerms Logical. If \code{TRUE}, then parallel coefficient terms
#' will be included in the model. \code{parallelTerms} and \code{nonparallelTerms}
#' cannot both be \code{FALSE}.
#' @param nonparallelTerms Logical. if \code{TRUE}, then nonparallel coefficient terms
#' will be included in the model. \code{parallelTerms} and \code{nonparallelTerms}
#' cannot both be \code{FALSE}.
#' @param parallelPenaltyFactor Nonnegative numeric value equal to one by
#' default. The penalty on all parallel terms is scaled by this factor (as well
#' as variable-specific \code{penaltyFactors}). Only used if
#' \code{parallelTerms=TRUE}.
#' @param lambdaVals An optional user-specified lambda sequence (vector). If \code{NULL},
#' a sequence will be generated based on \code{nLambda} and \code{lambdaMinRatio}.
#' In this case, the maximum lambda is the smallest value that sets all penalized
#' coefficients to zero, and the minimum lambda is the maximum value multiplied
#' by the factor \code{lambdaMinRatio}.
#' @param nLambda Positive integer. The number of lambda values in the solution path.
#' Only used if \code{lambdaVals=NULL}.
#' @param lambdaMinRatio A factor greater than zero and less than one. Only used
#' if \code{lambdaVals=NULL}.
#' @param includeLambda0 Logical. If \code{TRUE}, then zero is added to the end
#' of the sequence of \code{lambdaVals}. This is not done by default because
#' it can significantly increase computational time. An unpenalized saturated model
#' may have infinite coefficient solutions, in which case the fitting algorithm
#' will still terminate when the relative change in log-likelihood becomes small.
#' Only used if \code{lambdaVals=NULL}.
#' @param alphaMin \code{max(alpha, alphaMin)} is used to calculate the starting
#' lambda value when \code{lambdaVals=NULL}. In this case, the default lambda
#' sequence begins with the smallest lambda value such that all penalized
#' coefficients are set to zero (i.e. the value where the first penalized
#' coefficient enters the solution path). The purpose of this argument is to
#' help select a starting value for the lambda sequence when \code{alpha = 0},
#' because otherwise it would be infinite. Note that \code{alphaMin} is only
#' used to determine the default lamba sequence and that the model is always fit
#' using \code{alpha} to calculate the penalty.
#' @param pMin Value greater than zero, but much less than one. During the optimization
#' routine, the Fisher information is calculated using fitted probabilities. For
#' this calculation, fitted probabilities are capped below by this value to prevent
#' numerical instability.
#' @param stopThresh In the relative log-likelihood change between successive
#' lambda values falls below this threshold, then the last model fit is used for all
#' remaining lambda.
#' @param threshOut Convergence threshold for the coordinate descent outer loop.
#' The optimization routine terminates when the relative change in the
#' penalized log-likelihood between successive iterations falls below this threshold.
#' It is recommended to set \code{theshOut} equal to \code{threshIn}.
#' @param threshIn Convergence threshold for the coordinate descent inner loop. Each
#' iteration consists of a single loop through each coefficient. The inner
#' loop terminates when the relative change in the penalized approximate
#' log-likelihood between successive iterations falls below this threshold.
#' It is recommended to set \code{theshOut} equal to \code{threshIn}.
#' @param maxiterOut Maximum number of outer loop iterations.
#' @param maxiterIn Maximum number of inner loop iterations.
#' @param printIter Logical. If \code{TRUE}, the optimization routine progress is
#' printed to the terminal.
#' @param printBeta Logical. If \code{TRUE}, coefficient estimates are printed
#' after each coordinate descent outer loop iteration.
#' @param warn Logical. If \code{TRUE}, the following warning message is displayed
#' when fitting a cumulative probability model with \code{nonparallelTerms=TRUE}
#' (i.e. nonparallel or semi-parallel model).
#' "Warning message: For out-of-sample data, the cumulative probability model
#' with nonparallelTerms=TRUE may predict cumulative probabilities that are not
#' monotone increasing."
#' The warning is displayed by default, but the user may wish to disable it.
#' @param keepTrainingData Logical. If \code{TRUE}, then \code{x} and \code{y}
#' are saved with the returned "ordinalNet" object. This allows
#' \code{predict.ordinalNet} to return fitted values for the training data
#' without passing a \code{newx} argument.

#' @details
#' The \code{ordinalNet} function fits regression models for a categorical response
#' variable with \eqn{K+1} levels. Conditional on the covariate vector \eqn{x_i}
#' (the \eqn{i^{th}} row of \code{x}), each observation has a vector of \eqn{K+1}
#' class probabilities \eqn{(p_{i1}, \ldots, p_{i(K+1)})}. These probabilities
#' sum to one, and can therefore be parametrized by \eqn{p_i = (p_{i1}, \ldots, p_{iK})}.
#' The probabilities are mapped to a set of \eqn{K} quantities
#' \eqn{\delta_i = (\delta_{i1}, \ldots, \delta_{iK}) \in (0, 1)^K}, which depends on the choice
#' of model \code{family}. The elementwise \code{link} function maps
#' \eqn{\delta_i} to a set of \eqn{K} linear predictors. Together, the \code{family}
#' and \code{link} specifiy a link function between \eqn{p_i} and \eqn{\eta_i}.
#'
#' \strong{\emph{Model families:}}
#'
#' Let \eqn{Y} denote the random response variable for a single observation,
#' conditional on the covariates values of the observation. The random variable
#' \eqn{Y} is discrete with support \{\eqn{1, \ldots, K+1}\}. The following model
#' families are defined according to these mappings between the class
#' probabilities and the values \eqn{\delta_1, \ldots, \delta_K}:
#' \describe{
#'   \item{Cumulative probability}{\eqn{\delta_j = P(Y \le j)}}
#'   \item{Reverse cumulative probability}{\eqn{\delta_j = P(Y \ge j + 1)}}
#'   \item{Stopping ratio}{\eqn{\delta_j = P(Y = j | Y \ge j)}}
#'   \item{Reverse stopping ratio}{\eqn{\delta_j = P(Y=j + 1 | Y \le j + 1)}}
#'   \item{Continuation ratio}{\eqn{\delta_j = P(Y > j | Y \ge j)}}
#'   \item{Reverse continuation ratio}{\eqn{\delta_j = P(Y < j | Y \le j)}}
#'   \item{Adjacent category}{\eqn{\delta_j = P(Y = j + 1 | j \le Y \le j+1)}}
#'   \item{Reverse adjacent category}{\eqn{\delta_j = P(Y = j | j \le Y \le j+1)}}
#' }
#'
#' \strong{\emph{Parallel, nonparallel, and semi-parallel model forms:}}
#'
#' Models within each of these families can take one of three forms, which have
#' different definitions for the linear predictor \eqn{\eta_i}. Suppose each
#' \eqn{x_i} has length \eqn{P}. Let \eqn{b} be a length \eqn{P} vector of
#' regression coefficients. Let \eqn{B} be a \eqn{P \times K} matrix of regression
#' coefficient. Let \eqn{b_0} be a vector of \eqn{K} intercept terms.
#' The three model forms are the following:
#' \describe{
#'   \item{Parallel}{\eqn{\eta_i = b_0 + b^T x_i} (\code{parallelTerms=TRUE}, \code{nonparallelTerms=FALSE})}
#'   \item{Nonparallel}{\eqn{\eta_i = b_0 + B^T x_i} (\code{parallelTerms=FALSE}, \code{nonparallelTerms=TRUE})}
#'   \item{Semi-parallel}{\eqn{\eta_i = b_0 + b^T x_i + B^T x_i} (\code{parallelTerms=TRUE}, \code{nonparallelTerms=TRUE})}
#' }
#' The parallel form has the defining property of ordinal models, which is that
#' a single linear combination \eqn{b^T x_i} shifts the cumulative class probabilities
#' \eqn{P(Y \le j)} in favor of either higher or lower categories. The linear predictors
#' are parallel because they only differ by their intercepts (\eqn{b_0}). The nonparallel form
#' is a more flexible model, and it does not shift the cumulative probabilities together.
#' The semi-parallel model is equivalent to the nonparallel model, but the
#' elastic net penalty shrinks the semi-parallel coefficients toward a common
#' value (i.e. the parallel model), as well as shrinking all coefficients toward zero.
#' The nonparallel model, on the other hand, simply shrinks all coefficients toward zero.
#' When the response categories are ordinal, any of the three model forms could
#' be applied. When the response categories are unordered, only the nonparallel
#' model is appropriate.
#'
#' \strong{\emph{Elastic net penalty:}}
#'
#' The elastic net penalty is defined for each model form as follows. \eqn{\lambda}
#' and \eqn{\alpha} are the usual elastic net tuning parameters, where \eqn{\lambda}
#' determines the degree to which coefficients are shrunk toward zero, and \eqn{\alpha}
#' specifies the amound of weight given to the L1 norm and squared L2 norm penalties.
#' Each covariate is allowed a unique penalty factor \eqn{c_j}, which is specified with the
#' \code{penaltyFactors} argument. By default \eqn{c_j = 1} for all \eqn{j}.
#' The semi-parallel model has a tuning parameter \eqn{\rho} which determines the degree to
#' which the parallel coefficients are penalized. Small values of \eqn{\rho} will
#' result in a fit closer to the parallel model, and large values of \eqn{\rho}
#' will result in a fit closer to the nonparallel model.
#' \describe{
#'   \item{Parallel}{\eqn{\lambda \sum_{j=1}^P c_j \{ \alpha |b_j| +
#'                        \frac{1}{2} (1-\alpha) b_j^2 \}}}
#'   \item{Nonparallel}{\eqn{\lambda \sum_{j=1}^P c_j \{ \sum_{k=1}^K \alpha |B_{jk}| +
#'                           \frac{1}{2} (1-\alpha) B_{jk}^2 \}}}
#'   \item{Semi-parallel}{\eqn{\lambda [ \rho \sum_{j=1}^P c_j \{ \alpha |b_j| +
#'                             \frac{1}{2} (1-\alpha) b_j^2 \} +
#'                             \sum_{j=1}^P c_j \{ \sum_{k=1}^K \alpha |B_{jk}| +
#'                             \frac{1}{2} (1-\alpha) B_{jk}^2 \}]}}
#' }
#'
#' \code{ordinalNet} minimizes the following objective function. Let \eqn{N} be
#' the number of observations, which is defined as the sum of the \code{y} elements
#' when \code{y} is a matrix.
#' \deqn{objective = -1/N*loglik + penalty}
#'
#' @return An object with S3 class "ordinalNet".  Model fit information can be accessed
#' through the \code{coef}, \code{predict}, and \code{summary} methods.
#' \describe{
#'   \item{coefs}{Matrix of coefficient estimates, with each row corresponding to a lambda value.
#'   (If covariates were scaled with \code{standardize=TRUE}, the coefficients are
#'   returned on the original scale).}
#'   \item{lambdaVals}{Sequence of lambda values. If user passed a sequence to the
#'   \code{lambdaVals}, then it is this sequence. If \code{lambdaVals} argument
#'   was \code{NULL}, then it is the sequence generated.}
#'   \item{loglik}{Log-likelihood of each model fit.}
#'   \item{nNonzero}{Number of nonzero coefficients of each model fit, including intercepts.}
#'   \item{aic}{AIC, defined as \code{-2*loglik + 2*nNonzero}.}
#'   \item{bic}{BIC, defined as \code{-2*loglik + log(N)*nNonzero}.}
#'   \item{devPct}{Percentage deviance explained, defined as \eqn{1 - loglik/loglik_0},
#'   where \eqn{loglik_0} is the log-likelihood of the null model.}
#'   \item{iterOut}{Number of coordinate descent outer loop iterations until
#'   convergence for each lambda value.}
#'   \item{iterIn}{Number of coordinate descent inner loop iterations on last outer loop
#'   for each lambda value.}
#'   \item{dif}{Relative improvement in objective function on last outer loop
#'   for each lambda value. Can be used to diagnose convergence issues. If \code{iterOut}
#'   reached \code{maxiterOut} and \code{dif} is large, then \code{maxiterOut} should
#'   be increased. If \code{dif} is negative, this means the objective did not improve
#'   between successive iterations. This usually only occurs when the model is
#'   saturated and/or close to convergence, so a small negative value is not of concern.
#'   (When this happens, the algorithm is terminated for the current lambda value,
#'   and the coefficient estimates from the previous outer loop iteration are returned.)}
#'   \item{nLev}{Number of response categories.}
#'   \item{nVar}{Number of covariates in \code{x}.}
#'   \item{xNames}{Covariate names.}
#'   \item{args}{List of arguments passed to the \code{ordinalNet} function.}
#' }
#'
#' @examples
#' # Simulate x as independent standard normal
#' # Simulate y|x from a parallel cumulative logit (proportional odds) model
#' set.seed(1)
#' n <- 50
#' intercepts <- c(-1, 1)
#' beta <- c(1, 1, 0, 0, 0)
#' ncat <- length(intercepts) + 1  # number of response categories
#' p <- length(beta)  # number of covariates
#' x <- matrix(rnorm(n*p), ncol=p)  # n x p covariate matrix
#' eta <- c(x %*% beta) + matrix(intercepts, nrow=n, ncol=ncat-1, byrow=TRUE)
#' invlogit <- function(x) 1 / (1+exp(-x))
#' cumprob <- t(apply(eta, 1, invlogit))
#' prob <- cbind(cumprob, 1) - cbind(0, cumprob)
#' yint <- apply(prob, 1, function(p) sample(1:ncat, size=1, prob=p))
#' y <- as.factor(yint)
#'
#' # Fit parallel cumulative logit model
#' fit1 <- ordinalNet.mine(x, y, family="cumulative", link="logit",
#'                    parallelTerms=TRUE, nonparallelTerms=FALSE)
#' summary(fit1)
#' coef(fit1)
#' coef(fit1, matrix=TRUE)
#' predict(fit1, type="response")
#' predict(fit1, type="class")
#'
#' # Fit nonparallel cumulative logit model
#' fit2 <- ordinalNet.mine(x, y, family="cumulative", link="logit",
#'                    parallelTerms=FALSE, nonparallelTerms=TRUE)
#' fit2
#' coef(fit2)
#' coef(fit2, matrix=TRUE)
#' predict(fit2, type="response")
#' predict(fit2, type="class")
#'
#' # Fit semi-parallel cumulative logit model (with both parallel and nonparallel terms)
#' fit3 <- ordinalNet.mine(x, y, family="cumulative", link="logit",
#'                    parallelTerms=TRUE, nonparallelTerms=TRUE)
#' fit3
#' coef(fit3)
#' coef(fit3, matrix=TRUE)
#' predict(fit3, type="response")
#' predict(fit3, type="class")
#'
#' @export


## Enhanced version, adding the following arguments:
##    * parallel.factor: which columns of the design matrix to be included in the parallel part of the model,
##    * nonparallel.factor: which columns of the design matrix to be included in the non-parallel part of the model,
##    * nonparallel.category: which category of the response to include in the non-parallel part of the model
##
## 'nonparallel.factor' allowed us to exclude the offensive/defensive margin parameters from the non-parallel part of the model,
## while 'nonparallel.category' allowed to explicitly fit a model with one non-proportional odds category.

ordinalNet.mine <- function(x, y, alpha=1, standardize=TRUE, penaltyFactors=NULL, positiveID=NULL,
                            family=c("cumulative", "sratio", "cratio", "acat"), reverse=FALSE,
                            link=c("logit", "probit", "cloglog", "cauchit"), customLink=NULL,
                            parallelTerms=TRUE, nonparallelTerms=FALSE, 
                            parallel.factor=rep(TRUE, ncol(x)),  
                            nonparallel.factor=rep(TRUE, ncol(x)), 
                            nonparallel.categ=rep(TRUE, nlevels(y)-1),
                            parallelPenaltyFactor=1,
                            lambdaVals=NULL, nLambda=20, lambdaMinRatio=0.01, includeLambda0=FALSE, alphaMin=0.01,
                            pMin=1e-8, stopThresh=1e-8, threshOut=1e-8, threshIn=1e-8, maxiterOut=100, maxiterIn=100,
                            printIter=FALSE, printBeta=FALSE, warn=TRUE, keepTrainingData=TRUE)
{
  
  if (!nonparallelTerms){
    nonparallel.factor <- rep(FALSE, ncol(x))
    nonparallel.categ <- rep(FALSE, nlevels(y)-1)
  } 
  
  
  
  family <- match.arg(family)
  link <- match.arg(link)
  args <- as.list(environment())  # list of arguments to return
  if (!keepTrainingData) args$x <- args$y <- NULL
  
  # Initial argument checks
  if (!is.matrix(x))
    stop("x should be a matrix.")
  if (any(is.na(x)))
    stop("x must not contain missing values.")
  if (any(abs(x) == Inf))
    stop("x must not contain infinite values.")
  if (!is.factor(y) && !is.matrix(y))
    stop("y should be a factor or matrix.")
  if (any(is.na(y)))
    stop("y must not contain missing values.")
  
  # Variable definitions
  yMat <- if (is.matrix(y)) y else yFactorToMatrix(y)
  wts <- attr(yMat, "wts") <- rowSums(yMat)
  wtsum <- attr(yMat, "wtsum") <- sum(wts)
  nVar <- ncol(x)
  nLev <- ncol(yMat)
  if (reverse) yMat <- yMat[, nLev:1]
  
  
  if ((nonparallelTerms) & (length(nonparallel.factor) != ncol(x))){
    stop("Length of 'nonparallel.factor' has to be equal to the number of columns in x.")
  }
  
  if ((nonparallelTerms) & (length(nonparallel.categ) != nLev-1)){
    stop("Length of 'nonparallel.categ' has to be equal to the number of response level minus 1.")
  }
  
  
  # Other argument checks
  if (nrow(x) != nrow(yMat))
    stop("x and y dimensions do not match.")
  if (alpha<0 || alpha>1)
    stop("alpha should be a number such that 0 <= alpha <= 1.")
  if (!is.null(penaltyFactors) && length(penaltyFactors)!=nVar)
    stop(paste0("penaltyFactors should be a numeric vector of length equal to the number ",
                "of variables in x. Set penaltyFactor=NULL to penalize each variable equally."))
  if (!is.null(penaltyFactors) && any(penaltyFactors < 0))
    stop("penaltyFactors values should be nonnegative.")
  if (!is.null(positiveID) && length(positiveID)!=nVar)
    stop(paste0("positiveID should be a logical vector of length equal to the number ",
                "of variables in x. Set positiveID=NULL for no positive restrictions."))
  if (!is.null(positiveID) && any(!is.logical(positiveID)))
    stop("positiveID values should be logical.")
  if (!is.null(lambdaVals) && any(lambdaVals<0))
    stop("lambdaVals values should be nonnegative.")
  if (is.null(lambdaVals) && nLambda < 1)
    stop("nLambda should be >= 1.")
  if (is.null(lambdaVals) && lambdaMinRatio<=0)
    stop("lambdaMinRatio should be strictly greater than zero.")
  if (alpha<alphaMin && alphaMin<=0)
    stop("alphaMin should be strictly greater than zero.")
  if (length(parallelPenaltyFactor) > 1)
    stop("parallelPenaltyFactor should be a single value.")
  if (parallelTerms && parallelPenaltyFactor<0)
    stop("parallelPenaltyFactor should be >= 0.")
  if (!parallelTerms && parallelPenaltyFactor!=1)
    warning("parallelPenaltyFactor is not used when parallelTerms=FALSE.")
  if (!parallelTerms && !nonparallelTerms)
    stop("parallelTerms and nonparallelTerms cannot both be FALSE.")
  if (warn && family=="cumulative" && nonparallelTerms) {
    warning(paste0("For out-of-sample data, the cumulative probability model with ",
                   "nonparallelTerms=TRUE may predict cumulative probabilities that are not ",
                   "monotone increasing."))
  }
  if (!is.null(customLink)) {
    link <- customLink
    message("customLink should be a list containing:\n
                \ \ $linkfun := vectorized link function with domain (0, 1)\n
                \ \ $linkinv := vectorized inverse link with domain (-Inf, Inf)\n
                \ \ $mu.eta  := vectorized Jacobian of linkinv\n
                The customLink argument is not checked, so user should be cautious
                using it.")
  }
  if ((parallelPenaltyFactor != 1) && !parallelTerms)
    warning("parallelPenaltyFactor is not used when parallelTerms = FALSE.")
  
  # Create linkfun
  linkfun <- makeLinkfun(family, link)
  
  # Center x and create xList; also scale x if standardize=TRUE
  xMeans <- colMeans(x)
  if (standardize) {
    xSD <- sqrt(rowSums(wts*(t(x)-xMeans)^2) / wtsum)
    xSD[xSD==0] <- 1
    xStd <- t((t(x) - xMeans) / xSD)
  } else {
    xStd <- t(t(x) - xMeans)
  }
  
  xList <- makexList.mine(xStd, nLev, parallelTerms, nonparallelTerms, nonparallel.factor, nonparallel.categ, 
                          parallel.factor=parallel.factor)  # ADDED RECENTLY
  
  # Augment penaltyFactors to include all model coefficients
  if (is.null(penaltyFactors)) penaltyFactors <- rep(1, nVar)
  penaltyFactorsIntercept <- rep(0, nLev-1)
  penaltyFactorsParallel <- if (parallelTerms) penaltyFactors[parallel.factor] * parallelPenaltyFactor else NULL
  penaltyFactorsNonparallel <- if(nonparallelTerms) rep(1, sum(nonparallel.factor)*sum(nonparallel.categ)) else NULL
  penaltyFactors <- c(penaltyFactorsIntercept, penaltyFactorsParallel, penaltyFactorsNonparallel)
  
  
  # Augment positiveID to include all model coefficients
  if (is.null(positiveID)) positiveID <- rep(FALSE, nVar)
  positiveID <- c(rep(FALSE, nLev-1), 
                  if (!parallelTerms) NULL else positiveID[parallel.factor], 
                  if (!nonparallelTerms) NULL else rep(positiveID[nonparallel.factor], sum(nonparallel.categ)))
  
  
  # Initialize coefficient values to intercept-only model
  yFreq <- colSums(yMat) / wtsum
  interceptStart <- linkfun$g(yFreq[-nLev])
  interceptStart <- pmin(100, pmax(-100, interceptStart))
  noninterceptStart <- rep(0, sum(parallel.factor)*parallelTerms + nonparallelTerms*(sum(nonparallel.factor)*sum(nonparallel.categ)))
  
  betaStart <- c(interceptStart, noninterceptStart)
  
  mirlsNetFit <- mirlsNet(xList, yMat, alpha, penaltyFactors, positiveID, linkfun, betaStart,
                          lambdaVals, nLambda, lambdaMinRatio, includeLambda0, alphaMin,
                          pMin, stopThresh, threshOut, threshIn, maxiterOut, maxiterIn,
                          printIter, printBeta)
  betaHat <- mirlsNetFit$betaHat
  lambdaVals <- mirlsNetFit$lambdaVals
  loglik <- mirlsNetFit$loglik
  iterOut <- mirlsNetFit$iterOut
  iterIn <- mirlsNetFit$iterIn
  dif <- mirlsNetFit$dif
  
  # Change coefficient estimates back to original scale if standardize=TRUE
  intercepts0 <- betaHat[ , 1:(nLev-1), drop=FALSE]
  nonintercepts0 <- betaHat[ , -(1:(nLev-1)), drop=FALSE]
  unscaleFact <- if (standardize) xMeans / xSD else xMeans
  
  intAdjust <- matrix(0, nrow=nrow(betaHat), ncol=nLev-1)
  if (parallelTerms) intAdjust <- intAdjust +
    (nonintercepts0[ , 1:sum(parallel.factor), drop=FALSE] %*% unscaleFact[parallel.factor])[ , rep(1, nLev-1), drop=FALSE]
  
  
  if (nonparallelTerms){
    intAdjust[, nonparallel.categ] <- intAdjust[, nonparallel.categ] + sapply(1:sum(nonparallel.categ), function(i) {
      nonintercepts0[ , sum(parallel.factor)*parallelTerms + ((sum(nonparallel.factor)*(i-1) + 1):(sum(nonparallel.factor)*i)), drop=FALSE] %*% unscaleFact[nonparallel.factor]
    })
  }
  
  intercepts <- intercepts0 - intAdjust

  nonintercepts <- if (standardize) t(t(nonintercepts0) / c(if (parallelTerms) xSD[parallel.factor] else NULL, if (nonparallelTerms) rep(xSD[nonparallel.factor], sum(nonparallel.categ)) else NULL)) else nonintercepts0
  coefs <- cbind(intercepts, nonintercepts)
  
  # Create coefficient column names
  catOrder <- if (reverse) nLev:2 else 1:(nLev-1)
  interceptNames <- paste0("(Intercept):", catOrder)
  xNames <- if (is.null(colnames(x))) paste0("X", 1:nVar) else colnames(x)
  parallelNames <- nonparallelNames <- NULL
  if (parallelTerms) parallelNames <- xNames[parallel.factor]
  
  if (nonparallelTerms) nonparallelNames <- paste0(rep(xNames[nonparallel.factor], sum(nonparallel.categ)), ":", rep(catOrder[nonparallel.categ], each=sum(nonparallel.factor)))
  colnames(coefs) <- c(interceptNames, parallelNames, nonparallelNames)
  
  # Calculate approximate AIC, BIC, and %deviance
  nNonzero <- apply(coefs, 1, function(b) sum(b!=0))
  aic <- -2 * loglik + 2 * nNonzero
  bic <- -2 * loglik + log(wtsum) * nNonzero
  loglikNull <- getLoglikNull(yMat)
  devPct <- 1 - loglik / loglikNull
  
  fit <- list(coefs=coefs, lambdaVals=lambdaVals, loglik=loglik,
              nNonzero=nNonzero, aic=aic, bic=bic, devPct=devPct,
              iterOut=iterOut, iterIn=iterIn, dif=dif,
              nLev=nLev, nVar=nVar, xNames=xNames, args=args)
  class(fit) <- "ordinalNet"
  fit
}








#' Uses K-fold cross validation to obtain out-of-sample log-likelihood and
#' misclassification rates. Lambda is tuned within each cross validation fold.
#'
#' The data is divided into K folds. \code{ordinalNet} is fit \eqn{K} times, each time
#' leaving out one fold as a test set. For each of the \eqn{K} model fits, lambda
#' can be tuned by AIC or BIC, or cross validation. If cross validation is used,
#' the user can choose whether to user the best average out-of-sample log-likelihood,
#' misclassification rate, Brier score, or percentage of deviance explained.
#' The user can also choose the number of cross validation folds to use for tuning.
#' Once the model is tuned, the out of sample log-likelihood,
#' misclassification rate, Brier score, and percentage of deviance explained
#' are calculated on the held out test set.
#'
#' @param x Covariate matrix.
#' @param y Response variable. Can be a factor, ordered factor, or a matrix
#' where each row is a multinomial vector of counts. A weighted fit can be obtained
#' using the matrix option, since the row sums are essentially observation weights.
#' Non-integer matrix entries are allowed.
#' @param lambdaVals An optional user-specified lambda sequence (vector). If \code{NULL},
#' a sequence will be generated using the model fit to the full training data.
#' This default sequence is based on \code{nLambda} and \code{lambdaMinRatio},
#' which can be passed as additional arguments (otherwise \code{ordinalNet} default
#' values are used). The maximum lambda is the smallest value that sets all penalized
#' coefficients to zero, and the minimum lambda is the maximum value multiplied
#' by the factor \code{lambdaMinRatio}.
#' @param folds An optional list, where each element is a vector of row indices
#' corresponding to a different cross validation fold. Indices correspond to rows
#' of the \code{x} matrix. Each index number should be used in exactly one fold.
#' If \code{NULL}, the data will be randomly divided into equally-sized partitions.
#' It is recommended to call \code{set.seed} before calling \code{ordinalNetCV}
#' for reproducibility.
#' @param nFolds Numer of cross validation folds. Only used if \code{folds=NULL}.
#' @param nFoldsCV Number of cross validation folds used to tune lambda for each
#' training set (i.e. within each training fold). Only used of \code{tuneMethod} is
#' "cvLoglik", "cvMisclass", "cvBrier", or "cvDevPct.
#' @param tuneMethod Method used to tune lambda for each training set (ie. within
#' each training fold). The "cvLoglik", "cvMisclass", "cvBrier", and "cvDevPct"
#' methods use cross validation with \code{nFoldsCV} folds and select the
#' lambda value with the best average out-of-sample performance. The "aic" and "bic"
#' methods are less computationally intensive because they do not require the
#' model to be fit multiple times.
#' Note that for the methods that require cross validation, the fold splits are
#' determined randomly and cannot be specified by the user. The \code{set.seed()}
#' function should be called prior to \code{ordinalNetCV} for reproducibility.
#' @param printProgress Logical. If \code{TRUE} the fitting progress is printed
#' to the terminal.
#' @param warn Logical. If \code{TRUE}, the following warning message is displayed
#' when fitting a cumulative probability model with \code{nonparallelTerms=TRUE}
#' (i.e. nonparallel or semi-parallel model).
#' "Warning message: For out-of-sample data, the cumulative probability model
#' with nonparallelTerms=TRUE may predict cumulative probabilities that are not
#' monotone increasing."
#' The warning is displayed by default, but the user may wish to disable it.
#' @param ... Other arguments (besides \code{x}, \code{y}, \code{lambdaVals}, and \code{warn})
#' passed to \code{ordinalNet}.
#'
#' @details
#' \itemize{
#'   \item The fold partition splits can be passed by the user via the \code{folds}
#'   argument. By default, the data are randomly divided into equally-sized partitions.
#'   Note that if lambda is tuned by cross validation, the fold splits are
#'   determined randomly and cannot be specified by the user. The \code{set.seed}
#'   function should be called prior to \code{ordinalNetCV} for reproducibility.
#'   \item A sequence of lambda values can be passed by the user via the
#'   \code{lambdaVals} argument. By default, the sequence is generated by first
#'   fitting the model to the full data set (this sequence is determined by the
#'   \code{nLambda} and \code{lambdaMinRatio} arguments of \code{ordinalNet}).
#'   \item The \code{standardize} argument of \code{ordinalNet} can be modified through
#'   the additional arguments (...). If \code{standardize=TRUE}, then the data are scaled
#'   within each cross validation fold. If \code{standardize=TRUE} and lambda is tuned by
#'   cross validation, then the data are also scaled within each tuning sub-fold.
#'   This is done because scaling is part of the statistical procedure and should
#'   be repeated each time the procedure is applied.
#' }
#'
#' @return
#' An S3 object of class "ordinalNetCV", which contains the following:
#' \describe{
#'   \item{loglik}{Vector of out-of-sample log-likelihood values.
#'   Each value corresponds to a different fold.}
#'   \item{misclass}{Vector of out-of-sample misclassificaton rates.
#'   Each value corresponds to a different fold.}
#'   \item{brier}{Vector of out-of-sample Brier scores.
#'   Each value corresponds to a different fold.}
#'   \item{devPct}{Vector of out-of-sample percentages of deviance explained.
#'   Each value corresponds to a different fold.}
#'   \item{bestLambdaIndex}{The index of the value within the lambda sequence
#'   selected for each fold by the tuning method.}
#'   \item{lambdaVals}{The sequence of lambda values used for all cross validation folds.}
#'   \item{folds}{A list containing the index numbers of each fold.}
#'   \item{fit}{An object of class "ordinalNet", resulting from fitting
#'   \code{ordinalNet} to the entire dataset.}
#' }
#'
#' @examples
#' \dontrun{
#' # Simulate x as independent standard normal
#' # Simulate y|x from a parallel cumulative logit (proportional odds) model
#' set.seed(1)
#' n <- 50
#' intercepts <- c(-1, 1)
#' beta <- c(1, 1, 0, 0, 0)
#' ncat <- length(intercepts) + 1  # number of response categories
#' p <- length(beta)  # number of covariates
#' x <- matrix(rnorm(n*p), ncol=p)  # n x p covariate matrix
#' eta <- c(x %*% beta) + matrix(intercepts, nrow=n, ncol=ncat-1, byrow=TRUE)
#' invlogit <- function(x) 1 / (1+exp(-x))
#' cumprob <- t(apply(eta, 1, invlogit))
#' prob <- cbind(cumprob, 1) - cbind(0, cumprob)
#' yint <- apply(prob, 1, function(p) sample(1:ncat, size=1, prob=p))
#' y <- as.factor(yint)
#'
#' # Evaluate out-of-sample performance of the  cumulative logit model
#' # when lambda is tuned by cross validation (best average out-of-sample log-likelihood)
#' cv <- ordinalNetCV(x, y, tuneMethod="cvLoglik")
#' summary(cv)
#' }
#'
#' @export
ordinalNetCV.mine <- function(x, y, lambdaVals=NULL, folds=NULL, nFolds=5, nFoldsCV=5,
                              tuneMethod=c("cvLoglik", "cvMisclass", "cvBrier", "cvDevPct",
                                           "aic", "bic"),
                              printProgress=TRUE, warn=TRUE, ...)
{
  tuneMethod <- match.arg(tuneMethod)
  # cvID := indicator to use cross validation within folds
  cvID <- tuneMethod %in% c("cvLoglik", "cvMisclass", "cvBrier", "cvDevPct")
  # tuneMethod := argument passed to ordinalNetTune.mine
  if (tuneMethod == "cvLoglik") cvCriterion <- "loglik"
  if (tuneMethod == "cvMisclass") cvCriterion <- "misclass"
  if (tuneMethod == "cvBrier") cvCriterion <- "brier"
  if (tuneMethod == "cvDevPct") cvCriterion <- "devPct"
  
  # Argument checks
  if (is.matrix(y) && any(rowSums(y)!=1))
    warning(paste0("Data is split by row for cross validation, but note that ",
                   "y matrix rows have different weights. Be sure this is what you want."))
  if (!is.null(folds) && length(folds)<2)
    stop(paste0("\'folds\' should be a list of at least two vectors. ",
                "Each vector should contain indices of a cross validation fold. ",
                "Each index from 1:nrow(x) should be used exactly once."))
  if (!is.null(folds) && !setequal(unlist(folds), 1:nrow(x)))
    stop("\'folds\' should include each index from 1:nrow(x) exactly once.")
  
  yMat <- if (is.matrix(y)) y else yFactorToMatrix(y)  # for computing log-likelihood
  if (printProgress) cat("Fitting ordinalNet on full training data\n")
  fit <- ordinalNet.mine(x, y, lambdaVals=lambdaVals, warn=warn, ...)
  if (is.null(lambdaVals)) lambdaVals <- fit$lambdaVals
  
  if (is.null(folds))
  {
    n <- nrow(x)
    randIndex <- sample(n)
    folds <- split(randIndex, rep(1:nFolds, length.out=n))
  } else
  {
    nFolds <- length(folds)
  }
  
  nLambda <- length(lambdaVals)
  loglik <- misclass <- brier <- devPct <- bestLambdaIndex <- rep(NA, nFolds)
  names(loglik) <- names(misclass) <- names(brier) <- names(devPct) <-
    names(bestLambdaIndex) <- paste0("fold", 1:nFolds)
  for (i in 1:nFolds)
  {
    testFold <- folds[[i]]
    xTrain <- x[-testFold, , drop=FALSE]
    xTest <- x[testFold, , drop=FALSE]
    yTrain <- if (is.matrix(y)) y[-testFold, , drop=FALSE] else y[-testFold]
    yMatTest <- yMat[testFold, , drop=FALSE]
    if (printProgress) cat("Fitting ordinalNet on fold", i, "of", nFolds, '\n')
    
    if (cvID)
    {
      fitTrainCV <- ordinalNetTune.mine(xTrain, yTrain, lambdaVals=lambdaVals, folds=NULL,
                                        nFolds=5, printProgress=FALSE, warn=FALSE, ...)
      fitTrain <- fitTrainCV$fit
      if (cvCriterion %in% c("loglik", "devPct"))
        wm <- which.max
      if (cvCriterion %in% c("misclass", "brier"))
        wm <- which.min
      bestLambdaIndex[[i]] <- wm(rowMeans(fitTrainCV[[cvCriterion]]))
    } else  # tuneMethod is either "aic" or "bic"
    {
      fitTrain <- ordinalNet.mine(xTrain, yTrain, lambdaVals=lambdaVals, warn=FALSE, ...)
      bestLambdaIndex[[i]] <- which.min(fitTrain[[tuneMethod]])
    }
    
    pHatFull <- predict.ordinalNet.mine(fitTrain, newx=xTest, type="response", whichLambda=bestLambdaIndex[[i]])
    pHat <- pHatFull[, -ncol(pHatFull), drop=FALSE]
    loglik[i] <- getLoglik(pHat, yMatTest)
    misclass[i] <- getMisclass(pHat, yMatTest)
    brier[i] <- getBrier(pHat, yMatTest)
    loglikNull <- getLoglikNull(yMatTest)
    devPct[i] <- 1 - loglik[i] / loglikNull
  }
  
  if (printProgress) cat("Done\n")
  
  out <- list(loglik=loglik, misclass=misclass, brier=brier, devPct=devPct,
              bestLambdaIndex=bestLambdaIndex, lambdaVals=lambdaVals, folds=folds, fit=fit)
  class(out) <- "ordinalNetCV"
  out
}

#' Summary method for an "ordinalNetCV" object.
#'
#' Provides a data frame which summarizes the cross validation results, which
#' can be used as an estimate of the out-of-sample performance of a model tuned
#' by a particular method.
#'
#' @param object An "ordinalNetCV" S3 object
#' @param ... Not used. Additional summary arguments.
#'
#' @return A data frame containing a record for each cross validation fold.
#' Each record contains the following: lambda value, log-likelihood,
#' misclassification rate, Brier score, and percentage of deviance explained.
#'
#' @seealso
#' \code{\link{ordinalNetCV}}
#'
#' @examples
#' # See ordinalNetCV() documentation for examples.
#'
#' @export
summary.ordinalNetCV.mine <- function(object, ...)
{
  lambda <- object$lambdaVals[object$bestLambdaIndex]
  loglik <- object$loglik
  misclass <- object$misclass
  brier <- object$brier
  devPct <- object$devPct
  data.frame(lambda=lambda, loglik=loglik, misclass=misclass, brier=brier, devPct=devPct)
}

#' Print method for an "ordinalNetCV" object.
#'
#' Prints the data frame returned by the \code{summary.ordinalNetCV()} method.
#'
#' @param x An "ordinalNetCV" S3 object
#' @param ... Not used. Additional print arguments.
#'
#' @seealso
#' \code{\link{ordinalNetCV}}
#'
#' @examples
#' # See ordinalNetCV() documentation for examples.
#'
#' @export
print.ordinalNetCV.mine <- function(x, ...)
{
  cat("\nCross validation summary:\n\n")
  print(summary.ordinalNetCV(x))
  cat("\n")
  invisible(x)
}





#' Uses K-fold cross validation to obtain out-of-sample log-likelihood and
#' misclassification rates for a sequence of lambda values.
#'
#' The data is divided into K folds. \code{ordinalNet} is fit \eqn{K} times (\code{K=nFolds}),
#' each time leaving out one fold as a test set. The same sequence of lambda values is used
#' each time. The out-of-sample log-likelihood, misclassification rate, Brier score,
#' and percentage of deviance explained are obtained for each lambda value from
#' the held out test set. It is up to the user to determine how to tune the model
#' using this information.
#'
#' @param x Covariate matrix.
#' @param y Response variable. Can be a factor, ordered factor, or a matrix
#' where each row is a multinomial vector of counts. A weighted fit can be obtained
#' using the matrix option, since the row sums are essentially observation weights.
#' Non-integer matrix entries are allowed.
#' @param lambdaVals An optional user-specified lambda sequence (vector). If \code{NULL},
#' a sequence will be generated using the model fit to the full training data.
#' This default sequence is based on \code{nLambda} and \code{lambdaMinRatio},
#' which can be passed as additional arguments (otherwise \code{ordinalNet} default
#' values are used). The maximum lambda is the smallest value that sets all penalized
#' coefficients to zero, and the minimum lambda is the maximum value multiplied
#' by the factor \code{lambdaMinRatio}.
#' @param folds An optional list, where each element is a vector of row indices
#' corresponding to a different cross validation fold. Indices correspond to rows
#' of the \code{x} matrix. Each index number should be used in exactly one fold.
#' If \code{NULL}, the data will be randomly divided into equal-sized partitions.
#' It is recommended to use \code{set.seed} before calling this function to make
#' results reproducible.
#' @param nFolds Numer of cross validation folds. Only used if \code{folds=NULL}.
#' @param printProgress Logical. If \code{TRUE} the fitting progress is printed
#' to the terminal.
#' @param warn Logical. If \code{TRUE}, the following warning message is displayed
#' when fitting a cumulative probability model with \code{nonparallelTerms=TRUE}
#' (i.e. nonparallel or semi-parallel model).
#' "Warning message: For out-of-sample data, the cumulative probability model
#' with nonparallelTerms=TRUE may predict cumulative probabilities that are not
#' monotone increasing."
#' The warning is displayed by default, but the user may wish to disable it.
#' @param ... Other arguments (besides \code{x}, \code{y}, \code{lambdaVals}, and \code{warn})
#' passed to \code{ordinalNet}.
#'
#' @details
#' \itemize{
#'   \item The fold partition splits can be passed by the user via the \code{folds}
#'   argument. By default, the data are randomly divided into equally-sized partitions.
#'   The \code{set.seed} function should be called prior to \code{ordinalNetCV} for reproducibility.
#'   \item A sequence of lambda values can be passed by the user via the
#'   \code{lambdaVals} argument. By default, the sequence is generated by first
#'   fitting the model to the full data set (this sequence is determined by the
#'   \code{nLambda} and \code{lambdaMinRatio} arguments of \code{ordinalNet}).
#'   \item The \code{standardize} argument of \code{ordinalNet} can be modified through
#'   the additional arguments (...). If \code{standardize=TRUE}, then the data are scaled
#'   within each cross validation fold. This is done because scaling is part of
#'   the statistical procedure and should be repeated each time the procedure is applied.
#' }
#'
#' @return
#' An S3 object of class "ordinalNetTune", which contains the following:
#' \describe{
#'   \item{loglik}{Matrix of out-of-sample log-likelihood values.
#'   Each row corresponds to a lambda value, and each column corresponds to a fold.}
#'   \item{misclass}{Matrix of out-of-sample misclassificaton rates.
#'   Each row corresponds to a lambda value, and each column corresponds to a fold.}
#'   \item{brier}{Matrix of out-of-sample Brier scores. Each row corresponds
#'   to a lambda value, and each column corresponds to a fold.}
#'   \item{devPct}{Matrix of out-of-sample percentages of deviance explained.
#'   Each row corresponds to a lambda value, and each column corresponds to a fold.}
#'   \item{lambdaVals}{The sequence of lambda values used for all cross validation folds.}
#'   \item{folds}{A list containing the index numbers of each fold.}
#'   \item{fit}{An object of class "ordinalNet", resulting from fitting
#'   \code{ordinalNet} to the entire dataset.}
#' }
#'
#' @examples
#' \dontrun{
#' # Simulate x as independent standard normal
#' # Simulate y|x from a parallel cumulative logit (proportional odds) model
#' set.seed(1)
#' n <- 50
#' intercepts <- c(-1, 1)
#' beta <- c(1, 1, 0, 0, 0)
#' ncat <- length(intercepts) + 1  # number of response categories
#' p <- length(beta)  # number of covariates
#' x <- matrix(rnorm(n*p), ncol=p)  # n x p covariate matrix
#' eta <- c(x %*% beta) + matrix(intercepts, nrow=n, ncol=ncat-1, byrow=TRUE)
#' invlogit <- function(x) 1 / (1+exp(-x))
#' cumprob <- t(apply(eta, 1, invlogit))
#' prob <- cbind(cumprob, 1) - cbind(0, cumprob)
#' yint <- apply(prob, 1, function(p) sample(1:ncat, size=1, prob=p))
#' y <- as.factor(yint)
#'
#' # Fit parallel cumulative logit model; select lambda by cross validation
#' tunefit <- ordinalNetTune.mine(x, y)
#' summary(tunefit)
#' plot(tunefit)
#' bestLambdaIndex <- which.max(rowMeans(tunefit$loglik))
#' coef(tunefit$fit, whichLambda=bestLambdaIndex, matrix=TRUE)
#' predict(tunefit$fit, whichLambda=bestLambdaIndex)
#' }
#'
#' @export
ordinalNetTune.mine <- function(x, y, lambdaVals=NULL, folds=NULL, nFolds=5, printProgress=TRUE, warn=TRUE, ...)
{
  # Argument checks
  if (is.matrix(y) && any(rowSums(y)!=1))
    warning(paste0("Data is split by row for cross validation, but note that ",
                   "y matrix rows have different weights. Be sure this is what you want."))
  if (!is.null(folds) && length(folds)<2)
    stop(paste0("\'folds\' should be a list of at least two vectors. ",
                "Each vector should contain indices of a cross validation fold. ",
                "Each index from 1:nrow(x) should be used exactly once."))
  if (!is.null(folds) && !setequal(unlist(folds), 1:nrow(x)))
    stop("\'folds\' should include each index from 1:nrow(x) exactly once.")
  
  yMat <- if (is.matrix(y)) y else yFactorToMatrix(y)  # for computing log-likelihood
  if (printProgress) cat("Fitting ordinalNet on full training data\n")
  fit <- ordinalNet.mine(x, y, lambdaVals=lambdaVals, warn=warn, ...)
  if (is.null(lambdaVals)) lambdaVals <- fit$lambdaVals
  
  if (is.null(folds))
  {
    n <- nrow(x)
    randIndex <- sample(n)
    folds <- split(randIndex, rep(1:nFolds, length.out=n))
  } else
  {
    nFolds <- length(folds)
  }
  
  nLambda <- length(lambdaVals)
  loglik <- misclass <- brier <- devPct <- matrix(nrow=nLambda, ncol=nFolds)
  colnames(loglik) <- colnames(misclass) <- colnames(brier) <- colnames(devPct) <-
    paste0("fold", 1:nFolds)
  rownames(loglik) <- rownames(misclass) <- rownames(brier) <- rownames(devPct) <-
    paste0("lambda", 1:nLambda)
  for (i in 1:nFolds)
  {
    testFold <- folds[[i]]
    xTrain <- x[-testFold, , drop=FALSE]
    xTest <- x[testFold, , drop=FALSE]
    yTrain <- if (is.matrix(y)) y[-testFold, , drop=FALSE] else y[-testFold]
    yMatTest <- yMat[testFold, , drop=FALSE]
    if (printProgress) cat("Fitting ordinalNet on fold", i, "of", nFolds, '\n')
    fitTrain <- ordinalNet.mine(xTrain, yTrain, lambdaVals=lambdaVals, warn=FALSE, ...)
    nLambdaTrain <- length(fitTrain$lambdaVals)
    
    for (j in 1:nLambda)
    {
      if (j > nLambdaTrain) {
        loglik[j, i] <- NA
        misclass[j, i] <- NA
      } else {
        pHatFull <- predict.ordinalNet.mine(fitTrain, newx=xTest, type="response", whichLambda=j)
        pHat <- pHatFull[, -ncol(pHatFull), drop=FALSE]
        loglik[j, i] <- getLoglik(pHat, yMatTest)
        misclass[j, i] <- getMisclass(pHat, yMatTest)
        brier[j, i] <- getBrier(pHat, yMatTest)
        loglikNull <- getLoglikNull(yMatTest)
        devPct[j, i] <- 1 - loglik[j, i] / loglikNull
      }
    }
  }
  
  if (printProgress) cat("Done\n")
  
  out <- list(loglik=loglik, misclass=misclass, brier=brier, devPct=devPct,
              lambdaVals=lambdaVals, folds=folds, fit=fit)
  class(out) <- "ordinalNetTune"
  out
}

#' Summary method for an "ordinalNetTune" object.
#'
#' Provides a data frame which summarizes the cross validation results and may
#' be useful for selecting an appropriate value for the tuning parameter lambda.
#'
#' @param object An "ordinalNetTune" S3 object.
#' @param ... Not used. Additional summary arguments.
#'
#' @return A data frame containing a record for each lambda value in the solution
#' path. Each record contains the following: lambda value, average log-likelihood,
#' average misclassification rate, average Brier score, and average percentage
#' of deviance explained. Averages are taken across all cross validation folds.
#'
#' @seealso
#' \code{\link{ordinalNetTune}}
#'
#' @examples
#' # See ordinalNetTune.mine() documentation for examples.
#'
#' @export
summary.ordinalNetTune.mine <- function(object, ...)
{
  lambda <- object$lambdaVals
  loglik_avg <- unname(rowMeans(object$loglik))
  misclass_avg <- unname(rowMeans(object$misclass))
  brier_avg <- unname(rowMeans(object$brier))
  devPct_avg <- unname(rowMeans(object$devPct))
  data.frame(lambda=lambda, loglik_avg=loglik_avg, misclass_avg=misclass_avg,
             brier_avg=brier_avg, devPct_avg=devPct_avg)
}

#' Print method for an "ordinalNetTune" object.
#'
#' Prints the data frame returned by the \code{summary.ordinalNetTune.mine()} method.
#'
#' @param x An "ordinalNetTune" S3 object.
#' @param ... Not used. Additional print arguments.
#'
#' @seealso
#' \code{\link{ordinalNetTune}}
#'
#' @examples
#' # See ordinalNetTune.mine() documentation for examples.
#'
#' @export
print.ordinalNetTune.mine <- function(x, ...)
{
  cat("\nCross validation summary:\n\n")
  print(summary.ordinalNetTune.mine(x))
  cat("\n")
  invisible(x)
}

#' Plot method for "ordinalNetTune" object.
#'
#' Plots the average out-of-sample log-likelihood, misclassification rate,
#' Brier score, or percentage of deviance explained for each lambda value in the
#' solution path. The averae is taken over all cross validation folds.
#'
#' @param x An "ordinalNetTune" S3 object.
#' @param type Which performance measure to plot. Either "loglik", "misclass",
#' "brier", or "devPct".
#' @param ... Additional plot arguments.
#'
#' @seealso
#' \code{\link{ordinalNetTune}}
#'
#' @examples
#' # See ordinalNetTune.mine() documentation for examples.
#'
#'@export
plot.ordinalNetTune.mine <- function(x, type=c("loglik", "misclass", "brier", "devPct"), ...)
{
  type <- match.arg(type)
  y <- rowMeans(x[[type]])
  loglambda <- log(x$lambdaVals)
  if (type == "misclass")
    ylab <- "avg misclass rate"
  if (type == "loglik")
    ylab <- "avg loglik"
  if (type == "brier")
    ylab <- "avg brier score"
  if (type == "devPct")
    ylab <- "avg pct deviance explained"
  graphics::plot(y ~ loglambda, ylab=ylab, xlab="log(lambda)", ...)
}


