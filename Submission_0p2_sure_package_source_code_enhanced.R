#######
## MOST OF THIS CODE IS TAKEN FROM:
##    https://github.com/cran/ordinalNet/tree/master/R
## with modifications to the functions that end with ".mine" postfix.
##
## The most critical modifications are:
##    * creating the ".ordinalNet" versions of various functions to allow for
##      application of the surrogate residual approach to regularized fits
##    * adding the "which.categ" and "binarized=TRUE/FALSE" argument to the residuals-generating functions:
##        - "which.categ" allows user to specify which cumulative odds category 
##          they would like the residuals for.
##        - "binarized" indicator allows to pick whether to conduct full ordinal diagnostics,
##          or just the binarized diagnostics for the chosen cumulative odds category.
##      "which.categ" and "binarized=TRUE/FALSE are most useful for non-proportional odds models, 
##      where the fit differs from category to category.
#######



# various internal functions

################################################################################
# The following functions have been taken from truncdist, but have been modified
# to not throw warnings when vectors are passed to arguments a and b
################################################################################

#' @keywords internal
.rtrunc <- function (n, spec, a = -Inf, b = Inf, ...) {
  .qtrunc(runif(n, min = 0, max = 1), spec, a = a, b = b, ...)
}


#' @keywords internal
.qtrunc <- function (p, spec, a = -Inf, b = Inf, ...) {
  tt <- p
  G <- get(paste("p", spec, sep = ""), mode = "function")
  Gin <- get(paste("q", spec, sep = ""), mode = "function")
  G.a <- G(a, ...)
  G.b <- G(b, ...)
  pmin(pmax(a, Gin(G(a, ...) + p * (G(b, ...) - G(a, ...)), ...)), b)
}


#' @keywords internal
sim_trunc <- function(n, distribution, a, b, location = 0, scale = 1) {
  if (distribution == "norm") {
    .rtrunc(n, spec = distribution, a = a, b = b,
            mean = location, sd = scale)
  } else {
    .rtrunc(n, spec = distribution, a = a, b = b,
            location = location, scale = scale)
  }
}


################################################################################
# Gumbel distribution functions
################################################################################

# For log-log link

#' @keywords internal
pgumbel <- function(q, location = 0, scale = 1) {
  q <- (q - location) / scale
  exp(-exp(-q))
}


#' @keywords internal
qgumbel <- function(p, location = 0, scale = 1) {
  -scale * log(-log(p)) + location
}


#' @keywords internal
rgumbel <- function (n, location = 0, scale = 1) {
  qgumbel(runif(n, min = 0, max = 1), location = location, scale = scale)
}


# For complimentary log-log link

#' @keywords internal
pGumbel <- function(q, location = 0, scale = 1) {
  q <- (q - location) / scale
  1 - exp(-exp(q))
}


#' @keywords internal
qGumbel <- function(p, location = 0, scale = 1) {
  scale * log(-log(1 - p)) + location
}


#' @keywords internal
rGumbel <- function (n, location = 0, scale = 1) {
  qGumbel(runif(n, min = 0, max = 1), location = location, scale = scale)
}


################################################################################
# Generic function to extract truncation bounds for cumulative link models;
# these are used when sampling the surrogate values
################################################################################

#' @keywords internal
getBounds <- function(object, ...) {
  UseMethod("getBounds")
}


#' @keywords internal
getBounds.clm <- function(object, which.categ = 1, ...) {
  unname(
    c(-Inf, stats::coef(object)[seq_len(ncat(object) - 1)] -
        # stats::coef(object)[1L], Inf)
        stats::coef(object)[which.categ], Inf)
  )
}


#' @keywords internal
getBounds.glm <- function(object, ...) {
  y <- getResponseValues(object)
  c(ifelse(y == 0, yes = -Inf, no = 0), ifelse(y == 1, yes = Inf, no = 0))
}


#' @keywords internal
getBounds.lrm <- function(object, ...) {
  coefs <- -unname(stats::coef(object))
  c(-Inf, coefs[seq_len(ncat(object) - 1)] - coefs[1L], Inf)
}


#' @keywords internal
getBounds.orm <- function(object, ...) {
  coefs <- -unname(stats::coef(object))
  c(-Inf, coefs[seq_len(ncat(object) - 1)] - coefs[1L], Inf)
}


#' @keywords internal
getBounds.polr <- function(object, ...) {
  unname(
    c(-Inf, object$zeta - object$zeta[1L], Inf)
  )
}


#' @keywords internal
getBounds.vglm <- function(object, ...) {
  coefs <- if (object@misc$reverse) {
    -unname(stats::coef(object))
  } else {
    unname(stats::coef(object))
  }
  c(-Inf, coefs[seq_len(ncat(object) - 1)] - coefs[1L], Inf)
}


################################################################################
# Generic function for extracting the assumed cumulative distribution function
# from a cumulative link model
################################################################################

#' @keywords internal
getDistributionFunction <- function(object) {
  UseMethod("getDistributionFunction")
}


#' @keywords internal
getDistributionFunction.clm <- function(object) {
  switch(object$link,
         "logit" = plogis,
         "probit" = pnorm,
         "loglog" = pgumbel,
         "cloglog" = pGumbel,
         "cauchit" = pcauchy)
}


#' @keywords internal
getDistributionFunction.glm <- function(object) {
  switch(object$family$link,
         "logit" = plogis,
         "probit" = pnorm,
         # "loglog" = pgumbel,  # glm does not support this link function
         "cloglog" = pGumbel,
         "cauchit" = pcauchy)
}


#' @keywords internal
getDistributionFunction.lrm <- function(object) {
  plogis
}


#' @keywords internal
getDistributionFunction.orm <- function(object) {
  switch(object$family,
         "logistic" = plogis,
         "probit" = pnorm,
         "loglog" = pgumbel,
         "cloglog" = pGumbel,
         "cauchit" = pcauchy)
}


#' @keywords internal
getDistributionFunction.polr <- function(object) {
  switch(object$method,
         "logistic" = plogis,
         "probit" = pnorm,
         "loglog" = pgumbel,
         "cloglog" = pGumbel,
         "cauchit" = pcauchy)
}


#' @keywords internal
getDistributionFunction.vglm <- function(object) {
  switch(object@family@infos()$link,
         "logit" = plogis,
         "probit" = pnorm,
         "loglog" = pgumbel,
         "cloglog" = pGumbel,
         "cauchit" = pcauchy)
}


# ' @keywords internal
getDistributionFunction.ordinalNet <- function(object) {
  switch(object$args$link,
         "logit" = plogis,
         "probit" = pnorm,
         "loglog" = pgumbel,
         "cloglog" = pGumbel,
         "cauchit" = pcauchy)
}





################################################################################
# Generic function for extracting the name of the assumed distribution from a
# cumulative link model
################################################################################

#' @keywords internal
getDistributionName <- function(object) {
  UseMethod("getDistributionName")
}


#' @keywords internal
getDistributionName.clm <- function(object) {
  switch(object$link,
         "logit" = "logis",
         "probit" = "norm",
         "loglog" = "gumbel",
         "cloglog" = "Gumbel",
         "cauchit" = "cauchy")
}


#' @keywords internal
getDistributionName.glm <- function(object) {
  switch(object$family$link,
         "logit" = "logis",
         "probit" = "norm",
         # "loglog" = "gumbel",  # glm does not support this link function
         "cloglog" = "Gumbel",
         "cauchit" = "cauchy")
}


#' @keywords internal
getDistributionName.lrm <- function(object) {
  "logis"
}


#' @keywords internal
getDistributionName.orm <- function(object) {
  switch(object$family,
         "logistic" = "logis",
         "probit" = "norm",
         "loglog" = "gumbel",
         "cloglog" = "Gumbel",
         "cauchit" = "cauchy")
}


#' @keywords internal
getDistributionName.polr <- function(object) {
  switch(object$method,
         "logistic" = "logis",
         "probit" = "norm",
         "loglog" = "gumbel",
         "cloglog" = "Gumbel",
         "cauchit" = "cauchy")
}


#' @keywords internal
getDistributionName.vglm <- function(object) {
  switch(object@family@infos()$link,
         "logit" = "logis",
         "probit" = "norm",
         "loglog" = "gumbel",
         "cloglog" = "Gumbel",
         "cauchit" = "cauchy")
}


################################################################################
# Generic function for extracting the fitted probabilities from a cumulative
# link model
################################################################################

#' @keywords internal
getFittedProbs <- function(object) {
  UseMethod("getFittedProbs")
}


#' @keywords internal
getFittedProbs.clm <- function(object) {
  newdata <- stats::model.frame(object)
  vars <- as.character(attr(object$terms, "variables"))
  resp <- vars[1 + attr(object$terms, "response")]  # response name
  newdata <- newdata[!names(newdata) %in% resp]
  predict(object, newdata = newdata, type = "prob")$fit
}


#' @keywords internal
getFittedProbs.glm <- function(object) {
  prob <- object$fitted.values
  cbind(prob, 1 - prob)
}


#' @keywords internal
getFittedProbs.lrm <- function(object) {
  predict(object, type = "fitted.ind")
}


#' @keywords internal
getFittedProbs.orm <- function(object) {
  predict(object, type = "fitted.ind")
}


#' @keywords internal
getFittedProbs.polr <- function(object) {
  object$fitted.values
}


#' @keywords internal
getFittedProbs.vglm <- function(object) {
  object@fitted.values
}


################################################################################
# Generic function for extracting the fitted mean response from a cumulative
# link model
################################################################################

#' @keywords internal
getMeanResponse <- function(object) {  # for j = 1
  UseMethod("getMeanResponse")
}


#' @keywords internal
getMeanResponse.clm <- function(object) {
  # Have to do this the long way, for now! :(
  mf <- model.frame(object)
  if (!is.null(cl <- attr(object$terms, "dataClasses"))) {
    .checkMFClasses(cl, mf)
  }
  X <- model.matrix(object$terms, data = mf, contrasts = object$contrasts)
  if(sum(object$aliased$beta) > 0) {
    X <- X[, !c(FALSE, object$aliased$beta), drop = FALSE]
  }
  
  # drop(X[, -1L, drop = FALSE] %*% object$beta - object$alpha[1L])
  drop(sapply(object$alpha, function(x) X[, -1L, drop = FALSE] %*% object$beta - x))
}


#' @keywords internal
getMeanResponse.glm <- function(object) {
  object$linear.predictors
}


#' @keywords internal
getMeanResponse.lrm <- function(object) {
  # No negative sign since orm uses the reverse parameterization: Pr(Y >= j)
  predict(object, type = "lp", kint = 1L)
}


#' @keywords internal
getMeanResponse.orm <- function(object) {
  # No negative sign since orm uses the reverse parameterization: Pr(Y >= j)
  predict(object, type = "lp", kint = 1L)
}


#' @keywords internal
getMeanResponse.polr <- function(object) {
  # object$lp
  object$lp - object$zeta[1L]  # Xb - a1
}


#' @keywords internal
getMeanResponse.vglm <- function(object) {
  if (object@misc$reverse) {
    object@predictors[, 1L, drop = TRUE]
  } else {
    -object@predictors[, 1L, drop = TRUE]
  }
}


################################################################################
# Generic function for extracting the assumed quantile function from a
# cumulative link model
################################################################################

#' @keywords internal
getQuantileFunction <- function(object) {
  UseMethod("getQuantileFunction")
}


#' @keywords internal
getQuantileFunction.clm <- function(object) {
  switch(object$link,
         "logit" = qlogis,
         "probit" = qnorm,
         "loglog" = qgumbel,
         "cloglog" = qGumbel,
         "cauchit" = qcauchy)
}


#' @keywords internal
getQuantileFunction.glm <- function(object) {
  switch(object$family$link,
         "logit" = qlogis,
         "probit" = qnorm,
         "log" = qgumbel,
         "cloglog" = qGumbel,
         "cauchit" = qcauchy)
}


#' @keywords internal
getQuantileFunction.lrm <- function(object) {
  qlogis
}


#' @keywords internal
getQuantileFunction.orm <- function(object) {
  switch(object$family,
         "logistic" = qlogis,
         "probit" = qnorm,
         "loglog" = qgumbel,
         "cloglog" = qGumbel,
         "cauchit" = qcauchy)
}


#' @keywords internal
getQuantileFunction.polr <- function(object) {
  switch(object$method,
         "logistic" = qlogis,
         "probit" = qnorm,
         "loglog" = qgumbel,
         "cloglog" = qGumbel,
         "cauchit" = qcauchy)
}


#' @keywords internal
getQuantileFunction.vglm <- function(object) {
  switch(object@family@infos()$link,
         "logit" = qlogis,
         "probit" = qnorm,
         "loglog" = qgumbel,
         "cloglog" = qGumbel,
         "cauchit" = qcauchy)
}


################################################################################
# Generic function to extract the response values from a cumulative link or
# general model; returns an integer, not a factor!
################################################################################

#' @keywords internal
getResponseValues <- function(object, ...) {
  UseMethod("getResponseValues")
}


#' @keywords internal
getResponseValues.clm <- function(object, ...) {
  unname(as.integer(object$y))
}


#' @keywords internal
getResponseValues.glm <- function(object) {
  # FIXME: What about binomial models with matrix response, etc.?
  as.integer(as.factor(model.response(model.frame(object))))
}


#' @keywords internal
getResponseValues.lrm <- function(object) {
  as.integer(model.response(model.frame(object)))
}


#' @keywords internal
getResponseValues.orm <- function(object) {
  as.integer(model.response(model.frame(object)))
}


#' @keywords internal
getResponseValues.polr <- function(object, ...) {
  unname(as.integer(model.response(model.frame(object))))
}


#' @keywords internal
getResponseValues.vglm <- function(object, ...) {
  unname(apply(object@y, MARGIN = 1, FUN = function(x) which(x == 1)))
}


################################################################################
# Number of response categories
################################################################################

#' @keywords internal
ncat <- function(object) {
  UseMethod("ncat")
}


#' @keywords internal
ncat.clm <- function(object) {
  length(object$y.levels)
}


#' @keywords internal
ncat.glm <- function(object) {
  length(unique(getResponseValues(object)))
}


#' @keywords internal
ncat.lrm <- function(object) {
  object$non.slopes + 1
}


#' @keywords internal
ncat.orm <- function(object) {
  object$non.slopes + 1
}


#' @keywords internal
ncat.polr <- function(object) {
  length(object$lev)
}


#' @keywords internal
ncat.vglm <- function(object) {
  length(attributes(object)$extra$colnames.y)
}


################################################################################
# Surrogate and residual workhorse functions
################################################################################

#' @keywords internal
generate_surrogate <- function(object, method = c("latent", "jitter"),
                               jitter.scale = c("response", "probability"),
                               boot_id = NULL,
                               which.categ = 1) {
  
  # Match arguments
  method <- match.arg(method)
  
  # Generate surrogate response values
  s <- if (method == "latent") {  # latent variable approach
    
    # Get distribution name (for sampling)
    distribution <- getDistributionName(object)  # distribution name
    
    # Simulate surrogate response values from the appropriate truncated
    # distribution
    if (distribution %in% c("norm", "logis", "cauchy", "gumbel", "Gumbel")) {
      y <- getResponseValues(object)
      if (is.null(boot_id)) {
        boot_id <- seq_along(y)
      }
      mean_response <- getMeanResponse(object)  # mean response values
      if (!inherits(object, what = "lrm") && inherits(object, what = "glm")) {
        sim_trunc(n = length(y), distribution = distribution,
                  # {0, 1} -> {1, 2}
                  a = ifelse(y[boot_id] == 1, yes = -Inf, no = 0),
                  b = ifelse(y[boot_id] == 2, yes =  Inf, no = 0),
                  location = mean_response[boot_id + (which.categ - 1)*length(boot_id)], scale = 1)  # surrogate values
      } else {
        trunc_bounds <- getBounds(object, which.categ = which.categ)  # truncation bounds
        sim_trunc(n = length(y), distribution = distribution,
                  a = trunc_bounds[y[boot_id]],
                  b = trunc_bounds[y[boot_id] + 1L],
                  location = mean_response[boot_id + (which.categ - 1)*length(boot_id)], scale = 1)  # surrogate values
      }
    } else {
      stop("Distribution not supported.", call. = FALSE)
    }
    
  } else {  # jittering approach
    
    # Determine scale for jittering
    jitter.scale <- match.arg(jitter.scale)
    y <- getResponseValues(object)
    if (is.null(boot_id)) {
      boot_id <- seq_along(y)
    }
    y <- y[boot_id]
    prob <- getFittedProbs(object)[boot_id, ]
    if (jitter.scale == "response") {  # jittering on the response scale
      j <- seq_len(ncol(prob))
      jmat <- matrix(rep(j, times = nrow(prob)), ncol = ncol(prob), byrow = TRUE)
      runif(length(y), min = y, max = y + 1)
    } else {  # jittering on the probability scale
      if (getDistributionName(object) != "logis") {
        stop("Jittering on the probability scale is currently only supported",
             " for logit-type models.", call. = FALSE)
      }
      .min <- pbinom(y - 2, size = 1, prob = prob[, 1L, drop = TRUE])  # F(y-1)
      .max <- pbinom(y - 1, size = 1, prob = prob[, 1L, drop = TRUE])  # F(y)
      runif(length(y), min = .min, max = .max)  # S|Y=y - E(S|X)
    }
    
  }
  
  # Return results
  s
  
}


generate_residuals <- function(object, method = c("latent", "jitter"),
                               jitter.scale = c("response", "probability"),
                               boot_id = NULL,
                               which.categ = 1,
                               binarized = FALSE) {
  
  # Match arguments
  method <- match.arg(method)
  
  # Generate surrogate response values
  r <- if (method == "latent") {  # latent variable approach
    
    # Get distribution name (for sampling)
    distribution <- getDistributionName(object)  # distribution name
    
    # Simulate surrogate response values from the appropriate truncated
    # distribution
    if (distribution %in% c("norm", "logis", "cauchy", "gumbel", "Gumbel")) {
      y <- getResponseValues(object)
      if (is.null(boot_id)) {
        boot_id <- seq_along(y)
      }
      mean_response <- getMeanResponse(object)  # mean response values
      
      s <- if (binarized == FALSE){
        if (!inherits(object, what = "lrm") &&
            inherits(object, what = "glm")) {
          sim_trunc(n = length(y), distribution = distribution,
                    # {0, 1} -> {1, 2}
                    a = ifelse(y[boot_id] == 1, yes = -Inf, no = 0),
                    b = ifelse(y[boot_id] == 2, yes =  Inf, no = 0),
                    location = mean_response[boot_id  + (which.categ - 1)*length(y)], scale = 1)  # surrogate values
        } else {
          trunc_bounds <- getBounds(object, which.categ = which.categ)  # truncation bounds
          sim_trunc(n = length(y), distribution = distribution,
                    a = trunc_bounds[y[boot_id]],
                    b = trunc_bounds[y[boot_id] + 1L],
                    location = mean_response[boot_id + (which.categ - 1)*length(y)], scale = 1)  # surrogate values
        }
      } else {
        # Making y into a binary: either "<= specified category" or "> specified category"
        y <- ifelse(y <= which.categ, 1, 2)
        sim_trunc(n = length(y), distribution = distribution,
                  a = ifelse(y[boot_id] == 1, yes = -Inf, no = 0),
                  b = ifelse(y[boot_id] == 2, yes =  Inf, no = 0),
                  location = mean_response[boot_id  + (which.categ - 1)*length(y)], scale = 1)
      }
    }  else {
      stop("Distribution not supported.", call. = FALSE)
    }
    
    s - mean_response[boot_id + (which.categ - 1)*length(y)]  # surrogate residuals
  } else {  # jittering approach
    jitter.scale <- match.arg(jitter.scale)
    y <- getResponseValues(object)
    if (is.null(boot_id)) {
      boot_id <- seq_along(y)
    }
    y <- y[boot_id]
    prob <- getFittedProbs(object)[boot_id, ]
    if (jitter.scale == "response") {  # jittering on the response scale
      j <- seq_len(ncol(prob))
      jmat <- matrix(rep(j, times = nrow(prob)), ncol = ncol(prob), byrow = TRUE)
      runif(length(y), min = y, max = y + 1) - rowSums((jmat + 0.5) * prob)
    } else {  # jittering on the probability scale
      if (getDistributionName(object) != "logis") {
        stop("Jittering on the probability scale is currently only supported",
             " for logit-type models.", call. = FALSE)
      }
      .min <- pbinom(y - 2, size = 1, prob = prob[, 1L, drop = TRUE])  # F(y-1)
      .max <- pbinom(y - 1, size = 1, prob = prob[, 1L, drop = TRUE])  # F(y)
      runif(length(y), min = .min, max = .max) - 0.5  # S|Y=y - E(S|X)
    }
  }
  
  # Return results
  r
  
}
















#' Residual plots
#'
#' Residual-based diagnostic plots for cumulative link and general regression
#' models using \code{\link[ggplot2]{ggplot2}} graphics.
#'
#' @param object An object of class \code{\link[ordinal]{clm}},
#' \code{\link[stats]{glm}}, \code{\link[rms]{lrm}}, \code{\link[rms]{orm}},
#' \code{\link[MASS]{polr}}, or \code{\link[VGAM]{vglm}}.
#'
#' @param what Character string specifying what to plot. Default is \code{"qq"}
#' which produces a quantile-quantile plots of the residuals.
#'
#' @param x A vector giving the covariate values to use for residual-by-
#' covariate plots (i.e., when \code{what = "covariate"}).
#'
#' @param fit The fitted model from which the residuals were extracted. (Only
#' required if \code{what = "fitted"} and \code{object} inherits from class
#' \code{"resid"}.)
#'
#' @param distribution Function that computes the quantiles for the reference
#' distribution to use in the quantile-quantile plot. Default is \code{qnorm}
#' which is only appropriate for models using a probit link function. When
#' \code{jitter.scale = "probability"}, the reference distribution is always
#' U(-0.5, 0.5). (Only
#' required if \code{object} inherits from class \code{"resid"}.)
#'
#' @param ncol Integer specifying the number of columns to use for the plot
#' layout (if requesting multiple plots). Default is \code{NULL}.
#'
#' @param alpha A single values in the interval [0, 1] controlling the opacity
#' alpha of the plotted points. Only used when \code{nsim} > 1.
#'
#' @param xlab Character string giving the text to use for the x-axis label in
#' residual-by-covariate plots. Default is \code{NULL}.
#'
#' @param color Character string or integer specifying what color to use for the
#' points in the residual vs fitted value/covariate plot.
#' Default is \code{"black"}.
#'
#' @param shape Integer or single character specifying a symbol to be used for
#' plotting the points in the residual vs fitted value/covariate plot.
#'
#' @param size Numeric value specifying the size to use for the points in the
#' residual vs fitted value/covariate plot.
#'
#' @param qqpoint.color Character string or integer specifying what color to use
#' for the points in the quantile-quantile plot.
#'
#' @param qqpoint.shape Integer or single character specifying a symbol to be
#' used for plotting the points in the quantile-quantile plot.
#'
#' @param qqpoint.size Numeric value specifying the size to use for the points
#' in the quantile-quantile plot.
#'
#' @param qqline.color Character string or integer specifying what color to use
#' for the points in the quantile-quantile plot.
#'
#' @param qqline.linetype Integer or character string (e.g., \code{"dashed"})
#' specifying the type of line to use in the quantile-quantile plot.
#'
#' @param qqline.size Numeric value specifying the thickness of the line in the
#' quantile-quantile plot.
#'
#' @param smooth Logical indicating whether or not too add a nonparametric
#' smooth to certain plots. Default is \code{TRUE}.
#'
#' @param smooth.color Character string or integer specifying what color to use
#' for the nonparametric smooth.
#'
#' @param smooth.linetype Integer or character string (e.g., \code{"dashed"})
#' specifying the type of line to use for the nonparametric smooth.
#'
#' @param smooth.size Numeric value specifying the thickness of the line for the
#' nonparametric smooth.
#'
#' @param fill Character string or integer specifying the color to use to fill
#' the boxplots for residual-by-covariate plots when \code{x} is of class
#' \code{"factor"}. Default is \code{NULL} which colors the boxplots according
#' to the factor levels.
#'
#' @param ... Additional optional arguments to be passed onto
#' \code{\link[sure]{resids}}.
#'
#' @return A \code{"ggplot"} object.
#'
#' @rdname autoplot.resid
#'
#' @export
#'
#' @examples
#' # See ?resids for an example
#' ?resids
autoplot.resid <- function(
    object,
    what = c("qq", "fitted", "covariate"),
    which.categ = 1,
    binarized = FALSE,
    x = NULL,
    fit = NULL,
    distribution = qnorm,
    ncol = NULL,
    alpha = 1,
    xlab = NULL,
    color = "#444444",
    shape = 19,
    size = 2,
    qqpoint.color = "#444444",
    qqpoint.shape = 19,
    qqpoint.size = 2,
    qqline.color = "#888888",
    qqline.linetype = "dashed",
    qqline.size = 1,
    smooth = TRUE,
    smooth.color = "red",
    smooth.linetype = 1,
    smooth.size = 1,
    fill = NULL,
    ...
) {
  
  
  # What type of plot to produce
  what <- match.arg(what, several.ok = TRUE)
  
  # Figure out number of plots and layout
  np <- length(what)
  if (is.null(ncol)) {
    ncol <- length(what)
  }
  
  # Check that fitted mean response values are available
  if ("fitted" %in% what) {
    if (is.null(fit)) {
      stop("Cannot extract mean response. Please supply the original fitted",
           " model object via the `fit` argument.")
    }
    y.vec <- getResponseValues(fit)
    mr <- getMeanResponse(fit)[(which.categ - 1)*length(y.vec) + 1:length(y.vec)]
  }
  
  # Check that covariate values are supplied
  if ("covariate"  %in% what) {
    if (is.null(x)) {
      stop("No covariate to plot. Please supply a vector of covariate values",
           " via the `x` argument")
    }
    if (is.null(xlab)) {
      # xlab <- getColumnName(x)
      xlab <- deparse(substitute(x))
    }
  }
  
  # Deal with bootstrap replicates
  if (is.null(attr(object, "boot_reps"))) {
    nsim <- 1
    res <- object
    if ("qq" %in% what) {
      res.med <- object
    }
  } else {
    res.mat <- attr(object, "boot_reps")
    res <- as.numeric(as.vector(res.mat))
    nsim <- ncol(res.mat)
    if ("qq" %in% what) {
      res.med <- apply(apply(res.mat, MARGIN = 2, FUN = sort,
                             decreasing = FALSE), MARGIN = 1, FUN = median)
    }
    if ("fitted" %in% what) {
      mr <- mr[as.vector(attr(object, "boot_id"))]
    }
    if ("covariate" %in% what) {
      x <- x[as.vector(attr(object, "boot_id"))]
    }
  }
  
  # Quantile-quantile
  p1 <- if ("qq" %in% what) {
    if (!is.null(attr(object, "jitter.scale"))) {
      if (attr(object, "jitter.scale") == "response") {
        stop("Q-Q plots are not available for jittering on the response scale.")
      }
    }
    distribution <- match.fun(distribution)
    xvals <- distribution(ppoints(length(res.med)))[order(order(res.med))]
    qqline.y <- quantile(res.med, probs = c(0.25, 0.75),
                         names = FALSE, na.rm = TRUE)
    qqline.x <- distribution(c(0.25, 0.75))
    slope <- diff(qqline.y) / diff(qqline.x)
    int <- qqline.y[1L] - slope * qqline.x[1L]
    ggplot(data.frame(x = xvals, y = res.med), aes_string(x = "x", y = "y")) +
      geom_point(color = qqpoint.color, shape = qqpoint.shape,
                 size = qqpoint.size) +
      geom_abline(slope = slope, intercept = int, color = qqline.color,
                  linetype = qqline.linetype, size = qqline.size) +
      labs(x = "Theoretical quantile", y = "Sample quantile")
  } else {
    NULL
  }
  
  # Residual vs fitted value
  p2 <- if ("fitted" %in% what) {
    p <- ggplot(data.frame("x" = mr, "y" = res), aes_string(x = "x", y = "y")) +
      geom_point(color = color, shape = shape, size = size, alpha = alpha) +
      labs(x = "Fitted value", y = "Surrogate residual")
    if (smooth) {
      p <- p + geom_smooth(color = smooth.color, linetype = smooth.linetype,
                           size = smooth.size, se = FALSE)
    }
    p
  } else {
    NULL
  }
  
  # Residual vs covariate
  p3 <- if ("covariate" %in% what) {
    p <- ggplot(data.frame("x" = x, "y" = res), aes_string(x = "x", y = "y"))
    if (is.factor(x)) {
      if (is.null(fill)) {
        p <- p + geom_boxplot(aes_string(fill = "x"), alpha = alpha) +
          guides(fill = FALSE)
      } else {
        p <- p + geom_boxplot()
      }
    } else {
      p <- p + geom_point(color = color, shape = shape, size = size,
                          alpha = alpha)
      if (smooth) {
        p <- p + geom_smooth(color = smooth.color, linetype = smooth.linetype,
                             size = smooth.size, se = FALSE, method = "gam", formula = y ~ s(x, bs = "cs", k=3))
      }
    }
    p + labs(x = xlab, y = "Surrogate residual")
  } else {
    NULL
  }
  
  # Return plot(s)
  if (length(what) == 1) {  # return a single plot
    if (what == "qq") {
      p1
    } else if (what == "fitted") {
      p2
    } else {
      p3
    }
  } else {  # return multiple plots
    plots <- list(p1, p2, p3)
    grid.arrange(grobs = plots[!unlist(lapply(plots, FUN = is.null))],
                 ncol = ncol)
  }
  
}


#' @rdname autoplot.resid
#'
#' @export
autoplot.clm <- function(
    object,
    what = c("qq", "fitted", "covariate"),
    which.categ = 1,
    binarized = FALSE,
    x = NULL,
    fit = NULL,
    distribution = qnorm,
    ncol = NULL,
    alpha = 1,
    xlab = NULL,
    color = "#444444",
    shape = 19,
    size = 2,
    qqpoint.color = "#444444",
    qqpoint.shape = 19,
    qqpoint.size = 2,
    qqline.color = "#888888",
    qqline.linetype = "dashed",
    qqline.size = 1,
    smooth = TRUE,
    smooth.color = "red",
    smooth.linetype = 1,
    smooth.size = 1,
    fill = NULL,
    ...
) {
  
  # Compute residuals
  res <- resids(object, which.categ = which.categ, binarized = binarized, ...)
  
  # Quantile function to use for Q-Q plots
  qfun <- if (is.null(attr(res, "jitter.scale"))) {
    getQuantileFunction(object)
  } else {
    if (what == "qq" && attr(res, "jitter.scale") == "response") {
      stop("Quantile-quantile plots are not appropriate for residuals ",
           "obtained by jittering on the response scale.")
    }
    function(p) qunif(p, min = -0.5, max = 0.5)
  }
  
  # Default x-axis label
  if (is.null(xlab)) {
    xlab <- deparse(substitute(x))
  }
  
  # Call the default method
  autoplot.resid(
    res, what = what, which.categ = which.categ, x = x, distribution = qfun, fit = object, ncol = ncol,
    alpha = alpha, xlab = xlab, color = color, shape = shape, size = size,
    qqpoint.color = qqpoint.color, qqpoint.shape = qqpoint.shape,
    qqpoint.size = qqpoint.size, qqline.color = qqline.color,
    qqline.linetype = qqline.linetype, qqline.size = qqline.size,
    smooth = smooth, smooth.color = smooth.color,
    smooth.linetype = smooth.linetype, smooth.size = smooth.size, fill = fill
  )
  
}


#' @rdname autoplot.resid
#'
#' @export
autoplot.glm <- autoplot.clm


#' @rdname autoplot.resid
#'
#' @export
autoplot.lrm <- autoplot.clm


#' @rdname autoplot.resid
#'
#' @export
autoplot.orm <- autoplot.clm


#' @rdname autoplot.resid
#'
#' @export
autoplot.polr <- autoplot.clm


#' @rdname autoplot.resid
#'
#' @export
autoplot.vglm <- autoplot.clm




#' Simulated quadratic data
#'
#' Data simulated from a probit model with a quadratic trend. The data are
#' described in Example 2 of Liu and Zhang (2017).
#'
#' @docType data
#'
#' @keywords datasets
#'
#' @format A data frame with 2000 rows and 2 variables.
#' \itemize{
#'   \item \code{x} The predictor variable.
#'   \item \code{y} The response variable; an ordered factor.
#' }
#'
#' @references
#' Liu, Dungang and Zhang, Heping. Residuals and Diagnostics for Ordinal
#' Regression Models: A Surrogate Approach.
#' \emph{Journal of the American Statistical Association} (accepted).
#'
#' @name df1
#'
#' @usage
#' data(df1)
#'
#' @examples
#' head(df1)
NULL


#' Simulated heteroscedastic data
#'
#' Data simulated from a probit model with heteroscedasticity. The data are
#' described in Example 4 of Liu and Zhang (2017).
#'
#' @docType data
#'
#' @keywords datasets
#'
#' @format A data frame with 2000 rows and 2 variables.
#' \itemize{
#'   \item \code{x} The predictor variable.
#'   \item \code{y} The response variable; an ordered factor.
#' }
#'
#' @references
#' Liu, Dungang and Zhang, Heping. Residuals and Diagnostics for Ordinal
#' Regression Models: A Surrogate Approach.
#' \emph{Journal of the American Statistical Association} (accepted).
#'
#' @name df2
#'
#' @usage
#' data(df2)
#'
#' @examples
#' head(df2)
NULL


#' Simulated Gumbel data
#'
#' Data simulated from a log-log model with a quadratic trend. The data are
#' described in Example 3 of Liu and Zhang (2017).
#'
#' @docType data
#'
#' @keywords datasets
#'
#' @format A data frame with 2000 rows and 2 variables.
#' \itemize{
#'   \item \code{x} The predictor variable.
#'   \item \code{y} The response variable; an ordered factor.
#' }
#'
#' @references
#' Liu, Dungang and Zhang, Heping. Residuals and Diagnostics for Ordinal
#' Regression Models: A Surrogate Approach.
#' \emph{Journal of the American Statistical Association} (accepted).
#'
#' @name df3
#'
#' @usage
#' data(df3)
#'
#' @examples
#' head(df3)
NULL


#' Simulated proportionality data
#'
#' Data simulated from from two separate probit models. The data are described
#' in Example 5 of Liu and Zhang (2017).
#'
#' @docType data
#'
#' @keywords datasets
#'
#' @format A data frame with 4000 rows and 2 variables.
#' \itemize{
#'   \item \code{x} The predictor variable.
#'   \item \code{y} The response variable; an ordered factor.
#' }
#'
#' @references
#' Liu, Dungang and Zhang, Heping. Residuals and Diagnostics for Ordinal
#' Regression Models: A Surrogate Approach.
#' \emph{Journal of the American Statistical Association} (accepted).
#'
#' @name df4
#'
#' @usage
#' data(df4)
#'
#' @examples
#' head(df4)
NULL


#' Simulated interaction data
#'
#' Data simulated from from an ordered probit model with an interaction effect.
#'
#' @docType data
#'
#' @keywords datasets
#'
#' @format A data frame with 2000 rows and 3 variables.
#' \itemize{
#'   \item \code{x1} A continuous predictor variable.
#'   \item \code{x2} A factor with two levels: \code{"Control"} and
#'   \code{"Treatment"}.
#'   \item \code{y} The response variable; an ordered factor.
#' }
#'
#' @references
#' Liu, Dungang and Zhang, Heping. Residuals and Diagnostics for Ordinal
#' Regression Models: A Surrogate Approach.
#' \emph{Journal of the American Statistical Association} (accepted).
#'
#' @name df5
#'
#' @usage
#' data(df5)
#'
#' @examples
#' head(df5)
NULL













#' Goodness-of-Fit Simulation
#'
#' Simulate p-values from a goodness-of-fit test.
#'
#' @param object An object of class \code{\link[ordinal]{clm}},
#' \code{\link[stats]{glm}}, \code{\link[rms]{lrm}}, \code{\link[rms]{orm}},
#' \code{\link[MASS]{polr}}, or \code{\link[VGAM]{vglm}}.
#'
#' @param nsim Integer specifying the number of bootstrap replicates to use.
#'
#' @param test Character string specifying which goodness-of-fit test to use.
#' Current options include: \code{"ks"} for the Kolmogorov-Smirnov test,
#' \code{"ad"} for the Anderson-Darling test, and \code{"cvm"} for the
#' Cramer-Von Mises test. Default is \code{"ks"}.
#'
#' @param ... Additional optional arguments. (Currently ignored.)
#'
#' @param x An object of class \code{"gof"}.
#'
#' @return A numeric vector of class \code{"gof", "numeric"} containing the
#' simulated p-values.
#'
#' @details
#' Under the null hypothesis, the distribution of the p-values should appear
#' uniformly distributed on the interval [0, 1]. This can be visually
#' investigated using the \code{plot} method. A 45 degree line is indicative of
#' a "good" fit.
#'
#' @rdname gof
#'
#' @export
#'
#' @examples
#' # See ?resids for an example
#' ?resids
gof <- function(object, nsim = 10, test = c("ks", "ad", "cvm"), ...) {
  if (nsim <- as.integer(nsim) < 2) {
    stop("nsim must be a postive integer >= 2")
  }
  UseMethod("gof")
}


#' @rdname gof
#'
#' @export
gof.default <- function(object, nsim = 10, which.categ = 1, binarized = FALSE, test = c("ks", "ad", "cvm"), ...) {
  res <- resids(object, which.categ = which.categ, binarized = binarized, nsim = nsim)

  test <- match.arg(test)
  pfun <- getDistributionFunction(object)
  sim_pvals(res, test = test, pfun = pfun)
}


#' @rdname gof
#'
#' @export
plot.gof <- function(x, ...) {
  graphics::plot(stats::ecdf(x), xlab = "p-value", xlim = c(0, 1), ...)
  graphics::abline(0, 1, lty = 2)
}


#' @keywords internal
sim_pvals <- function(res, test, pfun) {
  gof_test <- switch(test, "ks" = stats::ks.test, "ad" = goftest::ad.test,
                     "cvm" = goftest::cvm.test)
  pvals <- apply(attr(res, "boot_reps"), MARGIN = 2, FUN = function(x) {
    gof_test(x, pfun)$p.value
  })
  class(pvals) <- c("gof", "numeric")
  pvals
}


















#' Arrange multiple grobs on a page
#'
#' See \code{\link[gridExtra]{grid.arrange}} for more details.
#'
#' @name grid.arrange
#'
#' @rdname grid.arrange
#'
#' @keywords internal
#'
#' @export
#'
#' @importFrom gridExtra grid.arrange
#'
#' @usage grid.arrange(..., newpage = TRUE)
NULL

















#' Surrogate residuals
#'
#' Simulate surrogate response values for cumulative link regression models
#' using the latent method described in Liu and Zhang (2017).
#'
#' @param object An object of class \code{\link[ordinal]{clm}},
#' \code{\link[stats]{glm}}, \code{\link[rms]{lrm}}, \code{\link[rms]{orm}},
#' \code{\link[MASS]{polr}}, or \code{\link[VGAM]{vglm}}.
#'
#' @param nsim Integer specifying the number of bootstrap replicates to use.
#' Default is \code{1L} meaning no bootstrap samples.
#'
#' @param method Character string specifying which method to use to generate the
#' surrogate response values. Current options are \code{"latent"} and
#' \code{"jitter"}. Default is \code{"latent"}.
#'
#' @param jitter.scale Character string specifyint the scale on which to perform
#' the jittering whenever \code{method = "jitter"}. Current options are
#' \code{"response"} and \code{"probability"}. Default is \code{"response"}.
#'
#' @param ... Additional optional arguments. (Currently ignored.)
#'
#' @return A numeric vector of class \code{c("numeric", "surrogate")} containing
#' the simulated surrogate response values. Additionally, if \code{nsim} > 1,
#' then the result will contain the attributes:
#' \describe{
#'   \item{\code{boot_reps}}{A matrix  with \code{nsim} columns, one for each
#'   bootstrap replicate of the surrogate values. Note, these are random and do
#'   not correspond to the original ordering of the data;}
#'   \item{\code{boot_id}}{A matrix  with \code{nsim} columns. Each column
#'   contains the observation number each surrogate value corresponds to in
#'   \code{boot_reps}. (This is used for plotting purposes.)}
#' }
#'
#' @note
#' Surrogate response values require sampling from a continuous distribution;
#' consequently, the result will be different with every call to
#' \code{surrogate}. The internal functions used for sampling from truncated
#' distributions are based on modified versions of
#' \code{\link[truncdist]{rtrunc}} and \code{\link[truncdist]{qtrunc}}.
#'
#' For \code{"glm"} objects, only the \code{binomial()} family is supported.
#'
#' @references
#' Liu, Dungang and Zhang, Heping. Residuals and Diagnostics for Ordinal
#' Regression Models: A Surrogate Approach.
#' \emph{Journal of the American Statistical Association} (accepted). URL
#' http://www.tandfonline.com/doi/abs/10.1080/01621459.2017.1292915?journalCode=uasa20
#'
#' Nadarajah, Saralees and Kotz, Samuel. R Programs for Truncated Distributions.
#' \emph{Journal of Statistical Software, Code Snippet}, 16(2), 1-8, 2006. URL
#' https://www.jstatsoft.org/v016/c02.
#'
#' @export
#'
#' @examples
#' # Generate data from a quadratic probit model
#' set.seed(101)
#' n <- 2000
#' x <- runif(n, min = -3, max = 6)
#' z <- 10 + 3 * x - 1 * x^2 + rnorm(n)
#' y <- ifelse(z <= 0, yes = 0, no = 1)
#'
#' # Scatterplot matrix
#' pairs(~ x + y + z)
#'
#' # Setup for side-by-side plots
#' par(mfrow = c(1, 2))
#'
#' # Misspecified mean structure
#' fm1 <- glm(y ~ x, family = binomial(link = "probit"))
#' scatter.smooth(x, y = resids(fm1),
#'                main = "Misspecified model",
#'                ylab = "Surrogate residual",
#'                lpars = list(lwd = 3, col = "red2"))
#' abline(h = 0, lty = 2, col = "blue2")
#'
#' # Correctly specified mean structure
#' fm2 <- glm(y ~ x + I(x ^ 2), family = binomial(link = "probit"))
#' scatter.smooth(x, y = resids(fm2),
#'                main = "Correctly specified model",
#'                ylab = "Surrogate residual",
#'                lpars = list(lwd = 3, col = "red2"))
#' abline(h = 0, lty = 2, col = "blue2")
resids <- function(object, nsim = 1L, method = c("latent", "jitter"),
                   jitter.scale = c("response", "probability"), 
                   which.categ = 1, binarized = FALSE, ...) {
  
  # Match arguments
  method <- match.arg(method)
  jitter.scale = match.arg(jitter.scale)
  
  # Issue warning for jittering method
  if (method == "jitter") {
    warning("Jittering is an experimental feature, use at your own risk!",
            call. = FALSE)
  }
  
  # Generate surrogate response values
  r <- generate_residuals(object, method = method, jitter.scale = jitter.scale, which.categ = which.categ, binarized = binarized)
  
  # Multiple samples
  if (nsim > 1L) {  # bootstrap
    boot_reps <- boot_id <- matrix(nrow = if (class(object) == "ordinalNet") length(object$args$y) else nobs(object), ncol = nsim)
    for(i in seq_len(nsim)) {
      boot_id[, i] <- sample(if (class(object) == "ordinalNet") length(object$args$y) else nobs(object), replace = TRUE)
      boot_reps[, i] <-
        generate_residuals(object, method = method, jitter.scale = jitter.scale,
                           boot_id = boot_id[, i, drop = TRUE],
                           which.categ = which.categ,
                           binarized = binarized)
    }
    attr(r, "boot_reps") <- boot_reps
    attr(r, "boot_id") <- boot_id
  }
  
  # Return residuals
  class(r) <- c("numeric", "resid")
  r
  
}























#' sure: An R package for constructing surrogate-based residuals and diagnostics
#' for ordinal and general regression models.
#'
#' The \code{sure} package provides surrogate-based residuals for fitted ordinal
#' and general (e.g., binary) regression models of class
#' \code{\link[ordinal]{clm}}, \code{\link[stats]{glm}}, \code{\link[rms]{lrm}},
#' \code{\link[rms]{orm}}, \code{\link[MASS]{polr}}, or
#' \code{\link[VGAM]{vglm}}.
#'
#' The development version can be found on GitHub:
#' \url{https://github.com/AFIT-R/sure}. As of right now, \code{sure} exports the
#' following functions:
#' \itemize{
#'   \item{\code{resids}} - construct (surrogate-based) residuals;
#'   \item{\code{autoplot}} - plot diagnostics using
#'   \code{\link[ggplot2]{ggplot2}}-based graphics;
#'   \item{\code{gof}} - simulate p-values from a goodness-of-fit test.
#' }
#'
#' @references
#' Liu, Dungang and Zhang, Heping. Residuals and Diagnostics for Ordinal
#' Regression Models: A Surrogate Approach.
#' \emph{Journal of the American Statistical Association} (accepted).
#'
#' @importFrom ggplot2  aes_string geom_abline geom_boxplot geom_point
#'
#' @importFrom ggplot2 geom_smooth ggplot ggtitle guides labs xlab ylab
#'
#' @importFrom stats .checkMFClasses lowess median model.frame model.matrix
#'
#' @importFrom stats model.response nobs pbinom pcauchy plogis pnorm ppoints
#'
#' @importFrom stats predict qcauchy qlogis qnorm qqline qqplot qqnorm quantile
#'
#' @importFrom stats qunif runif
#'
#' @docType package
#'
#' @name sure
NULL



















#' Surrogate response
#'
#' Simulate surrogate response values for cumulative link regression models
#' using the latent method described in Liu and Zhang (2017).
#'
#' @param object An object of class \code{\link[ordinal]{clm}},
#' \code{\link[stats]{glm}} \code{\link[rms]{lrm}}, \code{\link[rms]{orm}},
#' \code{\link[MASS]{polr}}, or \code{\link[VGAM]{vglm}}.
#'
#' @param nsim Integer specifying the number of bootstrap replicates to use.
#' Default is \code{1L} meaning no bootstrap samples.
#'
#' @param method Character string specifying which method to use to generate the
#' surrogate response values. Current options are \code{"latent"} and
#' \code{"jitter"}. Default is \code{"latent"}.
#'
#' @param jitter.scale Character string specifyint the scale on which to perform
#' the jittering whenever \code{method = "jitter"}. Current options are
#' \code{"response"} and \code{"probability"}. Default is \code{"response"}.
#'
#' @param ... Additional optional arguments. (Currently ignored.)
#'
#' @return A numeric vector of class \code{c("numeric", "surrogate")} containing
#' the simulated surrogate response values. Additionally, if \code{nsim} > 1,
#' then the result will contain the attributes:
#' \describe{
#'   \item{\code{boot_reps}}{A matrix  with \code{nsim} columns, one for each
#'   bootstrap replicate of the surrogate values. Note, these are random and do
#'   not correspond to the original ordering of the data;}
#'   \item{\code{boot_id}}{A matrix  with \code{nsim} columns. Each column
#'   contains the observation number each surrogate value corresponds to in
#'   \code{boot_reps}. (This is used for plotting purposes.)}
#' }
#'
#' @note
#' Surrogate response values require sampling from a continuous distribution;
#' consequently, the result will be different with every call to
#' \code{surrogate}. The internal functions used for sampling from truncated
#' distributions are based on modified versions of
#' \code{\link[truncdist]{rtrunc}} and \code{\link[truncdist]{qtrunc}}.
#'
#' For \code{"glm"} objects, only the \code{binomial()} family is supported.
#'
#' @references
#' Liu, Dungang and Zhang, Heping. Residuals and Diagnostics for Ordinal
#' Regression Models: A Surrogate Approach.
#' \emph{Journal of the American Statistical Association} (accepted). URL
#' http://www.tandfonline.com/doi/abs/10.1080/01621459.2017.1292915?journalCode=uasa20
#'
#' Nadarajah, Saralees and Kotz, Samuel. R Programs for Truncated Distributions.
#' \emph{Journal of Statistical Software, Code Snippet}, 16(2), 1-8, 2006. URL
#' https://www.jstatsoft.org/v016/c02.
#'
#' @export
#'
#' @examples
#' # Generate data from a quadratic probit model
#' set.seed(101)
#' n <- 2000
#' x <- runif(n, min = -3, max = 6)
#' z <- 10 + 3*x - 1*x^2 + rnorm(n)
#' y <- ifelse(z <= 0, yes = 0, no = 1)
#'
#' # Scatterplot matrix
#' pairs(~ x + y + z)
#'
#' # Setup for side-by-side plots
#' par(mfrow = c(1, 2))
#'
#' # Misspecified mean structure
#' fm1 <- glm(y ~ x, family = binomial(link = "probit"))
#' s1 <- surrogate(fm1)
#' scatter.smooth(x, s1 - fm1$linear.predictors,
#'                main = "Misspecified model",
#'                ylab = "Surrogate residual",
#'                lpars = list(lwd = 3, col = "red2"))
#' abline(h = 0, lty = 2, col = "blue2")
#'
#' # Correctly specified mean structure
#' fm2 <- glm(y ~ x + I(x ^ 2), family = binomial(link = "probit"))
#' s2 <- surrogate(fm2)
#' scatter.smooth(x, s2 - fm2$linear.predictors,
#'                main = "Correctly specified model",
#'                ylab = "Surrogate residual",
#'                lpars = list(lwd = 3, col = "red2"))
#' abline(h = 0, lty = 2, col = "blue2")
surrogate <- function(object, nsim = 1L, method = c("latent", "jitter"),
                      jitter.scale = c("response", "probability"),
                      which.categ = 1, ...) {
  
  # Match arguments
  method <- match.arg(method)
  jitter.scale = match.arg(jitter.scale)
  
  # Issue warning for jittering method
  if (method == "jitter") {
    warning("Jittering is an experimental feature, use at your own risk!",
            call. = FALSE)
  }
  
  # Generate surrogate response values
  s <- generate_surrogate(object, method = method, jitter.scale = jitter.scale, which.categ = which.categ)
  
  # Multiple (re)samples
  if (nsim > 1L) {  # bootstrap
    boot_reps <- boot_id <- matrix(nrow = if (class(object) == "ordinalNet") length(object$args$y) else nobs(object), ncol = nsim)
    for(i in seq_len(nsim)) {
      boot_id[, i] <- sample(if (class(object) == "ordinalNet") length(object$args$y) else nobs(object), replace = TRUE)
      boot_reps[, i] <-
        generate_surrogate(object, method = method, jitter.scale = jitter.scale,
                           boot_id = boot_id[, i, drop = TRUE],
                           which.categ = which.categ)
    }
    attr(s, "boot_reps") <- boot_reps
    attr(s, "boot_id") <- boot_id
  }
  
  # Return residuals
  class(s) <- c("numeric", "surrogate")
  s
  
}




########
## All the newly-added "ordinalNet"-specific function
########

getDistributionName.ordinalNet <- function(object) {
  switch(object$args$link,
         "logit" = "logis",
         "probit" = "norm",
         "loglog" = "gumbel",
         "cloglog" = "Gumbel",
         "cauchit" = "cauchy")
}

getResponseValues.ordinalNet <- function(object, ...) {
  unname(as.integer(object$args$y))
}


getMeanResponse.ordinalNet <- function(object) {
  if (object$args$reverse){
    predict.ordinalNet.mine(object, criteria="bic", type="link")[,c((object$nLev - 1):1)]
  } else {
    - predict.ordinalNet.mine(object, criteria="bic", type="link")
  }
}


ncat.ordinalNet <- function(object) {
  object$nLev
}


getBounds.ordinalNet <- function(object, which.categ=1, ...) {
  if (object$args$reverse) {
    coefs <- -unname(coef.ordinalNet.mine(object, criteria = "bic"))
    unname(c(-Inf, rev(coefs[seq_len(ncat(object) - 1)] - coefs[ncat(object) - which.categ]), Inf))
  } else {
    coefs <- unname(coef.ordinalNet.mine(object, criteria = "bic"))
    unname(c(-Inf, coefs[seq_len(ncat(object) - 1)] - coefs[which.categ], Inf))
  }
}


#' @keywords internal
getFittedProbs.ordinalNet <- function(object) {
  predict.ordinalNet.mine(object, criteria="bic")
}


#' @keywords internal
getQuantileFunction.ordinalNet <- function(object) {
  switch(object$args$link,
         "logit" = qlogis,
         "probit" = qnorm,
         "loglog" = qgumbel,
         "cloglog" = qGumbel,
         "cauchit" = qcauchy)
}


#' @rdname autoplot.resid
#'
#' @export
autoplot.ordinalNet <- autoplot.clm
