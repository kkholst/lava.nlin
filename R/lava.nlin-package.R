

#' Non-linear structural equation models.
#' 
#' Estimation of parameters in non-linear structural equation models via MCEM,
#' Laplace Approximation or Importance Sampling
#' 
#' @name lava.nlin-package
#' @aliases lava.nlin-package lava.nlin
#' @import lava
#' @import numDeriv
#' @importFrom survival Surv
#' @importFrom graphics matplot plot
#' @importFrom stats lm as.formula na.omit coef vcov sd var pnorm
#'     runif rexp residuals printCoefmat model.frame model.extract
#'     model.matrix terms window
#' @importFrom utils tail
#' @importFrom coda as.mcmc
#' @importFrom MASS eqscplot
#' @importFrom Rcpp evalCpp
#' @useDynLib lava.nlin
#' @docType package
#' @author Klaus K. Holst <k.k.holst@@biostat.ku.dk>
#' @keywords package
NULL


#' For internal use
#' 
#' for internal use
#'
#' @name Eval
#' @aliases Eval Lapl Mstep_nsem1 Mstep_nsem1b Mstep_weibullmm as.mcmc.StEM
#' coef.StEM hessian.weibull finalize lap logl.weibull merge.StEM obj.weibull
#' plot.StEM print.StEM restart restart.StEM score.weibull sim.StEM sim.phwreg
#' weibull.maximize weibullmm window.StEM
#' @author Klaus K. Holst
#' @keywords utilities
NULL


