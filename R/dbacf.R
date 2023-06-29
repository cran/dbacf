#' Difference-based (auto)covariance/correlation function estimation
#' 
#' Computes \emph{and by default plots} the (auto)covariance/correlation function
#' estimate without pre-estimating the underlying \emph{piecewise constant signal} 
#' of the observations. To that end, a class of second-order 
#' \emph{difference-based estimators} is implemented according to Eqs.(2.5)-(2.6)
#' of \cite{Tecuapetla-Gómez and Munk (2017)}. By default, this function computes
#' a subclass of estimates with minimal bias according to Eqs.(2.12)-(2.14) of the 
#' aforementioned paper.
#' 
#' @param    data numeric vector or a univariate object of class
#'                \code{\link[stats]{ts}} of length at least \code{2(m + 1)}.
#' @param       m integer scalar giving the underlying level of dependency.
#' @param       d numeric vector giving the weights used in difference-based
#'                estimation method. Only pertinent when \code{order=second}.
#'                If missing, the weights \code{d} are calculated according 
#'                to Eqs.(2.12)-(2.14) of \cite{Tecuapetla-Gómez and Munk (2017)}.
#'                When a single value \eqn{d^\ast}{d*} is specified, 
#'                \code{d = rep(}\eqn{d^\ast}{d*}\code{, m + 1)}.              
#' @param    type character string specifying whether covariance (default) 
#'                or correlation must be computed.
#' @param   order character specifying whether a \code{first} (default)
#'                or a \code{second} difference-based estimate should be employed.
#' @param    plot logical. If \code{TRUE} (default) the acf is plotted.
#' @param   \dots further arguments passed to \code{\link{plot.dbacf}}.
#'
#' @note Although the theoretical properties of the methods implemented 
#' in this function were derived for change point regression with stationary 
#' \emph{Gaussian} \eqn{m}-dependent errors, these methods have proven robust against 
#' non-normality of the errors and as efficient as other methods in which 
#' pre-estimation of an underlying smooth signal is required. For further 
#' details see Section 6 of \cite{Tecuapetla-Gómez and Munk (2017)}.
#' 
#' The first-order difference-based estimator was implemented following Eqs.(4)-(5)
#' of \cite{Levine and Tecuapetla-Gómez (2023)}. For the robustness of this estimator
#' see Section 4 of the just mentioned paper.
#'  
#' @export                             
#'  
#' @examples
#' ma2 <- arima.sim(n = 50, model = list(ma = c(0.4, -0.4), order = c(0, 0, 2)), 
#'                  sd = 0.25)
#' dbacf(data=ma2, m = 2)
#' dbacf(data=ma2, m = 2, order="first")
#'  
#' @seealso \code{\link[stats]{acf}}, \code{\link[dbacf]{plot.dbacf}}
# 
#' @return An object of class "dbacf" containing:
#' \item{acf}{numeric vector of length \code{m + 1} giving estimated
#'            (auto)covariance-correlation.}
#' \item{m}{integer giving underlying level of dependency.}
#' \item{d}{numeric vector containing the weights used to estimate acf.}
#' \item{acfType}{string indicating whether \code{covariance} or 
#'                \code{correlation} has been computed.}
#' \item{n}{integer giving \code{length(data)}.}
#' \item{series}{string with name of variable \code{data}.}
#' 
#' @references Tecuapetla-Gómez, I and Munk, A. (2017). \emph{Autocovariance
#' estimation in regression with a discontinuous signal and \eqn{m}-dependent errors: A 
#' difference-based approach}. Scandinavian Journal of Statistics, \bold{44(2)}, 346--368.
#' 
#' @references Levine, M. and Tecuapetla-Gómez, I. (2023). \emph{Autocovariance 
#' function estimation via difference schemes for a semiparametric change point model
#' with \eqn{m}-dependent errors}. Submitted.
#'
dbacf <- function(data, m, d, type = c("covariance", "correlation"), 
                  order = c("second", "first"), plot = TRUE, 
                  ...) {
  
  if (!is.numeric(data)) 
    stop("'data' must be numeric")
  
  if (length(data) > 2^31 - 1) {
    stop("length of 'data' must be smaller than 2^31 - 1")
  }  
  sampleT <- length(data)
  
  series <- deparse(substitute(data))
  
  if (missing(m)) {
    stop("'m' must be specified")
  } else {
    if ( !is.wholenumber(m) || m < 0)
      stop("m must be a positive integer")
  }
  
  m <- as.integer(m + 2 * .Machine$double.eps ^ 0.5)        
  
  if (length(data) <= 2L * (m + 1L))
    stop("length of data must be greater than 2 * (m + 1)")
  
  type <- match.arg(type)
  
  order <- match.arg(order)
  
  if(order == "first"){
    output <- dbacfFirstOrder(data=data, m=m)
  } else {
    output <- dbacfSecondOrder(data=data, m=m, d=d)
  }
  
  if (type == "correlation"){
    acf <- output$acf / output$acf[1]
  } else {
    acf <- output$acf
  }
  
  dbacf.out <- structure(list(acf = acf, m = m,
                              d = output$d, acfType = type, 
                              n = sampleT, series = series), 
                         class = "dbacf")  
  
  if (plot) {
    plot.dbacf(dbacf.out, ...)
    invisible(dbacf.out)
  }
  else dbacf.out
}
#
#' Plot autocovariance and autocorrelation functions
#'
#' This function returns the plot method for objects of class "dbacf".
#'  
#' @param           x an object of class "dbacf".
#' @param        type what type of plot should be drawn. For possible types see
#'                    \code{\link[graphics]{plot}}.
#' @param        xlab the x label of the plot.
#' @param        ylab the y label of the plot.
#' @param        xlim numeric vector of length 2 giving the \code{x} coordinates
#'                    range.
#' @param        main an overall title for the plot.
#' @param ltyZeroLine type of line used to draw horizontal line passing at 0.
#' @param colZeroLine string indicating color of horizontal line passing at 0.
#' @param       \dots extra arguments to be passed to plot.
#'  
#' @rdname plot.dbacf
#' @method plot dbacf
#' @export
#' 
#' @note \code{\link[dbacf]{dbacf}} documents the structure of objects of class "dbacf".
#' 
#' @return No return value
#' 
#' @importFrom graphics abline
#' 
#' @seealso \code{\link[stats]{acf}}, \code{\link[dbacf]{dbacf}}. 
#' 
plot.dbacf <- function(x, type = "h", xlab = "Lag", 
                       ylab = paste("ACF", ifelse(x$acfType == "covariance",
                                                  "(cov)", " ")),                        
                       xlim = c(0, x$m + 1), main = paste("Series", x$series),
                       ltyZeroLine = 3, colZeroLine = "blue", ...) {    
  plot(0:x$m, x$acf, type = type, xlab = xlab, 
       ylab = ylab, xlim = xlim, main = main, ...)
  abline(h = 0, lty = ltyZeroLine, col = colZeroLine)
}   
