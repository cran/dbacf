#' Robust dbacf in change point regression with AR(1) errors
#' 
#' In the context of change point regression with a stationary AR(1) error process, this function
#' estimates the autoregressive coefficient along with the autocovariance/correlation function
#' as a function of given lags.
#' 
#' @param data numeric vector or a univariate object of class \code{\link[stats]{ts}}.
#' @param type character string specifying whether covariance (default) or correlation must be computed.
#' @param lags numeric giving the number of lags to compute.
#' 
#' @export
#' 
#' @examples
#' ar1 <- arima.sim(n = 50, model = list(ar = c(0.5), order = c(1, 0, 0)), 
#'                  sd = 0.25)
#' dbacf_AR1(ar1, type="correlation", lags=10)
#' 
#' @importFrom stats median
#' 
#' @return An object of class "dbacf" containing:
#' \itemize{
#'    \item \code{acf} numeric vector of length \code{lags + 1} giving estimated (auto)covariance/correlation function
#'    \item \code{rho} numeric, estimate of autoregressive coefficient
#'    \item \code{acfType} string indicating whether \code{covariance} or \code{correlation} has been computed
#'    \item \code{n} integer giving \code{length(data)}
#' }
#' 
#' @references Chakar, S. and Lebarbier, E. and Lévy-Leduc, C. and Robin, S. (2017). \emph{A robust
#' approach for estimating change-points in the mean of an AR(1) process}, Bernoulli, \bold{23(2)},
#' 1408-1447
#' 
dbacf_AR1 <- function(data, type=c("covariance", "correlation"), lags){
  if (!is.numeric(data)) 
    stop("'data' must be numeric")
  
  if (length(data) > 2^31 - 1) {
    stop("length of 'data' must be smaller than 2^31 - 1")
  }  
  
  type <- match.arg(type)
  
  if (missing(lags)){
    stop("'lags' must be provided")
  } 
  
  rho <- median((diff(data, lag = 2))^2)/median(diff(data)^2) - 1
  
  acf <- numeric(lags+1)
  acf[1] <- 1 / (1-rho^2)
  for( h in 1:lags ){
    acf[h+1] <- rho^(h) * acf[1]
  }
  
  if (type == "correlation")
    acf <- acf / acf[1]
  
  dbacf.out <- structure(list(acf = acf, rho = rho,
                              acfType = type, n = length(data)), 
                         class = "dbacf")
  
  dbacf.out
}