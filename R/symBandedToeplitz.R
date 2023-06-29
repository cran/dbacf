#' Creates a symmetric banded Toeplitz matrix
#'  
#' @param x numeric vector or an object of class \code{\link[dbacf]{dbacf}}.
#' @param n integer specifying number of columns (rows) of banded matrix. 
#'  
#' @export                             
#'  
#' @examples
#' alphas <- c(-2, 0.5, -4)
#' (true.acf <- ARMAacf(ma = alphas))
#' symBandedToeplitz(true.acf, n = 10)
#' 
#' @seealso \code{\link[dbacf]{nearPDToeplitz}}, \code{\link[Matrix]{bandSparse}}
# 
#' @return An \eqn{\code{n} \times \code{n}}{n x n} symmetric banded Toeplitz 
#' matrix whose entries in main band are given by object \code{x}.
#' 
symBandedToeplitz <- function(x, n) {
  if (is.numeric(x)) {
    ret <- .symBandedToeplitz(x=x, n=n)
  }
  
  if (inherits(x, "dbacf")){
    ret <- .symBandedToeplitz(x=x$acf, n=x$n)
  }
  # else {
  #   stop("'x' must be a numeric")
  # }
  # else {
  #   UseMethod("symBandedToeplitz")
  # }
  ret
}

# symBandedToeplitz.dbacf <- function(x, n = x$n) {
#   .symBandedToeplitz(x=x$acf, n=n)
# }

.symBandedToeplitz <- function(x, n) {
  x <- x[1:min(length(x), n)]
  bLis <- as.data.frame(matrix(x, nrow =  n, ncol = length(x), byrow = TRUE))
  as.matrix(Matrix::bandSparse(n, k = 0:(length(x) - 1), diagonals = bLis, 
                               symmetric = TRUE))
}
