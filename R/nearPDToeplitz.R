#' Computes the nearest positive definite Toeplitz matrix
#' 
#' 
#' Computes the nearest positive definite Toeplitz matrix to an initial approximation, 
#' typically a covariance (correlation) matrix of a stationary process. This function 
#' implements an \emph{alternating projection algorithm} that combines \cite{Grigoriadis et al. (1994)} 
#' and \cite{Higham (2002)}. For further details see Section 5 of \cite{Tecuapetla-Gómez and Munk (2017)}.
#'            
#' @param               matrix a symmetric matrix.
#' @param                 type string indicating whether the elements of the main 
#'                             diagonal must be all equal to 1 (\code{"correlation"}) 
#'                             or not (\code{"covariance"}).
#' @param       toleranceEigen defines relative positiveness of eigenvalues 
#'                             compared to largest one.
#' @param toleranceConvergence numeric indicating convergence tolerance for 
#'                              alternating projection algorithm.
#' @param      tolerancePosDef tolerance for forcing positive definiteness (in the 
#'                              final step) of alternating projection algorithm.
#' @param        maxIterations integer giving maximum number of iterations 
#'                             allowed in alternating projection algorithm; when
#'                             this number is exceeded without convergence                              
#'                             a warning is displayed and matrix computed in step
#'                             \code{maxIterations} of algorithm is returned.
#' @param              doEigen logical indicating whether finding the closest positive
#'                             definite matrix -through a eigen step- should be applied 
#'                             to the result of the alternating projection algorithm.
#'  
#' @export
#'    
#' @details This function is based on an alternating projection algorithm which 
#' involves projecting sequentially and iteratively the initial matrix into the 
#' set of symmetric positive definite and into the space of Toeplitz matrices, 
#' respectively. The iteration process will stop because either a criterion
#' of convergence is met or \code{maxIterations} has been exceeded (without convergence).
#' Criterion of convergence: if the Frobenius norm of the difference of the
#' projection matrices computed in the last two iterations of the algorithm
#' is smaller than \code{toleranceConvergence}, then the algorithm stops returning
#' the projection matrix computed in the last iteration. 
#'  
#' When projecting onto the set of symmetric positive definite matrices, \code{toleranceEigen}
#' controls the relative magnitude of any eigenvalue \eqn{\lambda_k}{\lambda_k} with
#' respect to the largest one \eqn{\lambda_1}{\lambda_1} and all eigenvalues \eqn{\lambda_k}{\lambda_k}
#' are treated as zero if 
#' \eqn{\lambda_k / \lambda_1 \leq \code{toleranceEigen}}{\lambda_k / \lambda_1 \le toleranceEigen}.
#'  
#' @examples 
#' # Higham (2002), p. 334
#' (mat <- matrix(c(1, 1, 0, 1, 1, 1, 0, 1, 1), byrow = TRUE, ncol = 3))
#' matProj <- matrix(c(1, 0.7607, 0.1573, 0.7607, 1, 0.7607, 0.1573, 0.7607, 1), 
#'                   byrow = TRUE, ncol = 3)
#' nrPDT.mat <- nearPDToeplitz(mat, type = "correlation")
#' stopifnot( identical(unname(matProj), unname(round(as.matrix(nrPDT.mat$projection), 
#'                      digits=4) ) ))
#' eigen(nrPDT.mat$projection)$values
#' 
#' # Toeplitz banded matrix near to the covariance matrix of 100 realizations
#' # of an MA(5) with following parameters:
#'
#' n <- 1e2
#' alphas <- c(-2, 0.5, -4, 0, 0.75)
#' (true.acf <- ARMAacf(ma = alphas))
#' alphasMat <- symBandedToeplitz(true.acf, n = n)
#' stopifnot( min(eigen(alphasMat)$values) > 0 ) # alphasMat is a positive definite matrix
#'
#' (l <- length(true.acf))
#' (acf.modified <- c(true.acf[-c(l - 1, l)], 0.25)) # modifying original acf
#' x <- acf.modified
#' acfMat <- symBandedToeplitz(x, n = n)
#' 
#' # no. of non positive eigenvalues of acfMat (6)
#' length( eigen(acfMat)$values[eigen(acfMat)$values < 0 ] )
#' # acfMat is a 100 x 100 symmetric banded Toeplitz matrix
#' acfMat[1:15, 1:30]
#'
#' system.time(nrPDT.acfMat <- nearPDToeplitz(acfMat, type = "correlation"))
#' y <- eigen(nrPDT.acfMat$projection)$values
#' # no. of non positive eigenvalues of nrPDT.acfMat
#' length( y[ y < 0 ] ) # none!
#'
#' @seealso \code{\link[Matrix]{nearPD}}, \code{\link[dbacf]{projectToeplitz}}, 
#'          \code{\link[dbacf]{symBandedToeplitz}}, \code{\link[sfsmisc]{posdefify}}
#'           
#' @return A list containing:
#' \item{projection}{a matrix, the computed symmetric positive definite Toeplitz matrix.}
#' \item{normF}{Frobenius norm of the difference between original and projection 
#'               matrix.}
#' \item{iterations}{number of iterations used for alternating projection algorithm.}
#' \item{relativeTolerance}{numeric giving relative error (in Frobenius norm) of 
#'                           final approximation with respect to original matrix.}
#' \item{converged}{logical indicating if alternating projection algorithm converged.}
#'  
#' @references Grigoriadis, K.M., Frazho, A., Skelton, R. (1994). 
#' \emph{Application of alternating convex projection methods for computation 
#' of positive Toeplitz matrices}, IEEE Transactions on signal processing 
#' \bold{42(7)}, 1873--1875.
#'  
#' @references N. Higham (2002). \emph{Computing the nearest correlation matrix - a
#' problem from finance}, Journal of Numerical Analysis \bold{22}, 329--343.
#' 
#' @references Tecuapetla-Gómez, I and Munk, A. (2017). \emph{Autocovariance
#' estimation in regression with a discontinuous signal and \eqn{m}-dependent errors: A 
#' difference-based approach}. Scandinavian Journal of Statistics, \bold{44(2)}, 346--368.
#' 
nearPDToeplitz <- function(matrix, type = c("covariance", "correlation"),
                           toleranceEigen = 1e-06, toleranceConvergence = 1e-06,
                           tolerancePosDef = 1e-06, maxIterations = 1e02,
                           doEigen = TRUE) {

  if( !isSymmetric(unname(matrix)) )
    stop("argument 'matrix' must be a symmetric matrix")
  
  DykstraCorrection <- matrix
  DykstraCorrection[] <- 0
  
  X <- matrix
  
  iter <- 0
  converged <- FALSE
  conv <- Inf
  
  type <- match.arg(type)
  
  while (iter < maxIterations && !converged) {
    Y <- X
    
    increment <- Y - DykstraCorrection
    eigenDecomposition <- eigen(increment, symmetric = TRUE)
    
    eigenVectors <- eigenDecomposition$vectors 
    eigenValues <- eigenDecomposition$values
    positiveEigenvalues <- eigenValues > toleranceEigen * eigenValues[1]
    if ( !any(positiveEigenvalues) ) 
      stop("Matrix seems negative semi-definite")
    eigenVectors <- eigenVectors[, positiveEigenvalues, drop = FALSE]
    X <- tcrossprod(eigenVectors * rep(eigenValues[positiveEigenvalues], 
                                       each = nrow(eigenVectors)), eigenVectors)
    
    DykstraCorrection <- X - increment
    X <- (X + t(X)) / 2
    
    if (type == "correlation") 
      diag(X) <- 1
    
    X <- projectToeplitz(X)
    
    conv <- Matrix::norm(Y - X, "F") / Matrix::norm(Y, "F")
    iter <- iter + 1
    
    converged <- (conv <= toleranceConvergence)
  }
  if (!converged) 
    warning(gettextf("'nearPDToeplitz()' did not converge in %d iterations", 
                     iter), domain = NA)
  
  n <- ncol(matrix)
  if (doEigen) {
    eigenDecomposition <- eigen(X, symmetric = TRUE)
    eigenValues <- eigenDecomposition$values
    epsilon <- tolerancePosDef * abs(eigenValues[1])
    if (eigenValues[n] < epsilon) {
      eigenValues[eigenValues < epsilon] <- epsilon
      
      eigenVectors <- eigenDecomposition$vectors
      diagX <- Matrix::diag(X)
      X <- eigenVectors %*% (eigenValues * t(eigenVectors))
      D <- sqrt(pmax(epsilon, diagX) / Matrix::diag(X))
      X[] <- D * X * rep(D, each = n)  
    }
    if ( type == "correlation" ) 
      diag(X) <- 1      
  }
  
  list(projection = X, normF = Matrix::norm(matrix - X, "F"), iterations = iter, 
       relativeTolerance = conv, converged = converged)      
}
