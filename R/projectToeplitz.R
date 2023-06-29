#' Projection onto the set of symmetric Toeplitz matrices
#' 
#' Computes the orthogonal projection onto the space of symmetric Toeplitz 
#' matrices as given in \cite{Grigoriadis et al. (1994)}.
#'                        
#' @param matrix a symmetric matrix. 
#'  
#' @export
#'  
#' @examples  
#' A <- matrix(c(2, 1, 1, 1, 2, 0, 1, 0, 0), byrow = 3, nrow = 3)
#' projectToeplitz(A)
#'  
#' @seealso \code{\link{nearPDToeplitz}}
#'  
#' @return The computed projection matrix.
#'          
#' @references Grigoriadis, K.M., Frazho, A., Skelton, R. (1994). 
#' \emph{Application of alternating convex projection methods for computation 
#' of positive Toeplitz matrices}, IEEE Transactions on signal processing 
#' \bold{42(7)}, 1873--1875
#' 
projectToeplitz <- function(matrix) {  
  if( !isSymmetric(unname(matrix)) )
    stop("argument 'matrix' must be a symmetric matrix")
    
  .projectToeplitz(matrix=matrix)
}

.projectToeplitz <- function(matrix) {  
  firstRow <- computeFirstRow(matrix = matrix)
  bLis <- as.data.frame(matrix(firstRow, nrow =  nrow(matrix),
                               ncol = length(firstRow), byrow = TRUE))
  as.matrix(Matrix::bandSparse(ncol(matrix), k = 0:(length(firstRow) - 1),
                               diagonals = bLis, symmetric = TRUE))
}
