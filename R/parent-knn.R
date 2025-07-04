#' Find the Nearest Neighbours Among Each Location's Parents
#'
#' For each row in \code{x}, find the k-nearest neighbours among the earlier
#' rows.
#'
#' @param x A matrix or vector in which each row/element represents a point.
#' @param n_parents The number of parent nearest neighbours to find.
#' @param leaf_size Leaf size for the KD-tree algorithm used in nearest
#' neighbour search; the algorithm may run faster for different values.
#'
#' @return A matrix where the ith row contains the indices of the parent nearest
#' neighbours of the ith input point. The first column by convention is equal to
#' i, so the matrix has `n_parents + 1` columns. Some entries can be `NA`.
#'
#' @examples
#' x <- matrix(rnorm(100), nrow = 50)
#' parents <- parent_knn(x, 5)
#'
#' @export
parent_knn <- function(x, n_parents, leaf_size = 40) {
  if (is.null(dim(x))) {
    x <- matrix(x, nrow = length(x))
  }
  .parent_knn(x, n_parents, leaf_size)
}
