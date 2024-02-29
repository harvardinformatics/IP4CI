#' @title  jaccard
#' @description calculate jaccard between two sets
#' @param a first set
#' @param b second set
#' @examples
#'jaccard(a,b)
#'@export
###############################################################################
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}
