#' For array-based DNA methylation, particularly summarized across regions,
#' we can do better (a lot better) than MAD.  Since we know that 
#'
#' max(SD(X_j)) if X_j ~ Beta(a, b) < max(SD(X_j)) if X_j ~ Bernoulli(a/(a+b))
#'
#' for X having a known mean and SD, hence solvable for a + b by MoM, define
#'
#' extremality = sd(X_j) / bernoulliSD(mean(X_j))
#'
#' This function selects the k most extremal rows of x and returns their values.
#'
#' @param     x   a matrix of beta values
#' @param     k   how many rows to return (500)
#' 
#' @return    a smaller matrix (usually) than x, ranked by extremality 
#' 
#' @import    matrixStats
#' 
#' @export
byExtremality <- function(x, k=500) {
  means <- rowMeans(x, na.rm=TRUE)
  bernoulliSd <- sqrt(means * (1 - means))
  actualSd <- rowSds(x, na.rm=TRUE)
  k <- min(nrow(x), k)
  x[rev(order(actualSd / bernoulliSd))[seq_len(k)], ]
}
