#' as the title suggests.  Original author Zhang (cf. Genetic Epi, 2016)
#' 
#' @import mgcv 
#' 
#' @export
createPredictors <- function(covariates=NULL, funcs, kz=30, kb=30,
                             smooth.cov=FALSE) {
  kb <- min(kz, kb)
  p <- ifelse(is.null(covariates), 0, dim(covariates)[2])
  if (is.matrix(funcs)) {
    n <- nrow(funcs)
    Funcs <- list(length = 1)
    Funcs[[1]] <- funcs
  } else {
    n <- nrow(funcs[[1]])
    Funcs <- funcs
  }
  N.Pred <- length(Funcs)
  t <- phi <- psi <- CJ <- list(length = N.Pred)

  # surely this can be parallelized?
  for (i in 1:N.Pred) {
    t[[i]] <- seq(0, 1, length = dim(Funcs[[i]])[2])
    N_obs <- length(t[[i]])
    meanFunc <- apply(Funcs[[i]], 2, mean, na.rm = TRUE)
    resd <- sapply(1:length(t[[i]]), function(u) Funcs[[i]][, u] - meanFunc[u])
    Funcs[[i]] <- resd
    num <- kb - 2
    qtiles <- seq(0, 1, length = num + 2)[-c(1, num + 2)]
    knots <- quantile(t[[i]], qtiles)
    phi[[i]] <- cbind(1, t[[i]], 
                      sapply(knots, function(k) 
                             ((t[[i]] - k > 0) * (t[[i]] - k))))
    CJ[[i]] <- Funcs[[i]] %*% phi[[i]]
  }
  X <- cbind(rep(1, n), covariates)
  for (i in 1:N.Pred) {
    X <- cbind(X, CJ[[i]])
  }
  return(X)
}
