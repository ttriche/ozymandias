#' Pull a 450k dataset from GEO.
#' 
#' Doesn't matter what the annotations are or any other bullshit. Just get it.
#' 
#' @param GSE     The dataset.
#' @param ...     Additional arguments to be passed to prepForReanalysis().
#'
#' @return a GenomicRatioSet.
#'
#' @export
fetch450k <- function(GSE, ...) {
  gset <- getGEO(GSE)[[1]]
  library(FDb.InfiniumMethylation.hg19)
  hm450 <- get450k()
  res <- GenomicRatioSet(gr=hm450[featureNames(gset)], 
                         pData=pData(gset), 
                         Beta=exprs(gset))
  res <- prepForReanalysis(res, ...)
  return(res)
}
