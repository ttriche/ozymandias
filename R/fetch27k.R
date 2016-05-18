#' Pull a 27k dataset from GEO.
#' 
#' Doesn't matter what the annotations are or any other bullshit. Just get it.
#' 
#' @param GSE   The dataset.
#' @param ...     Additional arguments to be passed to prepForReanalysis().
#'
#' @return a GenomicRatioSet.
#'
#' @export
fetch27k <- function(GSE, ...) {

  gset <- getGEO(GSE)[[1]]
  library(FDb.InfiniumMethylation.hg19)
  hm27 <- get27k()
  gset <- gset[intersect(featureNames(gset), names(hm27)), ] 
  res <- GenomicRatioSet(gr=hm27[featureNames(gset)], 
                         pData=pData(gset),
                         Beta=exprs(gset))
  res <- prepForReanalysis(res, ...)
  return(res)

}
