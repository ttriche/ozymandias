#' Pull a 27k dataset from GEO.
#' 
#' Doesn't matter what the annotations are or any other bullshit. Just get it.
#' 
#' @param GSE   The dataset.
#' @param IDATs   Reanalyze from IDAT files if possible? (TRUE) 
#' @param ...     Additional arguments to be passed to prepForReanalysis().
#'
#' @return a GenomicRatioSet.
#'
#' @export
fetch27k <- function(GSE, IDATs=TRUE, ...) {

  gset <- getGEO(GSE)[[1]]
  library(FDb.InfiniumMethylation.hg19)
  hm27 <- get27k()
  gset <- gset[intersect(featureNames(gset), names(hm27)), ] 
  res <- GenomicRatioSet(gr=hm27[featureNames(gset)], 
                         pData=pData(gset),
                         Beta=exprs(gset))
  if (IDATs) {
    message("Attempting to reanalyze from raw IDAT files...")
    res <- prepForReanalysis(res, ask=!IDATs, ...)
  }
  return(res)

}
