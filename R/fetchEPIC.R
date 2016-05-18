#' Pull an EPIC dataset from GEO.
#' 
#' Doesn't matter what the annotations are or any other bullshit. Just get it.
#' 
#' @param GSE     The dataset.
#' @param ...     Additional arguments to be passed to prepForReanalysis().
#'
#' @return a GenomicRatioSet.
#'
#' @export
fetchEPIC <- function(GSE, ...) {
  gset <- getGEO(GSE)[[1]]
  # FIXME: add hg38 support here 
  library(FDb.InfiniumMethylation.hg19)
  hmEPIC <- getEPIC()
  res <- GenomicRatioSet(gr=hmEPIC[featureNames(gset)], 
                         pData=pData(gset), 
                         Beta=exprs(gset))
  res <- prepForReanalysis(res, ...)
  return(res)
}
