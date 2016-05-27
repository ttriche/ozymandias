#' to update old grSets so I don't have to reload, e.g. cnAML
#' note: this can be a pig if the dataset is huge; in that case, just reload it
#'
#' @param grSet   the old grSet
#'
#' @return        a new, valid grSet 
#' 
#' @import minfi
#'
#' @export
old2new <- function(grSet) {
  res <- try(granges(grSet), silent=TRUE)
  if (class(res) != "try-error") {
    res <- try(grSet@elementMetadata, silent=TRUE)
  }
  if (class(res) == "try-error") {
    cdat <- colData(grSet)
    anno <- grSet@annotation
    rdat <- try(grSet@rowData, silent=TRUE)
    if (class(rdat) == "try-error") rdat <- grSet@rowRanges
    for (name in c("Beta", "M", "CN")) {
      if (name %in% names(grSet@assays)) {
        dat <- grSet@assays[[name]]
        colnames(dat) <- rownames(cdat)
        assign(name, dat)
        rm(dat) 
      } else {
        assign(name, NULL)
      }
    }
    grSet <- GenomicRatioSet(gr=rdat, 
                             Beta=Beta, 
                             pData=cdat, 
                             CN=CN,
                             M=M,
                             annotation=anno)
  }
  return(grSet)
}
