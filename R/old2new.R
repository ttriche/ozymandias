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
  if (class(res) != "try-error") res <- try(grSet@elementMetadata, silent=TRUE)
  if (class(res) == "try-error") {
    cdat <- colData(grSet)
    Beta <- grSet@assays[["Beta"]]
    rdat <- try(grSet@rowData, silent=TRUE)
    if (class(rdat) == "try-error") rdat <- grSet@rowRanges
    anno <- grSet@annotation
    grSet <- GenomicRatioSet(gr=rdat, Beta=Beta, pData=cdat, annotation=anno)
  }
  return(grSet)
}
