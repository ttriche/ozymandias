#' Take the result of IDATsToRgSet and preprocess it (pvals, SNPs, etc.)
#' 
#' @param rgSet   an RGChannelSet 
#' @param pcutoff p-value cutoff to set a beta value to NA (0.01)
#' @param genome  genome build to annotate against if not using funnorm (hg19)
#' @param funnorm run preprocessFunnnorm() to get the grSet? (FALSE, use noob)
#' @param gr      an optional GenomicRanges to annotate the rows of the result
#'
#' @return  a GenomicRatioSet with metadata(grSet)$SNPs and p > 0.01 set to NA
#' 
#' @export
rgSetToGrSet <- function(rgSet, pcutoff=0.01, genome=c("hg19","hg38"),
                         funnorm=FALSE, gr=NULL) {

  if (funnorm) {
    grSet <- preprocessFunnorm(rgSet)
  } else { 
    # needs monkeypatching
    mSet <- ifelse(annotation(rgSet)["array"] == "IlluminaHumanMethylation27k",
                   preprocessNoob27k(rgSet),
                   preprocessNoob(rgSet))
    genome <- match.arg(genome)
    if (is.null(gr)) {
      library(paste0("FDb.InfiniumMethylation.", genome), character.only=TRUE)
      gr <- switch(annotation(rgSet)["array"],
                   "IlluminaHumanMethylation27k"=get27k(),
                   "IlluminaHumanMethylation450k"=get450k(),
                   "IlluminaHumanMethylationEPIC"=getEPIC())
    }
    grSet <- GenomicRatioSet(gr=gr[rownames(mset)], Beta=getBeta(mset), 
                             preprocessMethod=preprocessMethod(mset), 
                             CN=getCN(mset), pData=pData(mset),
                             annotation=annotation(mset))
  }

  metadata(grSet)$SNPs <- getSnpBeta(rgSet)
  assays(grSet)$pval <- matrix(NA_real_, ncol=ncol(grSet), nrow=nrow(grSet))
  assays(grSet)$pval <- detectionP(rgSet)[rownames(grSet), colnames(grSet)]
  is.na(assays(grSet)$Beta[which(assays(grSet)$pval > pcutoff)]) <- TRUE 
  return(grSet)

}
