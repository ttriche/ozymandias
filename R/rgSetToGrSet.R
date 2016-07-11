#' Take the result of IDATsToRgSet and preprocess it (pvals, SNPs, etc.)
#' 
#' @param rgSet     an RGChannelSet 
#' @param pcutoff   p-value cutoff to set a beta value to NA (0.01)
#' @param funnorm   run preprocessFunnnorm() to get the grSet? (FALSE, use noob)
#' @param dyeMethod How to fix dye bias: reference or single sample approach?
#' @param genome    genome build to annotate against if not using funnorm (hg19)
#' @param gr        an optional GenomicRanges to annotate the rows of the result
#'
#' @return  a GenomicRatioSet with metadata(grSet)$SNPs and p > 0.01 set to NA
#' 
#' @export
rgSetToGrSet <- function(rgSet, pcutoff=0.01, funnorm=FALSE, dyeMethod="single",
                         genome=c("hg19","hg38","hg18"), gr=NULL, ...) {

  platform <- annotation(rgSet)["array"]
  genome <- match.arg(genome)
  if (genome != "hg19") stop("Genomes besides hg19 are currently unsupported.")
  process <- ifelse(funnorm, "funnorm", 
                    ifelse(platform == "IlluminaHumanMethylation27k", "noob27k",
                           "noob")) # default to single-sample noob

  dm <- dyeMethod
  grSet <- switch(process,
                  funnorm=preprocessFunnorm(rgSet, ...),
                  noob27k=map27k(preprocessNoob27k(rgSet), genome=genome, ...),
                  noob=ratioConvert(mapToGenome(preprocessNoob(rgSet, ..., 
                                                               dyeMethod=dm))))
  metadata(grSet)$SNPs <- getSnpBeta(rgSet)
  assays(grSet)$pval <- matrix(NA_real_, ncol=ncol(grSet), nrow=nrow(grSet))
  assays(grSet)$pval <- detectionP(rgSet)[rownames(grSet), colnames(grSet)]
  is.na(assays(grSet)$Beta[which(assays(grSet)$pval > pcutoff)]) <- TRUE 
  return(grSet)
}
