#' plot SNP probes from DNA methylation arrays to avoid label swaps
#' 
#' @param x           SummarizedExperiment-like object with DNA methylation data
#' @param rotate      rotate the plot? (default is FALSE) 
#' 
#' @return            output from heatmap(...)
#'
plotSNPs <- function(x, rotate=FALSE, ...) {

  name <- as.character(match.call()["x"])
  if (class(x) == "MultiAssayExperiment") {
    tmp <- perSampleMetadata(x)$SNPs 
  } else if (is(x, "matrix")) {
    tmp <- x  
  } else {
    tmp <- metadata(x)$SNPs[, colnames(x)] 
  }

  if (nrow(tmp) < 2) {
    stop("Need something with SNP (rsXX) features... none found")
  }


  # fit w/mclust
  SNPfit <- Mclust(as.numeric(tmp), G=3)
  calls <- matrix(SNPfit$class - 1, ncol=ncol(tmp))
  rownames(calls) <- rownames(tmp) 
  colnames(calls) <- colnames(tmp) 

  # plot
  par(mfrow=c(1,1)) 
  if(rotate) calls <- t(calls)
  SNP <- c("blue","yellow","red")
  heading <- paste("SNPs for", ncol(x) , "samples in", name)
  heatmap(calls, col=SNP, scale="none", main=heading)

}
