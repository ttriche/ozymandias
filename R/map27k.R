#' wrapper for preprocessNoob27k
#' 
#' @param rgSet     a humanMethylation27 rgSet
#' @param genome    a genome assembly (currently only hg19 works properly)
#' @param ...       additional arguments to pass to 
#'
#' @return          a GenomicRatioSet
#'
#' @export
map27k <- function(mset, genome="hg19", ...) { 
  pkg <- paste0("FDb.InfiniumMethylation.", genome)
  message("Loading annotations from ", pkg, "...")
  library(pkg, character.only=TRUE)
  gr <- get27k()
  rset <- ratioConvert(mset[which(rownames(mset) %in% names(gr)), ])
  # essentially a poor man's mapToGenome until I add hm27 annotations
  GenomicRatioSet(gr=gr[rownames(rset)], Beta=getBeta(rset), 
                  preprocessMethod=preprocessMethod(rset),
                  CN=getCN(rset), pData=pData(rset),
                  annotation=annotation(rset))
}
