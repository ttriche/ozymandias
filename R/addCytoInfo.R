#' omnibus sub-chromosome splitting information (chr, arm, and band)
#'
#' @param x   a SummarizedExperiment-derived object or something like it
#' @param y   what kind of information to add (chr, arm, band)
#'
#' @return    x, but with y added to mcols(x)
#'
#' @import SummarizedExperiment
#' 
#' @export
addCytoInfo <- function(x, y=c("chr", "arm", "band")) {
  y <- match.arg(y) 
  if (y == "chr") mcols(x)$chr <- seqnames(x)
  if (!y %in% names(mcols(x))) {
    assembly <- unique(genome(x)) # need correct database
    if (assembly != "hg19") { # only hg19 is supported for now 
      message("Your ", class(x), " has ", paste(assembly, collapse=", "), ",")
      stop("but only hg19 is currently supported. Cannot proceed.")
    }
    if (!exists("hg19.cytobands")) {
      if (require("ozymandias")) data(hg19.cytobands, package="ozymandias")
      else load("hg19.cytobands.rda") 
    }
    olaps <- findOverlaps(x, hg19.cytobands) # refactor, use build.info
    mcols(x)[, y] <- mcols(hg19.cytobands)[subjectHits(olaps), y]
  }
  return(x)
}
