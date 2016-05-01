#' This might be useful for _all_ SummarizedExperiments; split x on mcols[[y]]
#' If the column 'y' isn't in mcols(x), it will be added by addCytoInfo(x, y)
#' 
#' @param x   a SummarizedExperiment-derived object or something like it
#' @param y   what kind of information to split on (chr, arm, band)
#'
#' @return    x, split by y
#'
#' @import SummarizedExperiment
#' 
#' @export
splitBy <- function(x, y=c("chr","chromosome", "arm", "band","cytoband")) {
  y <- match.arg(y) 
  y <- sub("chromosome", "chr", y)
  y <- sub("cytoband", "band", y)
  x <- addCytoInfo(x, y)
  split(x, mcols(x)[, y])
}

#' @describeIn splitBy split x by chromosome 
#' @export
byChr <- function(x) splitBy(x, "chr")

#' @describeIn splitBy split x by chromosome arm
#' @export
byArm <- function(x) splitBy(x, "arm")

#' @describeIn splitBy split x by cytoband 
#' @export
byCytoband <- function(x) splitBy(x, "band")
