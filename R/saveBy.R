#' save-by-subset functions for huge (n > 1000) datasets for e.g. corpcor(x)
#'
#' @param x     a SummarizedExperiment-derived object or something like it
#' @param y     how to split x before saving it into subsets via mcols(x)$y
#' @param name  a prefix for the saved RDS files, becomes "name_subset.rds"
#'
#' @export
saveBy <- function(x, y, name) {
  x <- addCytoInfo(x, y) # FIXME? 
  for (y_i in unique(mcols(x)[[y]])) {
    filename <- paste0(name, "_", y_i, ".rds")
    message("Saving ", y_i, " to ", filename, "...", appendLF=FALSE)
    saveRDS(splitBy(x, y)[[y_i]], file=filename)
    message("done.")
  }
}

#' @describeIn saveBy save x by chromosome
#' @export
saveByChr <- function(x, name) saveBy(x, y="chr", name=name)

#' @describeIn saveBy save x by chromosome arm 
#' @export
saveByArm <- function(x, name) saveBy(x, y="arm", name=name)

#' @describeIn saveBy save x by cytoband 
#' @export
saveByCytoband <- function(x, name) saveBy(x, y="band", name=name)
