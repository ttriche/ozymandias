#' as the name says, if there are IDATs attached to a grSet, rerun them.
#'
#' @param grSet   a grSet with $supplementary_file and $supplementary_file.1
#' @param name    the name of the resulting rgSet and grSet files
#'
#' @return        a reanalyzed and funnorm()'ed grSet
#' 
#' @export
#'
reanalyzeFromIDATs <- function(grSet, name) {
  grSet <- prepForReanalysis(grSet)
  if (!is.null(grSet$Basename)) {
    rgSetFile <- paste0(name, "_rgSet.rds")
    rgSet <- geoIDATsToRgSet(grSet) # don't bother ungzipping
    message(ncol(rgSet), " samples processed sucessfully.")
    message("Saving raw intensities to ", rgSetFile, "...")
    saveRDS(rgSet, rgSetFile)
    newGrSet <- processMeth(rgSet, name=name)
  } else { 
    message("Failed to prepare grSet for reanalysis!")
  }
  return(newGrSet)
}
