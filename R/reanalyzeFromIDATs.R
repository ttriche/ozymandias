#' As the name says, if there are IDATs attached to a grSet, rerun them.
#'
#' @param grSet     a grSet with $supplementary_file and $supplementary_file.1
#' @param name      the name of the resulting rgSet and grSet files
#'
#' @return          a reanalyzed and funnorm()'ed grSet
#' 
#' @export
#'
reanalyzeFromIDATs <- function(grSet, name, ...) {

  if (!is.null(grSet$Basename)) {
    rgSetFile <- paste0(name, "_rgSet.rds")
    rgSet <- geoIDATsToRgSet(grSet) # don't bother ungzipping
    message(ncol(rgSet), " samples processed sucessfully.")
    message("Saving raw intensities to ", rgSetFile, "...")
    saveRDS(rgSet, rgSetFile)

    message("Checking to make sure annotations are installed...")
    manifest <- paste0(annotation(rgSet)["array"], "manifest")
    if(!require(manifest, character.only=TRUE)) {
      library(BiocInstaller)
      biocLite(manifest)
    }
    if(!require(manifest, character.only=TRUE)) {
      stop("Cannot install ", manifest)
    }
    annot <- paste(annotation(rgSet), collapse=".")
    if(!require(annot, character.only=TRUE)) {
      library(BiocInstaller)
      biocLite(annot)
    }
    return(rgSetToGrSet(rgSet, ...))
  } else { 
    message("Failed to prepare grSet for reanalysis!")
    return(grSet)
  }

}
