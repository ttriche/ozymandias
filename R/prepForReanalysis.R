#' extract IDAT names, etc. for reprocessing from raw data
#' 
#' @param res       a SummarizedExperiment, GenomicRatioSet, or suchlike
#' @param ask       ask before downloading the supplemental data? (TRUE) 
#' @param ...       other arguments to (potentially) pass to rgSetToGrSet()
#' 
#' @return      the same thing but with tidied-up Basename, Sample_Name, etc.
#' 
#' @export
prepForReanalysis <- function(res, ask=TRUE, ...) {

  if (!is.null(res$supplementary_file)) {
    if (ask) { 
      dl <- toupper(substr(readline("Download IDAT files? [Y/N] "),1,1))
      if (toupper(dl) == "N") {
        return(res)
      }
    } 
    message("Attempting to fetch IDAT files for each sample...") 
    lapply(res$geo_accession, getGEOSuppFiles, makeDirectory=FALSE)
    message("Reprocessing IDATs to an RGChannelSet...")
    rgSet <- geoIDATsToRgSet(res)
    message("Reprocessing RGChannelSet to a GenomicRatioSet...")
    res <- rgSetToGrSet(rgSet, ...)
    message("Done.") 
  }
  return(res)

}
