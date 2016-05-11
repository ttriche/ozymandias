#' extract IDAT names, etc. for reprocessing from raw data
#' 
#' @param res   a SummarizedExperiment, GenomicRatioSet, or suchlike
#' @param ask   ask before downloading the supplemental data? (TRUE) 
#' @param ...   other arguments to (potentially) pass to getGEOSuppFiles()
#' 
#' @return      the same thing but with tidied-up Basename, Sample_Name, etc.
#' 
#' @export
prepForReanalysis <- function(res, ask=TRUE, ...) {

  if (!is.null(gse) && !is.null(res$supplementary_file)) {
    res$Sample_Name <- res$title
    res$Sample_Group <- "Experimental"
    res$Basename <- gsub(".gz","", 
                         gsub("_Grn.idat","", 
                              basename(res$supplementary_file)))
    if (ask) { 
      dl <- toupper(substr(readline("Download IDAT files? [Y/N] "),1,1))
      if (dl == "N") {
        return(res)
      }
    } 
    lapply(res$geo_accession, getGEOSuppFiles, makeDirectory=FALSE, ...)
    res <- reanalyzeFromIDATs(res, "reanalyzed")
    message("Your reanalyzed rgSet and GenomicRatioSet have been stored in")
    message("reanalyzed_rgSet.rds and reanalyzed_450k.rds respectively.")
    message("(If they were actually EPIC files, well, sorry about that.)")
  }
  return(res)

}
