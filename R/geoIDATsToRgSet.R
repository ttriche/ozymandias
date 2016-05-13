#' read in gzipped IDAT files corresponding to a GEO series with IDAT files 
#'
#' @param res       something with $supplementary_file and $supplementary_file.1
#' @param parallel  try to run multicore (default: FALSE)
#'
#' @return          an rgSet 
#' 
#' @export
#'
geoIDATsToRgSet <- function(res, parallel=FALSE) { 
  
  if (is(res, "SummarizedExperiment")) res <- colData(res) 
  covs <- names(res)

  cols <- c(Grn="supplementary_file", Red="supplementary_file.1")
  message("Checking for ", paste(paste0("$", cols), collapse=" and "), "...")
  if (!all(cols %in% covs)) stop("Supplementary files from GEO not found.")
  if (!"title" %in% covs) message("Sample titles from GEO not found.")

  if (all(grepl("Red", res$supplementary_file))) 
    cols["Red"] <- "supplementary_file"
  if (all(grepl("Grn", res$supplementary_file.1))) 
    cols["Grn"] <- "supplementary_file.1"
  if(!all(file.exists(res[, cols["Grn"]]))) 
    stop("Could not find all Grn channel IDATs.")
  if(!all(file.exists(res[, cols["Red"]]))) 
    stop("Could not find all Red channel IDATs.")

  out <- IDATsToRgSet(res, 
                      idatGrn=cols["Grn"], 
                      idatRed=cols["Red"], 
                      parallel=parallel)

  if ("title" %in% names(res)) sampleNames(out) <- out$title
  return(out) 

}
