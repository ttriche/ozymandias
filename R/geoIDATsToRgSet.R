#' short-circuit minfi's irritating insistence on un-gzipped IDAT files
#'
#' @param res     something with $supplementary_file and $supplementary_file.1
#'
#' @return        an rgSet 
#'
#' @import illuminaio
#' 
#' @export
#'
geoIDATsToRgSet <- function(res) { 
   
  covs <- names(colData(res)) 
  message("Checking for $supplementary_file and $supplementary_file.1 ...")

  if (!c("supplementary_file", "supplementary_file.1" ) %in% covs) {
    stop("Supplementary files from GEO not found.")
  } else if (!"title" %in% covs) {
    stop("Sample titles from GEO not found.")
  }

  if (all(grepl("Grn", res$supplementary_file))) {
    G.files <- basename(res$supplementary_file)
    R.files <- basename(res$supplementary_file.1)
  } else if (all(grepl("Grn", res$supplementary_file.1))) {
    G.files <- basename(res$supplementary_file.1)
    R.files <- basename(res$supplementary_file)
  } else { 
    message("Cannot determine which supplementary files are Red channel")
    stop("and which are Grn channel. Cannot proceed.")
  }

  if(!all(file.exists(G.files))) stop("Could not find all Grn channel IDATs.")
  if(!all(file.exists(R.files))) stop("Could not find all Red channel IDATs.")

  names(G.files) <- names(R.files) <- res$title 
  getMeans <- function(x) readIDAT(x)$Quants[, "Mean"]
  G.signal <- do.call(cbind, lapply(G.files, getMeans))
  R.signal <- do.call(cbind, lapply(R.files, getMeans))
  out <- new("RGChannelSet", Red=R.signal, Green=G.signal)

  anIDAT <- readIDAT(G.files[1])
  featureNames(out) <- rownames(anIDAT$Quants)
  annot <- c(array=switch(anIDAT$ChipType,
    "BeadChip 12x1"="IlluminaHumanMethylation27k",
    "BeadChip 12x8"="IlluminaHumanMethylation450k",
    "BeadChip 8x5"="IlluminaHumanMethylationEPIC"))
  annot["annotation"] <- switch(annot["array"],
    "IlluminaHumanMethylation27k"="",
    "IlluminaHumanMethylation450k"=minfi:::.default.450k.annotation,
    "IlluminaHumanMethylationEPIC"=minfi:::.default.epic.annotation)
  annotation(out) <- annot
  return(out)  

}
