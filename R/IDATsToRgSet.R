#' Function to sidestep minfi's irritating insistence on un-gzipped IDAT files
#'
#' @param res       a data.frame with c(idatGrn, idatRed) OR Basename column(s)
#' @param idatGrn   the column with the green IDAT filenames, if any
#' @param idatRed   the column with the red IDAT filenames, if any 
#' @param parallel  try to run in parallel? (default: FALSE) 
#'
#' @return        an rgSet 
#'
#' @import        illuminaio
#' @import        minfi
#' 
#' @export
#'
IDATsToRgSet <- function(res, idatGrn=NULL, idatRed=NULL, parallel=FALSE) { 
   
  message("Checking for $",idatGrn," and $",idatRed," columns...", appendLF=F)

  if (is(res, "SummarizedExperiment")) res <- colData(res)
  covs <- names(res)
  idatCols <- c(idatGrn, idatRed)
  channels <- c("Grn", "Red") 
  names(channels) <- channels
  names(idatCols) <- channels 

  checkgzip <- function(f) ifelse(file.exists(f), f, paste0(f, ".gz"))
  getIDAT <- function(barcode, ch) checkgzip(paste0(barcode, "_", ch, ".idat"))
  if (!all(idatCols %in% covs)) {
    message("Checking for $Basename column...", appendLF=F)
    if (!"Basename" %in% covs) {
      stop("Missing. Cannot determine IDAT files from data, aborting load.")
    } else { 
      message("Found. Reconstructing per-channel IDAT filenames...", appendLF=F)
      for (ch in channels) res[, idatCols[ch]] <- getIDAT(res$Basename, ch)
    }
  } else {
    message("...found.")
  }
  if (!"Basename" %in% covs) {
    message("Adding Basename column...", appendLF=F)
    res[,"Basename"] <- sub("_(Red|Grn)\\.idat(\\.gz)?", "", res[, idatCols[1]])
    message("...done.") 
  }

  lost <- do.call(c, 
                  lapply(idatCols, 
                         function(i) setdiff(basename(res[, i]), list.files())))
  if (length(lost) > 0) stop("Missing files: ", lost)

  getMeans <- function(IDAT) readIDAT(IDAT)$Quants[, "Mean"]

  if (parallel) {
    signal <- mclapply(idatCols, 
                       function(i) 
                         do.call(cbind, lapply(basename(res[, i]), getMeans)))
  } else { 
    signal <- lapply(idatCols, 
                     function(i) 
                       do.call(cbind, lapply(basename(res[, i]), getMeans)))
  }
  out <- new("RGChannelSet", 
             Red=signal[["Red"]], 
             Green=signal[["Grn"]])
  sampleNames(out) <- res[,"Basename"]
  anIDAT <- readIDAT(basename(res[1, idatCols["Grn"]]))
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
