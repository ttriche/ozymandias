#' Parse a TCGA or TARGET barcode and translate the codebook entries.
#' Note that tissue source site is inconsistent between TCGA and TARGET
#' and is discarded by default to avoid confusion.
#'
#' @param   barcode     a TCGA or TARGET barcode for a sample/specimen
#' @param   sourcesite  retain and translate the tissue Source.Site? (FALSE)
#'
#' @return            a vector of named fragments
#'
#' @export
parseBarcode <- function(barcode, sourcesite=FALSE) { 
  parsed <- strsplit(strsplit(barcode, ".", fixed=TRUE)[[1]][1], 
                     "-", fixed=TRUE)[[1]][1:5]
  names(parsed) <- c("Project", "TSS", "Participant", 
                     "SampleVial", "PortionAnalyte")

  parsed["Tissue"] <- .getTissue(substr(parsed["SampleVial"], 1, 2))
  parsed["Analyte"] <- .getAnalyte(substr(parsed["PortionAnalyte"], 3, 3))
  parsed["Source.Site"] <- .getSourceSite(parsed["TSS"])

  fragments <- c("Project", "Participant", "Tissue", "Analyte")
  if (sourcesite == TRUE) fragments <- c(fragments, "Source.Site")
  return(parsed[fragments])
}

.getTissue <- function(tissue) {
  data(TCGA.tissues, package="ozymandias")
  TCGA.tissues[tissue, "Short.Letter.Code"]
}

.getAnalyte <- function(analyte) {
  data(TCGA.analytes, package="ozymandias")
  TCGA.analytes[analyte, "Definition"]
}

.getSourceSite <- function(TSS) {
  data(TCGA.sourcesites, package="ozymandias")
  TCGA.sourcesites[TSS, "Source.Site"]
}
