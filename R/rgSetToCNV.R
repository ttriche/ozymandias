#' Take the result of IDATsToRgSet, run noob, grind out CNV calls from conumee.
#' 
#' @param rgSet     an RGChannelSet 
#' @param cases     which samples to run, segment, and annotate?
#' @param controls  which samples are the control samples (for normalization)?
#' @param ...       optional arguments for conumee::CNV.create_anno function 
#'
#' @return          a list of CNV.analysis objects suitable for plotting etc. 
#'
#' @import conumee
#' 
#' @export
rgSetToCNV <- function(rgSet, cases, controls, ...) {

  if (annotation(rgSet)["array"] == "IlluminaHumanMethylation27k") {
    stop("Calling CNV/CNA on 27k is pointless... aborting.")
  } else {
    if (!exclude_regions %in% names(...)) data(exclude_regions,
                                               package="conumee")
    if (!detail_regions %in% names(...)) data(detail_regions, 
                                              package="ozymandias")
    if (!chrXY %in% names(...)) chrXY <- TRUE
    anno <- with(..., 
                 CNV.create_anno(chrXY=chrXY, 
                                 exclude_regions=exclude_regions,
                                 detail_regions=detail_regions))
  }

  cnvData <- CNV.load(preprocessNoob(rgSet))
  names(cases) <- sampleNames(rgSet)[cases] 
  CNA <- lapply(cases, function(case) {
    CNV.segment(CNV.bin(CNV.fit(cnvData[case], cnvData[controls], anno)))
  })
  attr(CNA, "SNPs") <- getSnpBeta(rgSet)[, cases]
  return(CNA)

}

# simplifies a recurrent task of mine 
.widen <- function(gr, bases=1e6) resize(gr, width(gr) + bases, fix="center")

# used in the creation of the AML-centric detail_regions in this package 
.makeDetailRegion <- function(geneRange, probeRange) {  # {{{

  stopifnot(length(geneRange) == 1)
  if(!"name" %in% names(mcols(geneRange))) geneRange$name <- names(geneRange)
  stopifnot(all(grepl("^cg", names(probeRange))))

  promRange <- promoters(geneRange) 
  promProbes <- countOverlaps(promRange, probeRange)[[1]]
  geneProbes <- countOverlaps(geneRange, probeRange)[[1]]
  totalProbes <- promProbes + geneProbes
  fullRange <- reduce(c(promRange, geneRange))
  thickRange <- ranges(.widen(fullRange))
  symbol <- geneRange$name

  detailRegion <- fullRange
  mcols(detailRegion)$name <- symbol 
  mcols(detailRegion)$thick <- thickRange
  mcols(detailRegion)$score <- totalProbes
  mcols(detailRegion)$probes_gene <- geneProbes
  mcols(detailRegion)$probes_promoter <- promProbes
  return(detailRegion)
 
} # }}}
