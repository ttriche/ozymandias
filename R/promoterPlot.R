#' simple function to grab all probes that overlap a promoter and plot 'em.
#' FIXME: bolt this onto coMET for prettier plotting and LD-blocking stuff.
#' 
#' @param grSet     the GenomicRatioSet from which to get the data
#' @param promoter  the promoter or other range to overlap 
#' @param ...       other arguments passed on to Heatmap()
#' 
#' @import ComplexHeatmap
#'
#' @export
promoterPlot <- function(grSet, promoter, ...) { 
  res <- sort(unique(subsetByOverlaps(grSet, promoter)))
  rownames(res) <- as(granges(res), "character") 
  Heatmap(t(getBeta(res)), cluster_columns=FALSE, col=jet, name="Methylation", 
          ...)
} 
