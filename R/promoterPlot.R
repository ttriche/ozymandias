#' simple function to grab all probes that overlap a promoter and plot 'em.
#' FIXME: bolt this onto coMET for prettier plotting and LD-blocking stuff.
#' 
#' @param grSet     the GenomicRatioSet from which to get the data
#' @param promoter  the promoter or other range to overlap 
#' @param ...       other arguments passed on to Heatmap()
#' 
#' @import ComplexHeatmap
#' @import impute 
#' 
#' @export
promoterPlot <- function(grSet, promoter, ...) { 
  res <- sort(unique(subsetByOverlaps(grSet, promoter)))
  rownames(res) <- as(granges(res), "character") 
  betas <- getBeta(res)
  if (any(is.na(betas))) betas <- impute.knn(betas)$data
  Heatmap(t(betas), cluster_columns=FALSE, col=jet, name="Methylation", 
          row_names_gp=gpar(fontsize=7), column_names_gp=gpar(fontsize=7), ...)
} 
