#' simple function to grab all probes that overlap a promoter and plot 'em.
#' FIXME: bolt this onto coMET for prettier plotting and LD-blocking stuff.
#' 
#' @param grSet     the GenomicRatioSet from which to get the data
#' @param promoter  the promoter or other range to overlap 
#' @param rotate    rotate the plot so that samples are rows? (TRUE) 
#' @param flank     flank the given region(s) upstream? (TRUE) 
#' @param ...       other arguments passed on to Heatmap()
#' 
#' @import ComplexHeatmap
#' @import impute 
#' 
#' @export
promoterPlot <- function(grSet, promoter, rotate=TRUE, flank=TRUE, ...) { 
  if (flank == TRUE) promoter <- promoters(promoter)
  res <- sort(unique(subsetByOverlaps(grSet, promoter)))
  rownames(res) <- as(granges(res), "character") 
  betas <- getBeta(res)
  if (any(is.na(betas))) betas <- impute.knn(betas)$data
  if (rotate == TRUE) {
    Heatmap(betas, cluster_rows=FALSE, col=jet, name="Methylation", 
      row_names_gp=gpar(fontsize=6), column_names_gp=gpar(fontsize=6), ...)
  } else { 
    Heatmap(t(betas), cluster_columns=FALSE, col=jet, name="Methylation", 
      row_names_gp=gpar(fontsize=6), column_names_gp=gpar(fontsize=6), ...)
  }
} 
