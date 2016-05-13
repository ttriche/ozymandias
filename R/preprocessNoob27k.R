#' avoid consulting an annotation package for a platform with only Type I probes
#' 
#' @param rgSet   an RGChannelSet
#' @param offset  offset term for intensities (default is 15)
#' @param dyeCorr ignored (no need on 27k)
#' @param verbose make noise? (TRUE) 
#' 
#' @import minfi 
#' 
#' @return        a MethylSet
#'
#' @export 
preprocessNoob27k <- function(rgSet, offset=15, dyeCorr=FALSE, verbose = TRUE) {

  manifest <- getManifest(rgSet)
  subverbose <- max(as.integer(verbose) - 1L, 0)

  ## Extraction of the out-of-band controls
  controls <- getOOB(rgSet) # works fine on 27k
  names(controls) <- c("Cy3", "Cy5")
  mset <- preprocessRaw(rgSet)
  meth <- getMeth(mset)
  unmeth <- getUnmeth(mset)

  if (any(meth<=0)) meth[which(meth<=0)] <- 1
  if (any(unmeth<=0)) unmeth[which(unmeth<=0)] <- 1

  G.probes <- getProbeInfo(manifest, type="I-Green")$Name
  cy3.probes <- which(rownames(mset) %in% G.probes)
  R.probes <- getProbeInfo(manifest, type="I-Red")$Name
  cy5.probes <- which(rownames(mset) %in% R.probes)

  dat <- list(Cy3 = list(M =  as.matrix(meth[cy3.probes,]), 
                         U =  as.matrix(unmeth[cy3.probes,])),
              Cy5 = list(M =  as.matrix(meth[cy5.probes,]), 
                         U =  as.matrix(unmeth[cy5.probes,])))

  rows <- lapply(dat, function(ch) {
    sapply(names(ch), function(nch) {
      nrow(ch[[nch]])
    })
  })
  last <- lapply(rows, cumsum)
  first <- lapply(names(last), function(nch) {
    last[[nch]] - rows[[nch]] + 1
  })
  names(first) <- names(last)

  estimates <- lapply(names(dat), function(nch) { 
    xf <- rbind(dat[[nch]][['M']], dat[[nch]][['U']])
    xs <- minfi:::normexp.get.xs(xf = xf, controls = controls[[nch]], 
                                 offset=offset, verbose = subverbose)
    names(xs[['params']]) <- paste(names(xs[['params']]), nch, sep='.')
    names(xs[['meta']]) <- paste(names(xs[['meta']]), nch, sep='.')
    xs
  })
  names(estimates) <- names(dat)

  cy3.M <- first[['Cy3']][['M']]:last[['Cy3']][['M']]
  meth[cy3.probes, ] <- estimates[['Cy3']][['xs']][cy3.M,]
  cy3.U <- first[['Cy3']][['U']]:last[['Cy3']][['U']]
  unmeth[cy3.probes,] <- estimates[['Cy3']][['xs']][cy3.U,]

  cy5.M <- first[['Cy5']][['M']]:last[['Cy5']][['M']]
  meth[cy5.probes,] <- estimates[['Cy5']][['xs']][cy5.M,]
  cy5.U <- first[['Cy5']][['U']]:last[['Cy5']][['U']]
  unmeth[cy5.probes,] <- estimates[['Cy5']][['xs']][cy5.U,]

  assayDataElement(mset, "Meth") <- meth
  assayDataElement(mset, "Unmeth") <- unmeth
  mset@preprocessMethod <- c(preprocessMethod(mset), 
                             mu.norm = sprintf("Noob, dyeCorr=%s", dyeCorr))
  return(mset)
}
