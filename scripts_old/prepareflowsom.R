erika1 <- function (files, colsToUse, nCells = 1e+06, FlowSOM.params = list(xdim = 15, 
                                                                  ydim = 15, nClus = 30, scale = FALSE), transformList = NULL, 
          verbose = FALSE, seed = NULL, ...) 
  

{
  if (verbose) 
    message("Aggregating files ... ")
  if (!is.null(seed)) 
    set.seed(seed)
  o <- capture.output(ff <- erika(files, 
                                                         nCells, channels = colsToUse, ...))
  if (!is.null(transformList)) 
    ff <- flowCore::transform(ff, transformList)
  FlowSOM.params <- c(FlowSOM.params, list(input = ff, colsToUse = colsToUse, 
                                           seed = seed))
  FlowSOM.params <- FlowSOM.params[unique(names(FlowSOM.params))]
  if (verbose) 
    message("Running the FlowSOM algorithm ... ")
  fsom <- do.call(FlowSOM::FlowSOM, FlowSOM.params)
  fsom
}
  