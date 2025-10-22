#FlowSOM::AggregateFlowFrames
erika <- function (fileNames, cTotal, channels = NULL, writeOutput = FALSE, 
          outputFile = "aggregate.fcs", keepOrder = FALSE, silent = FALSE, 
          sampleWithReplacement = FALSE, ...) 
  
{
  nFiles <- length(fileNames)
  cFile <- ceiling(cTotal/nFiles)
  flowFrame <- NULL
  fileMatrix <- NULL
  diffNumberChannels <- FALSE
  diffMarkers <- FALSE
  for (i in seq_len(nFiles)) {
    if (is(fileNames, "flowSet")) {
      file_name <- sampleNames(fileNames)[i]
      f <- fileNames[[i]]
    }
    else {
      file_name <- fileNames[i]
      if (!silent) {
        message("Reading ", file_name)
      }
      f <- flowCore::read.FCS(fileNames[i], ...)
    }
    if (isTRUE(sampleWithReplacement & (nrow(f) < cFile))) {
      ids <- c(seq_len(nrow(f)), sample(seq_len(nrow(f)), 
                                        cFile - nrow(f), replace = TRUE))
    }
    else {
      ids <- sample(seq_len(nrow(f)), min(nrow(f), cFile))
    }
    if (keepOrder) 
      ids <- sort(ids)
    file_ids <- rep(i, min(nrow(f), cFile))
    m <- cbind(file_ids, file_ids + stats::rnorm(length(file_ids), 
                                                 0, 0.1), ids)
    f <- f[ids, ]
    if (!all(channels %in% colnames(f))) {
      diffNumberChannels <- TRUE
      channelsToAdd <- channels[!channels %in% colnames(f)]
      extracols <- matrix(0, nrow = flowCore::nrow(f), 
                          ncol = length(channelsToAdd), dimnames = list(NULL, 
                                                                        channelsToAdd))
      f <- flowCore::fr_append_cols(f, extracols)
    }
    if (is.null(flowFrame)) {
      fileMatrix <- m
      if (is.null(channels)) {
        channels <- colnames(f)
        flowFrame <- f
      }
      else {
        channels <- GetChannels(f, channels)
        flowFrame <- f[, channels, drop = FALSE]
      }
      flowCore::keyword(flowFrame)[["$FIL"]] <- basename(outputFile)
      flowCore::keyword(flowFrame)[["FILENAME"]] <- basename(outputFile)
    }
    else {
      cols_f <- flowCore::colnames(f)
      cols_flowFrame <- flowCore::colnames(flowFrame)
      commonCols <- intersect(cols_f, cols_flowFrame)
      if (length(commonCols) == 0) 
        stop("No common channels between files")
      if (!diffNumberChannels && length(cols_flowFrame) != 
          length(commonCols)) {
        diffNumberChannels <- TRUE
      }
      if (!diffMarkers && any(!flowCore::markernames(f)[commonCols] %in% 
                              flowCore::markernames(flowFrame)[commonCols])) {
        diffMarkers <- TRUE
      }
      flowCore::exprs(flowFrame) <- rbind(flowCore::exprs(flowFrame)[, 
                                                                     commonCols, drop = FALSE], flowCore::exprs(f)[, 
                                                                                                                   commonCols, drop = FALSE])
      fileMatrix <- rbind(fileMatrix, m)
    }
  }
  colnames <- c("File", "File_scattered", "Original_ID")
  prev_agg <- length(grep("File[0-9]*$", colnames(flowFrame)))
  if (prev_agg > 0) {
    colnames[c(1, 2)] <- paste0(colnames[c(1, 2)], prev_agg + 
                                  1)
  }
  prev_ids <- length(grep("Original_ID[0-9]*$", colnames(flowFrame)))
  if (prev_ids > 0) {
    colnames[3] <- paste0(colnames[3], prev_ids + 1)
  }
  colnames(fileMatrix) <- colnames
  flowFrame <- flowCore::fr_append_cols(flowFrame, fileMatrix)
  if (diffNumberChannels) {
    warning("Files do not contain the same number of channels/markers. ", 
            "Zeros might have been imputed for missing values.")
  }
  if (diffMarkers) {
    warning("Files do not contain the same markers")
  }
  if (writeOutput) {
    flowCore::write.FCS(flowFrame, filename = outputFile)
  }
  return(flowFrame)
}