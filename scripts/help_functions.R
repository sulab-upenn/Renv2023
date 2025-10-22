#' Inverted versions of in
#'
#' @noRd
#' @examples
#' 1 %!in% 1:10
`%!in%` <- Negate(`%in%`)

#' Utility function for NULL coalescing
#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Check colname
#' @noRd
check_colname <- function(df_colnames, col_name, location = "metadata") {
  if (!is.null(col_name)) {
    if (col_name %!in% df_colnames) {
      stop("Column \"", col_name, "\" was not found in the ", location)
    }}
}

#' Check if directory exists, if not, make it
#' @noRd
check_make_dir <- function(dir.path) {
  if (!dir.exists(dir.path)) {dir.create(dir.path)}
}

#' Wrapper for missing packages
#'
#' @noRd
check_package <- function(package, repo = "CRAN", git_repo = "") {
  
  if (repo == "CRAN") {
    install_function <- "install.packages('"
  } else if (repo == "github") {
    install_function <- paste0("devtools::install_github('", git_repo, "/")
  } else if (repo == "Bioc") {
    install_function <- "BiocManager::install('"
  }
  
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(
      paste0("Package ", package, " is not installed.\n",
             "Please run: ", install_function, package, "')"))
  }
  requireNamespace(package, quietly = TRUE)
}

missing_package <- function(...) {
  check_package(...)
}

#'rename the column name with metadata
#'
#' @noRd
# rename_columns <- function(df, metadata) {
#   # Loop through each column name in the dataframe
#   for (i in 1:dim(metadata)[1]) {
#     col = metadata$Channel.name[i]
#     idx = grep(col, colnames(df))
#     col
#     idx
#     matching_channel <- metadata$Channel.name[grepl( col, metadata$Channel.name)]
#     matching_channel
#     if (length(idx) > 0) {
#       # Get the corresponding marker
#       marker <- metadata$markers[metadata$Channel.name == matching_channel]
#       # Rename the column with the marker
#       colnames(df)[idx] <- marker
#     } else {stop(paste(sep = " ", i,idx,col,"\n", "colnames of data/ and ORIGINAL MARKERS.csv are not matched"))}
#   }
#   return(df)
# }

rename_columns <- function(df, metadata) {
  # Loop through each column name in the dataframe
  for (i in 1:dim(metadata)[1]) {
    col = metadata[i,1]
    idx = grep(col, colnames(df))
    col
    idx
    matching_channel <- metadata[grepl( col, metadata[,1]),1]
    matching_channel
    if (length(idx) > 0) {
      # Get the corresponding marker
      marker <- metadata[metadata[,1] == matching_channel,2]
      # Rename the column with the marker
      colnames(df)[idx] <- marker
    } else {stop(paste(sep = " ", i,idx,col,"\n", "colnames of data/ and ORIGINAL MARKERS.csv are not matched"))}
  }
  return(df)
}


read.cytofFiles <- function (file.loc = getwd(), file.type = ".csv", nrows = NULL, 
                             do.embed.file.names = TRUE, header = TRUE) 
{
  if (!is.element("Spectre", installed.packages()[, 1])) 
    stop("Spectre is required but not installed")
  if (!is.element("data.table", installed.packages()[, 1])) 
    stop("data.table is required but not installed")
  require(Spectre)
  require(data.table)
  orig_wd <- getwd()
  if (!dir.exists(paste(orig_wd, file.loc, sep = "/")) & !dir.exists(file.loc)) {
    warning("We were not able to find the directory specified by file.loc. Are you sure that location exists?")
  }
  setwd(file.loc)
  wd <- getwd()
  if (length(list.files(path = wd, pattern = file.type)) == 
      0) {
    warning("We did not find any files in that directory, are you sure this is the right place?")
  }
  data.list = list()
  ncol.check = list()
  colName.check = list()
  nrow.check = list()
  if (file.type == ".csv") {
    file.names <- list.files(path = wd, pattern = file.type)
    for (file in file.names) {
      if (is.null(nrows)) {
        tempdata <- data.table::fread(file, check.names = FALSE, 
                                      header = header)
      }
      if (!is.null(nrows)) {
        message(paste0("Reading ", nrows, " rows (cells) per file"))
        tempdata <- data.table::fread(file, check.names = FALSE, 
                                      header = header, nrows = nrows)
      }
      file <- gsub(".csv", "", file)
      data.list[[file]] <- tempdata
    }
    rm(tempdata)
    msg <- "CSV files have been imported into a list"
  }
  if (file.type == ".fcs") {
    if (!is.element("flowCore", installed.packages()[, 1])) 
      stop("flowCore is required but not installed")
    require(flowCore)
    file.names <- list.files(path = wd, pattern = file.type)
    for (file in file.names) {
      if (is.null(nrows)) {
        x <- flowCore::read.FCS(file, transformation = FALSE,truncate_max_range = FALSE)
      }
      if (!is.null(nrows)) {
        message(paste0("Reading ", nrows, " rows (cells) per file"))
        x <- flowCore::read.FCS(file, transformation = FALSE, 
                                which.lines = nrows,truncate_max_range = FALSE)
      }
      nms <- vector()
      for (o in c(1:nrow(x@parameters@data))) {
        pr <- x@parameters@data$name[[o]]
        st <- x@parameters@data$desc[[o]]
        if (!is.na(st)) {
          nms <- c(nms, paste0(pr, "_", st))
        }
        else {
          nms <- c(nms, pr)
        }
      }
      tempdata <- exprs(x)
      tempdata <- tempdata[1:nrow(tempdata), 1:ncol(tempdata)]
      tempdata <- as.data.table(tempdata)
      names(tempdata) <- nms
      file <- gsub(".fcs", "", file)
      data.list[[file]] <- tempdata
    }
    rm(tempdata)
    msg <- "FCS files have been imported into a list"
  }
  if (do.embed.file.names == TRUE) {
    all.file.names <- c(names(data.list))
    all.file.names
    all.file.nums <- c(1:(length(data.list)))
    all.file.nums
    for (a in all.file.names) {
      data.list[[a]]$FileName <- a
    }
    for (i in all.file.nums) {
      data.list[[i]]$FileNo <- i
    }
  }
  setwd(orig_wd)
  return(data.list)
  message(msg)
}

make.cytofheatmap <- function (dat, sample.col, plot.cols, annot.cols = NULL, feature.annots = NULL, 
                               annotation_colors = NULL, file.name = paste0("Pheatmap by ", 
                                                                            sample.col, ".png"), plot.title = paste0(sample.col, 
                                                                                                                     " heatmap"), transpose = FALSE, normalise = TRUE, is.fold = FALSE, 
                               fold.range = NULL, dendrograms = "both", cutree_rows = 1, 
                               cutree_cols = 1, row.sep = c(), col.sep = c(), cell.size = 15, 
                               standard.colours = "BuPu", fold.colours = "Spectre", path = NULL) 
{
  if (!is.element("pheatmap", installed.packages()[, 1])) 
    stop("pheatmap is required but not installed")
  if (!is.element("RColorBrewer", installed.packages()[, 1])) 
    stop("RColorBrewer is required but not installed")
  if (!is.element("scales", installed.packages()[, 1])) 
    stop("scales is required but not installed")
  require(pheatmap)
  require(RColorBrewer)
  require(scales)
  
  #install.packages("circlize")
  if (standard.colours == "rev(RdBu)") {
    colour.palette <- (colorRampPalette(RColorBrewer::brewer.pal(9, 
                                                                 "RdBu"))(31))
    colour.palette <- rev(colour.palette)
  }
  if (standard.colours == "BuPu") {
    colour.palette <- (colorRampPalette(RColorBrewer::brewer.pal(9, 
                                                                 "BuPu"))(31))
  }
  dat <- as.data.frame(dat)
  heatmap.data <- dat
  rownames(heatmap.data) <- t(dat[sample.col])
  heatmap.data
  if (is.null(annot.cols) == FALSE) {
    annot <- heatmap.data[annot.cols]
    heatmap.data <- heatmap.data[plot.cols]
    heatmap.data
  }
  if (is.null(annot.cols) == TRUE) {
    annot <- NULL
    heatmap.data <- heatmap.data[plot.cols]
    heatmap.data
  }
  if (transpose == TRUE) {
    heatmap.data.t <- as.data.frame(t(heatmap.data))
    heatmap.data <- heatmap.data.t
  }
  if (normalise == TRUE) {
    if (is.fold == FALSE) {
      row.nam <- row.names(heatmap.data)
      col.nam <- names(heatmap.data)
      norm.fun <- function(x) {
        (x - min(x, na.rm = TRUE))/(max(x, na.rm = TRUE) - 
                                      min(x, na.rm = TRUE))
      }
      heatmap.data.norm <- as.data.frame(lapply(heatmap.data, 
                                                norm.fun))
      names(heatmap.data.norm) <- col.nam
      max(heatmap.data.norm)
      heatmap.data.norm <- as.matrix(heatmap.data.norm)
      heatmap.data <- heatmap.data.norm
      rownames(heatmap.data) <- row.nam
    }
  }
  heatmap.data <- as.matrix(heatmap.data)
  if (dendrograms == "none") {
    row.clustering <- FALSE
    col.clustering <- FALSE
  }
  if (dendrograms != "none") {
    hclustfunc <- function(x) hclust(x, method = "complete")
    distfunc <- function(x) dist(x, method = "euclidean")
    if (dendrograms == "both") {
      row.clustering <- TRUE
      col.clustering <- TRUE
    }
    if (dendrograms == "column") {
      row.clustering <- FALSE
      col.clustering <- TRUE
    }
    if (dendrograms == "row") {
      row.clustering <- TRUE
      col.clustering <- FALSE
    }
  }
  
  if (is.fold == FALSE) {
    map.colour <- colour.palette
    sym.key <- FALSE
    sym.breaks <- FALSE
    heatmap.data
    my.max <- function(x) ifelse(!all(is.na(x)), max(x, 
                                                     na.rm = T), NA)
    my.min <- function(x) ifelse(!all(is.na(x)), min(x, 
                                                     na.rm = T), NA)
    # my.breaks <- c(seq(my.min(heatmap.data), median(heatmap.data),length.out=15),seq(median(heatmap.data),my.max(heatmap.data),length.out=17))
    
    my.breaks <- c(seq(my.min(heatmap.data),quantile(heatmap.data,c(0.5)),length.out=15),
                   seq(quantile(heatmap.data,c(0.55)),my.max(heatmap.data),length=17))
  }
  scale.set <- "none"
  title.text <- plot.title
  if (is.null(path)) {
    flnm <- file.name
  }
  if (!is.null(path)) {
    flnm <- paste0(path, "/", file.name)
  }
  
  
  pheatmap::pheatmap(mat = as.matrix(heatmap.data), main = title.text, 
                     cellwidth = cell.size, cellheight = cell.size, cluster_rows = row.clustering, 
                     cluster_cols = col.clustering, breaks = my.breaks, cutree_rows = cutree_rows, 
                     cutree_cols = cutree_cols, gaps_row = row.sep, gaps_col = col.sep, 
                     annotation_row = annot, annotation_col = feature.annots, 
                     annotation_colors = annotation_colors, color = map.colour, 
                     filename = flnm)
  message(paste0("A pheatmap has been saved to your working directory", 
                 paste0(path, file.name)))
}

amshaw <- function (dat, cellular.cols, cluster.cols, batch.col, sample.col = NULL, 
                    dir = getwd(), xdim = 5, ydim = 5, meta.k = 10, seed = 42, 
                    mem.ctrl = TRUE) 
{
  if (!is.element("Spectre", installed.packages()[, 1])) 
    stop("Spectre is required but not installed")
  if (!is.element("data.table", installed.packages()[, 1])) 
    stop("data.table is required but not installed")
  if (!is.element("CytoNorm", installed.packages()[, 1])) 
    stop("CytoNorm is required but not installed")
  if (!is.element("flowCore", installed.packages()[, 1])) 
    stop("flowCore is required but not installed")
  if (!is.element("Biobase", installed.packages()[, 1])) 
    stop("Biobase is required but not installed")
  require(data.table)
  require(CytoNorm)
  require(Biobase)
  setwd(dir)
  starting.dir <- getwd()
  message("Working directory is '", starting.dir, "'")
  if (nrow(dat) == 0) {
    setwd(dir)
    stop("Error -- your 'dat' has no rows, please check")
  }
  if (meta.k == 2) {
    setwd(dir)
    stop("Error -- cannot create '2' metaclusters -- please set meta.k to '1' or a value >3")
  }
  all.cols <- unique(c(cellular.cols, cluster.cols))
  all.cols
  value <- dat[, all.cols, with = FALSE]
  if (isFALSE(all(sapply(value, is.numeric)))) {
    message("It appears that one column in your dataset is non numeric")
    print(sapply(value, is.numeric))
    setwd(dir)
    stop("do.asinh stopped")
  }
  rm(value)
  if (mem.ctrl == TRUE) {
    gc()
  }
  if (!is.null(sample.col)) {
    dat <- as.data.table(dat)
    dat <- dat[, c(sample.col, batch.col, all.cols), with = FALSE]
  }
  if (is.null(sample.col)) {
    dat <- as.data.table(dat)
    dat <- dat[, c(batch.col, all.cols), with = FALSE]
  }
  if (meta.k == 1) {
    meta.k <- 5
    one.clust <- TRUE
  }
  else {
    one.clust <- FALSE
  }
  setwd(starting.dir)
  unlink("tmp-cytonorm-fsom", recursive = TRUE)
  dir.create("tmp-cytonorm-fsom", showWarnings = FALSE)
  setwd("tmp-cytonorm-fsom")
  setwd(starting.dir)
  if (!is.null(sample.col)) {
    samps <- unique(dat[[sample.col]])
    samps.num <- c(1:length(samps))
    smp.tb <- cbind(samps, samps.num)
    x <- as.data.table(dat[[sample.col]])
    names(x) <- sample.col
    x <- do.add.cols(x, sample.col, smp.tb, "samps")
    x <- x$samps.num
    dat[[sample.col]] <- x
    dat[[sample.col]] <- as.numeric(dat[[sample.col]])
    rm(x)
    rm(samps)
    rm(samps.num)
  }
  message("Step 1/4 - Splitting files for use with original FlowSOM function")
  dat.list <- unique(dat[[batch.col]])
  dat.list
  setwd(starting.dir)
  dir.create("tmp-cytonorm-fsom", showWarnings = FALSE)
  setwd("tmp-cytonorm-fsom")
  for (i in c(1:length(dat.list))) {
    a <- dat.list[[i]]
    temp <- dat[dat[[batch.col]] == a, ]
    if (!is.null(sample.col)) {
      temp <- temp[, c(sample.col, all.cols), with = FALSE]
    }
    if (is.null(sample.col)) {
      temp <- temp[, c(all.cols), with = FALSE]
    }
    write.files(temp, file.prefix = a, write.csv = FALSE, 
                write.fcs = TRUE)
    rm(i)
    rm(a)
  }
  rm(dat.list)
  files <- list.files(getwd(), ".fcs")
  files
  file.nums <- c(1:length(files))
  setwd(starting.dir)
  if (mem.ctrl == TRUE) {
    gc()
  }
  message("Step 2/4 - Running FlowSOM")
  setwd(starting.dir)
  dir.create("tmp-cytonorm-fsom", showWarnings = FALSE)
  setwd("tmp-cytonorm-fsom")
  fsom <- erika1(files, colsToUse = cluster.cols, nCells = NULL, 
                 FlowSOM.params = list(xdim = xdim, ydim = ydim, nClus = meta.k, 
                                       scale = FALSE), seed = seed)
  setwd(starting.dir)
  unlink("tmp-cytonorm-fsom", recursive = TRUE)
  if (nrow(fsom$data) != nrow(dat)) {
    stop("Error - the numer of rows (cells) is different in the starting dataset and the FlowSOM prepared dataset")
  }
  if (mem.ctrl == TRUE) {
    gc()
  }
  message("Step 3/4 - Preparing FlowSOM object and results data.table")
  if (one.clust == TRUE) {
    length(fsom$metaclustering)
    fsom$metaclustering <- rep(1, length(fsom$metaclustering))
  }
  A <- fsom$data
  B <- fsom$map$mapping[, 1]
  C <- fsom$metaclustering[fsom$map$mapping[, 1]]
  dt <- as.data.table(A)
  dt <- cbind(A, prep.fsom.cluster = B, prep.fsom.metacluster = C)
  dt <- as.data.table(dt)
  files <- gsub(".fcs", "", files)
  res <- named.list(fsom, dt, cellular.cols, cluster.cols, 
                    files, file.nums)
  rm(fsom)
  rm(dt)
  if (mem.ctrl == TRUE) {
    gc()
  }
  message("Step 4/4 - FlowSOM preparation for alignment model complete")
  setwd(starting.dir)
  return(res)
}

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

make.z_norm_pheatmap = function (dat, 
                                 sample.col, plot.cols, 
                                 annot.cols = NULL, 
                                 feature.annots = NULL, 
                                 annotation_colors = NULL, 
                                 file.name = paste0("Pheatmap by ", sample.col, ".png"),
                                 plot.title = paste0(sample.col, " heatmap"), 
                                 transpose = FALSE, 
                                 normalise = FALSE, 
                                 is.fold = FALSE, 
                                 fold.range = NULL, dendrograms = "both", dendrograms.sort = FALSE, 
                                 cutree_rows = 1, cutree_cols = 1, row.sep = c(), col.sep = c(), 
                                 cell.size = 15, standard.colours = "BuPu", fold.colours = "Spectre", 
                                 path = NULL) {
  if (!is.element("pheatmap", installed.packages()[, 1])) 
    stop("pheatmap is required but not installed")
  if (!is.element("RColorBrewer", installed.packages()[, 1])) 
    stop("RColorBrewer is required but not installed")
  if (!is.element("scales", installed.packages()[, 1])) 
    stop("scales is required but not installed")
  require(pheatmap)
  require(RColorBrewer)
  require(scales)
  if (standard.colours == "BuPu") {
    colour.palette <- (colorRampPalette(RColorBrewer::brewer.pal(9, 
                                                                 "BuPu"))(31))
  }
  if (standard.colours == "RdYlBu") {
    colour.palette <- (colorRampPalette(RColorBrewer::brewer.pal(9, 
                                                                 "RdYlBu"))(31))
    colour.palette <- rev(colour.palette)
  }
  if (standard.colours == "rev(RdBu)") {
    colour.palette <- (colorRampPalette(RColorBrewer::brewer.pal(9, 
                                                                 "RdBu"))(31))
    colour.palette <- rev(colour.palette)
  }
  if (standard.colours == "Blues") {
    colour.palette <- (colorRampPalette(RColorBrewer::brewer.pal(9, 
                                                                 "Blues"))(31))
  }
  if (standard.colours == "Reds") {
    colour.palette <- (colorRampPalette(RColorBrewer::brewer.pal(9, 
                                                                 "Reds"))(31))
  }
  if (standard.colours == "Greys") {
    colour.palette <- (colorRampPalette(RColorBrewer::brewer.pal(9, 
                                                                 "Greys"))(31))
  }
  if (standard.colours == "YlGnBu") {
    colour.palette <- (colorRampPalette(RColorBrewer::brewer.pal(9, 
                                                                 "YlGnBu"))(31))
  }
  if (standard.colours == "viridis") {
    colour.palette <- colorRampPalette(c((scales::viridis_pal(option = "viridis"))(50)))
    colour.palette <- colour.palette(31)
  }
  if (standard.colours == "spectral") {
    spectral.list <- colorRampPalette(RColorBrewer::brewer.pal(11, 
                                                               "Spectral"))(50)
    spectral.list <- rev(spectral.list)
    colour.palette <- colorRampPalette(c(spectral.list))
    colour.palette <- colour.palette(31)
  }
  if (standard.colours == "magma") {
    colour.palette <- colorRampPalette(c((scales::viridis_pal(option = "magma"))(50)))
    colour.palette <- colour.palette(31)
  }
  if (standard.colours == "inferno") {
    colour.palette <- colorRampPalette(c((scales::viridis_pal(option = "inferno"))(50)))
    colour.palette <- colour.palette(31)
  }
  if (fold.colours == "Spectre") {
    fold.palette <- colorRampPalette(rev(c("#ffeda0", "#fed976", 
                                           "#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#bd0026", 
                                           "#800026", "black", "#023858", "#045a8d", "#0570b0", 
                                           "#3690c0", "#74a9cf", "#a6bddb", "#d0d1e6", "#ece7f2")))
    fold.palette <- fold.palette(31)
  }
  if (fold.colours == "BuPu") {
    fold.palette <- (colorRampPalette(RColorBrewer::brewer.pal(9, 
                                                               "BuPu"))(31))
  }
  if (fold.colours == "RdYlBu") {
    fold.palette <- (colorRampPalette(RColorBrewer::brewer.pal(9, 
                                                               "RdYlBu"))(31))
    fold.palette <- rev(fold.palette)
  }
  if (fold.colours == "rev(RdBu)") {
    fold.palette <- (colorRampPalette(RColorBrewer::brewer.pal(9, 
                                                               "RdBu"))(31))
    fold.palette <- rev(fold.palette)
  }
  if (fold.colours == "Blues") {
    fold.palette <- (colorRampPalette(RColorBrewer::brewer.pal(9, 
                                                               "Blues"))(31))
  }
  if (fold.colours == "Reds") {
    fold.palette <- (colorRampPalette(RColorBrewer::brewer.pal(9, 
                                                               "Reds"))(31))
  }
  if (fold.colours == "Greys") {
    fold.palette <- (colorRampPalette(RColorBrewer::brewer.pal(9, 
                                                               "Greys"))(31))
  }
  if (fold.colours == "YlGnBu") {
    fold.palette <- (colorRampPalette(RColorBrewer::brewer.pal(9, 
                                                               "YlGnBu"))(31))
  }
  if (fold.colours == "viridis") {
    fold.palette <- colorRampPalette(c((scales::viridis_pal(option = "viridis"))(50)))
    fold.palette <- fold.palette(31)
  }
  if (fold.colours == "spectral") {
    spectral.list <- colorRampPalette(RColorBrewer::brewer.pal(11, 
                                                               "Spectral"))(50)
    spectral.list <- rev(spectral.list)
    fold.palette <- colorRampPalette(c(spectral.list))
    fold.palette <- fold.palette(31)
  }
  if (fold.colours == "magma") {
    fold.palette <- colorRampPalette(c((scales::viridis_pal(option = "magma"))(50)))
    fold.palette <- fold.palette(31)
  }
  if (fold.colours == "inferno") {
    fold.palette <- colorRampPalette(c((scales::viridis_pal(option = "inferno"))(50)))
    fold.palette <- fold.palette(31)
  }
  dat <- as.data.frame(dat)
  heatmap.data <- dat
  rownames(heatmap.data) <- t(dat[sample.col])
  heatmap.data
  if (is.null(annot.cols) == FALSE) {
    annot <- heatmap.data[annot.cols]
    heatmap.data <- heatmap.data[plot.cols]
    heatmap.data
  }
  if (is.null(annot.cols) == TRUE) {
    annot <- NULL
    heatmap.data <- heatmap.data[plot.cols]
    heatmap.data
  }
  if (transpose == TRUE) {
    heatmap.data.t <- as.data.frame(t(heatmap.data))
    heatmap.data <- heatmap.data.t
  }
  if (normalise == TRUE) {
    if (is.fold == FALSE) {
      row.nam <- row.names(heatmap.data)
      col.nam <- names(heatmap.data)
      norm.fun <- function(x) {
        (x - min(x, na.rm = TRUE))/(max(x, na.rm = TRUE) - 
                                      min(x, na.rm = TRUE))
      }
      heatmap.data.norm <- as.data.frame(lapply(heatmap.data, 
                                                norm.fun))
      names(heatmap.data.norm) <- col.nam
      max(heatmap.data.norm)
      heatmap.data.norm <- as.matrix(heatmap.data.norm)
      heatmap.data <- heatmap.data.norm
      rownames(heatmap.data) <- row.nam
    }
  }
  heatmap.data <- as.matrix(heatmap.data)
  if (dendrograms == "none") {
    row.clustering <- FALSE
    col.clustering <- FALSE
  }
  if (dendrograms != "none") {
    hclustfunc <- function(x) hclust(x, method = "complete")
    distfunc <- function(x) dist(x, method = "euclidean")
    if (dendrograms == "both") {
      row.clustering <- TRUE
      col.clustering <- TRUE
      if (isTRUE(dendrograms.sort)) {
        row.clustering <- hclustfunc(distfunc(heatmap.data))
        col.clustering <- hclustfunc(distfunc(t(heatmap.data)))
        require(dendsort)
        sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
        row.clustering <- sort_hclust(row.clustering)
        col.clustering <- sort_hclust(col.clustering)
      }
    }
    if (dendrograms == "column") {
      row.clustering <- FALSE
      col.clustering <- TRUE
      if (isTRUE(dendrograms.sort)) {
        col.clustering <- hclustfunc(distfunc(t(heatmap.data)))
        require(dendsort)
        sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
        col.clustering <- sort_hclust(col.clustering)
      }
    }
    if (dendrograms == "row") {
      row.clustering <- TRUE
      col.clustering <- FALSE
      if (isTRUE(dendrograms.sort)) {
        row.clustering <- hclustfunc(distfunc(t(heatmap.data)))
        require(dendsort)
        sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
        row.clustering <- sort_hclust(row.clustering)
      }
    }
  }
  if (is.fold == TRUE) {
    map.colour <- fold.palette
    sym.key <- FALSE
    sym.breaks <- TRUE
    if (is.null(fold.range)) {
      fld.max <- max(heatmap.data, na.rm = TRUE)
      fld.min <- min(heatmap.data, na.rm = TRUE)
      if (fld.max == -fld.min) {
        fold.max.range <- fld.max
        fold.min.range <- fld.min
      }
      if (fld.max > -fld.min) {
        fold.max.range <- fld.max
        fold.min.range <- -fld.max
      }
      if (fld.max < -fld.min) {
        fold.max.range <- -fld.min
        fold.min.range <- fld.min
      }
    }
    if (!is.null(fold.range)) {
      fold.max.range <- fold.range[1]
      fold.min.range <- fold.range[2]
    }
    my.breaks <- seq(fold.min.range, fold.max.range, length.out = 32)
  }
  if (is.fold == FALSE) {
    map.colour <- colour.palette
    sym.key <- FALSE
    sym.breaks <- FALSE
    heatmap.data
    my.max <- function(x) ifelse(!all(is.na(x)), max(x, na.rm = T), 
                                 NA)
    my.min <- function(x) ifelse(!all(is.na(x)), min(x, na.rm = T), 
                                 NA)
    my.breaks <- seq(my.min(heatmap.data), my.max(heatmap.data), 
                     length.out = 32)
  }
  scale.set <- "none"
  title.text <- plot.title
  if (is.null(path)) {
    flnm <- file.name
  }
  if (!is.null(path)) {
    flnm <- paste0(path, "/", file.name)
  }
  pheatmap::pheatmap(mat = as.matrix(heatmap.data), main = title.text, 
                     cellwidth = cell.size, cellheight = cell.size, cluster_rows = row.clustering, 
                     cluster_cols = col.clustering, breaks = my.breaks, cutree_rows = cutree_rows, 
                     cutree_cols = cutree_cols, gaps_row = row.sep, gaps_col = col.sep, 
                     annotation_row = annot, annotation_col = feature.annots, 
                     annotation_colors = annotation_colors, color = map.colour, scale = "row",
                     filename = flnm)
  message(paste0("A pheatmap has been saved to your working directory", 
                 paste0(path, file.name)))
}

run_with_fallback <- function(expr_primary, expr_fallback, verbose = TRUE) {
  tryCatch(
    {
      if (verbose) message("Running primary function...")
      result <- eval(expr_primary, envir = parent.frame())
      if (verbose) message("Primary function succeeded.")
      return(result)
    },
    error = function(e) {
      message("Primary function failed: ", e$message)
      message("Running fallback function...")
      fallback_result <- eval(expr_fallback, envir = parent.frame())
      message("allback function completed.")
      return(fallback_result)
    }
  )
}