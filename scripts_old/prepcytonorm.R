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