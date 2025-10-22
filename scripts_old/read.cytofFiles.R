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
