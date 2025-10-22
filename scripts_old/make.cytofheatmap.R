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
