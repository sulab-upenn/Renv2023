run.spectre <- function (phenok,metaFile,markerFile, do.plot, do.summary, do.batchAlign, do.fcsExport,flowType,coFactor,plot.against,ref.ctrls)
  
  # this version does not allow clustering with a subset of the markers yet. aka single marker staining plotted=markers used for clustering
{
  if (!is.element("FastPG", installed.packages()[, 1]))
    stop("fastPG is required but not installed")
  if (!is.element("CytoNorm", installed.packages()[, 1])& do.batchAlign!= FALSE)
    stop("CytoNorm is required but not installed")
  ## spectre lib
  library(Spectre)
  Spectre::package.check()    # Check that all required packages are installed
  Spectre::package.load()     # Load required packages
  ## fastPG lib
  #BiocManager::install("sararselitsky/FastPG") #install FastPG- the fast-performing version of phenograph
  #library(devtools)
  #install_github('saeyslab/CytoNorm')
  library(FastPG)
  
  
  ### Set DT threads
  
  getDTthreads()
  
  ### Set primary directory
  
  dirname(rstudioapi::getActiveDocumentContext()$path)            # Finds the directory where this script is located
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Sets the working directory to where the script is located
  getwd()
  PrimaryDirectory <- getwd()
  
  ### Set input directory
  
  setwd(PrimaryDirectory)
  
  setwd("data")
  InputDirectory <- getwd()
  InputDirectory
  setwd(PrimaryDirectory)
  if (flowType == "cytof"){
    source("read.cytofFiles.R")
    source("make.cytofheatmap.R")
  }
  
  ### Set metadata directory
  
  setwd(PrimaryDirectory)
  setwd("metadata")
  MetaDirectory <- getwd()
  MetaDirectory
  setwd(PrimaryDirectory)
  
  
  ### Set output directory
  
  setwd(PrimaryDirectory)
  dir.create("Output 1 - data prep", showWarnings = FALSE)
  dir.create("Output 3 - clustering and DR_fastPG", showWarnings = FALSE)
  if (do.batchAlign!=FALSE){
    dir.create("Output 2 - batch alignment", showWarnings = FALSE)
    Output2Directory <- "Output 2 - batch alignment"
  }
  Output1Directory <- "Output 1 - data prep"
  Output3Directory <- "Output 3 - clustering and DR_fastPG"
  
  setwd(PrimaryDirectory)
  
  
  
  
  ### Import data
  setwd(InputDirectory)
  list.files(InputDirectory, ".fcs")
  
  
  
  if (flowType == "cytof") {
    
    data.list <- read.cytofFiles (file.loc = InputDirectory,
                                  file.type = ".fcs",
                                  do.embed.file.names = TRUE)
  } else { 
    data.list <- Spectre::read.files(file.loc = InputDirectory,
                                     file.type = ".fcs",
                                     do.embed.file.names = TRUE)
  }
  
  
  
  #if cytof files, truncate max range in flowcore read.fcs is set to false
  
  markers <- read.csv(markerFile)[,2]
  markerLength <- as.numeric(length((markers)))
  colnames<- append(markers, c("FileName","FileNo"))
  
  data.list <- lapply(data.list, setNames, colnames)
  check <- do.list.summary(data.list)
  
  ### Merge data
  
  cell.dat <- Spectre::do.merge.files(dat = data.list)
  cell.dat
  
  
  ##########################################################################################################
  #### Add metadata
  ##########################################################################################################
  
  setwd(MetaDirectory)
  
  ### Read in sample metadata
  
  meta.dat <- fread(metaFile)
  
  
  sample.info <- meta.dat[,c(1:6)] 
  
  
  ### Add sample metadata to primary data.table
  
  cell.dat <- do.add.cols(cell.dat, "FileName", sample.info, "FileName", rmv.ext = TRUE)
  
  
  ##########################################################################################################
  #### Arcsinh transformation
  ##########################################################################################################
  setwd(PrimaryDirectory)
  setwd(Output1Directory)
  
  
  ### Transformation settings
  
  
  ### Make plots of transformed columns from the subsampled data
  
  
  
  if (flowType == "cytof") {
    cell.dat <- do.asinh(cell.dat, use.cols = names(cell.dat)[1:markerLength], cofactor = 5)
  } else if (flowType == "aurora") {
    if (is.null(coFactor)){ cell.dat <- do.asinh(cell.dat, use.cols = names(cell.dat)[1:markerLength], cofactor = 2000)}
    else {cell.dat <- do.asinh(cell.dat, use.cols = names(cell.dat)[1:markerLength], cofactor = coFactor)}
    
  } else if (flowType == "flow") {
    if (is.null(coFactor)){ cell.dat <- do.asinh(cell.dat, use.cols = names(cell.dat)[1:markerLength], cofactor = 200)}
    else {cell.dat <- do.asinh(cell.dat, use.cols = names(cell.dat)[1:markerLength], cofactor = coFactor)}
  } else {
    stop("flowType must be flow, cytof or aurora!")
  }
  
  
  
  dir.create("Output 1.2 - transformed plots")
  setwd("Output 1.2 - transformed plots")
  transf.cols <- names(cell.dat)[grepl('_asinh', names(cell.dat))]
  which(names(cell.dat) == plot.against)
  
  for(i in transf.cols){
    make.colour.plot(cell.dat, i, col.min.threshold = 0, plot.against)
  }
  
  ##########################################################################
  #### Write data to disk
  ##########################################################################################################
  setwd(PrimaryDirectory)
  setwd(Output1Directory)
  dir.create("Output 1.1 - transformed data")
  setwd("Output 1.1 - transformed data")
  
  ### Write cellular data and analysis  preferences to disk
  
  fwrite(cell.dat, "cell.dat.csv") # data
  
  
  
  
  
  
  
  
  
  
  
  
  
  ### Save session info to disk
  setwd(PrimaryDirectory)
  setwd(Output1Directory)
  dir.create("Output - info", showWarnings = FALSE)
  setwd("Output - info")
  
  sink(file = "session_info.txt", append=TRUE, split=FALSE, type = c("output", "message"))
  session_info()
  sink()
  
  message("Data pre-processing completed")
  
  
  
  
  
  if (do.batchAlign==TRUE){
    
    sample.col <- 'Sample'
    group.col <- 'Group'
    batch.col <- 'Batch'  
    donor.col <- "Donor" # added specifically for our study
    cellular.cols <- names(cell.dat)[(markerLength+8):(markerLength*2+7)]
    cluster.cols <- cellular.cols
    
    as.matrix(unique(cell.dat[[sample.col]]))
    
    ref.dat <- do.filter(cell.dat, use.col = sample.col, values =ref.ctrls)
    
   # View(ref.dat)
    
    
    ### Check reference in ref.dat and cell.dat
    
    unique(ref.dat[[batch.col]])
    unique(cell.dat[[batch.col]])
    
    unique(ref.dat[[sample.col]])
    unique(cell.dat[[sample.col]])
    
   # View(ref.dat)
    
    ##########################################################################################################
    #### Initial (pre-alignment) plots
    ##########################################################################################################
    setwd(PrimaryDirectory)
    setwd(Output2Directory)
    dir.create("Output 2.1 - pre-Align plots")
    setwd("Output 2.1 - pre-Align plots")
    
    ### Pre-alignment UMAP
    
    rm(sub)
    #sub <- do.subsample(cell.dat, 10000) 
    sub <- subset(cell.dat,Group!="CD4 BC" ) # excluding the bulk CD4 which are used for batch control
    #sub <- cell.dat
    sub <- run.umap(sub, cellular.cols)
    
    
    ### Create plots
    
    make.colour.plot(sub, "UMAP_X", "UMAP_Y", batch.col, col.type = 'factor', filename = paste0('Batches.png'))
    make.colour.plot(sub, "UMAP_X", "UMAP_Y", group.col, col.type = 'factor', filename = paste0('Groups.png'))
    make.colour.plot(sub, "UMAP_X", "UMAP_Y", donor.col,col.type = 'factor', filename = paste0('Donors.png'))
    make.multi.plot(sub, "UMAP_X", "UMAP_Y", cellular.cols, figure.title = 'Celluar markers')
    #make.multi.plot(sub, "UMAP_X", "UMAP_Y", cluster.cols, figure.title = 'Clustering markers')
    #make.multi.plot(sub, "UMAP_X", "UMAP_Y", to.align, figure.title = 'Markers to align',save.each.plot = TRUE)
    # only plotted cellular.cols bc they should be the same
    
    ### Cleanup
    
    rm(sub)
    ### UMAP of reference samples
    
    #sub <- do.subsample(ref.dat, 10000)
    sub <- ref.dat
    sub <- run.umap(sub, cellular.cols)
    
    
    ### Create plots
    
    make.colour.plot(sub, "UMAP_X", "UMAP_Y", batch.col, col.type = 'factor', filename = paste0('ref_Batches.png'))
    make.colour.plot(sub, "UMAP_X", "UMAP_Y", group.col, col.type = 'factor', filename = paste0('ref_Groups.png'))
    make.colour.plot(sub, "UMAP_X", "UMAP_Y", donor.col, col.type = 'factor', filename = paste0('ref_Donors.png'))
    make.multi.plot(sub, "UMAP_X", "UMAP_Y", cellular.cols, figure.title = 'ref_Celluar markers')
    #make.multi.plot(sub, "UMAP_X", "UMAP_Y", cluster.cols, figure.title = 'Clustering markers')
    #make.multi.plot(sub, "UMAP_X", "UMAP_Y", to.align, figure.title = 'Markers to align',save.each.plot = TRUE)
    
    ### Cleanup
    
    rm(sub)
    
    
    ##########################################################################################################
    #### Fine alignment with CytoNorm
    ##########################################################################################################
    setwd(PrimaryDirectory)
    setwd(Output2Directory)
    dir.create("Output 2.2 - fine alignment")
    setwd("Output 2.2 - fine alignment")
    fine.dir <- "Output 2.2 - fine alignment"
    
  
    

    
    ref.dat <- do.filter(cell.dat, sample.col,ref.ctrls)
    unique(ref.dat[[sample.col]])
    unique(ref.dat[[batch.col]])
    
    #ref.dat
    ### Settings
    ### Examine clustering and generate plots
    
    
    
    dir.create("A - Pre-Align ref data")
    setwd("A - Pre-Align ref data")
    
    #View(cell.dat)
    
    
    
    ### Prep FlowSOM and perform clustering using ref.dat
    
    cytnrm <- amshaw(dat = ref.dat,
                            cellular.cols = cellular.cols,
                            cluster.cols = cluster.cols,
                            batch.col = batch.col,
                            sample.col=sample.col,
                            
                            xdim = 10,
                            ydim = 10,
                            meta.k = 5)  # 10 is the default for meta.k
    
    
    #cytnrm.sub <- do.subsample(cytnrm$dt, 10000)
    
    cytnrm.sub <- run.umap(cytnrm$dt, use.cols=cluster.cols)
    
    make.colour.plot(cytnrm.sub, "UMAP_X", "UMAP_Y", "File", col.type = 'factor', filename = "Reference - batches.png")
    
    make.colour.plot(cytnrm.sub, "UMAP_X", "UMAP_Y", "prep.fsom.metacluster", col.type = 'factor', add.label = TRUE, filename = "Reference data - metaclusters.png")
    
    
    #cytnrm$files
    #cytnrm$file.nums
    #see which file number corresponds to each batch 
    rm(cytnrm.sub)
    
    ### Train the alignment conversions in the 'align.model' object
    
    setwd(PrimaryDirectory)
    setwd(Output2Directory)
    setwd(fine.dir)
    
    cytnrm <- train.cytonorm( 
      model <- cytnrm,
      align.cols = cellular.cols,
    )
    
    
    
    ### Run cytonorm
    
    cell.dat <- run.cytonorm(dat = cell.dat,
                             model = cytnrm,
                             batch.col = batch.col,
    )
    
    aligned.cols <- paste0(cellular.cols, '_aligned')
    
    ### Examine results
    setwd(PrimaryDirectory)
    setwd(Output2Directory)
    setwd(fine.dir)
    dir.create("B - ref cytonorm results")
    setwd("B - ref cytonorm results")
    ref.sub <- do.filter(cell.dat, sample.col, ref.ctrls)
    
    ref.sub <- run.umap(ref.sub, use.cols = aligned.cols)
    make.colour.plot(ref.sub, 'UMAP_X', 'UMAP_Y', batch.col, 'factor')
    
    make.colour.plot(ref.sub, 'UMAP_X', 'UMAP_Y', 'Alignment_MC_aligned', 'factor', add.label = TRUE)
    
    make.colour.plot(ref.sub, 'UMAP_X', 'UMAP_Y', group.col, 'factor')
    make.colour.plot(ref.sub, 'UMAP_X', 'UMAP_Y', sample.col, 'factor')
    
    
    ###### if ref controls clustering look suitable plot the whole data
    setwd(PrimaryDirectory)
    setwd(Output2Directory)
    setwd(fine.dir)
    dir.create("C - all aligned results")
    setwd("C - all aligned results")
    aligned.sub <- subset(cell.dat,Group!="CD4 BC") # excluding the bulk CD4 which are used for batch control
    aligned.sub <- run.umap(aligned.sub, use.cols = aligned.cols)
    
    make.colour.plot(aligned.sub, "UMAP_X", "UMAP_Y", batch.col, col.type = 'factor', filename = "Aligned batches.png")
    make.colour.plot(aligned.sub, "UMAP_X", "UMAP_Y", group.col, col.type = 'factor', filename = "Aligned groups.png")
    make.colour.plot(aligned.sub, "UMAP_X", "UMAP_Y", donor.col,col.type = 'factor', filename = "Aligned Donors.png")
    make.colour.plot(aligned.sub, "UMAP_X", "UMAP_Y", sample.col,col.type = 'factor', filename = "Aligned Samples.png")
    make.colour.plot(aligned.sub, "UMAP_X", "UMAP_Y", "Alignment_MC_aligned", col.type = 'factor', add.label = TRUE, filename = "Aligned metaclusters.png")
    
    #make.multi.plot(aligned.sub, "UMAP_X", "UMAP_Y", paste0(cellular.cols, crs.append), figure.title = "Target - markers - coarse")
    make.multi.plot(aligned.sub, "UMAP_X", "UMAP_Y", aligned.cols, figure.title = "Aligned  markers - fine",save.each.plot = TRUE)
    
    #make.multi.plot(aligned.sub, "UMAP_X", "UMAP_Y", cellular.cols, figure.title = "Target - markers - raw")
    
    rm(aligned.sub)
    
    
    
    ### Finalsie data
    
    cell.dat <- subset(cell.dat,Group!="CD4 BC" ) # excluding the bulk CD4 which are used for batch control
    
    
    
    setwd(PrimaryDirectory)
    setwd(Output2Directory)
    setwd(fine.dir)
    dir.create("D - Fine aligned data")
    setwd("D - Fine aligned data")
    
    fwrite(cell.dat, "cell.dat_allAligned.csv")
    
    write.files(cell.dat, 
                file.prefix = "Fine_aligned", 
                divide.by = sample.col, 
                write.csv = FALSE, 
                write.fcs = TRUE)
    
    ##########################################################################################################
    #### Save session info
    ##########################################################################################################
    setwd(PrimaryDirectory)
    setwd(Output2Directory)
    dir.create("Output - info", showWarnings = FALSE)
    setwd("Output - info")
    
    
    ### save RData
    save.image("FineAligned.RData")
    ### Save session info to disk
    setwd(PrimaryDirectory)
    setwd(Output2Directory)
    dir.create("Output - info", showWarnings = FALSE)
    setwd("Output - info")
    
    sink(file = "session_info.txt", append=TRUE, split=FALSE, type = c("output", "message"))
    session_info()
    sink()
    
    
  }
  
  
  
  
  
  ## clustering
  ##########################################################################################################
  #### Setup preferences
  ##########################################################################################################
  setwd(PrimaryDirectory)
  setwd(Output3Directory)
  dir.create(paste0("Output 3.1 - clustered_", phenok))
  setwd(paste0("Output 3.1 - clustered_", phenok))
  
  ### Run clustering
  ### Sample preferences
  
  sample.col <- "Sample"
  group.col <- "Group"
  donor.col <- "Donor"
  
  
  ### Clustering preferences
  
  ## Cellular cols
  if (do.batchAlign != FALSE){
    cellular.cols <- names(cell.dat)[(markerLength*2+8):(markerLength*3+7)]
  } else{
    cellular.cols <- names(cell.dat)[(markerLength+8):(markerLength*2+7)]}
  
  
  ## Columns for clustering
  # we're using all of our markers for clustering
  
  
  clustering.cols <- cellular.cols
  
  
  
  ##########################################################################################################
  #### Run clustering and dimensionality reduction
  ##########################################################################################################
  
  if (do.batchAlign != FALSE){
    data_fastPG <- as.matrix(cell.dat[,(markerLength*2+8):(markerLength*3+7)])
  }else{
    data_fastPG <- as.matrix(cell.dat[,(markerLength+8):(markerLength*2+7)])}
  #install.packages("tictoc")
  library(tictoc)
  tic("fastPG")
  output_fastPG <- FastPG::fastCluster( data= data_fastPG, k= phenok, num_threads=10 )
  cell.dat$fastPG_Clusters<- output_fastPG[[2]]
  fwrite(cell.dat, "Clustered.csv")
  toc()
  
  message("Custering completed. Running dimension reduction for visualization next.")
  ### Run DR
  
  
  #cell.sub <- do.subsample(cell.dat, sub.targets)
  cell.sub <- cell.dat
  cell.sub <- run.umap(cell.sub, clustering.cols)
  
  fwrite(cell.sub, "RD.sub.csv")
  message("Dimension reduction completed. Plotting next.")
  ### Make expression heatmap, using default normalizing by range or cutomized normalizing by z trasform
  
  exp <- do.aggregate(cell.dat, cellular.cols, by = "fastPG_Clusters")
  if (flowType == "cytof") {
    make.cytofheatmap(exp,"fastPG_Clusters", plot.cols = cellular.cols,normalise=FALSE,standard.colours = "rev(RdBu)")
  } else {make.pheatmap(exp,"fastPG_Clusters", plot.cols = cellular.cols,normalise=TRUE)
  }
  
  
  #heatmap pdf
  if (flowType == "cytof") {
    make.cytofheatmap(exp,"fastPG_Clusters", plot.cols = cellular.cols,normalise=FALSE,standard.colours = "rev(RdBu)",file.name="heatmap.pdf")
  } else {make.pheatmap(exp,"fastPG_Clusters", plot.cols = cellular.cols,normalise=TRUE,file.name="heatmap.pdf")
  }
  
  #exp <- do.aggregate(cell.dat, cellular.cols, by = "fastPG_Clusters",func="mean")
  #make.cytofheatmap(exp,"fastPG_Clusters", plot.cols = cellular.cols,normalise=FALSE,standard.colours = "rev(RdBu)")
  
  #make.pheatmap(exp,"fastPG_Clusters", plot.cols = cellular.cols,normalise=FALSE)
  ### Make expression plots
  
  make.colour.plot(cell.sub, "UMAP_X", "UMAP_Y", "fastPG_Clusters", col.type = 'factor', add.label = TRUE)
  
  
  if (do.plot!=FALSE){
    marker <- as.numeric(unlist(cell.dat[[i]]))
    percentile <- ecdf(marker)(0)
    make.colour.plot(cell.sub, "UMAP_X", "UMAP_Y",col.min.threshold = percentile, group.col, col.type = 'factor')
    
    make.colour.plot(cell.sub, "UMAP_X", "UMAP_Y",col.min.threshold = percentile, sample.col, col.type = 'factor')
    make.colour.plot(cell.sub, "UMAP_X", "UMAP_Y",col.min.threshold = percentile, donor.col, col.type = 'factor')
    
    make.multi.plot(cell.sub, "UMAP_X", "UMAP_Y", "fastPG_Clusters", group.col, col.type = 'factor',save.each.plot = TRUE)
    make.multi.plot(cell.sub, "UMAP_X", "UMAP_Y",col.min.threshold = percentile, cellular.cols,save.each.plot = TRUE)
    make.multi.plot(cell.sub, "UMAP_X", "UMAP_Y", "fastPG_Clusters", sample.col, col.type = 'factor',save.each.plot = TRUE)
    
    for(i in cellular.cols){
      marker <- as.numeric(unlist(cell.dat[[i]]))
      percentile <- ecdf(marker)(0)
      #percentile1 <- ecdf(marker)(3)
      make.multi.plot(cell.sub, "UMAP_X", "UMAP_Y", i, group.col, col.min.threshold = percentile, figure.title = paste0('Multiplot - ', i, '.png'))
    }
    
    
    
    
    
    ##save pdf figures
    
    make.colour.plot(cell.sub, "UMAP_X", "UMAP_Y", "fastPG_Clusters", col.type = 'factor', add.label = FALSE,filename = "PG_cluster UMAP.pdf")
    
    make.colour.plot(cell.sub, "UMAP_X", "UMAP_Y", group.col, col.type = 'factor',filename = "PG_group UMAP.pdf")
    
    make.colour.plot(cell.sub, "UMAP_X", "UMAP_Y", sample.col, col.type = 'factor',filename = "PG_sample UMAP.pdf")
    
    make.colour.plot(cell.sub, "UMAP_X", "UMAP_Y", donor.col, col.type = 'factor',filename = "PG_donor UMAP.pdf")
    
    
    
    
    #save cellular plots in pdf
    
    for(i in cellular.cols){
      marker <- as.numeric(unlist(cell.dat[[i]]))
      percentile <- ecdf(marker)(0)
      make.colour.plot(cell.sub, "UMAP_X", "UMAP_Y", col.min.threshold = percentile, i,filename = paste0('intensity UMAP - ', i, '.pdf'))
      make.colour.plot(cell.sub, "UMAP_X", "UMAP_Y", col.min.threshold = percentile, i,filename = paste0('intensity UMAP - ', i, '.png'))
    }
    
    #save group clusters in pdf
    group.names <- unique(cell.dat$Group)
    for (j in 1:length(group.names)){
      Idx <-which(cell.sub$Group==group.names[j])
      cell.group <- cell.sub[Idx,]
      make.colour.plot(cell.group, "UMAP_X", "UMAP_Y","fastPG_Clusters",col.type = 'factor',filename = paste0('group PG clusters_', group.names[j], '.pdf'))
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    #skipped annotation
    
    group.names <- unique(cell.dat$Group)
    for (j in 1:length(group.names)){
      cell.group <- cell.sub[cell.sub$Group==group.names[j],]
      make.multi.plot(cell.group, "UMAP_X", "UMAP_Y","fastPG_Clusters",sample.col, col.type = 'factor',figure.title = paste0('Multiplot - group ', group.names[j], '.png'))
    }}
  
  
  ##########################################################################################################
  #### Write summary data
  ##########################################################################################################
  if (do.summary!=FALSE){
    
    dyn.cols <- cellular.cols
    
    
    ### Setup cell count data
    counts <- meta.dat[,c(sample.col, 'Cells per sample'), with = FALSE]
    
    
    ### Create summary tables
    
    sum.dat <- create.sumtable(dat = cell.dat,
                               sample.col = sample.col,
                               #pop.col = "Population",
                               pop.col = "fastPG_Clusters",
                               use.cols = dyn.cols,
                               annot.cols = c(group.col),
                               counts = counts
                               #perc.pos = perc.pos
    )
    
    
    ### Write summary data
    
    fwrite(sum.dat, paste0("sum.dat_fastPGk=",phenok,".csv"))
  }
  ## added: save clustering results cell.dat as fcs files
  
  
  
  
  #save R data
  save.image(paste0("Step3_fastPGk=",phenok,"_.RData"))
  #only export fine aligned/arcsinh transformed data columns
  if (do.fcsExport!=FALSE){
    colLength <- length(names(cell.sub))
    cell.export <- cell.sub[,(markerLength+3):colLength]
    
    save.image(paste0("Step3_fastPGk=",phenok,"_.RData"))
    setwd(PrimaryDirectory)
    setwd(Output3Directory)
    dir.create(paste0("Output 3.6 - fcs files_fastPGk=",phenok))
    setwd(paste0("Output 3.6 - fcs files_fastPGk=",phenok))
    write.files(cell.export,
                file.prefix = "Clustered",
                divide.by = NULL,
                write.csv = FALSE,
                write.fcs = TRUE)
    
    
    write.files(cell.export,
                file.prefix = "Clustered",
                divide.by = "fastPG_Clusters",
                write.csv = TRUE,
                write.fcs = TRUE)
    
    write.files(cell.export,
                file.prefix = "Clustered",
                divide.by = group.col,
                write.csv = FALSE,
                write.fcs = TRUE)
    
    write.files(cell.export,
                file.prefix = "Clustered",
                divide.by = sample.col,
                write.csv = FALSE,
                write.fcs = TRUE)
  
   
  }
}

