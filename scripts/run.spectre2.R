run.spectre2 <- function (phenok,
                               metaFile,
                               markerFile, 
                               meta_col = c('Sample', 'Group', 'Batch', "Donor"),
                               do.plot = TRUE, 
                               do.summary = TRUE, 
                               do.batchAlign=TRUE,
                               do.fcsExport = TRUE,
                               do.Rerun = FALSE,
                               flowType,
                               coFactor,
                               plot.against,
                               ref.ctrls){
  ## spectre lib
  message("loading packages")
  library(Spectre)
  Spectre::package.check()    # Check that all required packages are installed
  Spectre::package.load()     # Load required packages
  
  library(dplyr)
  library(FastPG)
  library(CytoNorm)
  
  message("load all packages")
  
  ### Set DT threads
  getDTthreads()
  threads = getDTthreads()
  message(paste0("Thread number is ", threads))
  
  ### Set primary directory
  PrimaryDirectory <- dirname(rstudioapi::getActiveDocumentContext()$path)
  setwd(PrimaryDirectory)
  
  message(getwd())
  
  InputDirectory <- file.path(PrimaryDirectory, "data")
  MetaDirectory <- file.path(PrimaryDirectory, "metadata")
  
  ### Set output directory
  Output1Directory <- file.path(PrimaryDirectory, "Output 1 - data prep")
  Output2Directory <- file.path(PrimaryDirectory, "Output 2 - batch alignment")
  Output2.2Directory <- file.path(Output2Directory, "Output 2.2 - fine alignment")
  Output3Directory <- file.path(PrimaryDirectory, "Output 3 - clustering and DR_fastPG")
  
  ### Import data
  message("data importing...")
  
  setwd(InputDirectory)
  list.files(InputDirectory, ".fcs")
  
  data.list <- read.cytofFiles(file.loc = InputDirectory,
                               file.type = ".fcs",
                               do.embed.file.names = TRUE)
  message("data imported")
  #if cytof files, truncate max range in flowcore read.fcs is set to false
  
  # markers with column Channel.name and markers please remove strange characters 
  markers <- read.csv(markerFile) #! match meta file 
  markers[,1] = stringr::str_sub(markers[,1], start = 1,end = -2)
  markerLength <- as.numeric(length((markers[,2])))
  colnames<- append(markers[,2], c("FileName","FileNo"))
  
  message("metadata imported")
  
  data.list = lapply(data.list, rename_columns,markers)
  
  #check <- do.list.summary(data.list2)
  ### Merge data
  
  cell.dat <- Spectre::do.merge.files(dat = data.list)
  
  message("data merged")
  
  ##############################################################################
  #### Add metadata
  ##############################################################################
  
  setwd(MetaDirectory)
  
  ### Read in sample metadata
  meta.dat <- fread(metaFile)
  
  ### Add sample metadata to primary data.table
  cell.dat <- do.add.cols(cell.dat, "FileName", meta.dat, "FileName", rmv.ext = TRUE)
  
  sample.col <- meta_col[1]
  group.col <- meta_col[2]
  batch.col <- meta_col[3]  
  donor.col <- meta_col[4]
  
  message("Added metadata on merged data")
  
  if (do.Rerun == FALSE){
    
    message("Rerun is FALSE. generate data")
    ############################################################################## 
    #### Arcsinh transformation
    ############################################################################## output 1 start
    setwd(PrimaryDirectory)  
    dir.create("Output 1 - data prep", showWarnings = FALSE)
    setwd(Output1Directory)
    
    ### Transformation settings
    ### Make plots of transformed columns from the subsampled data
    if (flowType == "cytof") {
      cell.dat <- do.asinh(cell.dat, use.cols = names(cell.dat)[1:markerLength], cofactor = 5)
    } else if (flowType == "aurora") {
      if (is.null(coFactor)){ 
        cell.dat <- do.asinh(cell.dat, use.cols = names(cell.dat)[1:markerLength], cofactor = 2000)
      }
      else {
        cell.dat <- do.asinh(cell.dat, use.cols = names(cell.dat)[1:markerLength], cofactor = coFactor)
      }
      
    } else if (flowType == "flow") {
      if (is.null(coFactor)){ 
        cell.dat <- do.asinh(cell.dat, use.cols = names(cell.dat)[1:markerLength], cofactor = 200)
      }
      else {
        cell.dat <- do.asinh(cell.dat, use.cols = names(cell.dat)[1:markerLength], cofactor = coFactor)
      }
    } else {
      stop("flowType must be flow, cytof or aurora!")
    }
    
    message("data normalized with asinh")
    
    dir.create("Output 1.2 - transformed plots", showWarnings = F)
    setwd("Output 1.2 - transformed plots")
    
    transf.cols <- names(cell.dat)[grepl('_asinh', names(cell.dat))]
    which(names(cell.dat) == plot.against)
    
    for(i in transf.cols){
      make.colour.plot(cell.dat, i, col.min.threshold = 0, plot.against)
    }
    
    message("Output 1.2 - transformed plots finished")
    
    ##############################################################################
    #### Write data to disk
    ##############################################################################
    setwd(Output1Directory)
    dir.create("Output 1.1 - transformed data", showWarnings = F)
    setwd("Output 1.1 - transformed data")
    
    ### Write cellular data and analysis  preferences to disk
    fwrite(cell.dat, "cell.dat.csv") # data
    
    ### Save session info to disk
    setwd(Output1Directory)
    dir.create("Output - info", showWarnings = FALSE)
    setwd("Output - info")
    
    capture.output(sessionInfo(), file = "session_info.txt")
    
    message("Data pre-processing completed")
    
    cellular.cols = transf.cols
    cluster.cols = cellular.cols
    
    ############################################################################## output 1 end
    ## batch correction  #########################################################
    ############################################################################## output 2 start
    
    if (do.batchAlign == TRUE){
      
      ##############################################################################
      #### Initial (pre-alignment) plots
      ##############################################################################
      setwd(PrimaryDirectory)
      dir.create(Output2Directory, showWarnings = FALSE)
      setwd(Output2Directory)
      
      ref.dat <- do.filter(cell.dat, use.col = sample.col, values =ref.ctrls)
      
      ### Check reference in ref.dat and cell.dat
      
      unique(ref.dat[[batch.col]])
      unique(cell.dat[[batch.col]])
      
      cat(unique(ref.dat[[sample.col]]) == ref.ctrls)
      cat(paste0("the size of reference is", dim(ref.dat)[1], " and ", dim(ref.dat[2])))
      dir.create("Output 2.1 - pre-Align plots", showWarnings = F)
      setwd("Output 2.1 - pre-Align plots")

      name.ref.group = cell.dat %>% filter(Sample %in% ref.ctrls ) %>% select(Group) %>% table %>% names %>% unique ## group name of ref
      message(paste0("reference group name is ", name.ref.group))
      ### Pre-alignment UMAP
      sub <- subset(cell.dat,Group!=name.ref.group) # excluding the bulk CD4 which are used for batch control
      sub <- run.umap(sub, cellular.cols)
      
      ### Create plots
      make.colour.plot(sub, "UMAP_X", "UMAP_Y", batch.col, col.type = 'factor', filename = paste0('Batches.png'))
      make.colour.plot(sub, "UMAP_X", "UMAP_Y", group.col, col.type = 'factor', filename = paste0('Groups.png'))
      make.colour.plot(sub, "UMAP_X", "UMAP_Y", donor.col,col.type = 'factor', filename = paste0('Donors.png'))
      
      make.multi.plot(sub, "UMAP_X", "UMAP_Y", cellular.cols, figure.title = 'Celluar markers')
      
      message("Pre-Align plots are generated")
      ### Cleanup
      rm(sub)
      
      ##############################################################################
      #### Fine alignment with CytoNorm w/o reference
      ##############################################################################
      Output2.2Directory <- file.path(Output2Directory, "Output 2.2 - fine alignment")
      dir.create(Output2.2Directory, showWarnings = F)
      setwd(Output2.2Directory)
      
      dir.create("A - Pre-Align ref data", showWarnings = F)
      setwd("A - Pre-Align ref data")
      
      cytnrm = run_with_fallback(expr_primary = 
        quote(prep.cytonorm(dat = ref.dat,
                            cellular.cols = cellular.cols,
                            cluster.cols = cluster.cols,
                            batch.col = batch.col,
                            sample.col=sample.col,
                            xdim = 10,
                            ydim = 10,
                            meta.k = 5)),
        
        expr_fallback = quote(amshaw(dat = ref.dat,
                          cellular.cols = cellular.cols,
                          cluster.cols = cluster.cols,
                          batch.col = batch.col,
                          sample.col=sample.col,
                          xdim = 10,
                          ydim = 10,
                          meta.k = 5))
        )
      
      cytnrm.sub <- cytnrm$dt
      
      cytnrm.sub <- run.umap(cytnrm$dt, use.cols=cluster.cols)
      
      make.colour.plot(cytnrm.sub, "UMAP_X", "UMAP_Y",
                       "File", col.type = 'factor', filename = "batches.png")
      
      make.colour.plot(cytnrm.sub, "UMAP_X", "UMAP_Y",
                       "prep.fsom.metacluster", col.type = 'factor',
                       add.label = TRUE, filename = "metaclusters.png")
      
      rm(cytnrm.sub)
      
      ### Train the alignment conversions in the 'align.model' object
      message("train cytonorm")
      cytnrm <- train.cytonorm(
        model <- cytnrm,
        align.cols = cellular.cols
      )
      
      saveRDS(object = cytnrm, file = "model.rds")
      
      message("run cytonorm")
      ### Run cytonorm
      cell.dat <- run.cytonorm(dat = cell.dat,
                                       model = cytnrm,
                                       batch.col = batch.col
      )
      
      # Add aligned columns suffix
      aligned.cols <- paste0(cellular.cols, "_aligned")
      
      
      ##############################################################################
      setwd(Output2.2Directory)
      dir.create("B - ref cytonorm results", showWarnings = F)
      setwd("B - ref cytonorm results")
      
      message("generate aligned UMAP plots")
      ref.sub <- do.filter(cell.dat, sample.col, ref.ctrls)
      
      ref.sub <- run.umap(ref.sub, use.cols = aligned.cols)
      make.colour.plot(ref.sub, 'UMAP_X', 'UMAP_Y', batch.col, 'factor')
      
      make.colour.plot(ref.sub, 'UMAP_X', 'UMAP_Y', 'Alignment_MC_aligned', 'factor', add.label = TRUE)
      
      make.colour.plot(ref.sub, 'UMAP_X', 'UMAP_Y', group.col, 'factor')
      make.colour.plot(ref.sub, 'UMAP_X', 'UMAP_Y', sample.col, 'factor')
      
      ###### if ref controls clustering look suitable plot the whole data
      setwd(Output2.2Directory)
      dir.create("C - all aligned results", showWarnings = F)
      setwd("C - all aligned results")

      aligned.sub <- subset(cell.dat,Group!=name.ref.group) # excluding the bulk CD4 which are used for batch control
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
      
      ### Finalsie data ############################################################
      
      cell.dat <- subset(cell.dat,Group!=name.ref.group ) # excluding the bulk CD4 which are used for batch control
      
      setwd(Output2.2Directory)
      dir.create("D - Fine aligned data")
      setwd("D - Fine aligned data")

      fwrite(cell.dat, "cell.dat_allAligned.csv")
      
      write.files(cell.dat, 
                  file.prefix = "Fine_aligned", 
                  divide.by = sample.col, 
                  write.csv = FALSE, 
                  write.fcs = TRUE)
      
      message("aligned data saved")
      
      ############################################################################
      #### Save session info
      ############################################################################
      setwd(Output2Directory)
      dir.create("Output - info", showWarnings = FALSE)
      setwd("Output - info")
      
      ### save RData
      save.image("FineAligned.RData")
      capture.output(sessionInfo(), file = "session_info.txt")
      
    }
    message("Rdata saved")
    
  } else { 
    message("load previous data")
    
    setwd(Output2.2Directory)
    setwd("D - Fine aligned data")
    cell.dat = read.table("cell.dat_allAligned.csv", sep = ",", header = T) 
    cell.dat = cell.dat %>% as.data.table()
    message("reload from previous data")
  }
  
  ########################################################################################################## from here rerun=T
  ## output 3 : clustering
  ##########################################################################################################
  #### Setup preferences
  ########################################################################################################## output3 start
  setwd(PrimaryDirectory)
  
  dir.create(Output3Directory, showWarnings = FALSE)
  setwd(Output3Directory)
  
  dir.create(paste0("Output 3.1 - clustered_", phenok), showWarnings = FALSE)
  setwd(paste0("Output 3.1 - clustered_", phenok))
  
  ### Clustering preferences
  ## Cellular cols asinh ()or asinh_aligned (batch)
  if (do.batchAlign == TRUE){ 
    cellular.cols <- names(cell.dat)[grepl('_asinh_aligned', names(cell.dat))] 
  } else { 
    cellular.cols <- transf.cols 
  } 
  
  ## Columns for clustering
  # we're using all of our markers for clustering
  clustering.cols <- cellular.cols
  
  ##########################################################################################################
  #### Run clustering and dimensionality reduction
  ##########################################################################################################
  
  data_fastPG <- as.matrix(cell.dat %>% dplyr::select(clustering.cols))
  
  output_fastPG <- FastPG::fastCluster( data= data_fastPG, k= phenok, num_threads=threads )
  cell.dat$fastPG_Clusters <- output_fastPG[[2]]
  fwrite(cell.dat, "Clustered.csv")
  
  message("Custering completed. Running dimension reduction for visualization next.")
  
  ### Run DR
  cell.sub <- cell.dat
  cell.sub <- run.umap(cell.sub, use.cols = clustering.cols)
  
  fwrite(cell.sub, "RD.sub.csv")
  message("Dimension reduction completed. Plotting next.")
  ### Make expression heatmap, using default normalizing by range or cutomized normalizing by z trasform
  
  exp <- do.aggregate(as.data.table(cell.dat), cellular.cols, by = "fastPG_Clusters")
  
  if (flowType == "cytof") { 
    make.cytofheatmap(exp,"fastPG_Clusters", plot.cols = cellular.cols,normalise=FALSE,standard.colours = "rev(RdBu)")
  } else {
    make.pheatmap(exp,"fastPG_Clusters", plot.cols = cellular.cols,normalise=T, standard.colours = "rev(RdBu)")
    pheatmap(as.matrix(exp),filename = "fastPG_Clusters_znorm.png", scale = "row", labels_row = exp$fastPG_Clusters, show_rownames = T)
  }
  message("heatmap generated")
  
  ### Make expression plots
  make.colour.plot(cell.sub, "UMAP_X", "UMAP_Y", "fastPG_Clusters", col.type = 'factor', add.label = TRUE)
  make.colour.plot(cell.sub, "UMAP_X", "UMAP_Y", batch.col, col.type = 'factor', add.label = TRUE)
  
  if (do.plot == TRUE){
    
    for(i in cellular.cols){
      
      marker <- as.numeric(unlist(cell.dat[[i]])) ### ??????
      percentile <- ecdf(marker)(0) ### ??????
      
      make.colour.plot(cell.sub, "UMAP_X", "UMAP_Y",col.min.threshold = percentile, group.col, col.type = 'factor')
      make.colour.plot(cell.sub, "UMAP_X", "UMAP_Y",col.min.threshold = percentile, sample.col, col.type = 'factor')
      make.colour.plot(cell.sub, "UMAP_X", "UMAP_Y",col.min.threshold = percentile, donor.col, col.type = 'factor')
      
      make.multi.plot(cell.sub, "UMAP_X", "UMAP_Y", "fastPG_Clusters", group.col, col.type = 'factor',save.each.plot = TRUE)
      make.multi.plot(cell.sub, "UMAP_X", "UMAP_Y",col.min.threshold = percentile, cellular.cols,save.each.plot = TRUE)
      make.multi.plot(cell.sub, "UMAP_X", "UMAP_Y", "fastPG_Clusters", sample.col, col.type = 'factor',save.each.plot = TRUE)
      
      make.colour.plot(cell.sub, "UMAP_X", "UMAP_Y", "fastPG_Clusters", col.type = 'factor', add.label = FALSE,
                       filename = "PG_cluster UMAP.pdf")
      make.colour.plot(cell.sub, "UMAP_X", "UMAP_Y", group.col, col.type = 'factor',filename = "PG_group UMAP.pdf")
      make.colour.plot(cell.sub, "UMAP_X", "UMAP_Y", sample.col, col.type = 'factor',filename = "PG_sample UMAP.pdf")
      make.colour.plot(cell.sub, "UMAP_X", "UMAP_Y", donor.col, col.type = 'factor',filename = "PG_donor UMAP.pdf")
      make.colour.plot(cell.sub, "UMAP_X", "UMAP_Y", batch.col, col.type = 'factor',filename = "PG_batch UMAP.pdf")
      
      make.multi.plot(cell.sub, "UMAP_X", "UMAP_Y", i, group.col, col.min.threshold = percentile, 
                      figure.title = paste0('Multiplot - ', i, '.png'))
      
      make.colour.plot(cell.sub, "UMAP_X", "UMAP_Y", col.min.threshold = percentile, i,
                       filename = paste0('intensity UMAP - ', i, '.pdf'))
      make.colour.plot(cell.sub, "UMAP_X", "UMAP_Y", col.min.threshold = percentile, i,
                       filename = paste0('intensity UMAP - ', i, '.png'))
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
    }
    
  }
  
  #save R data
  save.image(paste0("Step3_fastPGk=",phenok,"_.RData"))
  
  ##########################################################################################################
  #### only export fine aligned/arcsinh transformed data columns
  ##########################################################################################################
  if (do.fcsExport == TRUE){
    colLength <- length(names(cell.sub))
    cell.export <- cell.sub[,(markerLength+3):colLength]
    
    save.image(paste0("Step3_fastPGk=",phenok,"_.RData"))
    setwd(PrimaryDirectory)
    setwd(Output3Directory)
    dir.create(paste0("Output 3.6 - fcs files_fastPGk=",phenok), showWarnings = F)
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
  
  ##########################################################################################################
  #### Write summary data
  ##########################################################################################################
  if (do.summary == TRUE){
    
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
    fwrite(sum.dat, paste0("sum.dat_fastPGk=",phenok,".csv"))
  }
  setwd(PrimaryDirectory)
}
