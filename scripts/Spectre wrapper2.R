# modified 10.08.2025
#renv::init()
#renv::snapshot() 
## Flowcytometry analysis with Spectre package.
## data/ : include all .fcs files
## metadata/ : include sample.details.csv, with have sample, group, batch, and donor information.
## version: Spectre==1.3, flowCore==2.16, FlowSOM==2.13.11, dplyr==1.1.4, CytoNorm==2.0.9 stringr==1.5.2
## source

# Finds the directory where this script is located
directory <- dirname(rstudioapi::getActiveDocumentContext()$path) 
setwd(directory)  
print(paste("Current working directory:", getwd()))

source_dir = .libPaths()[1]
source(paste0(source_dir, "/../../../../../scripts/run.spectre2.R"))
source(paste0(source_dir, "/../../../../../scripts/help_functions.R"))
# source("run.spectre2.R")
# source("help_functions.R")

# User Customization ----------------------------------------------------
#define your parameters and file input location
run.spectre2(phenok=800, ################################################## adjustable variable
             metaFile="sample.details.csv",
             markerFile = "ORIGINAL MARKERS.csv", 
             meta_col = c('Sample', 'Group', 'Batch', "Donor"),
             do.plot=TRUE,
             do.summary=TRUE,
             do.batchAlign=TRUE,
             do.fcsExport=TRUE, 
             do.Rerun = FALSE, ### Rerun? set it TRUE! save your time!
             flowType = "aurora",
             coFactor =7000, ############################################## adjustable variable
             plot.against= " CD45RO_asinh",
             ref.ctrls=c("20250728_batch-control-10-G030-v11_neg",
                         "20250729_batch-control-10-G030-v11_neg",
                         "20250805_batch-control-10-G030-v11_neg",
                         "20250806_batch-control-10-G030-v11_neg",
                         "20250827_batch-control-10-G030-v11_neg",
                         "20250828_batch-control-10-G030-v11_neg"
                              ) ###################################'Sample' names of reference data in metaFile
)

# The defaults for cofactor: cytof = 5,    aurora = 2000,   flow = 200
# do.Rerun = FALSE as Default, if you rerun this analysis with different phenok and coFactor
# plot.against= "marker_asinh" ex) CD45RO
# phenok <- #put klevel 
# metaFile <- "sample.details.csv" #put file name
# meta_col <- column names of metadata Type names of Sample, Group, Batch, Donor column. CASE DEPENDENT!
# if not typed, Sample, Group, Batch, Donor are used by default
# markerFile <- "ORIGINAL MARKERS.csv" #put file n
# all the defaults for summary, batchAlign, fcsExport are TRUE
# do.plot is for plots in addition to basic heatmap and clusters

#renv::snapshot() 