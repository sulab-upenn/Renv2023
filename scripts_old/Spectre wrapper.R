# modified 5.23.2023
#renv::init()
#renv::snapshot() 

directory <- dirname(rstudioapi::getActiveDocumentContext()$path)       # Finds the directory where this script is located
setwd(directory)  
print(paste("Current working directory:", getwd()))

#Recalling all the modified versions of the functions alongside the main pipeline 09-Oct
source("run.spectre.exppositive.R")  # Contains function amshaw which is a redefined version of prep.cytonorm function
# The redefinition targets the error  'argument of length zero' in 'if (sampleWithReplacement & (nrow(f) < cFile))'
# function amshaw is loaded from file "prepcytonorm.R"
source("prepcytonorm.R") # Redefines the function prepareFlowSOM as erika1
source("prepareflowsom.R")  # Redefines the function AggregateFlowFrames as erika
source("flowaggregate.R") # changes the if statement to boolean type "if (isTRUE(sampleWithReplacement & (nrow(f) < cFile)))" in line 25

# User Customization ----------------------------------------------------

#define your parameters and file input location
run.spectre.exppositive(phenok=40,metaFile="sample.details.csv",markerFile = "ORIGINAL MARKERS.csv",
            do.plot=TRUE,do.summary=TRUE,do.batchAlign=TRUE,do.fcsExport=TRUE, flowType = "aurora",coFactor =2000,
            plot.against= "CD45RA_asinh",ref.ctrls=c("10-B027_CD3, CD4 _batch-1","10-B027_CD3, CD4 _batch-2", "10-B027_CD3, CD4 _batch-3"))
###    The defaults for cofactor: cytof = 5,    aurora = 2000,   flow = 200

#plot.against= "marker_asinh"
#phenok <- #put klevel 
#metaFile <- "sample.details.csv" #put file name
#markerFile <- "ORIGINAL MARKERS.csv" #put file n
#all the defaults are TRUE except for batch and sumamry table
#ref.ctrls need to be sample name as entered in sample.details.csv only when do.batchAlign is true
#ref.ctrls=c("5-0035_CD3, CD4 _citru batch1","5-0035_CD3, CD4 subset_batch2","5-0035_CD3, CD4 _citru batch2", "5-0035_CD3, CD4 subset_batch1")
#when there's batch alignment, CD4 grp name should better be CD4 BC 
#do.plot is for plots in addition to basic heatmap and clusters

#renv::snapshot() 
