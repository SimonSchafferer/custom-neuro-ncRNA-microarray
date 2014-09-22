#'  Perform an analysis on the Second 2.0 chip
#' @param analysisTemplate Analysis template (chip version)
#' @param gprDir directory of gpr file
#' @param analysisPath path were the html and quality control files should be written
#' @param normMethod normalization method one wants to use (quantile or vsn)
#' @param ndups number of duplicated spots on the chip
#' @param spacing spacing between the duplicated spots
#' @param minDetectLim minimum detection limit obtained from self-self hybridization experiment
#' @param subsetByManualControl if data should be subsetted by spike in control (Chip 1 -> no spike in controls)
#' @param subsetByPMvsMM if data should be subsetted by PM > MM condition
#' @param orderingVar By which variable should spot ordering occur
#' @param addInterestCand Candidate IDs that are added as labels in the plots and in the table
#' @return NULL
#' @export
performAnalysis2_0 = function(gprDir = file.path(getwd(), "../"), analysisPath = file.path( getwd() ), normMethod = "quantile", 
                              ndups = 4, spacing = 1,  minDetectLim = 0.15, subsetByManualControl = TRUE,  subsetByPMvsMM = TRUE, 
                              orderingVar = "Name",  addInterestCand = c("900000","20285","20285_1","20285_2" ), executeKnitR=TRUE  ){    
  require(knitr)
  data(analysisTemplate2.0)
  
  message("modifying template")
  analysisTemplate = modifyAnalysis( analysisTemplate2.0, gprDir=gprDir,analysisPath=analysisPath,normMethod=normMethod, ndups=ndups,
                                       spacing=spacing,minDetectLim=minDetectLim,subsetByManualControl=subsetByManualControl,
                                       subsetByPMvsMM=subsetByPMvsMM,orderingVar=orderingVar,addInterestCand=addInterestCand )
  writeLines(analysisTemplate, "analysis.Rmd")
  
  if(executeKnitR){
    message("execute template ...")
    knit2html(file.path(analysisPath,"analysis.Rmd"), envir = globalenv())     
  } else{
    message("Rmd file was not executed, custom execution by: ")
    message( paste0("knit2html('",file.path(analysisPath,"analysis.Rmd"),"', envir = globalenv())") )
  }
  message("... finished!")
}

#'  Perform an analysis on the First 1.0 chip
#' @param analysisTemplate Analysis template (chip version)
#' @param gprDir directory of gpr file
#' @param analysisPath path were the html and quality control files should be written
#' @param normMethod normalization method one wants to use (quantile or vsn)
#' @param ndups number of duplicated spots on the chip
#' @param spacing spacing between the duplicated spots
#' @param minDetectLim minimum detection limit obtained from self-self hybridization experiment
#' @param subsetByManualControl if data should be subsetted by spike in control (Chip 1 -> no spike in controls)
#' @param subsetByPMvsMM if data should be subsetted by PM > MM condition
#' @param orderingVar By which variable should spot ordering occur
#' @param addInterestCand Candidate IDs that are added as labels in the plots and in the table
#' @return NULL
#' @export
performAnalysis1_0 = function(  gprDir = file.path(getwd(), "../"), analysisPath = file.path( getwd() ), normMethod = "quantile", 
                     ndups = 8, spacing = 1,  minDetectLim = 0.15, subsetByManualControl = FALSE,  subsetByPMvsMM = TRUE, 
                     orderingVar = "Name",  addInterestCand = c("900000","20285","20285_1","20285_2" ), executeKnitR=TRUE  ){
  require(knitr)
  data( analysisTemplate1.0 )  
  message("modifying template")
  
   analysisTemplate = modifyAnalysis( analysisTemplate1.0, gprDir=gprDir,analysisPath=analysisPath,normMethod=normMethod, ndups=ndups,
                                      spacing=spacing,minDetectLim=minDetectLim,subsetByManualControl=subsetByManualControl,
                                      subsetByPMvsMM=subsetByPMvsMM,orderingVar=orderingVar,addInterestCand=addInterestCand )
   
   writeLines(analysisTemplate, "analysis.Rmd")
  
   if(executeKnitR){
     message("execute template ...")
     knit2html(file.path(analysisPath,"analysis.Rmd"), envir = globalenv())     
   } else{
     message("Rmd file was not executed, custom execution by: ")
     message( paste0("knit2html('",file.path(analysisPath,"analysis.Rmd"),"', envir = globalenv())") )
   }
  message("... finished!")
}

#'  Perform an analysis on the 2.1 chip
#' @param analysisTemplate Analysis template (chip version)
#' @param gprDir directory of gpr file
#' @param analysisPath path were the html and quality control files should be written
#' @param normMethod normalization method one wants to use (quantile or vsn)
#' @param ndups number of duplicated spots on the chip
#' @param spacing spacing between the duplicated spots
#' @param minDetectLim minimum detection limit obtained from self-self hybridization experiment
#' @param subsetByManualControl if data should be subsetted by spike in control (Chip 1 -> no spike in controls)
#' @param subsetByPMvsMM if data should be subsetted by PM > MM condition
#' @param orderingVar By which variable should spot ordering occur
#' @param addInterestCand Candidate IDs that are added as labels in the plots and in the table
#' @return NULL
#' @export
performAnalysis2_1 = function(gprDir = file.path(getwd(), "../"), analysisPath = file.path( getwd() ), normMethod = "quantile", 
                              ndups = 6, spacing = 1,  minDetectLim = 0.15, subsetByManualControl = TRUE,  subsetByPMvsMM = TRUE, 
                              orderingVar = "Name",  addInterestCand = c("900000","20285","20285_1","20285_2" ), executeKnitR=TRUE  ){    
  require(knitr)
  data(analysisTemplate2.1)
  
  message("modifying template")
  analysisTemplate = modifyAnalysis( analysisTemplate2.1, gprDir=gprDir,analysisPath=analysisPath,normMethod=normMethod, ndups=ndups,
                                     spacing=spacing,minDetectLim=minDetectLim,subsetByManualControl=subsetByManualControl,
                                     subsetByPMvsMM=subsetByPMvsMM,orderingVar=orderingVar,addInterestCand=addInterestCand )
  writeLines(analysisTemplate, "analysis.Rmd")
  
  if(executeKnitR){
    message("execute template ...")
    knit2html(file.path(analysisPath,"analysis.Rmd"), envir = globalenv())     
  } else{
    message("Rmd file was not executed, custom execution by: ")
    message( paste0("knit2html('",file.path(analysisPath,"analysis.Rmd"),"', envir = globalenv())") )
  }
  message("... finished!")
}

#'  Modifies the analysis Template
#' @param analysisTemplate Analysis template (chip version)
#' @param gprDir directory of gpr file
#' @param analysisPath path were the html and quality control files should be written
#' @param normMethod normalization method one wants to use (quantile or vsn)
#' @param ndups number of duplicated spots on the chip
#' @param spacing spacing between the duplicated spots
#' @param minDetectLim minimum detection limit obtained from self-self hybridization experiment
#' @param subsetByManualControl if data should be subsetted by spike in control (Chip 1 -> no spike in controls)
#' @param subsetByPMvsMM if data should be subsetted by PM > MM condition
#' @param orderingVar By which variable should spot ordering occur
#' @param addInterestCand Candidate IDs that are added as labels in the plots and in the table
#' @return NULL
#' @export
modifyAnalysis = function( analysisTemplate, gprDir=gprDir,analysisPath=analysisPath,normMethod=normMethod, ndups=ndups,
                spacing=spacing,minDetectLim=minDetectLim,subsetByManualControl=subsetByManualControl,
                subsetByPMvsMM=subsetByPMvsMM,orderingVar=orderingVar,addInterestCand=addInterestCand
){
  
  require(googleVis)
  require(ggplot2)
  require(reshape2)
  require(limma)
  require("geneplotter")
  require(rCharts)
  require(ggplot2)
  require(sqldf)
  require(genefilter)
  require("Biostrings")
  require("RCurl")
  require(arrayQualityMetrics)
  require("vsn")
  require("MmPalateMiRNA")
  require("latticeExtra")  
  
  analysisTemplate[grep("gprDir.*VARDEFINITION",analysisTemplate)] = paste0("gprDir = ", paste0("'",gprDir,"'"))
  analysisTemplate[grep("analysisPath.*VARDEFINITION",analysisTemplate)] = paste0("analysisPath = ",paste0("'",analysisPath,"'"))
  analysisTemplate[grep("normMethod.*VARDEFINITION",analysisTemplate)] = paste0("normMethod = ",paste0("'",normMethod,"'"))
  analysisTemplate[grep("ndups.*VARDEFINITION",analysisTemplate)] = paste0("ndups = ",ndups)
  analysisTemplate[grep("spacing.*VARDEFINITION",analysisTemplate)] = paste0("spacing = ",spacing)
  analysisTemplate[grep("minDetectLim.*VARDEFINITION",analysisTemplate)] = paste0("minDetectLim = ",minDetectLim)
  analysisTemplate[grep("subsetByManualControl.*VARDEFINITION",analysisTemplate)] = paste0("subsetByManualControl = ",subsetByManualControl)
  analysisTemplate[grep("subsetByPMvsMM.*VARDEFINITION",analysisTemplate)] = paste0("subsetByPMvsMM = ",subsetByPMvsMM)
  analysisTemplate[grep("orderingVar.*VARDEFINITION",analysisTemplate)] = paste0("orderingVar = ",paste0("'",orderingVar,"'"))
  analysisTemplate[grep("addInterestCand.*VARDEFINITION",analysisTemplate)] = paste0("addInterestCand = c(",paste0( paste0("'",addInterestCand,"'"), collapse=","), ")" )
  
  return(analysisTemplate)
}

# 
# library(CustomMicroarrayPackage)
# setwd("/media/Data/PHDStudies/ADpaper/testing/Analysis")
# performAnalysis1_0()


