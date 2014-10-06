#'  Loads a R object with a given name into a desired environment
#' @param name Name of the object you want to load
#' @param  filename Path to the file and name of the file that should be loaded
#' @param enir The environment in which the object should be loaded
#' @return Tthe name of the just created object
#' @export
loadLibraryAs = function( name, filename, envir=globalenv() ){
  nameTmp = load(filename)
  command = paste0("assign(\"",substitute(name),"\",",nameTmp,", envir = envir)")
  eval(parse(text=( command )))
  message(paste0(nameTmp, " was assigned to variable ", substitute(name)))
  return(name)
}

#' Calculates the mean or median value of the spots on a chip by a given Identifier (deprecated)
#' @param x An object of type EListRaw, RGList or EList
#' @param byType An identifier that is used for grouping the values e.g. ID
#' @param method The method name as character for spot averaging (currently 'mean' and 'median' are supported) this will then call colMeans or colMedians
#' @return  The averaged object
#' @export
calcMeanValueOfSpots = function(x, byType="Name", method="mean"){  
  require(limma)
  require(miscTools)
  
  if( !( is(x, "EListRaw") || is(x, "RGList") || is(x, "EList") ) ){
    stop("Only RGList or EListRaw/EList objects are supported! (Type Not Recognized)")
  }
  if( !is.null(x$R) && !is.null(x$G) ){
    signalType = c("G","R" )#currently no background signal supported
  } else if( !is.null(x$E)){
    signalType = "E" #currently no background signal supported
  } else{
    stop("Only RGList or EListRaw/EList objects are supported! (Channel Intensities missing)")
  }
  for( st in signalType){	
    x$genes[,byType] = as.factor(x$genes[,byType]) 
    xL = split( as.data.frame(x[[st]]), x$genes[,byType] )
    xL = lapply( xL, function(t){
      if(tolower(method) == "mean"){
        t = colMeans(t, na.rm=TRUE)	
      } else if(tolower(method) == "median"){
        t = colMedians(t, na.rm=TRUE)
      } 
      else{
        stop("Unknown averaging method -> only mean and median are recognized!")
      }
      t = ifelse( t < 0, 0.5, t  )#Setting a minimum value of 0.5 if the background subtraction produces negative values (for logFC calculation)
      return(t)
    })
    xL = do.call(rbind,xL)
    x = x[!duplicated( x$genes[,byType] ),]
    xL = as.data.frame(xL)
    xL$SortIndex = rownames(xL)
    xL = merge( data.frame( "SortIndex"=as.character( x$genes[,byType] ) ), xL, by="SortIndex", sort=FALSE)
    x[[st]] = xL[,-1]
  }
  return(x)
}

#' Calculates the optimal table size for displaying a data.frame as googleVis table in the html file, when using markdown
#' @param data.frame that should represent the table to display
#' @param maxRowsToDisplay Maximum number of rows that should be displayed before paging occurs
#' @param charMultiplicationFactorInPixels Size of a char in the browser in pixels
#' @param tableHeightInPixels Height of the table rows in pixels
#' @return GvisTable object that can be called with plot() in a code chunk
#' @export
createGVisTableFromDF = function( tDF, maxRowsToDisplay=20, charMultiplicationFactorInPixels = 8, tableHeightInPixels = 32  ){
  require(googleVis)
  maxRowsToDisplay = maxRowsToDisplay
  charMultiplicationFactorInPixels = charMultiplicationFactorInPixels
  if( dim(tDF)[1] < maxRowsToDisplay ){
    height = ( dim(tDF)[1] + 1 )*tableHeightInPixels
  }else {
    height = maxRowsToDisplay*tableHeightInPixels 
  }
  maxRowChars = max( c( max( sapply( colnames(tDF), nchar) ), max( apply( tDF, 1, function( x ){ max(nchar(x))}  ) ) ) )
  width = maxRowChars * charMultiplicationFactorInPixels * dim(tDF)[2]
  resL = list( "height"=height, "width"=width)
  gvt = gvisTable(tDF,options=list( page='enable', pageSize=maxRowsToDisplay, height=height, width=width ))
  return( gvt )
}

#' This function creates a ggplot object where MA plots of all Arrays are plotted, present in the MA object
#' @param MA object of type MA
#' @param subset Integer vector of the column indices defining the subset
#' @param colorVar the factor for coloring the dots in the plot (could also be NULL) (this vector needs to be present in the data.frame MA$genes) 
#' @param ncol Number of columns that is used for facetting the plot
#' @param smooth Adds a LOESS curve to the plot
#' @export
plotMAofArrays = function( MA, subset=NULL, colorVar = "probe.type", ncol=2, smooth = FALSE ){
  require(reshape2)
  require(ggplot2)
  if(!is.null(subset)){
    MA = MA[ , subset ]
  }
  dat = melt(MA$M, id.vars=colnames(MA$M)) #reshaping the logFC values
  colnames(dat) = c("Index", "Array", "logFC")
  tmp = melt(MA$A, id.vars=colnames(MA$A)) #reshaping the AverageExpression Values
  dat$AveExpr = tmp$value
  rm(tmp)
  if(!is.null(colorVar)){
    dat$colorVar = rep( MA$genes[,colorVar], dim(MA)[2] )
  }
  
  g = ggplot( dat, aes(x=AveExpr, y=logFC, color=colorVar) ) + geom_point( size = 1 ) + facet_wrap( ~ Array, ncol = ncol) + 
    labs(title="MAPlot of All Arrays", x="Intensity (A)", y="log2 Fold Change (M)") 
  if( smooth ){
    g = g + geom_smooth(aes(x = AveExpr, y = logFC), data=dat, method="loess")
  }
  return(g)
}

#' A plot is created as a combination of Boxplot and Density plot of the array intensities 
#' @param RG Object of type RGList
#' @return Already performs plotting
#' @export
plotBoxDensityFromRG = function( RG ) {
  par(mfrow=c(2,2))
  plotformula = RG$G~col(RG$G)
  boxplot(plotformula, outline=FALSE, col="forestgreen", xlab="arrays", ylab=expression(log[2]~G), main="boxplot G Values")
  #  Distributions of G values on the different arrays
  multidensity(plotformula,main="densities G Values", xlab=expression(log[2]~G))
  
  plotformula = RG$R~col(RG$R)
  boxplot(plotformula, outline=FALSE, col="red", xlab="arrays", ylab=expression(log[2]~R), main="boxplot R Values")
  #	Distributions of G values on the different arrays
  multidensity(plotformula,main="densities R Values", xlab=expression(log[2]~R))
}

#' Extracts the Net Intensity Values by the identfier: Net intensity (mean) {Cy5} and Net intensity (mean) {Cy3}
#' @param RG Object of type RGList
#' @param deleteOther If true, the values defined under RG$other are set to NULL
#' @return RG object with Net Intensity Values as R and G
#' @export
extractNetIntensityValues = function( RG, deleteOther = TRUE ){
  RG$R = RG$other$`Net intensity (mean) {Cy5}`
  RG$G = RG$other$`Net intensity (mean) {Cy3}`
  RG$Rb = NULL
  RG$Gb = NULL
  if(deleteOther){
    RG$other = NULL
  }
  return(RG)
}

getDefaultDesign = function( MA ){
  
}

#' Normalizes a RG object by a given method (eather vsn or quantile normalization)
#' When choosing quantile normalization an offset of 0.1 for zero values is defined (methods used: printtiploess for within and Aquantile for between arrays))
#' @param RG Object of type RGList
#' @param method Choose between 'quantile' or 'vsn'
#' @return Normalized RG object
#' @export
normalizeTwoColorArray = function(RG, method ="quantile" ){
  require(vsn)
  require(limma)
  if( method == "vsn" ){
    message("VSN Normalization .... ")
    RG = normalizeVSN(RG)
  } else if(method == "quantile"){
    message("Quantile Normalization .... ")
    #Quantile Normalization with printtiploess within array
    offsetF = function(y){
      ifelse(y <= 0, 0.1, y)
    } 
    RG$R=apply(RG$R,2,offsetF);
    RG$G=apply(RG$G,2,offsetF); 
    RG = normalizeWithinArrays(RG, method = "printtiploess")
    RG = normalizeBetweenArrays(RG, method="Aquantile")    
  }
  return(RG)
}

#' Order by a given name
#' @param x Object that is two dimensional (x[ordering,] is called)
#' @param by The variable to order from
#' @return An ordered object
#' @export
orderTwoColorArray = function( x, by = "Name" ){
  x[order(x$genes[,by]),]
}

#' Performs a differential expression analysis via limma
#' @param MA An object of type MAList
#' @param aveDups The number of duplicated spots in the MA object
#' @param design Design of the experiment (if null, a dye swap design is considered (1,-1,1,-1,...))
#' @param block Blocking vector
#' @param blockpairs For creation of a default blocking vector (default creation would be 1,1,2,2,...)
#' @param ndups Number of duplicated spots
#' @param spacing If spacing is not defined a sorted MA object is considered
#' @return List containing a resulting Table from eBayes, the model fit and the fit after eBayes 
#' @export
performDiffExpAnalysis = function(MA, aveDups=FALSE, design=NULL, block=NULL, blockpairs = 2, ndups = 4, spacing=1 ){
  require(limma)
  if(is.null(design)){
    design = rep(c(1,-1), length(MA$targets[,1])/2)
    message("Design is considered as Dye Swap Design...")
    print(design)
  }
  
  if(is.null(block)){
    blockn = length(MA$targets[,1])/blockpairs
    block = rep( 1:blockn, blockpairs )
    block = block[order(block)]
    message("Design is considered as biological pairs of dye swaps...")
    print(block)
  }
  
  if(aveDups){
    message("Calculating the differential Expression by averaging the duplicated spots - Have you ordered the data?")
    message("By doing this, the technical replicates can be blocked in the linear model")
    message( paste0("Comparing Cy5 versus Cy3: ", MA$targets$Cy5[1], " versus ", MA$targets$Cy3[1]) )
    
    MA = avedups(MA, ndups=ndups, spacing=spacing)
    design = design
    corCons = duplicateCorrelation(MA, design = design, block=block)$consensus
    fit <- lmFit(MA, design = design,cor=corCons, block=block)
    fitE = eBayes(fit)
    MAres = topTable(fitE, number = dim(MA)[1])  
    message(paste0( "Number of significant candidates: ",dim(MAres[MAres$adj.P.Val < 0.05,])[1] ))
    return(list( "ResultTable"=MAres, "fit"=fit, "fitEBayes"=fitE )) 
    
  } 
  else {
    message("Calculating the differential Expression by adjusting for correlation of duplicated spots in the linear model")
    message("The technical replicates are treated as biological replicates - Please check if correlation of Technical replicates is low!")
    message(paste0("Comparing Cy5 versus Cy3: ", MA$targets$Cy5[1], " versus ", MA$targets$Cy3[1]))
    
    ##########################################################
    #  Seems reasonable taking the dye swap into account and consider the duplicated spot correlation instead of the correlation between replicates
    #	Correlation between replicates reaches ~0.12 and duplicated spot correlation ~0.5
    #	Accounting for the dye effect since it could be gene specific (gc content...)
    ##########################################################
    #duplicateCorrelation fpr the dye effect should be quite negative -> therefore dye effect cannot be ignored!
    
    dyeSwapCorr = duplicateCorrelation(MA,block=block)$consensus #should be quite negative -> therefore dye effect cannot be ignored!
    dyeSwapAdjTechnCorr <- duplicateCorrelation(MA,  design =design, block=block)$consensus
    corfitDups = duplicateCorrelation(MA, ndups = ndups, spacing=spacing, design = design)$consensus
    
    message( paste0( "Correlation between Dye Swap Blocks: ", dyeSwapCorr, " \nCorrelation between technical replicates dye swap adjusted: ", dyeSwapAdjTechnCorr, " \nCorrelation between duplicated spots: ", corfitDups ) )
    message("Therefore the correlation between duplicated spots is taken into account when calculating the linear model fit!")
    fit <- lmFit(MA, cor = corfitDups, ndups = ndups, spacing=spacing, design = design)
    fitE = eBayes(fit)
    MAres = topTable(fitE, number = dim(MA)[1]/4)	
    message(paste0( "Number of significant candidates: ",dim(MAres[MAres$adj.P.Val < 0.05,])[1] ))
    return(list( "ResultTable"=MAres, "fit"=fit, "fitEBayes"=fitE)) 
    
  }
}

#' Adds the probe and chipType to the RG$genes table based on the database connection given
#' @param RG Object of type RGList
#' @param dbcon Database connection
#' @param oligoTable If no database connection is provided the oligoTable has to be provided containing the columns category,chipType,oligoTableID
#' @param anchorName Name of the Anchor spots
#' @param noSignal Name of the spots lacking a signal
#' @return The 'genes' table of the RG object provided
#' @export
addProbeTypeAndChipTypeToRGGenes = function( RG, dbcon, oligoTable, anchorName = "unknown", noSignal="noSig", refCol = "Name" ){
  #geneTable: RG$genes
  #oligoTable: from Database containing category and chipType
  
  if( missing(oligoTable) & !missing(dbcon) ){
    oligoTable = dbGetQuery(db,paste("SELECT * FROM oligo o JOIN chipType ct ON o.chipType_chipTypeTableID = ct.chipTypeTableID;"))
    oligoTable = oligoTable[,c("oligoTableID","category","chipType")] 
    oligoTable$oligoTableID = paste( "e", oligoTable$oligoTableID, sep="")
  } else{
    stop("Please provide database connection to oligoTable or oligoTable!")
  }
  
  if( sum( c("category","chipType","oligoTableID") %in% colnames(oligoTable) ) != length(colnames(oligoTable)) ){
    stop("Please provide column: category,chipType,oligoTableID in oligoTable")
  }
  if( is.null(RG$genes) ){
    stop("Please provide list object with a data.frame slot named 'genes'")
  }
  
  RG$genes$probe.type = c("")
  RG$genes$chipType = c("")
  RG$genes[RG$genes[,refCol] == anchorName,]$probe.type = "AnchorSpot"
  RG$genes[RG$genes[,refCol] == anchorName,]$chipType = "AnchorSpot"
  RG$genes[RG$genes[,refCol] == noSignal,]$probe.type = "NoSignal"
  RG$genes[RG$genes[,refCol] == noSignal,]$chipType = "NoSignal"
  
  genesTable = RG$genes
  
  tableid = oligoTable$oligoTableID
  category = oligoTable$category
  chipType = oligoTable$chipType
  
  genesID = genesTable$ID
  genesProbeType = genesTable$probe.type
  genesChipType = genesTable$chipType
  
  for( i in 1:length(tableid) ){
    pos = which(genesID == tableid[i])
    genesProbeType[pos] = category[i]
    genesChipType[pos] = chipType[i]
  }
  
  genesTable$probe.type = genesProbeType
  genesTable$chipType = genesChipType
  
  return(genesTable)
  
}

#' Reads and processes the gpr files of two channel arrays. Attention changes commas to dots of the names in the first chip
#' @param currDir directory were the file is read in
#' @param columns custom column Names
#' @param customGridNames ID are taken from gridName Table if set to TRUE
#' @param specificationFile Text file in the directory of the gpr files were the gpr file is specified
#' @param source Source of the Laser Scanner 
#' @param na.strings The Name of the spots without informatio or data (as defined by the scanning software)
#' @return List containing an annotated data frame a target data.frame and the RG object
#' @export
readAndProcessChipDataFirstChip <- function( currDir, columns=NULL, customGridNames=FALSE, specificationFile="samplesInfo.txt", source="genepix", na.strings="<NO DATA>" ){
  require(limma)
  if(is.null(columns)){
    columns <- columns <- list( G="Raw intensity (mean) {Cy3}", R = "Raw intensity (mean) {Cy5}", Gb = "Background (mean) {Cy3}", Rb = "Background (mean) {Cy5}" )
  }
  if(!exists(as.character(substitute(currDir)))){
    currDir <- getwd()
  }
  otherColumns = list( Gnet = "Net intensity (mean) {Cy3}", Rnet = "Net intensity (mean) {Cy5}", 
                       GNetSat = "Net intensity (sat.pr.) {Cy3}", RNetSat="Net intensity (sat.pr.) {Cy5}",
                       GBgMed = "Background (med) {Cy3}", RBgMed = "Background (med) {Cy5}", 
                       GNet1sd = "Net intensity (pr.1sd) {Cy3}", RNet1sd = "Net intensity (pr.1sd) {Cy5}" )
  adf = read.AnnotatedDataFrame(specificationFile, path=currDir)
  targets = pData(adf)
  targets$FileName = row.names(targets)
  RG = read.maimages(targets, path=currDir, source=source,  columns=columns, na.strings=na.strings,other.columns = otherColumns) #reads in all data files -> $R, $G, $Rb, $Gb,    $targets, $genes (Rb red background
  if(customGridNames){
    RG$genes$ID <- gridNameTable$ID
  }
  #changing commas to dots
  RG$genes$ID <- sub( ",",".",RG$genes$ID )
  
  controlSpots <- grep( "Cy", RG$genes$ID, ignore.case=TRUE )
  emptySpots <- which( RG$genes$ID == "" | RG$genes$ID == "unknown" )
  mmOligos <-  grep( "_MM", RG$genes$ID, ignore.case=TRUE )
  pmOligos <- setdiff(1:length(RG$genes$ID), c(controlSpots, emptySpots, mmOligos) )
  RG$genes$probe.type <- c("")
  RG$genes$probe.type[controlSpots] <- "Control"
  RG$genes$probe.type[emptySpots] <- "Empty"
  RG$genes$probe.type[mmOligos] <- "MM Probes"
  RG$genes$probe.type[pmOligos] <- "PM Probes"
  chipDataL <- list("adf"=adf, "targets"=targets, "RG"=RG)
  print(paste("Data succesfully read in and processed from:", currDir))
  return(chipDataL)
}

#' Reads and processes the gpr files of two channel arrays For The NRC5 gal file or Chip 2.1.
#' @param currDir directory were the file is read in
#' @param specificationFile Text file in the directory of the gpr files were the gpr file is specified
#' @param source Source of the Laser Scanner 
#' @param na.strings The Name of the spots without informatio or data (as defined by the scanning software)
#' @param columns custom column Names
#' @param otherColumns Additional columns
#' @param GALFile the GAL file of the chip
#' @return List containing an annotated data frame a target data.frame and the RG object
#' @export
readAndProcessChipData2_1 <- function( currDir, specificationFile="samplesInfo.txt", source="genepix", na.strings="<NO DATA>", columns, otherColumns, GALFile, ... ){  
  require(limma)
  #actually it is not a genepix software that creates the data but the input is modified to have the same composition as a genepix table
  if(missing(columns)){
    columns <- list( G="Raw intensity (mean) {Cy3}", R = "Raw intensity (mean) {Cy5}", Gb = "Background (mean) {Cy3}", Rb = "Background (mean) {Cy5}"     )#Cy3=A   Cy5=B  	
  }
  if(missing(otherColumns)){
    otherColumns = list( Gnet = "Net intensity (mean) {Cy3}", Rnet = "Net intensity (mean) {Cy5}", 
                         GNetSat = "Net intensity (sat.pr.) {Cy3}", RNetSat="Net intensity (sat.pr.) {Cy5}",
                         GBgMed = "Background (med) {Cy3}", RBgMed = "Background (med) {Cy5}", 
                         GNet1sd = "Net intensity (pr.1sd) {Cy3}", RNet1sd = "Net intensity (pr.1sd) {Cy5}" )
  }
  if(missing(GALFile)){
    warning("GAL file not specified! Supposing GPR annotation as valid\n")
  }
  
  if(!exists(as.character(substitute(currDir)))){
    currDir <- getwd()
  }
  
  adf = read.AnnotatedDataFrame(specificationFile, path=currDir)
  targets = pData(adf)
  targets$FileName = row.names(targets)
  
  RG = read.maimages(targets, path=currDir, source=source,  columns=columns, na.strings=na.strings, other.columns = otherColumns) #reads in all data files -> $R, $G, $Rb, $Gb, $targets, $genes (Rb red background
  
  bgspots = which(is.na(RG$genes$Block))
  anchorSpots = which(RG$genes$Name == "unknown")
  
  cat(paste(dim(RG)[1],"spots are on the chip with",dim(RG)[2]," Samples\n"))
  
  table( apply( RG$G,2, is.na) )
  table( apply( RG$R,2, is.na) )
  naVals = FALSE
  if( sum(rowSums( apply( RG$G,2, is.na) ) > 0) > 0 ){
    naVals = TRUE
    warning( "!!!!!!!!!!! There are NA value spots in the Green Channel !!!!!!!!!!!" )
  }
  
  if( sum(rowSums( apply( RG$R,2, is.na) ) > 0) > 0 ){
    naVals = TRUE
    warning( "!!!!!!!!!!! There are NA value spots in the Red Channel !!!!!!!!!!!" )
  }
  
  if( !naVals ){
    if( min(RG$G) == 0 | min(RG$R) == 0 )
      warning( "!!!!!!!!!!! There are zero value spots in the data !!!!!!!!!!!" )
  }
  
  cat(paste(length(bgspots),"Background spots are present on the chip... \n"))
  cat(paste(length(anchorSpots),"Anchor spots are present on the chip... \n"))
  
  if(!missing(GALFile)){	
    cat("\nComparing GAL file with GPR File...\n")
    RGGal = RG[-bgspots,]
    
    if( dim(GALFile)[1] == dim(RGGal$genes)[1] ){
      
      rg = RGGal$genes[,c("Block","Row","Column","ID")]
      rg = rg[ order( rg$Block,rg$Row,rg$Column,rg$ID ),]
      gal = GALFile[,c("Block","Row","Column","ID")]
      gal = gal[order( gal$Block,gal$Row,gal$Column,rg$ID ),]
      
      if( table(rg == gal)/dim(gal)[2] == dim(gal)[1]){
        cat("GAL File is the same as the GPR File\n")			
      }
      
    } else{
      cat("GAL File is different from the GPR File, taking GPR Annotation...\n")			
    }
  }
  
  RG$genes$ID = ifelse( RG$genes$ID == "","noSig",RG$genes$ID )
  RG$genes$Name = ifelse( RG$genes$Name == "","noSig",RG$genes$Name )
  
  chipDataL <- list("adf"=adf, "targets"=targets, "RG"=RG)
  
  print(paste("Data succesfully read in and processed from:", currDir))
  return(chipDataL)
}


#' Reads and processes the gpr files of two channel arrays. Attention changes commas to dots of the names in the first chip
#' @param currDir directory were the file is read in
#' @param specificationFile Text file in the directory of the gpr files were the gpr file is specified
#' @param source Source of the Laser Scanner 
#' @param na.strings The Name of the spots without informatio or data (as defined by the scanning software)
#' @param columns custom column Names
#' @param otherColumns Additional columns
#' @param GALFile the GAL file of the chip
#' @return List containing an annotated data frame a target data.frame and the RG object
#' @export
readAndProcessChipData <- function( currDir, specificationFile="samplesInfo.txt", source="genepix", na.strings="<NO DATA>", columns, otherColumns, GALFile, ... ){  
  require(limma)
  #actually it is not a genepix software that creates the data but the input is modified to have the same composition as a genepix table
  if(missing(columns)){
    columns <- list( G="Raw intensity (mean) {Cy3}", R = "Raw intensity (mean) {Cy5}", Gb = "Background (mean) {Cy3}", Rb = "Background (mean) {Cy5}" )#Cy3=A 	Cy5=B		
  }
  if(missing(otherColumns)){
    otherColumns = list( Gnet = "Net intensity (mean) {Cy3}", Rnet = "Net intensity (mean) {Cy5}", 
                         GNetSat = "Net intensity (sat.pr.) {Cy3}", RNetSat="Net intensity (sat.pr.) {Cy5}",
                         GBgMed = "Background (med) {Cy3}", RBgMed = "Background (med) {Cy5}", 
                         GNet1sd = "Net intensity (pr.1sd) {Cy3}", RNet1sd = "Net intensity (pr.1sd) {Cy5}" )
  }
  if(missing(GALFile)){
    warning("GAL file not specified! Supposing GPR annotation as valid\n")
  }
  
  if(!exists(as.character(substitute(currDir)))){
    currDir <- getwd()
  }
  
  adf = read.AnnotatedDataFrame(specificationFile, path=currDir)
  targets = pData(adf)
  targets$FileName = row.names(targets)
  
  RG = read.maimages(targets, path=currDir, source=source,  columns=columns, na.strings=na.strings, other.columns = otherColumns) #reads in all data files -> $R, $G, $Rb, $Gb, $targets, $genes (Rb red background
  
  bgspots = which(is.na(RG$genes$Block))
  anchorSpots = which(RG$genes$Name == "unknown")
  
  cat(paste(dim(RG)[1],"spots are on the chip with",dim(RG)[2]," Samples\n"))
  #	cat(paste("Deleting",length(bgspots),"the Background spots (Block equals NA) ... \n"))
  #	if( length(bgspots) != 0){ RG = RG[-bgspots,] } else{ RG = RG}
  #	if( length(anchorSpots) != 0){ RG = RG[-anchorSpots,] } else{ RG = RG}
  #	cat(paste(dim(RG)[1],"spots are left on the chip with",dim(RG)[2]," Samples\n"))
  
  table( apply( RG$G,2, is.na) )
  table( apply( RG$R,2, is.na) )
  naVals = FALSE
  if( sum(rowSums( apply( RG$G,2, is.na) ) > 0) > 0 ){
    naVals = TRUE
    warning( "!!!!!!!!!!! There are NA value spots in the Green Channel !!!!!!!!!!!" )
  }
  
  if( sum(rowSums( apply( RG$R,2, is.na) ) > 0) > 0 ){
    naVals = TRUE
    warning( "!!!!!!!!!!! There are NA value spots in the Red Channel !!!!!!!!!!!" )
  }
  
  if( !naVals ){
    if( min(RG$G) == 0 | min(RG$R) == 0 )
      warning( "!!!!!!!!!!! There are zero value spots in the data !!!!!!!!!!!" )
  }
  
  cat(paste(length(bgspots),"Background spots are present on the chip... \n"))
  cat(paste(length(anchorSpots),"Anchor spots are present on the chip... \n"))
  
  if(!missing(GALFile)){	
    cat("\nComparing GAL file with GPR File...\n")
    RGGal = RG[-bgspots,]
    
    if( dim(GALFile)[1] == dim(RGGal$genes)[1] ){
      
      rg = RGGal$genes[,c("Block","Row","Column","ID")]
      rg = rg[ order( rg$Block,rg$Row,rg$Column,rg$ID ),]
      gal = GALFile[,c("Block","Row","Column","ID")]
      gal = gal[order( gal$Block,gal$Row,gal$Column,rg$ID ),]
      
      if( table(rg == gal)/dim(gal)[2] == dim(gal)[1]){
        cat("GAL File is the same as the GPR File\n")			
      }
      
    } else{
      cat("GAL File is different from the GPR File, taking GPR Annotation...\n")			
    }
  }
  
  RG$genes$ID = ifelse( RG$genes$ID == "","noSig",RG$genes$ID )
  RG$genes$Name = ifelse( RG$genes$Name == "","noSig",RG$genes$Name )
  
  chipDataL <- list("adf"=adf, "targets"=targets, "RG"=RG)
  
  print(paste("Data succesfully read in and processed from:", currDir))
  return(chipDataL)
}

#' Filtering of spot values by the condition PM > MM and spot intensity > spike in spot intensity
#' @param MA MAList object as defined by limma
#' @param normMethod Optional normalization method name that attached to the log file name (e.g. _vsn_)
#' @param subsetByManualControl Boolean, data may be subsetted by the spike in controls
#' @param subsetByPMvsMMBoolean, data may be subsetted by restriction PM > MM
#' @param ndups Number of duplicated spots
#' @param spacing Placement of the duplicated spots (in default mode it is assumed that spots are sorted)
#' @param weights Parameter for avedups function
#' @param design Design of the experiment is considered when caluclating restrictions (subsetByManualControl,subsetByPMvsMM)
#' @param writeLog Writes a log file of the filtering steps
#' @return Filtered MA object
#' @export
filterSpotsBySpikeIn = function(MA, normMethod="", subsetByManualControl=FALSE, subsetByPMvsMM = FALSE,ndups=4, spacing=1,weights=NULL, design, writeLog=TRUE){
  
  if( dim(MA$A)[2] < 1){
    stop("MAList has no values!!!")
  }
  
  if( dim(MA$A)[2] == 1){
    message("MA List has only one dimension -> copying the A values in order to be able to run this function (this will not affect the outcome since a mean is calculated)")
    MA$A = cbind( MA$A, MA$A )
    colnames(MA$A) = c( colnames(MA$A)[1], paste0(colnames(MA$A)[1],"_fake") )
    MA$M = cbind( MA$M, MA$M )
    colnames(MA$M) = c( colnames(MA$M)[1], paste0(colnames(MA$M)[1],"_fake") )
  }
  
  logstr=c()
  log = function(logEntry, env=sys.frame(-1)){
    #warning logstr must be defined in parent environmnet
    cat(logEntry)
    env$logstr = c(env$logstr,logEntry)
  }
  
  log(paste0("\n ########## Filtering for Normalization: ", normMethod, " ############## \n"))
  
  MAmean = avedups(MA, ndups=ndups, spacing=spacing, weights=weights)
  
  ######################
  #  Comparing PM vs MM probes
  ######################
  
  MAmeanPMvsMM = MAmean
  MAmeanPMvsMM = MAmeanPMvsMM[MAmeanPMvsMM$genes$probe.type == "PM",]
  MAmeanPMvsMM$A = matrix( FALSE, nrow = dim(MAmeanPMvsMM$A)[1], ncol=dim(MAmeanPMvsMM$A)[2] )
  
  pmvsmm = lapply( which(MAmean$genes$probe.type == "PM"), function(i) {
    currPM = MAmean[i,]
    currPM$A = ifelse( is.na(currPM$A), 0, currPM$A )
    currMM = MAmean[which( paste(currPM$genes$Name,"MM",sep="_") == MAmean$genes$Name ),]
    currMM$A = ifelse( is.na(currMM$A), 0, currMM$A )
    #checking if the mmpair contains values -> else there is no MM pair -> therefore candidate cannot be filtered!
    if( dim(currMM)[1] != 0 ){
      return( as.vector(currPM$A) > as.vector(currMM$A) )
    } else{
      return( rep(TRUE, length(currPM$A)) )
    }
  })
  pmvsmm = do.call(rbind, pmvsmm)  
  MAmeanPMvsMM$A = pmvsmm
  
  pmHMMtable = apply(MAmeanPMvsMM$A,2,table)	
  colnames(pmHMMtable) = colnames( MAmeanPMvsMM$M )
  pmHMMtable = apply(pmHMMtable, 2, as.character)
  pmHMMtable = rbind(pmHMMtable, apply( pmHMMtable, 2, function(x){ return(  paste0(round( (1-as.numeric(x[1])/as.numeric(x[2]))*100, digits=2),"%") ) }) )
  rownames(pmHMMtable) = c("FALSE", "TRUE", "PercentPassed")
  if(writeLog){
    write.table( pmHMMtable, file=paste0("pmHMMtable_",normMethod,".csv"),sep="\t", col.names=TRUE, row.names=TRUE )		
  }
  
  if( missing(design) ){
    pmHigherMM = rowSums( MAmeanPMvsMM$A ) == dim(MAmeanPMvsMM$A)[2]
    log( paste0("Design is not considered\n") )
  } else{
    log( paste0( "pmHigherMM by design: Therefore candidates need to be higher in eather WTvsTG or TGvsWT\n" ) )
    
    designBool = as.logical( ifelse(design == -1, 0, design) )
    pmHigherMMwt = rowSums( MAmeanPMvsMM[,designBool]$A ) == dim(MAmeanPMvsMM[,designBool]$A)[2]
    pmHigherMMtg = rowSums( MAmeanPMvsMM[,!designBool]$A ) == dim(MAmeanPMvsMM[,!designBool]$A)[2]
    pmHigherMM = rowSums( cbind( pmHigherMMwt, pmHigherMMtg) ) > 0 
    
  }
  
  log( paste0("\nThere are: ",sum(pmHigherMM)," PM probes higher than MM probes (" , round(sum(pmHigherMM)/length(pmHigherMM)*100,digits=2) , "%) from a total of ",length(pmHigherMM)," PM probes \n") )
  
  ######################
  #	Comparing PM versus Spike in probes
  ######################
  
  maxMin = colMeans( MAmean[MAmean$genes$chipType == "SpikeInControl",]$A, na.rm=TRUE )
  maxMin = ifelse( is.na(maxMin), 0, maxMin )
  
  MAmeanPMvsSpike = MAmean
  MAmeanPMvsSpike = MAmeanPMvsSpike[MAmeanPMvsSpike$genes$probe.type == "PM",]
  MAmeanPMvsSpike$A = matrix( FALSE, nrow = dim(MAmeanPMvsSpike$A)[1], ncol=dim(MAmeanPMvsSpike$A)[2] )
  
  pmvsSpikeIn = lapply( which(MAmean$genes$probe.type == "PM"), function(i) {
    currPM = MAmean[i,]
    currPM$A = ifelse( is.na(currPM$A), 0, currPM$A )
    return( as.vector(currPM$A) > maxMin )
  })
  
  pmvsSpikeIn = do.call(rbind, pmvsSpikeIn)	
  MAmeanPMvsSpike$A = pmvsSpikeIn
  
  pmSpikeIntable = apply(MAmeanPMvsSpike$A,2,table)	
  colnames(pmSpikeIntable) = colnames( MAmeanPMvsSpike$M )
  pmSpikeIntable = apply(pmSpikeIntable, 2, as.character)
  pmSpikeIntable = rbind(pmSpikeIntable, apply( pmSpikeIntable, 2, function(x){ return(  paste0(round((1-as.numeric(x[1])/as.numeric(x[2]))*100, digits=2),"%") ) }) )
  rownames(pmSpikeIntable) = c("FALSE", "TRUE", "PercentPassed")
  if(writeLog){
    write.table( pmSpikeIntable, file=paste0("pmSpikeIntable_",normMethod,".csv"),sep="\t", col.names=TRUE, row.names=TRUE )
  }
  
  if( missing(design) ){
    pmHigherSpikeIn = rowSums( MAmeanPMvsSpike$A ) == dim(MAmeanPMvsSpike$A)[2]		
    log( paste0("Design is not considered\n") )
  } else{
    log( paste0( "spikeIn by design: Therefore candidates need to be higher in eather WTvsTG or TGvsWT\n" ) )
    
    designBool = as.logical( ifelse(design == -1, 0, design) )
    pmHigherSpikeInwt = rowSums( MAmeanPMvsSpike[,designBool]$A ) == dim(MAmeanPMvsSpike[,designBool]$A)[2]
    pmHigherSpikeIntg = rowSums( MAmeanPMvsSpike[,!designBool]$A ) == dim(MAmeanPMvsSpike[,!designBool]$A)[2]
    pmHigherSpikeIn = rowSums( cbind( pmHigherSpikeInwt, pmHigherSpikeIntg) ) > 0 
    
  }
  
  log( paste0("There are: ",sum(pmHigherSpikeIn)," PM probes higher than Spike in probes (" , round(sum(pmHigherSpikeIn)/length(pmHigherSpikeIn)*100,digits=2) , "%) from a total of ",length(pmHigherSpikeIn)," PM probes \n") )
  
  ######################
  #	Subsetting the PM probes!
  ######################
  
  MAmean = MAmean[which(MAmean$genes$probe.type == "PM"),]
  toKeep = c()
  toKeepPM = c()
  toKeepSpike = c()
  if( subsetByPMvsMM ){
    #     toKeep = c(toKeep,which(pmHigherMM)) 
    toKeepPM = which(pmHigherMM)
  }else{
    log("\n##### subsetting by PM higher MM skipped! ##### \n")
  }
  if( subsetByManualControl ){
    #     toKeep = c(toKeep, which(pmHigherSpikeIn) ) 
    toKeepSpike = which(pmHigherSpikeIn) 
  }else{
    log("\n\n##### subsetting by manual controls skipped! ##### \n")
  }
  
  if( subsetByPMvsMM & subsetByManualControl ){
    toKeep = intersect( toKeepPM, toKeepSpike )
  } else if ( subsetByPMvsMM & !subsetByManualControl ){
    toKeep = toKeepPM 
  } else if( !subsetByPMvsMM & subsetByManualControl ){
    toKeep = toKeepSpike 
  } else{
    toKeep = c()
  }
  
  #   toKeep = unique(toKeep)
  
  if( length(toKeep) != 0 ){
    before = dim(MAmean)[1]
    MAmean = MAmean[toKeep,] 
    log( paste0("\nRemoving ", before-(length(toKeep)-8), " of all PM spots (",before,") that did not exceed the filtering threshold! \n") )
  }
  
  # removing the spike in controls (have been set to true, since they do not have a MM pair!)
  spikeInIndex = which(MAmean$genes$chipType != "SpikeInControl")
  if(length(spikeInIndex) != 0){
    MAmean = MAmean[spikeInIndex,]	
  }
  # removing the spike in controls (have been set to true, since they do not have a MM pair!)
  MAsub = MA[MA$genes$probe.type == "PM",]
  
  #Number of oligos separated in chip types before subsetting by filtering
  numberChipType = table(as.character(MAsub[which(MAsub$genes$chipType != "SpikeInControl"),]$genes$chipType))/ndups
  
  MAsub = MAsub[ which( MAsub$genes$ID %in% MAmean$genes$ID ),]
  
  log( paste0( "\nThis results in ", dim(MAmean)[1], " oligos, or ", dim(MAsub)[1], " spots! ", dim(MAmean)[1]*ndups == dim(MAsub)[1] ) )
  
  resTable = table(MAmean$genes$chipType)
  log( paste0( "\nThis results in: ", names(resTable), ": ", as.vector(resTable), " oligos (", round(as.vector(resTable)/as.vector(numberChipType)*100, digits=2) ,"%)") )
  if(writeLog){
    writeLines(logstr, con=file.path(getwd(),paste0("filterLog_",normMethod)))
    message("log file written")
  }
  return( MAsub )
}


