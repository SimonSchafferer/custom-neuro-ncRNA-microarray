custom-neuro-ncRNA-microarray
=============================
This package creates analysis templates for the analysis of a custom microarray and therefore it is intended for internal use. 

#Installation
```R
library(devtools)
install_github("SimonSchafferer/custom-neuro-ncRNA-microarray")
```

##Examples

The gpr files have to be in the specified directory (PATHTOGPRFILES) and hast to contain a samplesInfo.txt file that should contain the following entries: 

|filename	 |Cy3|Cy5|
|---------------|----|---|
|0302-1574si.gpr|WT01|TG01|
|0302-1578si.gpr|TG01|TG01|
| ... | ... | ... |

The executeKnitR option in each command specifies if the code should be run and the html files should be generated, or if only the code is printed for modification.

###NRC3_01 (Chip 1.0)
```R
library(CustomMicroarrayPackage)
performAnalysis1_0(gprDir="PATHTOGPRFILES", executeKnitR=FALSE)
```

###NRC4_01 (Chip 2.0)
```R
library(CustomMicroarrayPackage)
performAnalysis2_0(gprDir="PATHTOGPRFILES", executeKnitR=FALSE)
```

###NRC5_01 (Chip 2.1)
```R
library(CustomMicroarrayPackage)
performAnalysis2_1(gprDir="PATHTOGPRFILES", executeKnitR=FALSE)
```

##SQLite Database for annotation
The microarray annotations and informations about the oligonucleotides can be accessed through a sqlite database. 

```R
library(sqldf)
db <- dbConnect(SQLite(), dbname=system.file("data/ncrna.chip2.annotation.db",package="CustomMicroarrayPackage"))

#inspect tables
dbListTables(db)

#Fetch data from the database
#For example load annotation data for some differentially expressed candidates: 
diffExpCandidatesExample = c("e787",e56","e597")

additionalData = dbGetQuery(db, paste(
"SELECT *  
FROM oligoIDMapping oidm
JOIN oligo o ON oidm.oligo_oligoTableID = o.oligoTableID
JOIN secondary seco ON seco.secondaryTableID = o.secondary_secondaryTableID
LEFT OUTER JOIN contig c ON c.contigTableID = oidm.contig_contigTableID
LEFT OUTER JOIN annotationSummary ans ON ans.contigAnnotation_contig_contigTableID = c.contigTableID
JOIN geneAnnot ga ON ga.contigAnnotation_contig_contigTableID = c.contigTableID
LEFT OUTER JOIN predictionAdditional preda ON preda.oligo_oligoTableID = o.oligoTableID
WHERE o.oligoTableID IN( ",paste(  sub("e","",diffExpCandidatesExample) , collapse=","),");",sep=""))

```

## UML
![alt tag](https://github.com/SimonSchafferer/custom-neuro-ncRNA-microarray/blob/master/data/erdiagram.png)


##Publication of the microarray
###Abstract
We have generated a novel, neuro-specific ncRNA microarray, covering 1472 ncRNA species, to investigate their expression in different mouse models for central nervous system diseases. Thereby, we analyzed ncRNA expression in two mouse models with impaired calcium channel activity, implicated in Epilepsy or Parkinson's disease, respectively, as well as in a mouse model mimicking pathophysiological aspects of Alzheimer's disease. We identified well over a hundred differentially expressed ncRNAs, either from known classes of ncRNAs, such as miRNAs or snoRNAs or which represented entirely novel ncRNA species. Several differentially expressed ncRNAs in the calcium channel mouse models were assigned as miRNAs and target genes involved in calcium signaling, thus suggesting feedback regulation of miRNAs by calcium signaling. In the Alzheimer mouse model, we identified two snoRNAs, whose expression was deregulated prior to amyloid plaque formation. Interestingly, the presence of snoRNAs could be detected in cerebral spine fluid samples in humans, thus potentially serving as early diagnostic markers for Alzheimer's disease. In addition to known ncRNAs species, we also identified 63 differentially expressed, entirely novel ncRNA candidates, located in intronic or intergenic regions of the mouse genome, genomic locations, which previously have been shown to harbor the majority of functional ncRNAs.

<a href="http://www.ncbi.nlm.nih.gov/pubmed/25344396">Generation of a neuro-specific microarray reveals novel differentially expressed noncoding RNAs in mouse models for neurodegenerative diseases.</a>



