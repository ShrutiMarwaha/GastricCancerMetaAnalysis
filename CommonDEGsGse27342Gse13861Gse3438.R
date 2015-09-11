## i) to find commmon set of genes differentially expressed among different gastric cancer microarray studies
## ii) submit the commmon genes to cmap
## load the output (list of differentially expressed genes) of the microarray studies
load("/Users/shruti/Dropbox/SHRUTIM/Microarray/Microarray_R_Scripts/GastricCancer/GSE27342/final/2014/Gse27342DEGFinal.rda")
load("/Users/shruti/Dropbox/SHRUTIM/Microarray/Microarray_R_Scripts/GastricCancer/GSE13861/rough/DEG_gse13861.rda")
load(file="/Users/shruti/Dropbox/SHRUTIM/Microarray/Microarray_R_Scripts/GastricCancer/GSE3438/final/DEG_gse3438.rda")

setwd("/Users/shruti/Dropbox/SHRUTIM/Microarray/Microarray_R_Scripts/CommonGene/rough/3GastricCancer")

dim(Gse27342DEGFinal)
complete <- complete.cases(Gse27342DEGFinal)
Gse27342DEGFinal <- Gse27342DEGFinal[complete,]
Gse27342DEGFinal <- subset(Gse27342DEGFinal,GeneID!="",)
gse27342_up <- unique(subset(Gse27342DEGFinal,logFC>0)$GeneID)
gse27342_down <- unique(subset(Gse27342DEGFinal,logFC<0)$GeneID)

dim(DEG_gse13861)
complete <- complete.cases(DEG_gse13861)
DEG_gse13861 <- DEG_gse13861[complete,]
DEG_gse13861 <- subset(DEG_gse13861,Entrez_Gene_ID!="",)
gse13861_up <- unique(subset(DEG_gse13861,logFC>0)$Entrez_Gene_ID)
gse13861_down <- unique(subset(DEG_gse13861,logFC<0)$Entrez_Gene_ID)

dim(DEG_gse3438)
complete <- complete.cases(DEG_gse3438)
DEG_gse3438 <- DEG_gse3438[complete,]
DEG_gse3438 <- subset(DEG_gse3438,Gene.ID!="",)
grep("/",DEG_gse3438$Gene.ID)
# if GeneID contains more than 1 id, remove additional ones
GeneIds <- sapply(DEG_gse3438$Gene.ID,function(x) {
  sub("/.+","",x)
})
DEG_gse3438$Gene.ID <- GeneIds
gse3438_up <- unique(subset(DEG_gse3438,logFC>0)$Gene.ID)
gse3438_down <- unique(subset(DEG_gse3438,logFC<0)$Gene.ID)

library(limma)
VennIntersection <- function(set1,set2,set3)
{
  #total set
  universe <- sort(unique(c(set1,set2,set3)))
  
  # Generate a matrix, with the sets in columns and possible letters on rows
  Counts <- matrix(0, nrow=length(universe), ncol=3)
  # Populate the said matrix
  for (i in 1:length(universe)) 
  {
    Counts[i,1] <- universe[i] %in% set1
    Counts[i,2] <- universe[i] %in% set2
    Counts[i,3] <- universe[i] %in% set3
  }
  # Name the columns with the sample names
  colnames(Counts) <- c("GSE13861","GSE27342","GSE3438")
  rownames(Counts) <- universe
  # Specify the colors for the sets
  cols<-c("Red", "Green", "Blue")
  vennDiagram(vennCounts(Counts), circle.col=cols)
}

VennIntersection(gse13861_up,gse27342_up,gse3438_up)
VennIntersection(gse13861_down,gse27342_down,gse3438_down)
CommonUp_Gse3438_Gse13861_Gse27342 <- Reduce(intersect, list(gse3438_up,gse13861_up,gse27342_up))
CommonDown_Gse3438_Gse13861_Gse27342 <- Reduce(intersect, list(gse3438_down,gse13861_down,gse27342_down))

write.table(CommonUp_Gse3438_Gse13861_Gse27342,"CommonUp_Gse3438_Gse13861_Gse27342.txt",quote=F,col.names=F,row.names=F)
save(CommonUp_Gse3438_Gse13861_Gse27342,file="CommonUp_Gse3438_Gse13861_Gse27342.rda")
write.table(CommonDown_Gse3438_Gse13861_Gse27342,"CommonDown_Gse3438_Gse13861_Gse27342.txt",quote=F,col.names=F,row.names=F)
save(CommonDown_Gse3438_Gse13861_Gse27342,file="CommonDown_Gse3438_Gse13861_Gse27342.rda")

####################################################################
###### Cmap
####################################################################
## use biomaRt package to convert entrez ids to hgu1331a afffy probeset ids
#source("http://bioconductor.org/biocLite.R")
#biocLite("biomaRt")
library("biomaRt")
listMarts()
## select biomart database and dataset
database <- useMart("ensembl")
grep("sapiens",listDatasets(database)$description,)
listDatasets(database)[32,]
#data <- useDataset("hsapiens_gene_ensembl",mart=database)
## connect to a specified biomart database
data <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
## list the filters available
filters <- listFilters(data)
filters
test <- grep("hg.?u.?133.?a",filters$description,ignore.case=T)
filters[test,]
## attribites are values that you are interested in to retrieve
attributes <- listAttributes(data)
attributes
grep("entrez",attributes$description,ignore.case=T)
grep("hg.?u.?133.?a",attributes$description,ignore.case=T)
#attributes[c(55,56,100,101),]
attributes[c(59,60,108,109),]
#######################################################################################

CommonUpProbesetIDs <- getBM(attributes=c('entrezgene','affy_hg_u133a'), filters = 'entrezgene', values = CommonUp_Gse3438_Gse13861_Gse27342, mart = data)
CommonUpProbesetIDs <- subset(CommonUpProbesetIDs,affy_hg_u133a!="")
CommonUpCmapInput <- CommonUpProbesetIDs$affy_hg_u133a
write.table(CommonUpCmapInput,file="CommonUpCmapInput.txt.grp",sep="\t",quote=F,col.names=F,row.names=F)

CommonDownProbesetIDs <- getBM(attributes=c('entrezgene','affy_hg_u133a'), filters = 'entrezgene', values = CommonDown_Gse3438_Gse13861_Gse27342, mart = data)
CommonDownProbesetIDs <- subset(CommonDownProbesetIDs,affy_hg_u133a!="")
CommonDownCmapInput <- CommonDownProbesetIDs$affy_hg_u133a
write.table(CommonDownCmapInput,file="CommonDownCmapInput.txt.grp",sep="\t",quote=F,col.names=F,row.names=F)
