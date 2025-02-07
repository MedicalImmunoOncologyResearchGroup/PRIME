
setwd('./')
####################################
##Reading data of each cohort 
## visualizing the data using PCA plots to see whether batch effects are present or not
library(mixOmics)
library(factoextra)

## reading in TPM values of all cohorts
riazTPM <- read.table('./TPM_RIAZ_VALUES_FOUND_IN_CLINICAL_DATA.csv', sep = '\t', header = TRUE)
gideTPM <- read.table('./TPM_GIDE.csv', sep = '\t', header = TRUE)
hugoTPM <- read.table('./TPM_HUGO.csv', sep = '\t', header = TRUE)
liuTPM <- read.table('./Liu_TPM_Star_RSEM.csv', sep = '\t', header = TRUE)

row.names(riazTPM) <- riazTPM$Gene
row.names(gideTPM) <- gideTPM$Genes
row.names(hugoTPM) <- hugoTPM$Genes
row.names(liuTPM) <- liuTPM$Genes

## plotting PCA plots separately with grouping
riazGrouping <- gsub('_.*', '', gsub('Pt[0-9]+_', '', colnames(riazTPM)))[-1]
riazTPM <- riazTPM[rowSums(riazTPM[, 2:ncol(riazTPM)]) > 0,]
riazPCA <- pca(t(riazTPM[, 2:ncol(riazTPM)]), center = TRUE, scale = TRUE)

hugoClinical <- read_xls('./Clinical Outcome .xls')
hugoGrouping <- hugoClinical1$irRECIST[c(1:20, 20, 21:27)]
hugoGrouping <- gsub(' ', '_', hugoGrouping)
hugoTPM <- hugoTPM[rowSums(hugoTPM[, 2:ncol(hugoTPM)]) > 0,]
hugoPCA <- pca(t(hugoTPM[, 2:ncol(hugoTPM)]), center = TRUE, scale = TRUE)

gideGrouping <- gsub('_[0-9]+', '', colnames(gideTPM))[-1]
gideTPM <- gideTPM[rowSums(gideTPM[, 2:ncol(gideTPM)]) > 0,]
gidePCA <- pca(t(gideTPM[, 2:ncol(gideTPM)]), center = TRUE, scale = TRUE)

liuClinical <- read_xlsx('./41591_2019_654_MOESM4_ESM.xlsx', sheet = 'Supplemental Table 1', skip = 2)
liuClinical <- liuClinical[1:144, ]
liuTPMUpdated <-liuTPM[, colnames(liuTPM) %in% c('Genes', intersect(colnames(liuTPM), liuClinical$...1))] 
liuTPMUpdated <- liuTPMUpdated[, order(colnames(liuTPMUpdated))]
liuClinical <- liuClinical[liuClinical$...1 %in% intersect(colnames(liuTPM), liuClinical$...1),]
liuGrouping1 <- liuClinical$Tx
liuTPMUpdated <- liuTPMUpdated[rowSums(liuTPMUpdated[, 2:ncol(liuTPMUpdated)]) > 0,]
liuPCA <- pca(t(liuTPMUpdated[, 2:ncol(liuTPMUpdated)]), center = TRUE, scale = TRUE)

pdf('./PCA_Seperately.pdf')
par(mfrow=c(2,2))
plotIndiv(hugoPCA, group = hugoGrouping, legend = TRUE,  pch = c(16, 17, 18), col.per.group = c('red', 'blue', 'green'))
plotIndiv(riazPCA, group = riazGrouping, legend = TRUE,  pch = c(16,17), col.per.group = c('red', 'blue'))
plotIndiv(gidePCA, group = gideGrouping, legend = TRUE,  pch = c(15,16,17, 18), col.per.group = c('red', 'blue', 'green', 'black'))
plotIndiv(liuPCA, group = liuGrouping1, legend = TRUE,  pch = c(16,17), col.per.group = c('red', 'blue'))
dev.off()


## combining 3 data sets (Riaz, Gide and Hugo)

## reading in the data sets of different cohorts

intersect_all <- function(a,b,...){
  Reduce(intersect, list(a,b,...))
}

genes <- intersect_all(riazTPM$Gene, gideTPM$Genes, hugoTPM$Genes, liuTPM$Genes)

##Patient 61 was excluded from the Liu cohort
liuTPMCG <- liuTPMUpdated[liuTPMUpdated$Genes %in% genes, c(2:101, 103:ncol(liuTPMUpdated))]
riazTPMCG <- riazTPM[riazTPM$Gene %in% genes, 2: ncol(riazTPM)]
gideTPMCG <- gideTPM[gideTPM$Genes %in% genes, 2:ncol(gideTPM)]
hugoTPMCG <- hugoTPM[hugoTPM$Genes %in% genes, 2:ncol(hugoTPM)]

combinedTPM <- data.frame(liuTPMCG, riazTPMCG, gideTPMCG,hugoTPMCG)
combinedTPM <- combinedTPM[rowSums(combinedTPM[]) > 0, ]
liuGrouping1 <- liuGrouping1[-102]

combGroup1 <- c(rep('Liu', 120), rep('Riaz', 101), rep('Gide', 91), rep('Hugo', 28))
combGroup2 <- c(liuGrouping1, riazGrouping, gideGrouping, hugoGrouping)
pcaComb<- prcomp(t(combinedTPM), center = TRUE, scale = TRUE)
fviz_eig(pcaComb)
pdf('./PCA_Combined_Data.pdf')
fviz_pca_ind(pcaComb, geom = 'point', label="none", habillage= combGroup1) 
dev.off()

write.table(combinedTPM, file = './CombinedTPM.csv', sep = '\t', row.names = FALSE, quote = FALSE)

##############################
## Performing Combat analysis on the combined dataset
library(sva)
combTPMBatch <- ComBat(combinedTPM, batch = combGroup1)

pcaCombBatchCor <- prcomp(t(combTPMBatch), center = TRUE, scale. = TRUE)
pdf('./PCA_BATCH_EFFECT_CORRECTED.pdf')


fviz_pca_ind(pcaCombBatchCor, geom = 'point', label="none", habillage= combGroup1)
dev.off()
combTPMBatch <- cbind('Genes' =row.names(combTPMBatch), combTPMBatch)
write.table(combTPMBatch, file = './CombinedTPM_Batch_Corrected.csv', sep = '\t', row.names = FALSE, quote = FALSE)

write.table(liuTPMBC, file = './Liu_Batch_Corrected.csv', sep = '\t', row.names = FALSE, quote = FALSE)
write.table(riazTPMBC, file = './Riaz_Batch_Corrected.csv', sep = '\t', row.names = FALSE, quote = FALSE)
write.table(gideTMPBC, file = './GideTPM_Batch_Corrected.csv', sep = '\t', row.names = FALSE, quote = FALSE)
write.table(hugoTMPBC, file = './HugoTPM_Batch_Corrected.csv', sep = '\t', row.names = FALSE, quote = FALSE)





