# Set working directory
setwd('/Users/arfamehmood/Documents/Projects_Saara/Project_CD74/')
## libraries
library(dplyr)

#######################################
## Preparing count matrix for Gide Data
## set the paths to the Kallisto folder 
samples <- list.dirs('../Kallisto')[-1]
samplesName <- list.dirs('../Kallisto', full.names = FALSE)[-1]

## getting the TPM values for each sample and combing them into a matrix
for(i in 1:length(samples)){
  
  tpms <- read.table(paste0(samples[i], '/abundance.tsv'), header=TRUE, sep = '\t')
  tpms <- tpms[order(tpms$target_id), ]
  if(i==1){
    gideTPM <- tpms[, c(1, 5)]
  }
  else{
    gideTPM <- data.frame(gideTPM, tpms$tpm)  
  }
}
colnames(gideTPM) <- c('target_id', samplesName)


## reading gene and transcript mapping file and merging transcripts and gene name

txToGenes <- read.table('./transcripts_to_genes.txt', sep = '\t' )
colnames(txToGenes) <- c('target_id', 'Genes', 'GeneSymbol')
gideTPM <- merge(gideTPM, txToGenes, by = 'target_id')

## converting transcript TPMs to gene TPMS and summing the values
gideTPMGene <- gideTPM[, c(2:92, 94)] %>% group_by( GeneSymbol) %>% summarise(across(everything() , list(sum)))

## changing the ENA accession numbers of samples to the biological names using metadata file
gideMetadata <- read.table('./filereport_read_run_PRJEB23709_tsv.txt', sep= '\t', header = TRUE)
sn <- gsub('_R1.*', '', gideMetadata$submitted_ftp)
sn <- gsub('.*/', '', sn)
sn <- data.frame(gideMetadata$run_accession, sn)
sn <- sn[order(sn$gideMetadata.run_accession), ]
gideTPMGene <- gideTPMGene[, order(colnames(gideTPMGene))]
colnames(gideTPMGene) <- c(sn$sn, 'Genes')
gideTPMGene <- gideTPMGene[, c(92, 1:91)]
## saving the matrix
write.table(gideTPMGene, file = './Gide_TPM.csv',sep = '\t', row.names = FALSE, quote = FALSE)


###########################################
## preparing count matrix for Hugo data set


samples <- list.dirs('OtherCohorts/Hugo/Kallisto')[-1]
samplesName <- list.dirs('OtherCohorts/Hugo/Kallisto', full.names = FALSE)[-1]

## getting the TPM values for each sample and combing them into a matrix
for(i in 1:length(samples)){
  
  tpms <- read.table(paste0(samples[i], '/abundance.tsv'), header=TRUE, sep = '\t')
  tpms <- tpms[order(tpms$target_id), ]
  if(i==1){
    TPM <- tpms[, c(1, 5)]
  }
  else{
    TPM <- data.frame(TPM, tpms$tpm)  
  }
}
colnames(TPM) <- c('target_id', samplesName)

## reading transcript gene mapping file
txToGenes <- read.table('./transcripts_to_genes.txt', sep = '\t' )
colnames(txToGenes) <- c('target_id', 'Genes', 'GeneSymbol')
TPM <- merge(TPM, txToGenes, by = 'target_id')

## conserting transcript TPMs to gene TPMS and summing the values
TPMGene <- TPM[, c(2:29, 31)] %>% group_by( GeneSymbol) %>% summarise(across(everything() , list(sum)))

## changing the ENA accession numbers of samples to the biological names using metadata file
samplesHugo <- read.table('./SraRunInfo.csv', sep = ',', header = TRUE)
colnames(TPMGene) <- c('Genes', samplesHugo$SampleName)

## saving the matrix
write.table(TPMGene, file = './Hugo_TPM.csv',sep = '\t', row.names = FALSE, quote = FALSE)









############################################
## Preparing count matrix for Riaz Data
samples <- list.dirs('OtherCohorts/Riaz/Kallisto')[-1]
samplesName <- list.dirs('OtherCohorts/Riaz/Kallisto', full.names = FALSE)[-1]

## getting the TPM values for each sample and combing them into a matrix
for(i in 1:length(samples)){
  
  tpms <- read.table(paste0(samples[i], '/abundance.tsv'), header=TRUE, sep = '\t')
  tpms <- tpms[order(tpms$target_id), ]
  if(i==1){
    TPM <- tpms[, c(1, 5)]
  }
  else{
    TPM <- data.frame(TPM, tpms$tpm)  
  }
}
colnames(TPM) <- c('target_id', samplesName)


## reading gene and transcript mapping file and merging transcrip and gene name

txToGenes <- read.table('OtherCohorts/Kallisto_index/homo_sapiens/transcripts_to_genes.txt', sep = '\t' )
colnames(txToGenes) <- c('target_id', 'Genes', 'GeneSymbol')
TPM <- merge(TPM, txToGenes, by = 'target_id')

## conserting transcript TPMs to gene TPMS and summing the values
TPMGene <- TPM[, c(2:110, 112)] %>% group_by( GeneSymbol) %>% summarise(across(everything() , list(sum)))
colnames(TPMGene) <- gsub('_.*', '', colnames(TPMGene))
#samplesRiaz <- read.table('/Users/arfamehmood/Downloads/SraRunInfo (1).csv', sep = ',', header = TRUE)

## getting metadata to set the samplenames
riazInfo1 <- read.table('OtherCohorts/Riaz/sra_result.csv', header = TRUE, sep = ',')
riazInfo2 <- read.table('OtherCohorts/Riaz/SraRunInfo (5).csv', header = TRUE, sep = ',')
colnames(riazInfo1)[1] <- 'Experiment'
riazInfo <- merge(riazInfo1, riazInfo2, by = 'Experiment')
sn <- gsub('.*: ', '', riazInfo$Experiment.Title)
sn <- gsub('_.*', '', sn)
sn <- sn[order(sn$riazInfo.Run),]
write.table(sn, file = 'OtherCohorts/Riaz/Riaz_SampleInfo.csv',sep = '\t', row.names = FALSE, quote = FALSE)




## saving the matrix
write.table(TPMGene, file = 'OtherCohorts/Riaz/Riaz_TPM.csv',sep = '\t', row.names = FALSE, quote = FALSE)

