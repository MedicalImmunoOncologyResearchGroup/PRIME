## Reference for tximport and use of DESeq2 with such data:
## Importing transcript abundance with tximport Michael I. Love, Charlotte Soneson, Mark D. Robinson, 2023-10-24

## Set working directory
setwd('/Users/oipulk/Documents/prime_data')

## Load libraries
library(rhdf5)
library(tximport)
library(tximportData)
library("readxl")

###### Gide cohort ##########

## Get paths and sample names 
sampleNames <- list.dirs('raw_data/Gide/Kallisto', full.names = FALSE)[-1]

## Load table that connects transcript id's to gene id's and gene names
txToGenes <- read.table('raw_data/Kallisto_index/homo_sapiens/transcripts_to_genes.txt', sep = '\t' )
colnames(txToGenes) <- c('target_id', 'Genes', 'GeneSymbol')

## Define files and load 
files <- file.path('/Users/oipulk/Documents/prime_data/raw_data/Gide', "Kallisto", sampleNames, "abundance.tsv")
txi.Gide <- tximport(files, type = "kallisto", tx2gene = txToGenes, ignoreAfterBar = TRUE)
colnames(txi.Gide$counts)<-sampleNames
colnames(txi.Gide$length)<-sampleNames

## Fetch the patient id's and construct a table connecting them to file names
gideFileReport <- read.table('raw_data/Gide/filereport_read_run_PRJEB23709_tsv.txt', sep= '\t', header = TRUE)
sn <- gsub('_R1.*', '', gideFileReport$submitted_ftp)
sn <- gsub('.*/', '', sn)
sn <- data.frame(gideFileReport$run_accession, sn)
sn <- sn[order(sn$gideFileReport.run_accession), ]

## Rename count matrix columns by patient id's
txi.Gide$length <- txi.Gide$length[, order(colnames(txi.Gide$counts))]
txi.Gide$counts <- txi.Gide$counts[, order(colnames(txi.Gide$counts))]
colnames(txi.Gide$counts) <- sn$sn
colnames(txi.Gide$length) <- sn$sn

## Extract pretherapy samples and remove '_PRE' from the patient id's
pretherapySamples <- grep('_PRE', colnames(txi.Gide$counts), value=TRUE)

txi.Gide$length<-txi.Gide$length[,grep('_PRE', colnames(txi.Gide$counts))]
txi.Gide$counts<-txi.Gide$counts[,grep('_PRE', colnames(txi.Gide$counts))]

colnames(txi.Gide$length) <- gsub('_PRE', '', colnames(txi.Gide$length))
txi.Gide$length <- txi.Gide$length[,order(colnames(txi.Gide$length))]

colnames(txi.Gide$counts) <- gsub('_PRE', '', colnames(txi.Gide$counts))
txi.Gide$counts <- txi.Gide$counts[,order(colnames(txi.Gide$counts))]

## Load patient metadata for DEA
metadataMono <- read.table("/Users/oipulk/Documents/prime_data/raw_data/Gide/Gide_clinical_mono.tsv", sep='\t', header=TRUE)
metadataMono <- unique(metadataMono)
rownames(metadataMono) <- 1:nrow(metadataMono)

metadataCombination <- read.table("/Users/oipulk/Documents/prime_data/raw_data/Gide/Gide_clinical_combination.tsv", sep='\t', header=TRUE)
metadataCombination <- unique(metadataCombination)
rownames(metadataCombination) <- 1:nrow(metadataCombination)

## Rename patients to reflect therapy and combine metadata
metadataMono$Patient <- paste('PD1_',metadataMono$Patient, sep="")
metadataCombination$Patient <- paste('ipiPD1_',metadataCombination$Patient, sep="")

metadata <- rbind(metadataMono,metadataCombination)
metadata <- metadata[order(metadata$Patient),]
rownames(metadata) <- 1:nrow(metadata)

## Include only patients with sequencing data and clinical data
commonPatients <- intersect(metadata$Patient,colnames(txi.Gide$counts))
txi.Gide$counts <- txi.Gide$counts[,commonPatients]
txi.Gide$length <- txi.Gide$length[,commonPatients]
metadata <- metadata[metadata$Patient%in%commonPatients,]


## Choose the therapy type (mono or combi or both) for DEA
all_therapies <- 1:length(metadata$Patient)
combination_therapies = grep('ipi',metadata$Patient)
monotherapies = setdiff(all_therapies, combination_therapies)

chosen_therapies <- all_therapies
## chosen_therapies <- monotherapies
## chosen_therapies <- combination_therapies

txi.Gide$counts <- txi.Gide$counts[,chosen_therapies]
txi.Gide$length <- txi.Gide$length[,chosen_therapies]
metadata <- metadata[chosen_therapies,]

sampleTable_Gide <- data.frame( batch =factor(rep('GIDE', length(factor(metadata$Response)))), gender = factor(metadata$Gender), condition = factor(metadata$Response) )
rownames(sampleTable_Gide) <- colnames(txi.Gide$counts)


###### Hugo cohort ##########

## Get paths and sample names
sampleNames <- list.dirs('raw_data/Hugo/Kallisto', full.names = FALSE)[-1]

## Load table that connects transcript id's to gene id's and gene names
txToGenes <- read.table('raw_data/Kallisto_index/homo_sapiens/transcripts_to_genes.txt', sep = '\t' )
colnames(txToGenes) <- c('target_id', 'Genes', 'GeneSymbol')

## Define files and load 
files <- file.path('/Users/oipulk/Documents/prime_data/raw_data/Hugo', "Kallisto", sampleNames, "abundance.tsv")
txi.Hugo <- tximport(files, type = "kallisto", tx2gene = txToGenes, ignoreAfterBar = TRUE)
colnames(txi.Hugo$counts)<-sampleNames
colnames(txi.Hugo$length)<-sampleNames

## Load patient metadata for DEA
metadata <- read.table("/Users/oipulk/Documents/prime_data/raw_data/Hugo/Clinical_Outcome.tsv", sep='\t', header=TRUE)

## Select only pre-treatment samples with RNA Seq data
metadata<-metadata[metadata$Biopsy.Time=="pre-treatment",]
metadata<-metadata[!is.na(metadata$SRA.Run.ID..tumor.RNA),]

patients_with_RNAseq_meta <- intersect(metadata$SRA.Run.ID..tumor.RNA,colnames(txi.Hugo$length))

patients_with_RNAseq_meta <- patients_with_RNAseq_meta[!(metadata$irRECIST[metadata$SRA.Run.ID..tumor.RNA%in%patients_with_RNAseq_meta]=="")]

metadata <- metadata[metadata$SRA.Run.ID..tumor.RNA%in%patients_with_RNAseq_meta,]
txi.Hugo$length <- txi.Hugo$length[,patients_with_RNAseq_meta]
txi.Hugo$counts <- txi.Hugo$counts[,patients_with_RNAseq_meta]

sampleTable_Hugo <- data.frame( batch =factor(rep('HUGO', length(factor(metadata$irRECIST)))), gender = factor(metadata$Gender),  condition = factor(as.integer(as.logical((metadata$irRECIST=="Complete Response")|(metadata$irRECIST=="Partial Response")))))

rownames(sampleTable_Hugo) <- colnames(txi.Hugo$counts)



###### Riaz cohort ##########

## Get paths and sample names
sampleNames <- list.dirs('raw_data/Riaz/Kallisto', full.names = FALSE)[-1]

## Load table that connects transcript id's to gene id's and gene names
txToGenes <- read.table('raw_data/Kallisto_index/homo_sapiens/transcripts_to_genes.txt', sep = '\t' )
colnames(txToGenes) <- c('target_id', 'Genes', 'GeneSymbol')

## Define files and load 
files <- file.path('/Users/oipulk/Documents/prime_data/raw_data/Riaz', "Kallisto", sampleNames, "abundance.tsv")
txi.Riaz <- tximport(files, type = "kallisto", tx2gene = txToGenes, ignoreAfterBar = TRUE)
colnames(txi.Riaz$counts)<-sampleNames
colnames(txi.Riaz$length)<-sampleNames

## Load patient metadata for DEA
metadata <- read.table("/Users/oipulk/Documents/prime_data/raw_data/Riaz/Clinical_characteristics_mmc2.tsv", sep='\t', header=TRUE)

## Load files that connects patient id's to sample names
riazInfo1 <- read.table('raw_data/Riaz/sra_result.csv', header = TRUE, sep = ',')
riazInfo2 <- read.table('raw_data/Riaz/SraRunInfo.csv', header = TRUE, sep = ',')
colnames(riazInfo1)[1] <- 'Experiment'
riazInfo <- merge(riazInfo1, riazInfo2, by = 'Experiment')
sn <- gsub('.*: ', '', riazInfo$Experiment.Title)

##sn <- gsub('_.*', '', sn)
## sn <- sn[order(sn$riazInfo.Run),]

preTherapySamples <- sn%in%grep('_Pre', sn, value=TRUE)

txi.Riaz$counts <- txi.Riaz$counts[,colnames(txi.Riaz$counts)%in%riazInfo$Run[preTherapySamples]]
txi.Riaz$length <- txi.Riaz$length[,colnames(txi.Riaz$length)%in%riazInfo$Run[preTherapySamples]]

## Label patients by Patient id's Pt1 etc.
sn <- gsub('_.*', '',sn[preTherapySamples])
colnames(txi.Riaz$counts) <- sn
colnames(txi.Riaz$length) <- sn

sn <- sn[order(sn)]
txi.Riaz$counts <- txi.Riaz$counts[,order(colnames(txi.Riaz$counts))]
txi.Riaz$length <- txi.Riaz$length[,order(colnames(txi.Riaz$length))]

metadata <- metadata[metadata$Patient%in%sn,]

## Acral and UVEAL/OCULAR subtypes have zero response. Let's leave them out
mask <- (metadata$Subtype!='ACRAL')& (metadata$Subtype!='OCULAR/UVEAL')
txi.Riaz$counts <- txi.Riaz$counts[, mask]
txi.Riaz$length <- txi.Riaz$length[, mask]
metadata <- metadata[mask,]

sampleTable_Riaz <- data.frame( batch =factor(rep( 'RIAZ', length(factor(metadata$Response)))), gender =factor(rep( 'UNKNOWN', length(factor(metadata$Response)))), condition = factor(as.integer(as.logical((metadata$Response=="CR")|(metadata$Response=="PR")))))

rownames(sampleTable_Riaz) <- colnames(txi.Riaz$counts)

############ Stack data ###################

sampleTable = rbind(sampleTable_Gide, sampleTable_Hugo, sampleTable_Riaz)

txi.kallisto.tsv <- txi.Gide
txi.kallisto.tsv$counts <-  cbind(txi.Gide$counts, txi.Hugo$counts, txi.Riaz$counts )
txi.kallisto.tsv$length <- cbind(txi.Gide$length, txi.Hugo$length, txi.Riaz$length )
txi.kallisto.tsv$abundance <- NULL

############ Start DEA #######################

library(DESeq2)

dds <- DESeqDataSetFromTximport(txi.kallisto.tsv, sampleTable,  ~ batch + condition)

## Filter low count rows
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

## Run DESeq
dds <- DESeq(dds)
res <- results(dds)

## Order rows by p-value
resOrdered <- res[order(res$pvalue),]

## Rename rows by common gene names
genesToSymbols <- txToGenes[c('Genes','GeneSymbol')]
rownames(resOrdered) <- genesToSymbols$GeneSymbol[match(rownames(resOrdered),genesToSymbols$Genes)]

## Export results
write.csv(resOrdered, "DESeq2_results.csv", row.names=TRUE)
