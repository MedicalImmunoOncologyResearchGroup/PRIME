# Perform Survival analysis (Kaplan Meier) for the Transcript per million (batch 
#corrected using combatseq) from 4 different ICT cohorts. The survival analysis is performed
# by selecting the quantile having a minimum p-value from quartile 0.25 to 0.75.
# input: gene expression data, survival/clinical data
# output: Kaplan Meier plots, KM Statistics (p-value, selected quantile)

# loading the libraries
library(survival)
library(survminer)
source('code/KMTesting.R')
# set working directory
setwd('./PRIME')

# loading the data sets'
ictGeneExp <- read.table('./ICT/CombinedTPM_Batch_Corrected.csv', sep = '\t', header = TRUE)
row.names(ictGeneExp) <- ictGeneExp$Genes
ictGeneExp$Genes <- NULL
ictClinical <- read.table('./ICT/ClinicalData.txt', sep = '\t', header = TRUE)

# Selecting the PRE samples from ICT cohorts
ictClinical$SurvivalTimeDays <- gsub(',', '.', ictClinical$SurvivalTimeDays)
## removing EDT and On samples and renaming the patients similarly as clinical data
ictGeneExp <- ictGeneExp[, -grep('_On', colnames(ictGeneExp))]
ictGeneExp <- ictGeneExp[, -grep('_EDT', colnames(ictGeneExp))]

samples <- gsub('_Pre_.*', '', colnames(ictGeneExp))
samples <- gsub('_PRE', '', samples)
samples <- gsub('.*_', '', samples)
colnames(ictGeneExp) <- samples
ictGeneExp <- ictGeneExp[, colnames(ictGeneExp) %in% ictClinical$sample]
colnames(ictClinical)[3:4] <- c('OS.time', 'OS')

# Setting up the variables
quan <- NULL
ratio <- NULL
pval <- NULL
q <- NULL
gene <- NULL
folderPath <- './ICT/KM'
fileName <- 'KM_ICT_AllGenes.pdf'
pdf(file = paste0(folderPath, '/', fileName))
for(i in 1:10){
  #nrow(ictGeneExp)){
  tryCatch({
    #print(nsGenesHighCor[i, 1])
    geneExp <- data.frame('sample' = colnames(ictGeneExp)[2:ncol(ictGeneExp)], 
                          'Exp' = t(ictGeneExp[i, 2: ncol(ictGeneExp)]))
    
    colnames(geneExp)[2] <- 'Exp'
    
    ictClinical <- merge(ictClinical, geneExp, by = 'sample')
    ########getting p-values with different thresholds
    
    quanPval <- testingKMForDifferentValues(ictClinical)
    # checking if all the values in a vector are NAs
    if(all(is.na(quanPval$Pval))){
      print('TRUE')
      ictClinical <- ictClinical[, - (5:ncol(ictClinical))]
      #rm(km.fit)
      #j <- c(j, i)
    }
    print(i)
    minPval <- quanPval[quanPval$Pval %in% min(quanPval$Pval, na.rm = TRUE), ]
    try(ictClinical <- transform(ictClinical, 
                                 'exp' = ifelse(ictClinical$Exp >= quantile(ictClinical$Exp, as.numeric(minPval$quan)), 'high', 'low')), silent = TRUE)
    # fitting the data using the survfit function 
    
    km.fit <- survfit(Surv(time = OS.time, event = OS) ~ exp, data = ictClinical)
    
    # getting the p-value
    pVal <- surv_pvalue(km.fit, data = ictClinical)[4]
    # ploting the fit and adding the legend
    plot(km.fit, col = c('red', 'blue'), xlab = 'Days', ylab = 'Survival%',
         main = row.names(ictGeneExp)[i], mark.time = TRUE,  lwd = 2.5, cex.lab =1.5)
    legend('topright', legend = paste0(levels(as.factor(ictClinical$exp)) , '  n = '  , ifelse(levels(as.factor(ictClinical$exp)) == 'high', length(which(ictClinical$exp %in% 'high')), length(which(ictClinical$exp %in% 'low') ))),
           col = c('red', 'blue'), lty = 1, lwd = 2,  box.lty = 0.0, cex = 1.5)
    text(700, 0.1, pVal )
    text(500, 0.15, minPval$quan)
    text(500, 0.2, ifelse((summary(km.fit)$table[1, 'median'] / summary(km.fit)$table[2, 'median']) >= 1, 'good', 'poor') )
    # saving the results
    pval <- c(pval, surv_pvalue(km.fit, data = ictClinical)[1, 2])
    print(pval)
    quan <- c(quan, unlist(survdiff(Surv(time = OS.time, event = OS) ~ exp, data = ictClinical)[5]))
    ratio <- c(ratio, summary(km.fit)$table[1, 'median'] / summary(km.fit)$table[2, 'median'])
    q <- c(q, as.numeric(minPval$quan))
    gene <- c(gene, row.names(ictGeneExp)[i])
    # removing the gene information from sckmSurvivalData
    ictClinical <- ictClinical[, - (5:ncol(ictClinical))]
    print(head(ictClinical))
    rm(km.fit)
  }, error=function(cond) {
    message(cond)
  })
}

resStatICT <- data.frame('Genes' = gene, 'ChiSq' = quan ,
                         'P-value' = pval, 'Ratio' = ratio, 
                         'BestQuantile' = q, 'adjPval' = p.adjust(pval), 
                         row.names = gene)


