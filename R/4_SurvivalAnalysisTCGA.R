# Perform Survival analysis (Kaplan Meier Plots) for the gene expression values (FPKM_UQ)
# from Melanoma cancer downloaded from TCGA (Xena). The survival analysis is performed
# by selecting the quantile having a minimum p-value from quartile 0.25 to 0.75.
# input: gene expression data, survival data
# output: Kaplan Meier plots, KM Statistics (p-value, selected quantile)

#libraries
library(survival)
library(survminer)
source('code/KMTesting.R')

# set working directory
#setwd('./PRIME')

# loading the data sets
load('raw_data/TCGA/skcmSurvival.RData')
load('raw_data/TCGA/geneExpressionDataFPKMUQ.RData')
dim(geneExpFPKMUQ)
dim(sckmSurvivalData)
# updating the data sets
row.names(geneExpFPKMUQ) <- geneExpFPKMUQ$Gene
sckmSurvivalData$sample <- gsub('-', '.', sckmSurvivalData$sample)

# selecting the same samples in gene expression and survival data sets
samples <- intersect(colnames(geneExpFPKMUQ), sckmSurvivalData$sample)
sckmSurvivalData <- sckmSurvivalData[sckmSurvivalData$sample %in% samples, ]
geneExpFPKMUQ <- geneExpFPKMUQ[, which(colnames(geneExpFPKMUQ) %in% samples)]

# making the order of the samples in both dataset same
sckmSurvivalData <- sckmSurvivalData[order(sckmSurvivalData$sample), ]
geneExpFPKMUQ <- geneExpFPKMUQ[, order(colnames(geneExpFPKMUQ))]
identical(sckmSurvivalData$sample, colnames(geneExpFPKMUQ))
geneExpFPKMUQ <- geneExpFPKMUQ[comGenes, ]
# running the KM testing for all the genes
quan <- NULL
ratio <- NULL
pval <- NULL
q <- NULL
gene <- NULL
folderPath <- '.z /TCGA/KM'

system(paste0('mkdir -p ', folderPath))

fileName <- 'KM_TCGA_All.pdf'
pdf(file = paste0(folderPath, '/', fileName))
for(i in 1:nrow(geneExpFPKMUQ)){
  tryCatch({
    print(i)
    geneExp <- data.frame('sample' = colnames(geneExpFPKMUQ)[2:ncol(geneExpFPKMUQ)], 
                          'Exp' = t(geneExpFPKMUQ[i, 2: ncol(geneExpFPKMUQ)]))
    
    colnames(geneExp)[2] <- 'Exp'
    
    sckmSurvivalData <- merge(sckmSurvivalData, geneExp, by = 'sample')
    
    ########getting p-values with different thresholds
    
    quanPval <- testingKMForDifferentValues(sckmSurvivalData)
    #print(quanPval)
    # checking if all the values in a vector are NAs
    if(all(is.na(quanPval$Pval))){
      print('TRUE')
      sckmSurvivalData <- sckmSurvivalData[, 1:4]
      #rm(km.fit)
      #j <- c(j, i)
    }
    minPval <- quanPval[quanPval$Pval %in% min(quanPval$Pval, na.rm = TRUE), ]
    try(sckmSurvivalData <- transform(sckmSurvivalData, 
                                      'exp' = ifelse(sckmSurvivalData$Exp >= quantile(sckmSurvivalData$Exp, as.numeric(minPval$quan)), 'high', 'low')), silent = TRUE)
    # fitting the data using the survfit function 
    #print(i)
    km.fit <- survfit(Surv(time = OS.time, event = OS) ~ exp, data = sckmSurvivalData)
    
    # getting the p-value
    #print(i)
    pVal <- surv_pvalue(km.fit, data = sckmSurvivalData)[4]
    # ploting the fit and adding the legend
    plot(km.fit, col = c('red', 'blue'), xlab = 'Days', ylab = 'Survival%',
         main = row.names(geneExpFPKMUQ)[i], mark.time = TRUE,  lwd = 2.5, cex.lab =1.5)
    legend('topright', legend = paste0(levels(as.factor(sckmSurvivalData$exp)) , '  n = '  , ifelse(levels(as.factor(sckmSurvivalData$exp)) == 'high', length(which(sckmSurvivalData$exp %in% 'high')), length(which(sckmSurvivalData$exp %in% 'low') ))),
           col = c('red', 'blue'), lty = 1, lwd = 2,  box.lty = 0.0, cex = 1.5)
    text(700, 0.1, pVal )
    text(500, 0.15, minPval$quan)
    text(500, 0.2, ifelse((summary(km.fit)$table[1, 'median'] / summary(km.fit)$table[2, 'median']) >= 1, 'good', 'poor') )
    # saving the results
    pval <- c(pval, surv_pvalue(km.fit, data = sckmSurvivalData)[1, 2])
    quan <- c(quan, unlist(survdiff(Surv(time = OS.time, event = OS) ~ exp, data = sckmSurvivalData)[5]))
    ratio <- c(ratio, summary(km.fit)$table[1, 'median'] / summary(km.fit)$table[2, 'median'])
    q <- c(q, as.numeric(minPval$quan))
    gene <- c(gene, row.names(geneExpFPKMUQ)[i])
    # removing the gene information from sckmSurvivalData
    sckmSurvivalData <- sckmSurvivalData[, 1:4]
    rm(km.fit)
  }, error=function(cond) {
    message(cond)
  })
}
dev.off()
resStatTCGA <- data.frame('Genes' = gene, 'ChiSq' = quan ,
                          'P-value' = pval, 'Ratio' = ratio, 
                          'BestQuantile' = q, 'adjPval' = p.adjust(pval), 
                          row.names = gene)
saveRDS(resStatTCGA, 'raw_data/TCGA/SurvivalStat.RData')
