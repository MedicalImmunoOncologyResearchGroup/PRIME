# Ploting the circos plot 

categories <- read.table('raw_data/categories.txt', sep = '\t', header = FALSE)
colnames(categories) <- c('T_Cell_Dysfunction', 'T_Cell_Acivation', 'Antigen_Presentation',
                          'Immune_Suppression', 'Immune_Activation', 'Autoimmunity', 'Epigenetic',
                          'Metabolism', 'Unassigned', 'Categories')
MyMerge <- function(x, y){
  df <- merge(x, y, by= "Genes", all.x= TRUE, all.y= TRUE)
  return(df)
}
# setting working directory
set.seed(100)
setwd('./PRIME/')
#setwd('/Users/arfamehmood/Documents/Projects_Saara')
#libraries
library(circlize)
library(reshape2)
library(dplyr)

## GSP GENES
# getting the GSP Genes sublocal cellurization data
gspLocal <- read.table('raw_data/GSP_Local.txt', header = TRUE, sep = '\t')

gspLocal$CYTOSKELETON[which(gspLocal$CYTOSOL == 'X')] <- 'X'
gspLocal$CYTOSOL <- NULL
gspLocal$ER[unique(c(which(gspLocal$GOLGI == 'X'), which(gspLocal$VESICLES == 'X')))] <- 'X'

gspLocal$GOLGI <- NULL
gspLocal$VESICLES <- NULL

## 
# res
tcgaStats <- read.table(file = 'processed_data/TCGA/TCGA_Sur_Sig_Spearman_cor.csv',
                        header = TRUE)
ictStats <- read.table(file = 'processed_data/ICT/ICT_Sur_Sig_Spear_cor.csv',
                       header = TRUE)
# Reading in TCGA data
#ictGeneExp <- read.table('Project_CD74/OtherCohorts/TPM/BatchCorrectedTPM/CombinedTPM_Batch_Corrected.csv',
#                         header = TRUE, row.names = 1)
load('raw_data/TCGA/geneExpressionDataFPKMUQ.RData')
row.names(geneExpFPKMUQ) <- geneExpFPKMUQ$Gene
colnames(geneExpFPKMUQ)[1] <- 'Genes'

# getting the good prognosis genes expression values
gspGenesExp <- geneExpFPKMUQ[gspLocal$GeneName, ]

gspGenesExp$Variance <- apply(gspGenesExp[, 2:ncol(gspGenesExp)], 1, var)

gspGenesExp$log2Var <- log2(gspGenesExp$Variance)

#### PSP GENES 
# getting the poor prognosis local sub cellularization data

pspLocal <- read.table('raw_data/PSP_Local.txt', header = TRUE, sep = '\t')

# combining the data of 
pspLocal$CYTOSKELETON[which(pspLocal$CYTOSOL == 'X')] <- 'X'
pspLocal$CYTOSOL <- NULL
pspLocal$ER[unique(c(which(pspLocal$GOLGI == 'X'), which(pspLocal$VESICLES == 'X')))] <- 'X'

pspLocal$GOLGI <- NULL
pspLocal$VESICLES <- NULL


pspGenesExp <- geneExpFPKMUQ[pspLocal$GeneName, ]

pspGenesExp$Variance <- apply(pspGenesExp[, 2:ncol(pspGenesExp)], 1, var)

pspGenesExp$log2Var <- log2(pspGenesExp$Variance)

# getting the categories
dirPath <- './'

## reading nanostring data
nsGenes <- read.table(paste0(dirPath, './allPanelsUniqueFinal.csv'), header = TRUE, sep = '\t', check.names = F)
colnames(nsGenes)[c(1:4, 209)] <- c('Genes', 'Naive', 'NaiveTcells', 'dTcellsfunctions', 'NF-B')

# getting the psp genes from nanostring gene panel

pspNet <- nsGenes[nsGenes$Genes %in% pspLocal$GeneName,]
pspStats <- tcgaStats[tcgaStats$Genes %in% pspNet$Genes, ]
pspStats <- merge(pspStats, pspGenesExp[, c(1, 474)], by =  'Genes')

pspLocal <- data.frame('Genes' = pspLocal$GeneName, ifelse(pspLocal[, 2:ncol(pspLocal)] == 'X', 1, 0))

pspNet <- data.frame('Genes' = pspNet$Genes , ifelse(pspNet[, 2:ncol(pspNet)] == '-', 0, 1))
merge(pspStats, pspNet, by = 'Genes')

pspFinal <- Reduce(MyMerge, list(pspStats, pspNet)) #pspLocal[pspLocal$Genes %in% pspNet$Genes, ]))
row.names(pspFinal) <- pspFinal$Genes
#pspFinalM  <- melt(pspFinal[, c(1,16:ncol(commonGenes))] , id.vars = 'Genes')


# Clustering of genes

pspFinalClust <- hclust(dist(as.matrix(pspFinal[, c(2, 4, 8, 10, 12)])))
plot(pspFinalClust)
row.names(pspFinal)[pspFinalClust$order]

# plotting circos plot whithout chord daigram 

pspFinalM  <- melt(pspFinal[, c(1,13:ncol(pspFinal))] , id.vars = 'Genes')
pspFinalM <- pspFinalM[pspFinalM$value != 0, ]

# arranging the genes according to the clustering
pspFinalM <- pspFinalM %>% arrange(factor(pspFinalM$Genes, 
                                          levels = row.names(pspFinal)[pspFinalClust$order]),
                                   desc(variable), desc(value))

factors = 1:2
circos.clear()
circos.par(cell.padding = c(0, 0, 0, 0)) 
circos.initialize(factors, xlim = c(0,1)) 
circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.05, bg.col = rand_color(2), bg.border = NA ) 


cdm_plt <- chordDiagram(pspFinalM, grid.col = NA, directional = TRUE, 
                        annotationTrack =  c('grid'), transparency = 0,
                        #annotationTrackHeight = 0.1#,
                        link.target.prop = FALSE, preAllocateTracks = list(track.height = 0.1))

circos.track(track.index = 1, panel.fun = function(x, y) {
  #print(get.cell.meta.data("sector.index"))
  si = get.cell.meta.data("sector.index")
  print(si)
  circos.text(CELL_META$xcenter, CELL_META$ylim[1] + 2, CELL_META$sector.index,
              facing = "clockwise", cex = 0.20 , niceFacing = TRUE, adj = c(0, 0.5))
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
}, bg.border = NA)
###############################################
## Making the circos plot
#changing the shape of the matrix

commonGenes <- merge(commonGenes, commonGenesExp[, c(1, 474, 475)], by = 'Genes')

colnames(nsGenes)[c(1:4, 209)] <- c('Genes', 'Naive', 'Naive T cells', 'dT cells functions', 'NF-B')
nsGenesCG <- nsGenes[nsGenes$Genes %in% commonGenes$Genes, ] 
nsGenesCG <- data.frame('Genes' = nsGenesCG$Genes , ifelse(nsGenesCG[, 2:ncol(nsGenesCG)] == '-', 0, 1))
#
commonGenes <- merge(commonGenes, nsGenesCG, by = 'Genes')
commonGenesM  <- melt(commonGenes[, c(1,16:ncol(commonGenes))] , id.vars = 'Genes')
commonGenesM <- commonGenesM[commonGenesM$value != 0, ]
commonGenesM <- merge(commonGenesM, commonGenes[, c(14, 15, 1:13)], by = 'Genes')






## marking the chord color
commonGenesM$color <- 'sienna3'
commonGenesM[commonGenesM$adjPval >= 0.01 &  commonGenesM$adjPval  <= 0.05, 'color'] <- 'lightgrey'
commonGenesM[(commonGenesM$adjPval > 0.001 &  commonGenesM$adjPval <= 0.01), 'color'] <- 'darkgrey'
commonGenesM[(commonGenesM$adjPval < 0.0001), 'color'] <- 'black'

chordCol <- commonGenesM$color
names(chordCol) <- commonGenesM$Genes

## marking the correlation color

commonGenesM$cd74color <- 'navy'
commonGenesM[(abs(commonGenesM$CD74Cor >= 0.5) &  abs(commonGenesM$CD74Cor) < 0.75 ), 'cd74color'] <- 'skyblue'
commonGenesM[(abs(commonGenesM$CD74Cor) >= 0.75), 'cd74color'] <- 'navy'

commonGenesM$cd8Acolor <- 'navy'
commonGenesM[(abs(commonGenesM$CD8ACor >= 0.5) &  abs(commonGenesM$CD8ACor) < 0.75 ), 'cd8Acolor'] <- 'skyblue'
commonGenesM[(abs(commonGenesM$CD8ACor) >= 0.75), 'cd8Acolor'] <- 'navy'



circos.clear()
cdm_plt <- chordDiagram(commonGenesM[1:100, 1:4], link.decreasing = FALSE, link.sort = TRUE,  
                        grid.col = NA, #directional = TRUE, 
                        col = chordCol,  
                        annotationTrack =  c('grid'), transparency = 0
                        #annotationTrackHeight = 0.1#,
                        #link.target.prop = FALSE, preAllocateTracks = list(track.height = 0.1)
)

cnm <- 19
cn8a <- 20
circos.track(track.index = 1, panel.fun = function(x, y) {
  #print(get.cell.meta.data("sector.index"))
  si = get.cell.meta.data("sector.index")
  print(si)
  circos.text(CELL_META$xcenter, CELL_META$ylim[1] + 2, CELL_META$sector.index,
              facing = "clockwise", cex = 0.20 , niceFacing = TRUE, adj = c(0, 0.5))
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  # print(ylim[1])
  # print(ylim[2])
  gn = commonGenesM[commonGenesM$variable %in% si, 1]
  if(length(gn) > 1){
    for(i in 1:length(gn)){
      print(paste0(gn[i],'and' ,si))
      con <- cdm_plt[which(cdm_plt$rn %in% gn[i] & cdm_plt$cn %in% si),]
      # print(con)
      # print(con[, 'x2'])
      # print(con[, 'x2']-1)
      circos.rect(con[, 'x2'], ylim[1], 
                  con[, 'x2'] - abs(con[, 'value2']), 0.75,#ylim[2],
                  col = commonGenesM[commonGenesM$variable %in% si & commonGenesM$Genes %in% gn[i], cnm] , border = NA)
      circos.rect(con[, 'x2'], 0.90,
                  con[, 'x2'] - abs(con[, 'value2']), 1.6,
                  col = commonGenesM[commonGenesM$variable %in% si & commonGenesM$Genes %in% gn[i], cn8a] , border = NA)
      
    }
  } else{
    breaks = seq(xlim[1], xlim[2], by = 0.1)
    n_breaks = length(breaks)
    circos.rect(breaks[-n_breaks], ylim[1],#rep(ylim[1], n_breaks - 1),
                breaks[-1], 0.75, #rep(ylim[2], n_breaks - 1),
                col = commonGenesM[commonGenesM$Genes %in% si, cnm] , border = NA)
    circos.rect(breaks[-n_breaks], 0.90,
                breaks[-1], 1.6,
                col = commonGenesM[commonGenesM$Genes %in% si, cn8a] , border = NA)
    
    
    if(length(gn) ==1){
      circos.rect(breaks[-n_breaks], ylim[1], #rep(ylim[1], n_breaks - 1),
                  breaks[-1], 0.75,#rep(ylim[2], n_breaks - 1),
                  col = commonGenesM[commonGenesM$Genes %in% gn[1], cnm] , border = NA)
      circos.rect(breaks[-n_breaks], 0.90,
                  breaks[-1], 1.6,
                  col = commonGenesM[commonGenesM$Genes %in% gn[1], cn8a] , border = NA)
      
    }
  }
  
}, bg.border = NA)


