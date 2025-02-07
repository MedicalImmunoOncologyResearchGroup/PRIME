# Set working directory 
setwd('C:/Users/arfam/Desktop/Arfa_Mehmood/MIORG/NanoStringPanels/')
library('readxl')
library('ggplot2')
library('dplyr')

## reading immunology panels
myeloidInnateImmunity <- read_excel('ImmunityPanel/LBL-10397-01_nCounter_Hs_Myeloid_Innate_Immunity_V2.xlsx', 
                                    sheet = 'Annotations', skip = 1)
myeloidInnateImmunity <- myeloidInnateImmunity[-nrow(myeloidInnateImmunity), ]

immunologyPanelGL <- read_excel('ImmunityPanel/LBL-10540-02_nCounter_NHP_Immunology_V2_Panel_Gene_List.xlsx', 
                                sheet = 'Annotations', skip = 1)
immunologyPanelGL <- immunologyPanelGL[-nrow(immunologyPanelGL), ]
immunologyPanelGL <- immunologyPanelGL[, - 2]

humanAutoImmunePanel <- read_excel('ImmunityPanel/LBL-10560-02_Human_AutoImmune_Profiling_Panel.xlsx', 
                                   sheet = 'Annotations', skip = 1)
humanAutoImmunePanel <- humanAutoImmunePanel[-nrow(humanAutoImmunePanel), ]
humanAutoImmunePanel <- humanAutoImmunePanel[, -2]

fibrosisConsortiumGL <- read_excel('ImmunityPanel/LBL-10606-01_nCounter_Fibrosis_Consortium_Panel_Gene_List.xlsx', 
                                   sheet = 'Annotations', skip = 1)
fibrosisConsortiumGL <- fibrosisConsortiumGL[-nrow(fibrosisConsortiumGL), ]

metabolicPathwayPanelGL<- read_excel('ImmunityPanel/LBL-10725-01_Metabolic_Pathways_Panel_Gene_List.xlsx', 
                                     sheet = 'Annotations', skip = 1)
metabolicPathwayPanelGL<- metabolicPathwayPanelGL[-nrow(metabolicPathwayPanelGL), ]
metabolicPathwayPanelGL <- metabolicPathwayPanelGL[, -2]

organTransplantPanel <- read_excel('ImmunityPanel/LBL-10743-01_Human_Organ_Transplant_Panel.xlsx', 
                                   sheet = 'Annotations', skip = 1)
organTransplantPanel<- organTransplantPanel[-nrow(organTransplantPanel), ]
organTransplantPanel <- organTransplantPanel[, -2]

hostResponseGL <- read_excel('ImmunityPanel/LBL-10805-01_nCounter_Host_Response_Gene_List.xlsx', 
                             sheet = 'Annotations', skip = 1)
hostResponseGL <- hostResponseGL[-nrow(hostResponseGL), ]
hostResponseGL <- hostResponseGL[, -2]

humanImmunologyGL <- read_excel('ImmunityPanel/LBL-C0269-03_nCounter-Human-Immunology-V2-Panel-Gene-List.xlsx', 
                                sheet = 'Annotations', skip = 1)
humanImmunologyGL <- humanImmunologyGL[-nrow(humanImmunologyGL), ]


humanFibrosisGL <- read_excel('ImmunityPanel/LD-00171-02-Human-Fibrosis-V2-Gene-List.xlsx', 
                              sheet = 'Annotations', skip = 1)
humanFibrosisGL <- humanFibrosisGL[-nrow(humanFibrosisGL), ]
humanFibrosisGL <- humanFibrosisGL[, -2]

tcrDiversityGL <- read_excel('ImmunityPanel/LD-00408-02-TCR-Diversity-Gene-List.xlsx', 
                             sheet = 'Annotations', skip = 1)
tcrDiversityGL <- tcrDiversityGL[-nrow(tcrDiversityGL), ]

## reading oncology panel
io360FunctionalAnnotation <- read_excel('OncologyPanel/LBL-10498-02_IO_360_Gene_List.xlsx', 
                                        sheet = 'Functional Annotations', skip = 1)
io360FunctionalAnnotation <- io360FunctionalAnnotation[-nrow(io360FunctionalAnnotation), ]
io360FunctionalAnnotation <- io360FunctionalAnnotation[, -2]
io360CancerImmunityCycle <- read_excel('OncologyPanel/LBL-10498-02_IO_360_Gene_List.xlsx', 
                                       sheet = 'Cancer-Immunity Cycle', skip = 1)
io360CancerImmunityCycle  <- io360CancerImmunityCycle[-nrow(io360CancerImmunityCycle), ]
io360CancerImmunityCycle <-io360CancerImmunityCycle[, -2] ## removing cell type column

oncoTCRDiversityGL <- read_excel('OncologyPanel/LD-00408-02-TCR-Diversity-Gene-List.xlsx', 
                                 sheet = 'Annotations', skip = 1)
oncoTCRDiversityGL <- oncoTCRDiversityGL[-nrow(oncoTCRDiversityGL), ]

## Reading BAP1 uveal file

BAP1 <- read_excel('path5384-sup-0005-tables1.xlsx', 
                   sheet = 'Combined Immune Categories', skip = 1)
BAP1 <- BAP1[1:168, -c(1, 3, 4, 70)]
colnames(BAP1)[1] <- 'Gene'

## Merging the immunology panels 

mergedPanels <- function(panelsName){
  mergedPanels <- Reduce(
    function(x, y, ...) merge(x, y, all = TRUE, ...),
    panelsName
  )
}

mergedPanelsByGene <- function(panelsName){
  mergedPanels <- Reduce(
    function(x, y, ...) merge(x, y, by = 'Gene', all = TRUE, ...),
    panelsName
  )
}

## Merging the immunity Panel, Oncology Panel and BAP1


allPanelsList <-list(myeloidInnateImmunity, immunologyPanelGL, humanAutoImmunePanel, 
                     fibrosisConsortiumGL, hostResponseGL, humanFibrosisGL,
                     humanImmunologyGL,metabolicPathwayPanelGL, 
                     organTransplantPanel, tcrDiversityGL, io360CancerImmunityCycle,
                     io360FunctionalAnnotation, oncoTCRDiversityGL, BAP1)
allPanels <- mergedPanels(allPanelsList) # merging the panels

## sorting the colnames of immunityPanels
allPanels <- allPanels %>% select(sort(names(allPanels)))
i <- which(colnames(allPanels) == 'Gene') ## getting the column number of Gene
allPanels <- allPanels[, c(i, 1:(i-1),(i+1):ncol(allPanels))] ##bringing the gene column to the start

## saving the results with NAs
write.table(allPanels, 'CombinedPanels/allPanelsRaw.csv', sep = '\t',
            row.names = FALSE, quote = F )


#allPanels[is.na(allPanels)] <- '-' # replacing NAs with - 7351 rows

allPanels <- allPanels %>% distinct() # removing duplicate rows 7351 rows
duplicateGenes <- unique(allPanels$Gene[duplicated(allPanels$Gene) ]) # 1609 rows
allPanelsDup <- allPanels[allPanels$Gene %in% duplicateGenes, ] # 5419 rows and 1593 genes
allPanelsUnique <- allPanels[allPanels$Gene %in% setdiff(allPanels$Gene, duplicateGenes), ] # 1753 genes

intersect(allPanelsUnique$Gene, unique(allPanelsDup$Gene))

#################################


checkNA <- function(x){
  ind <- which(!is.na(x))
  
  if(length(ind) > 0 & length(unique(x[ind])) == 1){
    #print(ind)
    rep(x[ind[1]], length(x))
  } else if (length(ind) == 0){
    #print(ind)
    c(rep(NA, length(x)))
  }else{
    #print(x)
    return(x)
  }
}

checkNAUpdated <- function(x){
  ind <- which(!is.na(x))
  
  if(length(ind) > 0){
    #print(ind)
    rep(x[ind[1]], length(x))
  } else if (length(ind) == 0){
    #print(ind)
    c(rep(NA, length(x)))
  }else{
    #print(x)
    return(x)
  }
}

getData <- function(geneCat){
  #print(geneCat)
  as.data.frame(apply(geneCat[, 2:ncol(geneCat)], 2, checkNAUpdated) )
  
}
allPanelsSort <- NULL
#aggregate(allPanels[1:20, ], by = list(allPanels$Gene[1:20]), FUN = function(allPanels) getData(allPanels))
for(gene in unique(allPanelsDup$Gene)){
  allPanelsSort <-  rbind(allPanelsSort, getData(allPanelsDup[allPanelsDup$Gene %in% gene,  ])  )
}
allPanelsSort$Gene <- allPanelsDup$Gene

allPanelsSortUnique <- allPanelsSort %>% distinct()
duplicateGenes <- unique(allPanelsSortUnique$Gene[duplicated(allPanelsSortUnique$Gene) ])
allPanelsSortUnique <- allPanelsSortUnique[allPanelsSortUnique$Gene %in% setdiff(allPanelsSortUnique$Gene, duplicateGenes), ]
#allPanelsSortUnique <- allPanelsSortUnique[, c(ncol(allPanelsSortUnique), 1:(ncol(allPanelsSortUnique)-1) )]

## making the colnames in the same order and making the 2 dataframes to 1


allPanelsUnique <- allPanelsUnique[, order(colnames(allPanelsUnique))]
allPanelsSortUnique <- allPanelsSortUnique[, order(colnames(allPanelsSortUnique))]
allPanelsUnique <- allPanelsUnique[ ,c(2, 1, 3:ncol(allPanelsUnique))]
allPanelsFinal <- rbind(allPanelsUnique, allPanelsSortUnique)
i <- which(colnames(allPanelsFinal) %in% 'Gene')
allPanelsFinal <- allPanelsFinal[, c (i, 1: (i-1), (i+1) : ncol(allPanelsFinal))]
allPanelsFinal[is.na(allPanelsFinal)] <- '-' 
## Saving unique genes
write.table(allPanelsFinal, 'CombinedPanels/allPanelsUniqueRawFinal3346.csv', sep = '\t',
            row.names = FALSE, quote = F)



write.table(allPanelsFinal, 'CombinedPanels/allPanelsUniqueFinal.csv', sep = '\t',
            row.names = FALSE, quote = F)

allPanelsFinal <- read.table('CombinedPanels/allPanelsUniqueFinal.csv', header = TRUE,sep = '\t', fill = TRUE)

## sorting duplicated genes



allPanelsSortDup[is.na(allPanelsSortDup)] <- '-' 
write.table(allPanelsSortDup[, c(337, 1:336)], 'CombinedPanels/allPanelsDupFinal.csv', sep = '\t',
            row.names = FALSE, quote = F)
