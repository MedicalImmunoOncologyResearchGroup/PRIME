## making scripts for Gide Data to perform Kallisto alignment.
## getting the sample files
samples <- list.dirs('./Data/fastq', full.names = FALSE)[-1]
path <- './Gide/Data'
path0 <- './Gide/Scripts'
cmdHead <- '#! /bin/bash'
## for Gide data
for(sn in samples){
  print(sn)
  
  cmdKal <- paste0(cmdHead, '\nmkdir ', path, './Kallisto/', sn, '\nconda run kallisto quant -t 3 -i ', 
                   path, 'Kallisto_index/homo_sapiens/transcriptome.idx -o ', path, './Kallisto/', sn, ' ',  
                   path, 'Gide/fastq/', sn, '/', sn, '_1.fastq.gz ', path, 'Gide/fastq/', sn, '/', sn, '_2.fastq.gz' )
  cat(cmdKal, file = paste0(path0, '/', sn,  '.sh'))
  
}

######################################
## Hugo Dataset
path1 <- './Hugo/Scripts'
## for Hugo dataset
samplesHugo <- read.table('../../Downloads/SraRunInfo.csv', sep = ',', header = TRUE)

for(sn in samplesHugo$Run){
  print(sn)
  
  cmdKal <- paste0(cmdHead, '\nmkdir ', path, '/Hugo/Kallisto/', sn, '\n/Users/arfamehmood/miniconda3/bin/conda run kallisto quant -t 3 -i ', 
                   path, 'Kallisto_index/homo_sapiens/transcriptome.idx -o ', path, 'Hugo/Kallisto/', sn, ' ',  
                   path, 'Hugo/fastq/', sn, '_1.fastq.gz ', path, 'Hugo/fastq/', sn, '_2.fastq.gz' )
  cat(cmdKal, file = paste0(path1, '/', sn,  '.sh'))
  
}


####################################
## RIAZ DATASET SCRIPTS 

path1 <- '.`/Riaz/'
samplesRiaz <- read.table('/Users/arfamehmood/Downloads/SraRunInfo (1).csv', sep = ',', header = TRUE)
cmdKal <- NULL
## creatingscript to download all SRA for Riaz dataset
for(i in 81:nrow(samplesRiaz)){
  #print(sn)
  
  cmdKal <- c(cmdKal, paste0('curl --output Downloads/',samplesRiaz$Run[i] ,' https://sra-pub-run-odp.s3.amazonaws.com/sra/',
                             samplesRiaz$Run[i], '/', samplesRiaz$Run[i], '\n' ))
  
  
}
cat(cmdHead, '\n', cmdKal, file = paste0(path1, '/Download81.sh'))


## script to convert SRA files to fastq files and zipping them
## removing the SRA files
path <- '/Users/arfamehmood/Documents/Projects_Saara/Project_CD74/OtherCohorts/Riaz/fastq/'

cmdKal <- NULL
## creatingscript to download all SRA for Riaz dataset
for(i in 81:nrow(samplesRiaz)){
  
  cmdKal <- c(cmdKal, paste0('conda run fasterq-dump --outdir ', path,  
  ' --split-3 ./',  samplesRiaz$Run[i], '\n', 
  'gzip ', path, samplesRiaz$Run[i], '_1.fastq', '\ngzip ', path, samplesRiaz$Run[i], '_2.fastq',
  '\nrm ./', samplesRiaz$Run[i], '\n'))

}
cat(cmdHead, '\n', cmdKal, file = paste0(path1, '/ConvertingSRAToFastq81.sh'))


## writing scripts for Kallisto tool

path0 <- './Riaz/Scripts'
cmdHead <- '#! /bin/bash'
## for Riaz data
for(sn in samplesRiaz$Run){
  print(sn)
  
  cmdKal <- paste0(cmdHead, '\nmkdir ', path, '/Riaz/Kallisto/', sn, '\nconda run kallisto quant -t 3 -i ', 
                   path, 'Kallisto_index/homo_sapiens/transcriptome.idx -o ', path, 'Riaz/Kallisto/', sn, ' ',  
                   path, 'Riaz/fastq/',  '/', sn, '_1.fastq.gz ', path, 'Riaz/fastq/', sn, '_2.fastq.gz' )
  cat(cmdKal, file = paste0(path0, '/', sn,  '.sh'))
  
}

