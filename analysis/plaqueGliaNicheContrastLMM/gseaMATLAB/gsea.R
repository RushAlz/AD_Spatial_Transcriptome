# Perform fgsea to each GWAS
# http://bioconductor.org/packages/devel/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html
rm(list = ls())

library(fgsea)
library(data.table)
library(stringr)

# Parameters
args = commandArgs(trailingOnly=TRUE)
gsPath = args[1]
gs = args[2]
tfolder = args[3]

# Load geneset
geneset <- gmtPathways(paste(gsPath,gs,sep=""))

# Apply GSEA
ranks <- read.table(paste(tfolder,"rankedPval.rnk",sep=""),header=TRUE,colClasses=c("character","numeric"))
ranks <- setNames(ranks$pS,ranks$geneS)

fgseaRes <- fgsea(geneset,ranks,scoreType="std") #,minSize=20,maxSize=1000)
fwrite(fgseaRes,file=paste(tfolder,str_remove(gs,"gmt"),"gsea.std.txt",sep=""),sep="\t",sep2=c("", " ", ""))

fgseaRes <- fgsea(geneset,ranks,scoreType="pos") #,minSize=20,maxSize=1000)
fwrite(fgseaRes,file=paste(tfolder,str_remove(gs,"gmt"),"gsea.pos.txt",sep=""),sep="\t",sep2=c("", " ", ""))

fgseaRes <- fgsea(geneset,ranks,scoreType="neg") #,minSize=20,maxSize=1000)
fwrite(fgseaRes,file=paste(tfolder,str_remove(gs,"gmt"),"gsea.neg.txt",sep=""),sep="\t",sep2=c("", " ", ""))

