
# demonstration using existing dataset.

if(FALSE){
  BiocManager::install("HSMMSingleCell")
  
  ??HSMMSingleCell
  
  ?HSMMSingleCell
  
  install.packages("monocle")
  BiocManager::install("monocle")
}
library(HSMMSingleCell)

library(monocle)

library(rgsepd)
library(AnnotationDbi)
library(org.Hs.eg.db)

x <- load_HSMM()
#rownames(x@assayData$exprs )
ensg = rownames(x@assayData$exprs )

x@phenoData@data$Media
table(x@phenoData@data$State, x@phenoData@data$Hours)

fixed_names = sapply(strsplit(ensg, ".", fixed=T), function(x) x[1])
lookedup <- select(org.Hs.eg.db,keytype="ENSEMBL",
                   keys=fixed_names, multiVals="first",
                   columns=c('REFSEQ'))

lookedup <- subset(lookedup, !is.na(REFSEQ))
lookedup <- subset(lookedup, !duplicated(ENSEMBL))
lookedup <- subset(lookedup, !duplicated(REFSEQ))

# THESE ARE FPKM from cufflinks 
fc = x@assayData$exprs
rownames(fc) <- fixed_names
fc <- fc[lookedup$ENSEMBL,]
#nrow(fc)
#fc[10:11,1:3]
#x@assayData$exprs[sapply(c("ENSG00000001167","ENSG00000001460"), function(i) grep(i, ensg) ),1:3]
rownames(fc) <- lookedup$REFSEQ
ncol(fc)
# 271 samples !


x@phenoData$Well

SM = data.frame(Sample = colnames(fc),
                Condition = sprintf("H%02d",as.numeric(as.character(x@phenoData$Hours))),
                SHORTNAME = paste( sapply(strsplit(colnames(fc), "_", fixed=T), function(x) x[1]),
                sapply(strsplit(colnames(fc), "_", fixed=T), function(x) x[3]), sep="_" ) ) 

head(SM)

# want to reduce the number of samples to some similar ones.
# my rgsepd system is failing when the groups are so diverse
# these timepoints are super varied
# some of these samples are outliers - different paths vs the group.


# select the middle-most??
pick_samples <- function(CONDITION){
  sC = t(log(fc[, SM$Condition==CONDITION] +0.1) )
  pox = princomp(t(sC))
  midsection <- function(x, thresh=0.2){
    j=which( x > quantile(x,thresh) & x < quantile(x,1-thresh))
    j
  }
  i=intersect(midsection(pox$loadings[,'Comp.1']),
              midsection(pox$loadings[,'Comp.2']))
  normal_samples = rownames(sC)[i]
  which(SM$Sample %in% normal_samples)
}


CUTDOWN=FALSE
if(CUTDOWN){
# then we'll subset for speed please.
#i = unlist(lapply( unique(SM$Condition), pick_samples))
i=which(SM$Sample %in% rownames(x@phenoData@data)[x@phenoData@data$State==2])
fc = fc[,i]
SM = SM[i,]
}

G <- GSEPD_INIT(Output_Folder="OUT",
                finalCounts=round(fc*1000),
                sampleMeta=SM,
                COLORS=c("green","black","red"))

#loosen statistics to ensure this example gets some gene sets 'enriched'
#G$LIMIT$LFC=0.1;
#G$LIMIT$HARD=FALSE
#G$LIMIT$GO_PVAL <- 0.5
#G$LIMIT$PADJ <- 0.1
#G$LIMIT$Seg_P <- 0.1
#G$LIMIT$baseMean <- 2

G$LIMIT$goseq_padj_method <- "none"

G <- GSEPD_ChangeConditions( G, c("H24","H48")) #set testing groups first!    
G <- GSEPD_Process( G ) #have to have processed results to plot them




