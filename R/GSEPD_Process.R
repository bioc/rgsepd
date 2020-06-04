
#a simple interface to DESeq, 
#assuming you have a finalCounts matrix with 'samples' columns
#part of the rgsepd-0.3.3 package KStamm June 2014

GSEPD_Process <- function(GSEPD){
  
  if(is.null(GSEPD$Conditions))
    stop("Must specify conditions to test with GSEPD_ChangeConditions before GSEPD_Process")
  
  #this is a fine time to clean up the ID converter
  #since you've provided all relevant transcript IDs.
  EnsureCleanIDs(rownames(GSEPD$finalCounts))
  
  GSEPD$Conditions <- sort(GSEPD$Conditions)
  sampleMeta <- GSEPD$sampleMeta
  conds <- GSEPD$Conditions
  
  #this would be a good place to make sure "conds" are in sampleMeta$Conditions
  if(!(all(conds %in% sampleMeta$Condition))) 
    stop("specified conditions not found in sample metadata table")
  
  samples <- sampleMeta$Sample[sampleMeta$Condition %in% conds]
  samples_conds <-  sampleMeta$Condition[sampleMeta$Condition %in% conds]
  
  #Count the subjects of each type and affix xN to their types.
  tconds=table(samples_conds)
  conds=paste(samples_conds,tconds[samples_conds],sep=GSEPD$C2T_Delimiter)
  conds <- factor( conds)
  C2T = sort(  unique(as.character(conds))  )  
  GSEPD$C2T<-C2T
  
  GSEPD <- GSEPD_CheckCounts(GSEPD)
  if(!file.exists(DESEQ_CountsFile(GSEPD))) {
    if(!GSEPD$QUIET)Message_Generate(DESEQ_CountsFile(GSEPD), TRUE)
    write.csv( GSEPD$normCounts, DESEQ_CountsFile(GSEPD))
  }else{
    #file already exists
    if(!GSEPD$QUIET)Message_Generate(DESEQ_CountsFile(GSEPD), FALSE)
  }
  
  if(!file.exists(GSEPD_HMA_File(GSEPD))){ #the final HMA is not done, make it from
    if(!file.exists(GSEPD_MFile(GSEPD))){ #if the M file is not done, make it
      if(!file.exists(GSEPD_GO2File(GSEPD))){ #if the GO2 file is not done, make it
        if(!file.exists(DESEQ_AFile(GSEPD))){ #if the A file is not done, make it
          if(!file.exists(DESEQ_RFile(GSEPD))){ #if the RES file is not done, make it
            ProcessDESEQ(GSEPD)
          }
          AnnotateTable(GSEPD)
          GSEPD_DEGHeatmap(GSEPD) ; #as soon as we made the .Annote file
        }
        suppressWarnings({ #dont worry about GOSeq's complaints. 
          AnnotateTable.GO(GSEPD)
        })
      }
      FileMerge(GSEPD)
    }
    
    tryCatch({ #sometimes the PCA fails...
      GSEPD_PCA_Plot(GSEPD) ## make another picture
      }, warning = function(w) {
        warning(w)
      }, error = function(e) {
        warning(e) #but it's no big deal if a picture isn't made.
      })
    
    GSEPD_ProjectionProcessor(GSEPD) #and do a lot of work.
  }
  MergeResultsTables(GSEPD)
  if(!GSEPD$QUIET) message("All Done!")
  #finally return
  GSEPD
}

MergeResultsTables <- function(G){
#finally, to match the resulting HMA we want to pick up exciting GO terms
  if(!G$QUIET)Message_Generate(GSEPD_HMACSV_File(G))
  #we'll produce the filtered file because .Alpha and .Beta already contain everything. 
  # don't really need extra filters for .Beta linearity, just making the shortlist final output
  GO<-read.csv(GSEPD_GOFile(G),as.is=TRUE,header=TRUE,row.names=1)
  SEG <- read.csv(GSEPD_Seg_File(G),as.is=TRUE,header=TRUE,row.names=1) #these row names are go terms
  
  Limit_GO  <- ifelse(G$LIMIT$HARD , G$LIMIT$GO_PVAL , max(c(median(GO$over_represented_pvalue) , G$LIMIT$GO_PVAL)) )
  Limit_Seg <- ifelse(G$LIMIT$HARD , G$LIMIT$Seg_P   , max(c(median(SEG$Segregation_PV) , G$LIMIT$Seg_P)))
  
  OM<-merge( subset(GO, GO$over_represented_pvalue <= Limit_GO),
             subset(SEG,SEG$Segregation_PV <= Limit_Seg), by.x="category", by.y="row.names")
  write.csv(OM,GSEPD_HMACSV_File(G));  
  
}

ProcessDESEQ <- function(G, 
                         deseq.contrast = "NULL",
                         deseq.lfcShrink = "normal"){
  if(missing(deseq.contrast))
    deseq.contrast=c("Condition",G$Conditions[1],G$Conditions[2])
  
  dds <- DESeq(GSEPD_Export_DESeq(G))
  res <- results(dds,
                 contrast=deseq.contrast,
                 lfcThreshold = as.numeric(G$LIMIT$LFC)[1])
  res <- lfcShrink(dds, 
                   type=deseq.lfcShrink,
                   contrast=deseq.contrast,
                   res = res)
  
#  message("trying to drop genes with zero baseMean");
#  res = subset(res,res$baseMean > 0.01)
#  message(paste("Filtered entries remaining: ",length(res$id)));
  
  #in DESeq1 we had columns  
  #  id  baseMean	baseMeanA	baseMeanB	foldChange	log2FoldChange	pval	padj
  #now we have columns 
  #   baseMean log2FoldChange      lfcSE       stat        pvalue          padj
  #so to mash it into the old format:
  res$pvalue[is.na(res$pvalue)] <- 1;#DESeq2 gave me NA's which broke later filters
  res$padj[is.na(res$padj)] <- 1; #so we'll mark them clearly not significant.
  res$id<-rownames(res); res$pval<-res$pvalue;res$foldChange<-2^res$log2FoldChange
  #it's a little more work to generate the C2T columns as baseMeanA/B.
  Columns.A <- which(colnames(G$finalCounts) %in% 
                       subset(G$sampleMeta, G$sampleMeta$Condition==G$Conditions[1])$Sample)
  Columns.B <- which(colnames(G$finalCounts) %in% 
                       subset(G$sampleMeta, G$sampleMeta$Condition==G$Conditions[2])$Sample)

  if(length(Columns.A)>1){
    res$baseMeanA <- apply( G$normCounts[,Columns.A],1,mean) #but normCounts are logspace, so this mean is biased low...
  }else{
    res$baseMeanA <- G$normCounts[,Columns.A] # mean of one element
  }
  if(length(Columns.B)>1){
    res$baseMeanB <- apply( G$normCounts[,Columns.B],1,mean)
  }else{
    res$baseMeanB <- G$normCounts[,Columns.B] # mean of one element
  }

  outfilename=paste(G$Output_Folder,"/DESEQ.Volcano.",G$C2T[1],".",G$C2T[2],".png",sep="")
  if(!G$QUIET)Message_Generate(outfilename)
  png(filename=outfilename); # volcano with 30K points was too slow as PDF
  plotDE( res )
  dev.off();
  
  if(!G$QUIET)Message_Generate(DESEQ_RFile(G))
  write.csv(res,file=DESEQ_RFile(G));
}

#run the above, iterated by possible pairings
GSEPD_ProcessAll <- function(G){
  X<-as.character(unique(G$sampleMeta$Condition))
  Xn<-length(X)
  message(sprintf("Preparing to run %d conditions as %d pairs",Xn,Xn*(Xn-1)/2))
  for(i in 1:(Xn-1))
    for(j in (i+1):Xn){
      G<-GSEPD_ChangeConditions(G,c(X[i],X[j]))
      G<-GSEPD_Process(G)
    }
  G
}

