

#later when I want 'global' variables in the rgsepd code, theyll be stored under pkg_globals
pkg_globals <- new.env();

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(sprintf("Loading R/GSEPD %s",packageDescription("rgsepd")$Version))
  packageStartupMessage("Building human gene name caches")
  
  #found the syntax for some kind of local/global variables here
  # http://stackoverflow.com/questions/12598242/global-variables-in-packages-in-r
  #assign("pkg_globals", new.env(), envir=parent.env(environment()))

  RefSeqs <- grep("^N[MR]_",
                  AnnotationDbi::keys(org.Hs.eg.db,keytype="REFSEQ"),
                  value=TRUE)
  lookedup <- select(org.Hs.eg.db,keytype="REFSEQ",
                     keys=RefSeqs, multiVals="first",
                     columns=c('SYMBOL','ENTREZID'))
  assign("R2H", hash(RefSeqs,lookedup$SYMBOL), pkg_globals)
  assign("R2E", hash(RefSeqs,lookedup$ENTREZID), pkg_globals)
  assign("H2R", hash(lookedup$SYMBOL,RefSeqs), pkg_globals)
  Egs <- AnnotationDbi::keys(org.Hs.eg.db,keytype="ENTREZID")
  lookedup <- select(org.Hs.eg.db,keytype="ENTREZID",
                     keys=Egs, multiVals="first",
                     columns=c('SYMBOL'))
  assign("H2E", hash(lookedup$SYMBOL,Egs), pkg_globals)
  assign("E2H", hash(Egs,lookedup$SYMBOL), pkg_globals)
         
}



## function definitions for GSEPD
## will run DESEQ and will run GOSeq
##then merge results, and process into GeneSetEnrichmentProfileDisplay

DESEQ_RFile <- function(G) 
  paste(G$Output_Folder,"/DESEQ.RES.",G$C2T[1],".",G$C2T[2],".csv",sep="")
DESEQ_AFile <- function(G) 
  paste(G$Output_Folder,"/DESEQ.RES.",G$C2T[1],".",G$C2T[2],".Annote.csv",sep="")
DESEQ_AFFile <- function(G) 
  paste(G$Output_Folder,"/DESEQ.RES.",G$C2T[1],".",G$C2T[2],".Annote_Filter.csv",sep="")
GSEPD_GOFile <- function(G) 
  paste(G$Output_Folder,"/GOSEQ.RES.",G$C2T[1],".",G$C2T[2],".GO.csv",sep="")
GSEPD_GOUpFile <- function(G) 
  paste(G$Output_Folder,"/GOSEQ.RES.",G$C2T[1],".",G$C2T[2],".GO-UP.csv",sep="")
GSEPD_GODownFile <- function(G) 
  paste(G$Output_Folder,"/GOSEQ.RES.",G$C2T[1],".",G$C2T[2],".GO-DOWN.csv",sep="")
GSEPD_GO2File <- function(G) 
  paste(G$Output_Folder,"/GSEPD.RES." , G$C2T[1],".",G$C2T[2],".GO2.csv",sep="")
GSEPD_MFile <- function(G) 
  paste(G$Output_Folder,"/GSEPD.RES." , G$C2T[1],".",G$C2T[2],".MERGE.csv",sep="")
GSEPD_HMA_File <- function(G) 
  paste(G$Output_Folder,"/GSEPD.HMA." , G$C2T[1],".",G$C2T[2],".pdf",sep="")
DESEQ_CountsFile <- function(G) 
  paste(G$Output_Folder,"/DESEQ.counts." , G$C2T[1],".",G$C2T[2],".csv",sep="")
GSEPD_Seg_File <- function(G) 
  paste(G$Output_Folder,"/GSEPD.Segregation_P." , G$C2T[1],".",G$C2T[2],".csv",sep="")
GSEPD_HMACSV_File <- function(G) 
  paste(G$Output_Folder,"/GSEPD.HMA." , G$C2T[1],".",G$C2T[2],".csv",sep="")
#new for 2016, let's calculate the distance to each centroid. 
GSEPD_HMG_File <- function(G) 
  paste(G$Output_Folder,"/GSEPD.HMG." , G$C2T[1],".",G$C2T[2],".pdf",sep="")
GSEPD_HMG1CSV_File <- function(G) 
  paste(G$Output_Folder,"/GSEPD.HMG1." , G$C2T[1],".",G$C2T[2],".csv",sep="")
GSEPD_HMG2CSV_File <- function(G) 
  paste(G$Output_Folder,"/GSEPD.HMG2." , G$C2T[1],".",G$C2T[2],".csv",sep="")


FileMerge<-function(GSEPD){
  
  if(!GSEPD$QUIET)Message_Generate(GSEPD_MFile(GSEPD))
  
  data.GO <- read.csv(GSEPD_GOFile(GSEPD),as.is=TRUE,header=TRUE)
  sGO1<-subset(data.GO[,c("category","over_represented_padj","Term")] , 
                   data.GO$over_represented_padj < (GSEPD$LIMIT$GO_PVAL * 1.25));
  sGO1$GOSEQ_DEG_Type <- rep("Mix", nrow(sGO1)) ; rm(data.GO)
  data.GO <- read.csv(GSEPD_GOUpFile(GSEPD),as.is=TRUE,header=TRUE)
  sGO2<-subset(data.GO[,c("category","over_represented_padj","Term")] , 
               data.GO$over_represented_padj < (GSEPD$LIMIT$GO_PVAL * 1.25));
  sGO2$GOSEQ_DEG_Type <- rep("Up", nrow(sGO2)) ; rm(data.GO)
  data.GO <- read.csv(GSEPD_GODownFile(GSEPD),as.is=TRUE,header=TRUE)
  sGO3<-subset(data.GO[,c("category","over_represented_padj","Term")] , 
               data.GO$over_represented_padj < (GSEPD$LIMIT$GO_PVAL * 1.25));
  sGO3$GOSEQ_DEG_Type <- rep("Down", nrow(sGO3)) ; rm(data.GO)
  
  sdata.GO <- rbind(sGO3,sGO2,sGO1) ; rm(sGO1,sGO2,sGO3)
  sdata.GO <- subset(sdata.GO, !duplicated(sdata.GO$category)) 
  
  data.GO2 <- read.csv(GSEPD_GO2File(GSEPD),as.is=TRUE,header=TRUE)
  
  tcats <- table(data.GO2$category)
  cats<-unique(data.GO2$category)
  cats=setdiff(cats,names(tcats)[tcats<GSEPD$MinGenesInSet | tcats>GSEPD$MaxGenesInSet])
  cats <- unique(c(cats,GSEPD$Force_GO_Include))

  mdata<-merge(sdata.GO, data.GO2, by.x="category",by.y="category")
  mdata<-subset(mdata, mdata$category %in% cats,
                select=c("category", "over_represented_padj.x","Term.x","GOSEQ_DEG_Type","HGNC", "LOG2.X.Y.","REFSEQ",
                         GSEPD$C2T[1],GSEPD$C2T[2], "PVAL","PADJ") )
  write.csv(mdata,GSEPD_MFile(GSEPD))
}


Message_Generate <- function(filename,overwrite=TRUE){
  if(file.exists(filename)){
    if(overwrite)
      message(sprintf("Generating %s, overwriting previous existing version.",filename))
    else
      message(sprintf("Would generate %s, but one is already present.",filename))
  }else{
    message(sprintf("Generating %s",filename))
  }
}

#Restrict elements of x to be between two numbers. 
cap<-function(x,L,H){x[x>H]=H;x[x<L]=L;x}

#Volcano Plot:
plotDE <- function( res ) {
  plot( 
    cap(res$log2FoldChange,-7,7), 
    0-log(res$pval,base=10), 
    pch=20, cex=.33, 
    col = ifelse( res$padj < .10, "red", "black" ),
    main=sprintf("Volcano for %d Genes",nrow(res)),
    xlab="Log2 FoldChange",
    ylab="Log10 p-value",
    sub="red points are significant with FDR < 10%")
}


#res is the table prepopulated with $id and other fields NA
LocalDB_Lookup <- function(res) 
  merge(res , build_Annotation_Matrix(res$id),
        by.x="id", by.y="REFSEQ", all.x=TRUE)


biomaRt_Lookup <- function(res) {  

  #now subset the table to those we didn't know yet.
  gns=unique( subset(res,is.na(res$ENTREZ))$id)
  #listMarts(host="www.ensembl.org")
  ensembl = useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org")
  ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
  attrs=listAttributes(ensembl)
  aa=c("refseq_mrna" , "chromosome_name" , "start_position" ,  "end_position" ,
       "entrezgene_id" ,      "hgnc_symbol",      "description" );
  
  #well of those that start with NM we want to query refseq_mrna
  #and for those that start withj NR we want to query refseq_ncrna
  gnsNM=gns[grep("NM_",gns)]
  if(length(gnsNM)>0){
    message(paste("Querying for",length(aa),"attributes from biomaRt on",length(gnsNM)," NM numbers"));
    bmt_NM=getBM(attributes = aa, filters = "refseq_mrna",values=gns[grep("NM_",gns)], mart = ensembl)
    message(paste("We got back",length(bmt_NM$refseq_mrna),"entries, which will be added to the local database."));
    if(nrow(bmt_NM)>0) {
    for(i in 1:length(bmt_NM$refseq_mrna)){
      rightrow=which(res$id==bmt_NM$refseq_mrna[i])[1];
      res$ENTREZ[rightrow] <- as.character(bmt_NM$entrezgene_id[i])
      res$HGNC[rightrow] <- as.character(bmt_NM$hgnc_symbol[i])
      res$LOCUS[rightrow] <- paste(bmt_NM$chromosome_name[i], ":",bmt_NM$start_position[i],"-",bmt_NM$end_position[i],sep="");
      res$DESC[rightrow] <- bmt_NM$description[i]
    } }
  }
  gnsNR=gns[grep("NR_",gns)]
  if(length(gnsNR)>0){
    message(paste("Querying for",length(aa),"attributes from biomaRt on",length(gnsNR)," NR numbers"));
    aa[1]="refseq_ncrna";
    bmt_NR=getBM(attributes = aa, filters = "refseq_ncrna",values=gns[grep("NR_",gns)], mart = ensembl)
    message(paste("We got back",length(bmt_NR$refseq_ncrna),"entries, which will be added to the local database."));
    if(nrow(bmt_NR)>0) {
    for(i in 1:length(bmt_NR$refseq_ncrna)){
      rightrow=which(res$id==bmt_NR$refseq_ncrna[i])[1];
      res$ENTREZ[rightrow] <- as.character(bmt_NR$entrezgene[i])
      res$HGNC[rightrow] <- as.character(bmt_NR$hgnc_symbol[i])
      res$LOCUS[rightrow] <- paste(bmt_NR$chromosome_name[i], ":",bmt_NR$start_position[i],"-",bmt_NR$end_position[i],sep="");
      res$DESC[rightrow] <- bmt_NR$description[i]
    } }
  }
  res
}



AnnotateTable <- function(G){
  C2T<-G$C2T
  infile=DESEQ_RFile(G)
  outfile=DESEQ_AFile(G)
  
  res=read.csv(infile,header=TRUE,as.is=TRUE)
  #res has got columns like
  #id	baseMean	baseMeanA	baseMeanB	foldChange	log2FoldChange	pval	padj
  # print("dumping rows with PADJ>0.90")
  # res=subset(res,padj <= 0.90)
  message("dumping rows with raw PVAL>0.990 OR baseMean < 1 OR on excludes list")
  res=subset(res,res$pval <= 0.990 & res$baseMean >= 1 & !(res$id %in% G$EXCLUDES))
  res$id = sub(" +$", "", res$id)
  message(paste("Rows remaining:",length(res$id)));
  if(length(unique(res$id)) < 1){
    warning("AnnotateTable : no rows significant enough to annotate.  Quitting table re-gen.");
    return();
  }
  
  message("Converting identifiers with local DB");
  res<-LocalDB_Lookup(res);
  #and the rest
  message("Converting identifiers with biomart");
  #fix strings?
  res$ENTREZ <- as.character(res$ENTREZ)
  res$HGNC <- as.character(res$HGNC)
  
   res = tryCatch({
     biomaRt_Lookup(res)
   }, warning = function(w) {
     warning(w)
     res
   }, error = function(e) {
     warning(e) #if bioMart is down, just keep moving.
     res
   })
  res<-subset(res,select=c("id","baseMean","baseMeanA","baseMeanB","foldChange","log2FoldChange","lfcSE","pval","padj","HGNC","ENTREZ"))
  columns = c("REFSEQ","baseMean",C2T[1],C2T[2],"(X/Y)","LOG2(X/Y)","lfcSE","PVAL","PADJ","HGNC","ENTREZ");
  if(!G$QUIET)Message_Generate(outfile)
  write.table(res,file=outfile,row.names=FALSE,sep=",",col.names=columns);
  ## and a post-filter for shortlisting
  write.table(subset(res, res$padj <= G$LIMIT$PADJ &
                          abs(res$log2FoldChange) >= G$LIMIT$LFC &
                          res$baseMean >= G$LIMIT$baseMean )
              ,file=DESEQ_AFFile(G),row.names=FALSE,sep=",",col.names=columns);
  
}

AnnotateTable.GO <- function(G){
  C2T<-G$C2T
  infile=DESEQ_AFile(G)
  outfile=GSEPD_GOFile(G)
  outfile2=GSEPD_GO2File(G)
  plotFile=paste(G$Output_Folder,"/GOSEQ.PWF.",C2T[1],".",C2T[2],".pdf",sep="")
  
  data=read.csv(infile,as.is=TRUE,header=TRUE)
  data=subset(data,!duplicated(data$ENTREZ))
  sdata=subset(data,data$ENTREZ != "NA")
  
  #default hard filtering:
  DEG=ifelse(  sdata$PADJ < G$LIMIT$PADJ  &
               sdata$baseMean >= G$LIMIT$baseMean &
               abs(sdata$LOG2.X.Y.) >= G$LIMIT$LFC, 1,0)
  
  #but if we're not set to HARD mode, and we got only a handful of genes:
  if( !(G$LIMIT$HARD) && sum(DEG)<(0.075*length(DEG))){
    ToLevel=round((0.075*nrow(sdata)))
    PThresh = sort(sdata$PVAL,decreasing=FALSE)[ToLevel]
    DEG=ifelse(sdata$PVAL < PThresh , 1,0)
    message(sprintf("Not many genes found differentially expressed, (re)moving the filters to raw p=%f so we can use %d genes with GOSEQ.",PThresh,sum(DEG)))
  }
  if(sum(DEG)<2){
    stop("I can't do pathways or GO terms when so few genes pass significance. Try changing your G$LIMIT parameters.")
  }
  
  DEG.up = ifelse(DEG & (sdata$LOG2.X.Y. > 0) , 1,0)
  DEG.dn = ifelse(DEG & (sdata$LOG2.X.Y. < 0) , 1,0)
  
  names(DEG) <- sdata$ENTREZ;
  names(DEG.up) <- sdata$ENTREZ;
  names(DEG.dn) <- sdata$ENTREZ;
  
  pdf(plotFile)
  pwf <- nullp(DEG, G$GOSEQ$genome, G$GOSEQ$system)
  #rownames(pwf)<-sdata$REFSEQ # we're calling goseq() with gene system 'refGene'
  #error. 
  GO.wall <- goseq(pwf, G$GOSEQ$genome,
   id=G$GOSEQ$system, use_genes_without_cat=G$GOSEQ$use_genes_without_cat)
  title(main="GOSeq pwf Allgenes")
  if(sum(DEG.up)>1) {
   pwf <- nullp(DEG.up, G$GOSEQ$genome, G$GOSEQ$system)
   #rownames(pwf)<-sdata$REFSEQ # we're calling goseq() with gene system 'refGene'
   GO.wall.Up <- goseq(pwf, G$GOSEQ$genome,
    id=G$GOSEQ$system, use_genes_without_cat=G$GOSEQ$use_genes_without_cat)
   title(main="GOSeq pwf Upreg")
  }else{
   GO.wall.Up <- GO.wall[FALSE,]
  }
  if(sum(DEG.dn)>1){
    pwf <- nullp(DEG.dn, G$GOSEQ$genome, G$GOSEQ$system)
    #rownames(pwf)<-sdata$REFSEQ # we're calling goseq() with gene system 'refGene'
    GO.wall.Dn <- goseq(pwf, genome=G$GOSEQ$genome,
    id=G$GOSEQ$system, use_genes_without_cat=G$GOSEQ$use_genes_without_cat)
    title(main="GOSeq pwf Downregulated genes.")
   }else{
    GO.wall.Dn <- GO.wall[FALSE,]
   }
  dev.off();
  
  ## here's a good place to send a warning if some whitelisted GO terms are not available
  igo <- setdiff(G$Force_GO_Include, "GO:0000000")
  Missing_GO_terms <- igo[!(igo%in%GO.wall$category)]
  if(length(Missing_GO_terms)>0) # allow one for the default GO:00000
    warning(sprintf("Some %d Whitelisted GO terms from Force_GO_Include.csv not found in database: %s",
     length(Missing_GO_terms),
     paste(Missing_GO_terms, collapse=" ")))
  
  FilterGOSEQ <- function(GOR, PLIM) {
    # do a Benjamini-Hochberg correction of p-values.
    GOR$over_represented_padj = p.adjust(GOR$over_represented_pvalue, method="BH")
    GOR$under_represented_padj = p.adjust(GOR$under_represented_pvalue, method="BH")
    #filter reported set to something with a pvalue:
    GOR$minp <- apply(cbind( 
      GOR$over_represented_padj , GOR$under_represented_padj),
                      1,min)
    GOR <- subset(GOR,GOR$minp <= PLIM | GOR$category %in% igo)
    gots = select(GO.db,keys=as.character(GOR$category),
                  columns=c("GOID","TERM","ONTOLOGY","DEFINITION"),
                  keytype="GOID")
    colnames(gots) <- c("GOID","Term","Ontology","Definition")
    merge(GOR, gots, by.x="category",by.y="GOID", all.x=TRUE)
  }
  GO.wall    <- FilterGOSEQ(GO.wall   , PLIM=G$LIMIT$GO_PVAL*1.5)
  GO.wall.Up <- FilterGOSEQ(GO.wall.Up, PLIM=G$LIMIT$GO_PVAL*1.5)
  GO.wall.Dn <- FilterGOSEQ(GO.wall.Dn, PLIM=G$LIMIT$GO_PVAL*1.5)
  
  Check_Missing_ForcedGO(GO.wall,igo)
  
  if(!G$QUIET)Message_Generate(outfile)
  write.csv(GO.wall   ,outfile);
  write.csv(GO.wall.Up,GSEPD_GOUpFile(G));
  write.csv(GO.wall.Dn,GSEPD_GODownFile(G));
  
  message("Written GO categories, now reverse mapping");
  test.cats = c("GO:CC", "GO:BP", "GO:MF")
  gene2cat = getgo(rownames(pwf), G$GOSEQ$genome,G$GOSEQ$system, fetch.cats = test.cats)
  names(gene2cat) = rownames(pwf)
  cat2gene = reverseMap(gene2cat)
  gene2cat = reverseMap(cat2gene)
  
  GO.wall <- rbind(GO.wall.Dn, GO.wall.Up, GO.wall) ; rm(GO.wall.Dn, GO.wall.Up)
  GO.wall <- subset(GO.wall, !(duplicated(GO.wall$category)))
  #this will pick whoever was first, probably wall.Dn, in the case of duplications.

  #okay awesome, let's make a table of those genes in GO.wall with over_represented_pvalue and genes in sdata 
  OM=matrix(nrow=0,ncol=2); colnames(OM)<-c("category","ENTREZ");
  AllCategories <- unique( unlist(GO.wall$category) )
  for(gocat in AllCategories){
    GID=cat2gene[[gocat]] #as vector of strings, but numeric entrez ids
    GID=GID[GID %in% sdata$ENTREZ]
    if(length(GID)>0)
      OM=rbind(OM, cbind(rep(gocat,length(GID)), GID))
  }
  
  merged=merge(OM,sdata,by="ENTREZ")
  #now making the GO2 file out of just the single GOseq. 
  merged <- merge(merged, GO.wall, by="category", all.x=TRUE)
  #count genes per category:
  Repetition <- table(merged$category)
  merged$GO_GCount <- Repetition[merged$category];
  if(!G$QUIET)Message_Generate(outfile2)
  keepCol=c("REFSEQ","HGNC",G$C2T[1],G$C2T[2], "LOG2.X.Y.",  "PVAL",
            "PADJ",	"category",	"over_represented_pvalue",	"under_represented_pvalue",
            "over_represented_padj",	"under_represented_padj",
            "numDEInCat",	"numInCat",	"Term",	"Ontology")
 #keepCol[ !(keepCol %in% colnames(merged))]
  merged=subset(merged,merged$numInCat<800 |  merged$category %in% igo, select=keepCol)
  merged[,c(G$C2T[1],G$C2T[2], "LOG2.X.Y.")] = signif( merged[,c(G$C2T[1],G$C2T[2], "LOG2.X.Y.")], digits=3 )
  write.csv(merged,outfile2)
}

Check_Missing_ForcedGO <- function(go2, include_list){
  m <- include_list[  !(include_list%in% go2$category) ]
  if(length(m)>0)
    warning(sprintf("Trying to include some GO terms but they arent found by GO.db : %s",paste(m,collapse=" ")))
}


GSEPD_Export_DESeq<-function(G){
  colData<-G$sampleMeta
  rownames(colData)<-G$sampleMeta$Sample
  colData$Condition <- as.factor(colData$Condition)
  colData$Condition <- relevel(colData$Condition, G$Conditions[1])
  dds <- DESeqDataSetFromMatrix(countData = G$finalCounts,
                                colData = colData,
                                design = ~ Condition)
  return(dds)
}


GSEPD_CheckCounts <- function(GSEPD) {
  if(is.null(GSEPD$normCounts)) {
    #rld <- rlog()
    vsd <- DESeq2::varianceStabilizingTransformation(GSEPD_Export_DESeq(GSEPD), 
                                                     blind=GSEPD$vstBlind)
    #rlogMat <- assay(rld)
    vstMat <- assay(vsd)
    GSEPD$normCounts <- vstMat
  }
GSEPD
}

