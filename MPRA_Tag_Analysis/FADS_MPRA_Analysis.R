rm(list=ls())
stringsAsFactors = FALSE
library(DESeq2)
library(ggplot2)
library(tidyr)
library(splitstackshape)
library(rtracklayer)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)

args = commandArgs(trailingOnly=TRUE)

countfile<-args[1]
wd<-dirname(dirname(countfile))
file_prefix=args[2]
file_date=format(Sys.Date(), "%Y%m%d")
setwd(paste0(wd,"/count_analysis/"))
dir.create("plots")
dir.create("results")

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

message(countfile)
message(file_prefix)
message(wd)

inputFullAttributesFileName=paste0(wd,"/MPRA_Tag_Analysis/OL13_FADS_Encode.attributes")

countsFileName=paste0(countfile)
inputConditionsFileName=paste0(wd,"/MPRA_Tag_Analysis/OL13_FADS_K562.cond")

tag_counts<-read.table(countsFileName,header=TRUE)
fullattributeData=read.table(inputFullAttributesFileName, header=TRUE, row.names= NULL, check.names=FALSE,skipNul=TRUE)

if(!("haplotype" %in% colnames(fullattributeData)))
{
  fullattributeData$haplotype<-"NA"
}

counts<-read.table(file=countsFileName,header=TRUE)

counts_oligo<-aggregate(. ~Oligo, data=tag_counts[,c(-1,-3,-4,-5,-6)], FUN=sum)

countData<-counts_oligo[,-1]
rownames(countData)<-counts_oligo[,1]
countDataRowNames=rownames(countData)

colData=read.table(inputConditionsFileName, header=FALSE, row.names=1)


colnames(colData)[1]="condition"
colData[,1]=factor(colData[,1])

colData$condition=relevel(colData$condition, "DNA")
rowNames=rownames(countData)
colnames(countData)<-row.names(colData)

#Just doing this so the DF can be reordered....
colData$DNA<-0
colData[colData$condition=="DNA",]$DNA<-1 


################# initial processing of files #############



colData$condition=relevel(colData$condition, "DNA")
rowNames=rownames(countData)
colnames(countData)<-row.names(colData)

dds <- DESeqDataSetFromMatrix(countData=countData, colData=colData, design=~condition)
dds$condition<-relevel(dds$condition, "DNA")
dds_results <- DESeq(dds,fitType='local')

orig_dds_results<-dds_results

#############
###### Normalization
#############

dds_results<-orig_dds_results

###Code for neg (all or ctrl) normalization
#correct on neg control
nonexpA<-fullattributeData[fullattributeData$project=="negCtrl",]$ID
#correct from all
#nonexpA<-fullattributeData$ID

#for (celltype in c("K562") {
for (celltype in levels(colData$condition)) {
  if(celltype=="DNA") next
  message(celltype)
  output_tmpA = results(dds_results, contrast=c("condition",celltype,"DNA"))

  ##Code for summit shift normalization
  summit<-which.max(density(output_tmpA$log2FoldChange)$y)
  log_offset<-2^(density(output_tmpA$log2FoldChange)$x[summit])
  sizeFactors(dds_results)[which(colData$condition==celltype)]<-sizeFactors(dds_results)[which(colData$condition==celltype)]*(log_offset)
}

dds_results<-estimateDispersions(dds_results,fitType='local')
dds_results<-nbinomWaldTest(dds_results)

for (celltype in levels(colData$condition)) {
  if(celltype=="DNA") next
  
  output_tmpA = results(orig_dds_results, contrast=c("condition",celltype,"DNA"))
  outputA = results(dds_results, contrast=c("condition",celltype,"DNA"))
  
  pdf(paste0("plots/Normalized_FC_Density_",celltype,".pdf"),width=10,height=10)
  plot(density(output_tmpA[nonexpA,]$log2FoldChange,na.rm=TRUE),xlim=c(-3,3),col="grey",main=paste0("Normalization - ",celltype))
  lines(density(output_tmpA$log2FoldChange,na.rm=TRUE),xlim=c(-3,3),col="black")
  lines(density(outputA$log2FoldChange,na.rm=TRUE),xlim=c(-3,3),col="red")
  lines(density(outputA[nonexpA,]$log2FoldChange,na.rm=TRUE),xlim=c(-3,3),col="salmon")
  text(1.5,1.4,adj=c(0,0),labels="All - baseline",col="black")
  text(1.5,1.35,adj=c(0,0),labels="All - corrected",col="red")
  text(1.5,1.3,adj=c(0,0),labels="NegCtrl - baseline",col="grey")
  text(1.5,1.25,adj=c(0,0),labels="NegCtrl - corrected",col="salmon")
  abline(v=0)
  dev.off()
}



#############
#############
#############

counts <- counts(dds_results, normalized=TRUE)

#Expand IDs that denote duplicate oligos in count/DESeq results
expand_dups<-function(table)
{
  orig_table<-table
  table<-cbind(rownames(table), table)
  colnames(table)[1]<-"Row.names"
  dups<-table[grep("\\(.*\\)$",table$Row.names),]
  dups$Row.names<-gsub("^\\((.*)\\)$","\\1",dups$Row.names)
  final_table<-table[-(grep("\\(.*\\)$",table$Row.names)),] 
  final_table<-final_table[!(is.na(final_table$Row.names)),]
  if(nrow(dups)>0) {
    for(i in 1:nrow(dups))
    {
      dup_id<-unlist(strsplit(dups$Row.names[i],";"))
      dup_exp<-dups[rep(i, length(dup_id)), ]
      dup_exp$Row.names<-dup_id
      final_table<-rbind(dup_exp,final_table)
    }
    rownames(final_table)<-final_table[,1]
    final_table<-final_table[,-1]
  } else {
    final_table<-orig_table
  }
  return(final_table)
}

### Function to perform TTest on individual cell types
CellSpecific_Ttest<-function(tmp_attributeData, counts, func_output, Ctrl.Mean, Exp.Mean, ctrl_cols, exp_cols) {
  
  snp_data<-subset(tmp_attributeData,allele=="ref" | allele=="alt")
  snp_data$comb<-paste(snp_data$SNP,"_",snp_data$window,"_",snp_data$strand,"_",snp_data$haplotype,sep="")
  tmp_ct<-as.data.frame(table(snp_data$comb))
  
  snp_data_pairs<-snp_data[snp_data$comb %in% tmp_ct[tmp_ct$Freq==2,]$Var1,]
  
  snp_data_rejected<-snp_data[snp_data$comb %in% tmp_ct[tmp_ct$Freq!=2,]$Var1,]
  
  snp_data_ctdata_pairs<-merge(snp_data_pairs,counts,by.x="ID",by.y="row.names",all.x=TRUE,no.dups=FALSE)
  snp_data_ctdata_pairs<-snp_data_ctdata_pairs[order(snp_data_ctdata_pairs$SNP,snp_data_ctdata_pairs$window,snp_data_ctdata_pairs$strand,snp_data_ctdata_pairs$haplotype,snp_data_ctdata_pairs$allele),]
  snp_data_expdata_pairs<-merge(snp_data_pairs,func_output,by.x="ID",by.y="row.names",all.x=TRUE,no.dups=FALSE)
  snp_data_expdata_pairs<-snp_data_expdata_pairs[order(snp_data_expdata_pairs$SNP,snp_data_expdata_pairs$window,snp_data_expdata_pairs$strand,snp_data_expdata_pairs$haplotype,snp_data_expdata_pairs$allele),]
  
  
  ###This is confusing.... but switched the odd/even order here since I am sorting on alt/ref. If ref/alt becomes A/B switch the order here instead of updating the code below. Regardless if this is wrong it will display the A allele in output file correctly.
  evens <- seq(1, nrow(snp_data_pairs), by=2)
  odds <- seq(2, nrow(snp_data_pairs), by=2)
  
  out <- cbind(
    snp_data_expdata_pairs[odds,c(1,2,3,4,5,6,7,9)],
    within(data.frame(
      A.Ctrl.Mean=snp_data_expdata_pairs[odds, "Ctrl.Mean"],
      A.Exp.Mean=snp_data_expdata_pairs[odds, "Exp.Mean"],
      A.log2FC=snp_data_expdata_pairs[odds, "log2FoldChange"],
      A.log2FC_SE=-log10(snp_data_expdata_pairs[odds, "lfcSE"]),
      A.logP=-log10(snp_data_expdata_pairs[odds, "pvalue"]),
      A.logPadj_BH=-log10(snp_data_expdata_pairs[odds, "padj"]),    #BF Correction
      A.logPadj_BF=-log10(snp_data_expdata_pairs[odds, "pvalue"]*(nrow(snp_data_expdata_pairs)/2)),    #BF Correction
      B.Ctrl.Mean=snp_data_expdata_pairs[evens, "Ctrl.Mean"],
      B.Exp.Mean=snp_data_expdata_pairs[evens, "Exp.Mean"],
      B.log2FC=snp_data_expdata_pairs[evens, "log2FoldChange"],
      B.log2FC_SE=-log10(snp_data_expdata_pairs[evens, "lfcSE"]),
      B.logP=-log10(snp_data_expdata_pairs[evens, "pvalue"]),
      B.logPadj_BH=-log10(snp_data_expdata_pairs[evens, "padj"]),    #BF Correction
      B.logPadj_BF=-log10(snp_data_expdata_pairs[evens, "pvalue"]*(nrow(snp_data_expdata_pairs)/2))),{   #BF Correction
        A.logP[is.na(A.logP)] <- 0
        A.logP[A.logP == Inf] <- max(A.logP[is.finite(A.logP)])
        A.logPadj_BH[A.logPadj_BH < 0]<-0
        A.logPadj_BH[A.logPadj_BH == Inf] <- max(A.logPadj_BH[is.finite(A.logPadj_BH)])
        A.logPadj_BF[A.logPadj_BF < 0]<-0
        A.logPadj_BF[A.logPadj_BF == Inf] <- max(A.logPadj_BF[is.finite(A.logPadj_BF)])
        B.logP[is.na(B.logP)] <- 0
        B.logP[B.logP == Inf] <- max(B.logP[is.finite(B.logP)])
        B.logPadj_BH[B.logPadj_BH < 0]<-0
        B.logPadj_BH[B.logPadj_BH == Inf] <- max(B.logPadj_BH[is.finite(B.logPadj_BH)])
        B.logPadj_BF[B.logPadj_BF < 0]<-0
        B.logPadj_BF[B.logPadj_BF == Inf] <- max(B.logPadj_BF[is.finite(B.logPadj_BF)])
      }))
  
  # Don't try to do the t test for ones with all zeros.
  ignore_idx <- which(rowMeans(snp_data_ctdata_pairs[odds,ctrl_cols]) < 10 | rowMeans(snp_data_ctdata_pairs[odds, exp_cols]) < 10 |
                        rowMeans(snp_data_ctdata_pairs[evens, ctrl_cols]) < 10 | rowMeans(snp_data_ctdata_pairs[evens, exp_cols]) < 10 | 
                        is.na(rowMeans(snp_data_ctdata_pairs[odds,ctrl_cols]))  | is.na(rowMeans(snp_data_ctdata_pairs[odds, exp_cols])) |
                        is.na(rowMeans(snp_data_ctdata_pairs[evens, ctrl_cols])) | is.na(rowMeans(snp_data_ctdata_pairs[evens, exp_cols])) )
  
  # For the numerator, set zero values to 1 so that the log-ratio is defined.
  counts1 <- snp_data_ctdata_pairs
  counts1[counts1 == 0] <- 1
  
  # t test
  ratios.A <- log((counts1[odds, exp_cols]) / rowMeans(snp_data_ctdata_pairs[odds, ctrl_cols]))
  ratios.B <- log((counts1[evens, exp_cols]) / rowMeans(snp_data_ctdata_pairs[evens, ctrl_cols]))
  t.pvalue <- sapply(1:nrow(ratios.A), function(i) 
    if(i %in% ignore_idx) { NA } else{ t.test(as.numeric(ratios.A[i,]), as.numeric(ratios.B[i,]), var.equal=FALSE, paired=TRUE)$p.value})
  t.stat <- sapply(1:nrow(ratios.A), function(i) 
    if(i %in% ignore_idx) { NA } else{ t.test(as.numeric(ratios.A[i,]), as.numeric(ratios.B[i,]), var.equal=FALSE, paired=TRUE)$statistic})
  
  out$LogSkew <- out$B.log2FC - out$A.log2FC
  out$Skew.logP <- ifelse(is.na(t.pvalue), 0, -log10(t.pvalue))
  
  OE_threshold <- -log10(.01)
  is_OE <- out$A.logPadj_BF >= OE_threshold | out$B.logPadj_BF >= OE_threshold
  out$Skew.logFDR <- rep(NA, dim(out)[1])
  q_idx <- intersect(which(is_OE), which(!is.na(t.pvalue)))
  out$Skew.logFDR[q_idx] <- -log10(p.adjust(t.pvalue[q_idx],method="BH"))
  
  
  return(out)
}

###############
#### Loop through cell types


round_df <- function(x, digits) {
    # round all numeric variables
    # x: data frame 
    # digits: number of digits to round
    numeric_columns <- sapply(x, mode) == 'numeric'
    x[numeric_columns] <-  round(x[numeric_columns], digits)
    x
}


full_output<-list()
full_output_var<-list()

for (celltype in levels(colData$condition)) {
  if(celltype=="DNA") next
  message(celltype)
  
  outputA = results(dds_results, contrast=c("condition",celltype,"DNA"))
  
  ctrl_cols<-row.names(colData[colData$condition=="DNA",])
  exp_cols<-row.names(colData[colData$condition==celltype,])
  
  Ctrl.Mean=rowMeans(counts[, colnames(counts) %in% ctrl_cols])
  Exp.Mean=rowMeans(counts[, colnames(counts) %in% exp_cols])
  output_2<-cbind(Ctrl.Mean,Exp.Mean,outputA[,-1])
  output_undup<-expand_dups(output_2)
  
  full_outputA<-merge(fullattributeData, as.matrix(output_undup),by.x="ID",by.y="row.names",all.x=TRUE,no.dups=FALSE)
  full_output[[celltype]]<-full_outputA
 
  full_outputA$log2FoldChange<-round_df(full_outputA$log2FoldChange,7)
  full_outputA$Ctrl.Mean<-round_df(full_outputA$Ctrl.Mean,5)
  full_outputA$Exp.Mean<-round_df(full_outputA$Exp.Mean,5)
  full_outputA$lfcSE<-round_df(full_outputA$lfcSE,7)
 
  write.table(full_outputA,paste0("results/",file_prefix,"_",celltype,"_",file_date,".out"), row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
  
  outA<-CellSpecific_Ttest(fullattributeData, counts, output_undup, Ctrl.Mean, Exp.Mean, ctrl_cols, exp_cols) 
  full_output_var[[celltype]]<-outA
  
  outA$A.log2FC<-round_df(outA$A.log2FC,7)
  outA$A.Ctrl.Mean<-round_df(outA$A.Ctrl.Mean,5)
  outA$A.Exp.Mean<-round_df(outA$A.Exp.Mean,5)
  outA$A.log2FC_SE<-round_df(outA$A.log2FC_SE,7)
  outA$A.logP<-round_df(outA$A.logP,7)
  outA$A.logPadj_BH<-round_df(outA$A.logPadj_BH,7)
  outA$A.logPadj_BF<-round_df(outA$A.logPadj_BF,7)


  outA$B.log2FC<-round_df(outA$B.log2FC,7)
  outA$B.Ctrl.Mean<-round_df(outA$B.Ctrl.Mean,5)
  outA$B.Exp.Mean<-round_df(outA$B.Exp.Mean,5)
  outA$B.log2FC_SE<-round_df(outA$B.log2FC_SE,7)
  outA$B.logP<-round_df(outA$B.logP,7)
  outA$B.logPadj_BH<-round_df(outA$B.logPadj_BH,7)
  outA$B.logPadj_BF<-round_df(outA$B.logPadj_BF,7)
  
  outA$LogSkew<-round_df(outA$LogSkew,7)
  outA$Skew.logP<-round_df(outA$Skew.logP,7)
  outA$Skew.logFDR<-round_df(outA$Skew.logFDR,7)
  
  write.table(outA,paste0("results/",file_prefix,"_",celltype,"_emVAR_",file_date,".out"), row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
  
  full_bed_outputA<-merge(fullattributeData, as.matrix(output_undup),by.x="ID",by.y="row.names",all.x=TRUE,no.dups=FALSE)
  printbed<-full_bed_outputA[,c("chr","start","stop","ID","strand","log2FoldChange","Ctrl.Mean","Exp.Mean","pvalue","padj","lfcSE","cigar","md-tag","project")]      
  printbed$score<-"."
  printbed<-printbed[,c("chr","start","stop","ID","score","strand","log2FoldChange","Ctrl.Mean","Exp.Mean","pvalue","padj","lfcSE","cigar","md-tag","project")]      
  colnames(printbed)<-c("chr","start","stop","id","score","strand","log2fc","input-count","output-count","log10pval","log10fdr","lfc-se","cigar","md-tag","project")
  printbed$strand[printbed$strand=="fwd"]="+"
  printbed$strand[printbed$strand=="rev"]="-"
  printbed$log10pval=-log10(printbed$log10pval)
  printbed$log10fdr=-log10(printbed$log10fdr)
  
  printbed$log2fc<-round_df(printbed$log2fc,7)
  printbed[,"input-count"]<-round_df(printbed[,"input-count"],5)
  printbed[,"output-count"]<-round_df(printbed[,"output-count"],5)
  printbed[,"lfc-se"]<-round_df(printbed[,"lfc-se"],7)
  printbed$log10pval<-round_df(printbed$log10pval,7)
  printbed$log10fdr<-round_df(printbed$log10fdr,7)



  write.table(printbed,paste0("results/",file_prefix,"_",celltype,"_",file_date,".bed"),row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)

  full_bed_outputA<-merge(fullattributeData, as.matrix(output_undup),by.x="ID",by.y="row.names",all.x=TRUE,no.dups=FALSE)
  printbed<-full_bed_outputA[,c("chr","start","stop","ID","strand","log2FoldChange","Ctrl.Mean","Exp.Mean","pvalue","padj","lfcSE","cigar","md-tag","project")]      
  printbed$score<-"."
  printbed<-printbed[,c("chr","start","stop","ID","score","strand","log2FoldChange","Ctrl.Mean","Exp.Mean","pvalue","padj","lfcSE","cigar","md-tag","project")]      
  colnames(printbed)<-c("chr","start","stop","id","score","strand","log2fc","input-count","output-count","log10pval","log10fdr","lfc-se","cigar","md-tag","project")
  printbed$strand[printbed$strand=="fwd"]="+"
  printbed$strand[printbed$strand=="rev"]="-"
  printbed$log10pval=-log10(printbed$log10pval)
  printbed$log10fdr=-log10(printbed$log10fdr)

  printbed$log2fc<-round_df(printbed$log2fc,7)
  printbed[,"input-count"]<-round_df(printbed[,"input-count"],5)
  printbed[,"output-count"]<-round_df(printbed[,"output-count"],5)
  printbed[,"lfc-se"]<-round_df(printbed[,"lfc-se"],7)
  printbed$log10pval<-round_df(printbed$log10pval,7)
  printbed$log10fdr<-round_df(printbed$log10fdr,7)


  write.table(printbed,paste0("results/",file_prefix,"_",celltype,"_",file_date,".bed"),row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)

}


#### End cell type analysis
###############
###########
##########################
##########################
#########
###### General QC Plots 
#########
counts <- counts(dds_results, normalized=FALSE)
write.table(counts,paste0("results/",file_prefix,"_Counts_",file_date,".out"), row.names=TRUE,col.names=NA,sep="\t",quote=FALSE)
counts <- counts(dds_results, normalized=TRUE)
write.table(counts,paste0("results/",file_prefix,"_NormCounts_",file_date,".out"), row.names=TRUE,col.names=NA,sep="\t",quote=FALSE)


panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
panel.lm <- function (x, y,  pch = par("pch"), col.lm = "red",  ...) 
{   
  ymin <- min(y)
  ymax <- max(y)
  xmin <- min(x)
  xmax <- max(x)
  ylim <- c(min(ymin,xmin),max(ymax,xmax))
  xlim <- ylim
  points(x, y, pch = pch,ylim = ylim, xlim= xlim,col=rgb(144,144,144,75,maxColorValue=255),...)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)) 
    abline(lm(y[ok]~ x[ok]), 
           col = col.lm, ...)
}

panel.nlm <- function (x, y,  pch = par("pch"), col.lm = "red",  ...) 
{   
  ymin <- min(y)
  ymax <- max(y)
  xmin <- min(x)
  xmax <- max(x)
  ylim <- c(min(ymin,xmin),max(ymax,xmax))
  xlim <- ylim
  points(x, y, pch = pch,ylim = ylim, xlim= xlim,col=rgb(144,144,144,75,maxColorValue=255),...)
}

png(file="plots/Cor_mat_log.png",width=3000,height=3000) 
pairs(counts,upper.panel=panel.cor,lower.panel=panel.lm,log="xy",pch=16)
dev.off()
png(file="plots/Cor_mat.png",width=3000,height=3000) 
pairs(counts,upper.panel=panel.cor,lower.panel=panel.lm,pch=16)
dev.off()
png(file="plots/Cor_mat_2.png",width=3000,height=3000) 
pairs(counts,upper.panel=panel.cor,lower.panel=panel.lm,xlim=c(0,2000),ylim=c(0,2000),pch=16)
dev.off()


ggplot_mpra_scatter<-function(counts, sampleX, sampleY,xmax,ymax) {
  count_df<-as.data.frame(counts)
  ggplot_output<-ggplot(count_df, aes_string(sampleX,sampleY)) +
    theme_minimal() +
    geom_point(alpha = .3,size=1) +
    labs(x=sampleX,y=sampleY) +
    coord_fixed(ratio = 1,xlim=c(0,xmax), ylim=c(0,ymax)) +
    geom_abline(intercept = 0, slope = 1,linetype = 2, size=.75, color=rgb(255,140,0,150,maxColorValue=255))
  return(ggplot_output)
}

xmax<-quantile(counts,.99)
ymax<-quantile(counts,.99)

sampleX<-"plasmid_r1"
sampleY<-"plasmid_r2"
ggplot_output<-ggplot_mpra_scatter(counts, sampleX, sampleY,xmax,ymax)
ggsave("plots/plasmid_cor.png",ggplot_output,units="in",width=4,height=4,device="png")

sampleX<-"plasmid_r1"
sampleY<-"K562_r1"
ggplot_output<-ggplot_mpra_scatter(counts, sampleX, sampleY,xmax,ymax)
ggsave("plots/plasmid-K562_cor.png",ggplot_output,units="in",width=4,height=4,device="png")






plot_logFC<-function(full_output, sample) {
  exp.values<-full_output[full_output$Ctrl.Mean > 10 & !is.na(full_output$Ctrl.Mean),]
  exp.values$Exp.Mean[is.na(exp.values$Exp.Mean)]<-1
  exp.values$log2FoldChange[is.na(exp.values$log2FoldChange)]<-0
  exp.values$padj[is.na(exp.values$padj)]<-1
  exp.values$sig<-"Not Significant"
  exp.values$sig[exp.values$padj <= 0.00001]<-"Active"
  levels(exp.values$sig)<-c("Not Significant", "Active")
  exp.values$sig<-factor(exp.values$sig,levels=c("Not Significant", "Active"))
  
  tmp_plotA<-ggplot(exp.values,aes(x=Ctrl.Mean,y=log2FoldChange,color=sig)) +
    theme_bw() + theme(panel.grid.major = element_line(size = .25,colour = rgb(0,0,0,75,maxColorValue=255)), panel.grid.minor = element_blank()) +
    scale_colour_manual(values=c("Not Significant"=rgb(0,0,0,200,maxColorValue=255),"Active"=rgb(55,126,184,255,maxColorValue=255))) +
    geom_point(alpha = .3,size=1) +
    scale_x_log10() +
    #coord_cartesian(xlim = c(10, 1000),ylim = c(-1.5,7.5)) +
    labs(x="Normalized Tag Count - Plasmids",y=paste0(sample," Expression Fold Change log2(RNA/Plasmid)")) +
    theme(legend.position = c(.15, .90), 
          legend.key = element_blank(),
          legend.background = element_rect(color=rgb(0,0,0,150,maxColorValue=255), fill = "white", size = .5, linetype = "solid")) +
    guides(colour = guide_legend(override.aes = list(size=3,alpha=.7), title=NULL)) + 
    geom_abline(intercept = 0, slope = 0,linetype = 1, size=.75, color=rgb(255,140,0,150,maxColorValue=255))
  
  cbPalette <- c("#56B4E9", "#999999", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  tmp_plotB<-1
  tmp_plotB<-ggplot(exp.values,aes(x=Ctrl.Mean,y=log2FoldChange,color=project)) +
    theme_bw() + theme(panel.grid.major = element_line(size = .25,colour = rgb(0,0,0,75,maxColorValue=255)), panel.grid.minor = element_blank()) +
    scale_colour_manual(values=cbPalette) +
    geom_point(alpha = .2,size=1) +
    geom_point(data = subset(exp.values, project == 'negCtrl'), aes(x=Ctrl.Mean,y=log2FoldChange),alpha = .8,size=2) +
    geom_point(data = subset(exp.values, project == 'expCtrl'), aes(x=Ctrl.Mean,y=log2FoldChange),alpha = .8,size=2) +
    geom_point(data = subset(exp.values, project == 'emVarCtrl'), aes(x=Ctrl.Mean,y=log2FoldChange),alpha = .8,size=2) +
    scale_x_log10() +
    #coord_cartesian(xlim = c(10, 1000),ylim = c(-1.5,7.5)) +
    labs(x="Normalized Tag Count - Plasmids",y=paste0(sample," Expression Fold Change log2(RNA/Plasmid)")) +
    theme(legend.position = c(.15, .80), 
          legend.key = element_blank(),
          legend.background = element_rect(color=rgb(0,0,0,150,maxColorValue=255), fill = "white", size = .5, linetype = "solid")) +
    guides(colour = guide_legend(override.aes = list(size=3,alpha=.7), title=NULL)) + 
    geom_abline(intercept = 0, slope = 0,linetype = 1, size=.75, color=rgb(0,0,0,0,maxColorValue=255))
  
  return(list(tmp_plotA,tmp_plotB))
}

for (celltype in levels(colData$condition)) {
  if(celltype=="DNA" | celltype=="EXCLUDE_CELL" ) next
  message(celltype)
  output_tmp<-full_output[[celltype]]
  plot_list<-plot_logFC(output_tmp, celltype)
  ggsave(paste0("plots/logFC_",celltype,".pdf"),plot_list[[1]],units="in",width=8,height=6,device="pdf")
  ggsave(paste0("plots/logFC_",celltype,"_controls.pdf"),plot_list[[2]],units="in",width=8,height=6,device="pdf")
}

##############
#BED files 
##############
##############

############################
tmp<-as.data.frame(full_output[["K562"]])
tmp<-tmp[which(tmp$project=="Tiles"),] 
tmp$strand <- as.character(tmp$strand)
tmp[which(tmp$project=="Tiles" & tmp$strand=="fwd" ),]$strand<-"+"
tmp[which(tmp$project=="Tiles" & tmp$strand=="rev" ),]$strand<-"-"
tmp<-tmp[!is.na(tmp$log2FoldChange),]

gr <- GRanges(seqnames=Rle(tmp$chr),
              ranges = IRanges(tmp$start, end=tmp$stop),
              strand = Rle(strand(tmp$strand)),
              score = tmp$log2FoldChange)

  newStyle <- mapSeqlevels(seqlevels(gr),"UCSC")
  ucsc_gr <- renameSeqlevels(gr, newStyle)
  genome(ucsc_gr) <- "hg19"
  
  seqlengths(ucsc_gr) = seqlengths(Hsapiens)[11]
  
  ucsc_gr_pos<-ucsc_gr[ucsc_gr@strand=="+",]
  ucsc_gr_neg<-ucsc_gr[ucsc_gr@strand=="-",]
 
	paste0("results/",file_prefix,"_",celltype,"_",file_date,"_pos.bedgraph")
  
  export(ucsc_gr_pos,paste0("results/",file_prefix,"_",celltype,"_tiling_",file_date,"_pos.bedgraph"),format="bedGRaph")
  export(ucsc_gr_neg,paste0("results/",file_prefix,"_",celltype,"_tiling_",file_date,"_pos.bedgraph"),format="bedGRaph")
  
  
  ##########
  ### General Analysis
  
  tmp<-full_output[["K562"]]
  tmp$mergeName<-paste0(tmp$chr,":",tmp$start,"-",tmp$stop)
  tmp_fwd<-tmp[which(tmp$strand=="fwd" & tmp$project=="Tiles"),]
  tmp_rev<-tmp[which(tmp$strand=="rev" & tmp$project=="Tiles"),]
  
  merged_output<-merge(tmp_fwd,tmp_rev,by="mergeName")  
  merged_output<-merged_output[which(merged_output$Ctrl.Mean.x > 50 & merged_output$Ctrl.Mean.y > 50 ),]
  merged_output$sig<-"NS"
  merged_output[which(merged_output$padj.x < 0.000001),]$sig<-"Active fwd"
  merged_output[which(merged_output$padj.y < 0.000001),]$sig<-"Active rev"
  merged_output[which(merged_output$padj.x < 0.000001 & merged_output$padj.y < 0.000001),]$sig<-"Active both"
  
    
    #mid <- mean(merged_output_sig$Ctrl.Mean.x,na.rm = T)
  #mid <- mean(merged_output_sig$Ctrl.Mean.x,na.rm = T)
  
    tmp_plotB<-ggplot(merged_output,aes(x=log2FoldChange.x,y=log2FoldChange.y,color=sig)) +
    theme_bw() + theme(panel.grid.major = element_line(size = .25,colour = rgb(0,0,0,75,maxColorValue=255)), panel.grid.minor = element_blank()) +
    #scale_colour_gradientn(limits=c(0,200),colours = terrain.colors(10)) +
    scale_colour_brewer(palette = "Set1") +
    geom_point(alpha = .4,size=1) +
    coord_cartesian(xlim = c(-5, 10),ylim = c(-5,10)) +
    labs(x="log2(RNA/DNA) - FWD",y="log2(RNA/DNA) - REV",color="Rev Plasmid Count") +
    geom_abline(intercept = 0, slope = 1,linetype = 1, size=.75, color=rgb(0,0,0,0,maxColorValue=255))
    ggsave(paste0("plots/FADS_tile_strand.png"),tmp_plotB,units="in",width=12,height=12,device="png")
    
  
  