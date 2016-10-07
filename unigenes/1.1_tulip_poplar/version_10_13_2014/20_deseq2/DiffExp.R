#source("http://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
library("DESeq2")
library("ggplot2")
library("gplots")
library("grid")
library("VennDiagram")

#########
# Ozone #
#########

setwd("/lustre/projects/staton/projects/hardwoods/unigenes/1.1_tulip_poplar/version_10_13_2014/20_deseq2/htseq/")
sampleFiles <- grep(".txt.filtered",list.files(),value=TRUE)
sampleFiles <- c("14DAy-Co_S13.counts.txt.filtered",
                 "14DAy-80ppb_S14.counts.txt.filtered",
                 "14DAy-125ppb_S15.counts.txt.filtered",
                 "14DAy-225ppb_S16.counts.txt.filtered",
                 "28Day-CO_S9.counts.txt.filtered",
                 "28Day-80_S10.counts.txt.filtered",
                 "28Day-125_S11.counts.txt.filtered",
                 "28Day-225_S12.counts.txt.filtered",
                 "7HR-Co_S5.counts.txt.filtered",
                 "7HR-80ppb_S6.counts.txt.filtered",
                 "7HR-125ppb_S7.counts.txt.filtered",
                 "7HR-225PPB_S9.counts.txt.filtered")

time=c( rep("336",4),
        rep("672",4),
        rep("7",4))

treatment=c(rep( c("0","80","125","225"), 3))

sampleTable<-data.frame(
  sampleName=sampleFiles,
  filename=sampleFiles,
  time=time,
  treatment=treatment)

directory <- ("/lustre/projects/staton/projects/hardwoods/unigenes/1.1_tulip_poplar/version_10_13_2014/20_deseq2/htseq/")

#------------------------------------------
# transcripts responding to ozone treatment... ignoring interaction term
#------------------------------------------
dds_trmt <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,
                                     directory=directory,
                                     design= ~time + treatment)
                                     #design= ~time + treatment + time:treatment)

#dds$treatment<-relevel(dds$treatment, "control")
#dds$time<-relevel(dds$time, "7")
dds_trmt<-DESeq(dds_trmt, test = c("LRT"), reduced = ~time)
#dds_trmt <- DESeq(dds_trmt)

res_trmt<-results(dds_trmt)
resOrdered_trmt<-res_trmt[order(res_trmt$padj),]
resSig_trmt <- subset(resOrdered_trmt, padj < 0.05)
resSigUp <- subset(resSig_trmt, log2FoldChange > 0)
resSigDown <- subset(resSig_trmt, log2FoldChange < 0)
resSigUp
resSigDown

# write outputs
setwd("/lustre/projects/staton/projects/hardwoods/unigenes/1.1_tulip_poplar/version_10_13_2014/20_deseq2/")
write.csv(as.data.frame(resSigUp),
          file="Ozone.DE.up.csv")
write.csv(as.data.frame(resSigDown),
          file="Ozone.DE.dw.csv")

#--------
# quad treatment
#---------
dds_quad_trmt <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,
                                          directory=directory,
                                          design= ~treatment + time + I(treatment^2))
dds_quad_trmt<-DESeq(dds_quad_trmt, test = c("LRT"), reduced = ~ treatment + time)
res<-results(dds_quad_trmt)
resOrdered<-res[order(res$padj),]
resSigQuadTrmt <- subset(resOrdered, padj < 0.05)
resSigUp <- subset(resSigQuadTrmt, log2FoldChange > 0)
resSigDown <- subset(resSigQuadTrmt, log2FoldChange < 0)

#------------------------------------------
# heatmap experiment 1 - betas
#------------------------------------------
betas <- coef(dds_trmt)
colnames(betas)

library("grid")
library("pheatmap")
draw_colnames_45 <- function (coln, ...) {
  m = length(coln)
  x = (1:m)/m - 1/2/m
  grid.text(coln, x = x, y = unit(0.96, "npc"), vjust = .5, 
            hjust = 0, rot = 270, gp = gpar(...)) ## Was 'hjust=0' and 'rot=270'
}
## 'Overwrite' default draw_colnames with your own version 
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))
#topGenes <- head(order(resTC$pvalue),20)
#mat <- betas[resSig_trmt, -c(1,2)]


topGenes <- head(order(res_trmt$padj),350)
mat <- betas[topGenes, -c(1)]
thr <- 3 # threshold for plotting
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr
matOrder <- mat[,c(1,2,3,4,5)]
matOrder <- mat[,c(1,2,5,3,4)]

#withlabels
tiff("OzoneHeatmap_labeled.tiff", width = 4, height = 10, units = 'in', res = 100)
pheatmap(matOrder, breaks=seq(from=-thr, to=thr, length=101), show_rownames=FALSE,cluster_col=FALSE, cutree_rows=8, labels_col=c("14 days vs 7 hours","28 days vs 7 hours","80ppb vs normal atmospheric ozone","125ppb vs normal atmospheric ozone","225ppb vs normal atmospheric ozone"))
dev.off()

#withoutlabels
tiff("OzoneHeatmap_nolabel.tiff", width = 10, height = 10, units = 'in', res = 100)
pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101), show_rownames=FALSE, show_colnames=FALSE,cluster_col=FALSE, cutree_rows=8)
dev.off()

#white labels
tiff("OzoneHeatmap_whitelabel.tiff", width = 10, height = 10, units = 'in', res = 100)
pheatmap(matOrder, breaks=seq(from=-thr, to=thr, length=101), show_rownames=FALSE, labels_col=c("                                                ","      ","      ","      ","      "),cluster_col=FALSE, cutree_rows=8)
dev.off()



#------------------------------------------
# heatmap experiment 2 -vst
#------------------------------------------

library("RColorBrewer")
library("gplots")
select <- order(rowMeans(counts(dds_trmt,normalized=TRUE)),decreasing=TRUE)[1:350]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
vsd <- varianceStabilizingTransformation(dds_trmt, blind=TRUE)
heatmap.2(assay(vsd)[select,], col = hmcol,
          Rowv = TRUE, Colv = FALSE, scale="none",
          dendrogram="row", trace="none", margin=c(10, 6))
dev.copy(png,"DESeq2_heatmap3")
dev.off()





#------------------------------------------
# transcripts responding to ozone treatment differently across time
#------------------------------------------
dds_interaction<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,
                                            directory=directory,
                                            design= ~time + treatment + time:treatment)

dds_interaction<-DESeq(dds_interaction, test = c("LRT"), reduced = ~time + treatment)

res_interaction<-results(dds_interaction)
resOrdered_interaction<-res_interaction[order(res_interaction$padj),]
resSig_interaction <- subset(resOrdered_interaction, padj < 0.01)








#####################
# Unused Code
#####################

setwd("~/1_greenAsh/analysis/1_differentialExpression/counts-sorted/")
sampleFiles <- grep("ozone",list.files(),value=TRUE)
sampleFiles <- c("GA1-ozone-7hrs-control.txt",
                 "GA2-ozone-7hrs-80ppb.txt",
                 "GA3-ozone-7hrs-125ppb.txt",
                 "GA4-ozone-7hrs-225ppb.txt",
                 "GA5-ozone-14days-control.txt",
                 "GA6-ozone-14days-80ppb.txt",
                 "GA7-ozone-14days-125ppb.txt",
                 "GA8-ozone-14days-225ppb.txt",
                 "GA9-ozone-28days-control.txt",
                 "GA10-ozone-28days-80ppb.txt",
                 "GA11-ozone-28days-125ppb.txt",
                 "GA12-ozone-28days-225ppb.txt")


time=c( rep("7hr",4),
        rep("14d",4),
        rep("28d",4))

timeReg=c( rep(7,4),
           rep(14,4),
           rep(28,4))

treatment=rep( c("ctrl","80pb", "125ppb", "225ppb"), 3)
trmtReg=rep( c(0,80,125,225), 3)

sampleTable<-data.frame(
  sampleName=sampleFiles,
  filename=sampleFiles,
  time=time,
  treatment=treatment,
  timeReg=timeReg,
  trmtReg=trmtReg)

directory<-("~/1_greenAsh/analysis/1_differentialExpression/counts-sorted")
setwd("~/1_greenAsh/analysis/1_differentialExpression")


write.csv(as.data.frame(resSigUp),
          file="Ozone_QuadTrmt.Up.DE.csv")
write.csv(as.data.frame(resSigDown),
          file="Ozone_QuadTrmt.Dw.DE.csv")

# quad time
dds_quad_time<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,
                                          directory=directory,
                                          design= ~trmtReg + timeReg + I(timeReg^2))
dds_quad_time<-DESeq(dds_quad_time, test = c("LRT"), reduced = ~ trmtReg + timeReg)
res<-results(dds_quad_time)
resOrdered<-res[order(res$padj),]
resSigQuadTime <- subset(resOrdered, padj < 0.05)
resSigUp <- subset(resSigQuadTime, log2FoldChange > 0)
resSigDown <- subset(resSigQuadTime, log2FoldChange < 0)
write.csv(as.data.frame(resSigUp),
          file="Ozone_QuadTime.Up.DE.csv")
write.csv(as.data.frame(resSigDown),
          file="Ozone_QuadTime.Dw.DE.csv")

# time
dds_time<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,
                                     directory=directory,
                                     design= ~trmtReg + timeReg)
dds_time<-DESeq(dds_time, test = c("LRT"), reduced = ~ trmtReg)
res<-results(dds_time)
resOrdered<-res[order(res$padj),]
resSigTime <- subset(resOrdered, padj < 0.05)
resSigUp <- subset(resSigTime, log2FoldChange > 0)
resSigDown <- subset(resSigTime, log2FoldChange < 0)
v_reduceTime <- row.names(resSigTime)
write.csv(as.data.frame(resSigUp),
          file="Ozone_Time.Up.DE.csv")
write.csv(as.data.frame(resSigDown),
          file="Ozone_Time.Dw.DE.csv")

# treatment
dds_trmt<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,
                                     directory=directory,
                                     design= ~trmtReg + timeReg)

dds_trmt<-DESeq(dds_trmt, test = c("LRT"), reduced = ~ timeReg)

res<-results(dds_trmt)
resOrdered<-res[order(res$padj),]
resSigTrmt <- subset(resOrdered, padj < 0.05)
resSigUp <- subset(resSigTrmt, log2FoldChange > 0)
resSigDown <- subset(resSigTrmt, log2FoldChange < 0)
v_reduceTrmt <- row.names(resSigTrmt)
write.csv(as.data.frame(resSigUp),
          file="Ozone_Trmt.Up.DE.csv")
write.csv(as.data.frame(resSigDown),
          file="Ozone_Trmt.Dw.DE.csv")


# time FACTOR
dds_timeFac<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,
                                        directory=directory,
                                        design= ~treatment + time)
dds_timeFac<-DESeq(dds_timeFac, test = c("LRT"), reduced = ~ treatment)
res<-results(dds_timeFac)
resOrdered<-res[order(res$padj),]
resSigTimeFac <- subset(resOrdered, padj < 0.05)
resSigUp <- subset(resSigTimeFac, log2FoldChange > 0)
resSigDown <- subset(resSigTimeFac, log2FoldChange < 0)
v_time <- row.names(resSigTimeFac)
write.csv(as.data.frame(resSigUp),
          file="Ozone_TimeFac.Up.DE.csv")
write.csv(as.data.frame(resSigDown),
          file="Ozone_TimeFac.Dw.DE.csv")

# treatment FACTOR
dds_trmtFac<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,
                                        directory=directory,
                                        design= ~treatment + time)

dds_trmtFac<-DESeq(dds_trmtFac, test = c("LRT"), reduced = ~ time)
res<-results(dds_trmtFac)
resOrdered<-res[order(res$padj),]
resSigTrmtFac <- subset(resOrdered, padj < 0.05)
resSigUp <- subset(resSigTrmtFac, log2FoldChange > 0)
resSigDown <- subset(resSigTrmtFac, log2FoldChange < 0)
v_trmt<- row.names(resSigTrmtFac)
write.csv(as.data.frame(resSigUp),
          file="Ozone_TrmtFac.Up.DE.csv")
write.csv(as.data.frame(resSigDown),
          file="Ozone_TrmtFac.Dw.DE.csv")