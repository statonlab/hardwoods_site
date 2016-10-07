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

time=c( rep(336,4),
        rep(672,4),
        rep(7,4))

treatment=c(rep( c(0,80,125,225), 3))

sampleTable<-data.frame(
  sampleName=sampleFiles,
  filename=sampleFiles,
  time=time,
  treatment=treatment)

directory <- ("/lustre/projects/staton/projects/hardwoods/unigenes/1.1_tulip_poplar/version_10_13_2014/20_deseq2/htseq/")

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
resSiqQuadTrmt <- subset(resSigQuadTrmt, !is.na(padj))
resSigUp <- subset(resSigQuadTrmt, log2FoldChange > 0)
resSigDown <- subset(resSigQuadTrmt, log2FoldChange < 0)

# write outputs
setwd("/lustre/projects/staton/projects/hardwoods/unigenes/1.1_tulip_poplar/version_10_13_2014/20_deseq2/")
write.csv(as.data.frame(resSigUp),
          file="OzoneQuad.DE.up.csv")
write.csv(as.data.frame(resSigDown),
          file="OzoneQuad.DE.dw.csv")