#source("http://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
library("DESeq2")

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

time=c( rep("14day",4),
        rep("28day",4),
        rep("7hr",4))

treatment=c(rep( c("0ppb","80ppb","125ppb","ppb225"), 3))

sampleTable<-data.frame(
  sampleName=sampleFiles,
  filename=sampleFiles,
  time=time,
  treatment=treatment)

directory <- ("/lustre/projects/staton/projects/hardwoods/unigenes/1.1_tulip_poplar/version_10_13_2014/20_deseq2/htseq/")

#--------
# quad treatment
#---------
dds <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,
                                          directory=directory,
                                          design= ~treatment + time + time:treatment)
#dds<-DESeq(dds, test = c("LRT"), reduced = ~ time)
#dds<-DESeq(dds, test = c("LRT"), reduced = ~ treatment + time)
#dds<-DESeq(dds, test = c("Wald"))
dds <- DESeq(dds)

list <- resultsNames(dds)
mres<-sapply(list, function(x) results(dds, name=x))
sapply(mres, function(x) subset(x, padj < 0.05))

res<-results(dds)
resOrdered<-res[order(res$padj),]
resSig <- subset(resOrdered, padj < 0.05)
resSiq <- subset(resSig, !is.na(padj))
resSigUp <- subset(resSig, log2FoldChange > 0)
resSigDown <- subset(resSig, log2FoldChange < 0)]

# write outputs
setwd("/lustre/projects/staton/projects/hardwoods/unigenes/1.1_tulip_poplar/version_10_13_2014/20_deseq2/")
write.csv(as.data.frame(resSigUp),
          file="OzoneQuad.DE.up.csv")
write.csv(as.data.frame(resSigDown),
          file="OzoneQuad.DE.dw.csv")