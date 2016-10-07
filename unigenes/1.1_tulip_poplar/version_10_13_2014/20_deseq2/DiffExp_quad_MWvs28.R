#source("http://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
library("DESeq2")

#########
# Ozone #
#########
# compared 28 and 29 day experiments...

setwd("/lustre/projects/staton/projects/hardwoods/unigenes/1.1_tulip_poplar/version_10_13_2014/20_deseq2/htseq/")
sampleFiles <- grep(".txt.filtered",list.files(),value=TRUE)
sampleFiles <- c("28Day-CO_S9.counts.txt.filtered",
                 "28Day-80_S10.counts.txt.filtered",
                 "28Day-125_S11.counts.txt.filtered",
                 "28Day-225_S12.counts.txt.filtered",
                 "29DAY-80_S10.counts.txt.filtered",
                 "29Day-125_S15.counts.txt.filtered",
                 "29Day-225_S16.counts.txt.filtered",
                 "29Day-CO_S13.counts.txt.filtered")

wounding=c( rep("ctrl",4),
            rep("wounded",4))

time=c( rep("28day",4),
        rep("29day",4))

ozone=c(rep( c(0,80,125,225), 2))

sampleTable<-data.frame(
  sampleName=sampleFiles,
  filename=sampleFiles,
  wounding=wounding,
  ozone=ozone)

directory <- ("/lustre/projects/staton/projects/hardwoods/unigenes/1.1_tulip_poplar/version_10_13_2014/20_deseq2/htseq/")

#--------
# quad treatment
#---------
dds_quad_trmt <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,
                                          directory=directory,
                                          design= ~ozone + wounding)
dds_quad_trmt<-DESeq(dds_quad_trmt, test = c("LRT"), reduced = ~ ozone + wounding)
res<-results(dds_quad_trmt)
resOrdered<-res[order(res$padj),]
resSigQuadTrmt <- subset(resOrdered, padj < 0.05)
resSiqQuadTrmt <- subset(resSigQuadTrmt, !is.na(padj))
resSigUp <- subset(resSigQuadTrmt, log2FoldChange > 0)
resSigDown <- subset(resSigQuadTrmt, log2FoldChange < 0)
resSigUp
resSigDown

# write outputs
setwd("/lustre/projects/staton/projects/hardwoods/unigenes/1.1_tulip_poplar/version_10_13_2014/20_deseq2/")
write.csv(as.data.frame(resSigUp),
          file="OzoneQuadMWvs28.DE.up.csv")
write.csv(as.data.frame(resSigDown),
          file="OzoneQuadMWvs28.DE.dw.csv")