setwd("D:/Experiment (실험)/11.MDR2 KO/2.Colesevelam/Akkermansia/3차실험_20231013/RNA-seq")
library(tximport)
library(GenomicFeatures)
#txdb2 <- makeTxDbFromGFF(file="Mus_musculus.GRCm39.111.gff3.gz")
txdb2 <- makeTxDbFromGFF(file="gencode.vM34.annotation.gff3.gz")
transcripts(txdb2, columns=c("tx_id", "tx_name"))
genes(txdb2)
txdb2
columns(txdb2)


library(org.Mm.eg.db)
library(tximeta)
columns(org.Mm.eg.db)

### extract gene_name in '.gff3'
library(rtracklayer)
#aa <- import.gff3("Mus_musculus.GRCm39.111.gff3.gz") # 25s
aa2 <- import.gff3("gencode.vM34.annotation.gff3.gz")
gene.DF <- data.frame(aa2$gene_name,aa2$gene_id)
u.gene.DF<- unique(gene.DF)
head(u.gene.DF)

u.gene.DF[,2] <- gsub("\\..*", "", u.gene.DF[,2])

### 
k <- keys(txdb2, keytype = "TXNAME")
tx2gene <- select(txdb2, k, "GENEID", "TXNAME")
head(tx2gene)

#tx2gene에서 dot remove하기
tx2gene[,1] <- gsub("\\..*", "", tx2gene[,1])

###
setwd("D:/Experiment (실험)/11.MDR2 KO/2.Colesevelam/Akkermansia/3차실험_20231013/RNA-seq/Salmon")
dir <- getwd()
files = paste0(dir, "/",c("WV1","WV3","WV6","WA3","WA8","WA14","KV4", "KV6","KV7","KA6","KA7","KA10"), "/quant.sf")

#paste0("KO_Col_", 1:13) 1씩 증가

names(files) <- paste0( c("WV1","WV3","WV6","WA3","WA8","WA14","KV4", "KV6","KV7","KA6","KA7","KA10"))
file.exists(files)
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)

#write.csv(as.data.frame(u.gene.DF), file = "gene_name_20240401.csv")
#write.csv(as.data.frame(tx2gene), file = "gene_name_tx2gene_20240401.csv")
#write.csv(as.data.frame(txi.salmon$counts), file = "gene_count_20240401.csv")
#write.csv(as.data.frame(txi.salmon$abundance), file = "gene_abundance_20240401.csv")
#write.csv(as.data.frame(txi.salmon$countsFromAbundance), file = "gene_countsfromabundance_20240401.csv")


###deseq2 - using counts
setwd("D:/Experiment (실험)/11.MDR2 KO/2.Colesevelam/Akkermansia/3차실험_20231013/RNA-seq")
library(DESeq2)
sampleTable <- data.frame(
  condition = factor(c(rep("WT_VEH",3), rep("WT_AKK",3), rep("KO_VEH",3), rep("KO_AKK", 3))),
  genotype = factor(c(rep("WT",6), rep("KO",6))),
  treatment = factor(c(rep("VEH",3), rep("AKK",3), rep("VEH",3), rep("AKK",3)))
)
rownames(sampleTable) <- colnames(txi.salmon$counts)
dds <- DESeqDataSetFromTximport(txi.salmon, sampleTable, ~condition)
deseqresult <- DESeq(dds)
deseqresult

#condition1 WT_VEH vs KO_VEH
res01 <- results(deseqresult, contrast = c("condition","WT_VEH","KO_VEH"), alpha = 0.05)
res01
summary(res01)
plotMA(res01, ylim=c(-10,10))
#write.csv(as.data.frame(res01), file = "1.WT_VEHvsKO_VEH_results_count_20240401.csv")

#condition2 WT_AKK vs KO_AKK
res02 <- results(deseqresult, contrast = c("condition","WT_AKK","KO_AKK"), alpha = 0.05)
res02
summary(res02)
plotMA(res02, ylim=c(-10,10))
#write.csv(as.data.frame(res02), file = "2.WT_AKKvsKO_AKK_results_count_20240401.csv")

#condition3 WT_VEH vs WT_AKK
res03 <- results(deseqresult, contrast = c("condition","WT_VEH","WT_AKK"), alpha = 0.05)
res03
summary(res03)
plotMA(res03, ylim=c(-10,10))
#write.csv(as.data.frame(res03), file = "3.WT_VEHvsWT_AKK_results_count_20240401.csv")

#condition4 KO_VEH vs KO_AKK
res04 <- results(deseqresult, contrast = c("condition","KO_VEH","KO_AKK"), alpha = 0.05)
res04
summary(res04)
plotMA(res04, ylim=c(-10,10))
#write.csv(as.data.frame(res04), file = "4.KO_VEHvsKO_AKK_results_count_20240401.csv")

#####volcanoplot1 - WT-VEH vs KO-VEH (WT_VEH에 많으면 +, KO_VEH에 많으면 -)
res01genename <- read.csv("1.WT_VEHvsKO_VEH_results_count_20240401_genename.csv", stringsAsFactors = F)
cut_lfc <- 0.5
cut_pvalue <- 0.05
par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)
topT.res01 <- as.data.frame(res01genename)
# Adjusted P values
with(topT.res01, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", col='grey', 
                      cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~Q~value)))
with(subset(topT.res01, padj<cut_pvalue & log2FoldChange>cut_lfc), points(log2FoldChange, -log10(padj), pch=20, col='#E63943', cex=1.5))
with(subset(topT.res01, padj<cut_pvalue & log2FoldChange<(-cut_lfc)), points(log2FoldChange, -log10(padj), pch=20, col='#2A9D8F', cex=1.5))
## Add lines for FC and P-value cut-off
abline(v=0, col='black', lty=3, lwd=1.0)
abline(v=-cut_lfc, col='black', lty=4, lwd=2.0)
abline(v=cut_lfc, col='black', lty=4, lwd=2.0)
abline(h=-log10(max(topT.res01$padj[topT.res01$padj<cut_pvalue], na.rm=TRUE)), col='black', lty=4, lwd=2.0)

## labeling on volcano plot
#x.extrem.res01 <- na.omit(topT.res01[topT.res01$log2FoldChange < (-4) & -log10(topT.res01$padj) > 3, ])
#y.extrem.res01 <- na.omit(topT.res01[topT.res01$log2FoldChange < 0 & -log10(topT.res01$padj) > 20, ])
#
##with(x.extrem, text(log2FoldChange, -log10(padj), labels=X, pos=4, cex=0.8) )
#with(y.extrem.res01, text(log2FoldChange, -log10(padj), labels=X, pos=2, cex=0.7) )
#
#x.extrem2.res01 <- na.omit(topT.res01[topT.res01$log2FoldChange > 10 ,])
#y.extrem2.res01 <- na.omit(topT.res01[topT.res01$log2FoldChange > 0 & -log10(topT.res01$padj) > 20, ])
#
#with(x.extrem2.res01, text(log2FoldChange, -log10(padj), labels=X, pos=1, cex=0.7) )
#with(y.extrem2.res01, text(log2FoldChange, -log10(padj), labels=X, pos=c(4,3,2,1,1,1), cex=0.7) )


#####volcanoplot2 - WT-AKK vs KO-AKK
res02genename <- read.csv("2.WT_AKKvsKO_AKK_results_count_20240401_genename.csv", stringsAsFactors = F)
par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)
topT.res02 <- as.data.frame(res02genename)
# Adjusted P values
with(topT.res02, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", col='grey', 
                      cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~Q~value)))
with(subset(topT.res02, padj<cut_pvalue & log2FoldChange>cut_lfc), points(log2FoldChange, -log10(padj), pch=20, col='#4477AA', cex=1.5))
with(subset(topT.res02, padj<cut_pvalue & log2FoldChange<(-cut_lfc)), points(log2FoldChange, -log10(padj), pch=20, col='#DAA520', cex=1.5))
## Add lines for FC and P-value cut-off
abline(v=0, col='black', lty=3, lwd=1.0)
abline(v=-cut_lfc, col='black', lty=4, lwd=2.0)
abline(v=cut_lfc, col='black', lty=4, lwd=2.0)
abline(h=-log10(max(topT.res02$padj[topT.res02$padj<cut_pvalue], na.rm=TRUE)), col='black', lty=4, lwd=2.0)

## labeling on volcano plot
#x.extrem.res02 <- na.omit(topT.res02[topT.res02$log2FoldChange < (-4) & -log10(topT.res02$padj) > 2, ])
#y.extrem.res02 <- na.omit(topT.res02[topT.res02$log2FoldChange < 0 & -log10(topT.res02$padj) > 16, ])
#
#with(x.extrem.res02, text(log2FoldChange, -log10(padj), labels=X, pos=3, cex=0.7) )
#with(y.extrem.res02, text(log2FoldChange, -log10(padj), labels=X, pos=c(3,2,2,2,3), cex=0.7) )
#
#x.extrem2.res02 <- na.omit(topT.res02[topT.res02$log2FoldChange > 10 ,])
#y.extrem2.res02 <- na.omit(topT.res02[topT.res02$log2FoldChange > 0 & -log10(topT.res02$padj) > 16, ])
#
#with(x.extrem2.res02, text(log2FoldChange, -log10(padj), labels=X, pos=3, cex=0.7) )
#with(y.extrem2.res02, text(log2FoldChange, -log10(padj), labels=X, pos=c(3,3,1), cex=0.7) )


#####volcanoplot3 - WT-VEH vs WT-AKK
res03genename <- read.csv("3.WT_VEHvsWT_AKK_results_count_20240401_genename.csv", stringsAsFactors = F)
par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)
topT.res03 <- as.data.frame(res03genename)
# Adjusted P values
with(topT.res03, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", col='grey', 
                      cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~Q~value)))
with(subset(topT.res03, padj<cut_pvalue & log2FoldChange>cut_lfc), points(log2FoldChange, -log10(padj), pch=20, col='#E63943', cex=1.5))
with(subset(topT.res03, padj<cut_pvalue & log2FoldChange<(-cut_lfc)), points(log2FoldChange, -log10(padj), pch=20, col='#4477AA', cex=1.5))
## Add lines for FC and P-value cut-off
abline(v=0, col='black', lty=3, lwd=1.0)
abline(v=-cut_lfc, col='black', lty=4, lwd=2.0)
abline(v=cut_lfc, col='black', lty=4, lwd=2.0)
abline(h=-log10(max(topT.res03$padj[topT.res03$padj<cut_pvalue], na.rm=TRUE)), col='black', lty=4, lwd=2.0)

## labeling on volcano plot
#x.extrem.res03 <- na.omit(topT.res03[topT.res03$log2FoldChange < (-9) & -log10(topT.res03$padj) > 2, ])
#y.extrem.res03 <- na.omit(topT.res03[topT.res03$log2FoldChange < 0 & -log10(topT.res03$padj) > 80, ])
#
#with(x.extrem.res03, text(log2FoldChange, -log10(padj), labels=X, pos=c(3,2,3,1), cex=0.7) )
#with(y.extrem.res03, text(log2FoldChange, -log10(padj), labels=X, pos=c(3,3,1,1,3,3,4,3,1), cex=0.7) )
#
#x.extrem2.res03 <- na.omit(topT.res03[topT.res03$log2FoldChange > 7 ,])
#y.extrem2.res03 <- na.omit(topT.res03[topT.res03$log2FoldChange > 0 & -log10(topT.res03$padj) > 50, ])
#
#with(x.extrem2.res03, text(log2FoldChange, -log10(padj), labels=X, pos=3, cex=0.7) )
#with(y.extrem2.res03, text(log2FoldChange, -log10(padj), labels=X, pos=c(3,3,4), cex=0.7) )

#####volcanoplot4 - KO-VEH vs KO-AKK
res04genename <- read.csv("4.KO_VEHvsKO_AKK_results_count_20240401_genename.csv", stringsAsFactors = F)
par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)
topT.res04 <- as.data.frame(res04genename)
# Adjusted P value
with(topT.res04, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", col='grey', 
                      cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~Q~value)))
with(subset(topT.res04, padj<cut_pvalue & log2FoldChange>cut_lfc), points(log2FoldChange, -log10(padj),
                                                                          pch=20, col='#2A9D8F', cex=1.5))
with(subset(topT.res04, padj<cut_pvalue & log2FoldChange<(-cut_lfc)), points(log2FoldChange, -log10(padj), 
                                                                             pch=20, col='#DAA520', cex=1.5))
## Add lines for FC and P-value cut-off
abline(v=0, col='black', lty=3, lwd=1.0)
abline(v=-cut_lfc, col='black', lty=4, lwd=2.0)
abline(v=cut_lfc, col='black', lty=4, lwd=2.0)
abline(h=-log10(max(topT.res04$padj[topT.res04$padj<cut_pvalue], na.rm=TRUE)), col='black', lty=4, lwd=2.0)

## labeling on volcano plot
#x.extrem.res04 <- na.omit(topT.res04[topT.res04$log2FoldChange < (-8) ,])
#y.extrem.res04 <- na.omit(topT.res04[topT.res04$log2FoldChange < 0 & -log10(topT.res04$padj) > 47, ])
#
#with(x.extrem.res04, text(log2FoldChange, -log10(padj), labels=X, pos=c(3,2,4,4), cex=0.7) )
#with(y.extrem.res04, text(log2FoldChange, -log10(padj), labels=X, pos=c(2,1,4,2,3,4), cex=0.7) )

#sort(as.vector(na.omit((-log10(topT$padj))[ topT$log2FoldChange >0  &  -log10(topT$padj) > 20 ])))

#x.extrem2.res04 <- na.omit(topT.res04[topT.res04$log2FoldChange > 5 & topT.res04$padj < 0.05,])
#y.extrem2.res04 <- na.omit(topT.res04[topT.res04$log2FoldChange > 0 & -log10(topT.res04$padj) > 30, ])
#
#with(x.extrem2.res04, text(log2FoldChange, -log10(padj), labels=X, pos=3, cex=0.7) )
#with(y.extrem2.res04, text(log2FoldChange, -log10(padj), labels=X, pos=3, cex=0.7) )

genes_of_interest <- c("S100a9", "S100a8", "Lbp", "Mmp8", "Ccl9", "Ccl6", 
                       "Casp7", "Nlrc3", "Gps2", "Setd6", "Fanca", "Tlr5", "Nr1h4")

label_data <- subset(topT.res04, gene %in% genes_of_interest)

#with(label_data, text(log2FoldChange, -log10(padj), labels=gene, pos=c(4,4,4,4,2,4,2,2,2,2,4,4,2), cex=1.0, col="black"))
with(label_data, text(log2FoldChange, -log10(padj), labels=gene, pos=c(4,4,3,4,2,4,3,2,2,3,4,4,3), cex=0.6, col="black"))

# 라벨 및 선 추가
for (i in 1:nrow(label_data)) {
  gene <- label_data$gene[i]
  x <- label_data$log2FoldChange[i]
  y <- -log10(label_data$padj[i])
  text_pos <- ifelse(x > 0, 4, 2) # x 좌표에 따라 라벨 방향 결정
  text(x, y, labels=gene, pos=text_pos, cex=1.0, col="black")
  segments(x, y, x + ifelse(x > 0, 0.5, -0.5), y + 2, col="black", lty=2) # 점과 라벨 연결
}

#Extracting transformed values - vst, rlog, ntd
vsd <- vst(deseqresult, blind = "FALSE")
rld <- rlog(deseqresult, blind =FALSE)
ntd <- normTransform(deseqresult)
head(assay(vsd), 3)
head(assay(rld), 3)
head(assay(ntd), 3)

#Pca plot
colData(vsd)
par(mar=c(4,4,3,3))
plotPCA(vsd, intgroup = c("condition"))
plotPCA(rld, intgroup = c("condition"))
plotPCA(ntd, intgroup = c("condition"))

#Heatmap
library(pheatmap)

### matching gene name using all 4group (오류남)
head(u.gene.DF)
g.t.DF <- data.frame("aa2.gene_id" = rownames(assay(ntd)))
matching.g.name <- merge(g.t.DF,u.gene.DF, by= "aa2.gene_id")
dim(g.t.DF)
dim(matching.g.name)
matching.g.name

#heatmap using ntd
dds2<-estimateSizeFactors(dds)
heatmapselect <- order(rowMeans(counts(dds2,normalized=TRUE)),decreasing=TRUE)
heatmapdf <- as.data.frame(colData(dds)[,c("treatment","genotype")])
pheatmap(assay(ntd)[heatmapselect,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=heatmapdf)

#table(rownames(assay(ntd)) == matching.g.name$aa2.gene_id) # must have all true
#a.ntd2 <- assay(ntd) 
##write.csv(a.ntd2, file = "ntdform_genecount_20220727.csv")
#
#rownames(a.ntd2)<- matching.g.name$aa.gene_name
#table(matching.g.name$aa.gene_id == names(rowMeans(counts(deseqresult, normalized=TRUE))))  # must have all true
##select <- order(rowMeans(counts(deseqresult, normalized=TRUE)), decreasing=T)[1:200] #값 수정시 gene개수 늘어남
##df <- as.data.frame(colData(deseqresult)[,c("condition")])
##df <- as.data.frame(colData(deseqresult)[,c("treatment", "genotype")])
##pheatmap(a.ntd2[select,], cluster_rows=T, show_rownames=T, cluster_cols=T, annotation_col=df)
#
#select.all <- order(rowMeans(counts(deseqresult, normalized=TRUE)), decreasing=T)
#length(select.all)
#
#a.ntd3 <- a.ntd2[select.all, ]
#na.a.ntd3 <- na.omit(a.ntd3)
#pheatmap(na.a.ntd3[1:200,], cluster_rows=T, show_rownames=T, cluster_cols=T, annotation_col=df)
#
##heatmap using rld
#table(rownames(assay(rld)) == matching.g.name$aa.gene_id) # must have all true
#a.rld2 <- assay(rld)
##write.csv(a.rld2, file = "rldform_genecount_20220726.csv")
#rownames(a.rld2)<- matching.g.name$aa.gene_name
#table(matching.g.name$aa.gene_id == names(rowMeans(counts(deseqresult, normalized=TRUE)))) 
#select.all <- order(rowMeans(counts(deseqresult, normalized=TRUE)), decreasing=T)
#length(select.all)
#
#a.rld3 <- a.rld2[select.all, ]
#pheatmap(a.rld3[1:200,], cluster_rows=T, show_rownames=T, cluster_cols=T, annotation_col=df)
#
##heatmap using vst
#table(rownames(assay(vsd)) == matching.g.name$aa.gene_id) # must have all true
#a.vsd2 <- assay(vsd)
##write.csv(a.vsd2, file = "vstform_genecount_20220726.csv")
#rownames(a.vsd2)<- matching.g.name$aa.gene_name
#table(matching.g.name$aa.gene_id == names(rowMeans(counts(deseqresult, normalized=TRUE))))  # must have all true
#select.all <- order(rowMeans(counts(deseqresult, normalized=TRUE)), decreasing=T)
#length(select.all)
#
#a.vsd3 <- a.vsd2[select.all, ]
#pheatmap(a.vsd3[1:200,], cluster_rows=T, show_rownames=T, cluster_cols=T, annotation_col=df)

#annotation color list
#df <- as.data.frame(colData(dds)[,"condition"])
#anno_color=list(
#  condition=c("WT_VEH"="#E63943", "WT_AKK"="#4477AA", "KO_VEH"="#2A9D8F", "KO_AKK"="#DAA520"))



# heatmap
bileacidtarget.csv <- read.csv("Genelist/bile acid gene.csv", stringsAsFactors = F )
rownames(bileacidtarget.csv) <- bileacidtarget.csv$Bile.acid
bileacidtarget.csv <- bileacidtarget.csv[-1]

bileacidtarget.mat <- apply(bileacidtarget.csv, 1, scale)
rownames(bileacidtarget.mat) <- colnames(bileacidtarget.csv)
bileacidtarget.mat <- t(bileacidtarget.mat)
bileacidtarget.mat <- na.omit(bileacidtarget.mat)
pheatmap(bileacidtarget.mat, cluster_rows=T, show_rownames=T, cluster_cols=F, annotation_col=heatmapdf, annotation_colors = anno_color)
#pheatmap(bileacidtarget.mat, cluster_rows=T, show_rownames=T, cluster_cols=F, annotation_col=df, annotation_colors = anno_color)

#
is.na(inflamtarget.csv)

inflamtarget.csv <- read.csv("Genelist/Inflammation_gene.csv", header = TRUE, stringsAsFactors = F)
rownames(inflamtarget.csv) <- inflamtarget.csv$Inflammation.gene.name
inflamtarget.csv <- inflamtarget.csv[-1]


inflamtarget.mat <- apply(inflamtarget.csv, 1, scale)
rownames(inflamtarget.mat) <- colnames(inflamtarget.csv)
inflamtarget.mat <- t(inflamtarget.mat)
inflamtarget.mat <- na.omit(inflamtarget.mat)
pheatmap(inflamtarget.mat, cluster_rows=T, show_rownames=F, cluster_cols=F, annotation_col=heatmapdf)
pheatmap(inflamtarget.mat[1:350,], cluster_rows=T, show_rownames=F, cluster_cols=F, annotation_col=heatmapdf)

#######q<0.05
inflamtarget.csv <- read.csv("Genelist/Inflammation_gene_q0.05(kvka).csv", header = TRUE, stringsAsFactors = F)
rownames(inflamtarget.csv) <- inflamtarget.csv$Inflammation.gene.name
inflamtarget.csv <- inflamtarget.csv[-1]


inflamtarget.mat <- apply(inflamtarget.csv, 1, scale)
rownames(inflamtarget.mat) <- colnames(inflamtarget.csv)
inflamtarget.mat <- t(inflamtarget.mat)
inflamtarget.mat <- na.omit(inflamtarget.mat)
pheatmap(inflamtarget.mat, cluster_rows=T, show_rownames=T, cluster_cols=T, fontsize_row = 4)

#가로로
transposed_mat <- t(inflamtarget.mat)

# 히트맵 생성
pheatmap(
  transposed_mat, 
  cluster_rows = TRUE, 
  show_rownames = TRUE, 
  cluster_cols = TRUE, 
  fontsize_row = 8,  # sample 이름 글씨 크기
  fontsize_col = 8, # gene 이름 글씨 크기
  angle_col=90
)

# 관심 있는 gene list만 선택
selected_genes <- c("Nlrc3", "Gps2", "Setd6", "Fanca", "Tlr5", "Nr1h4", 
                    "S100a9", "S100a8", "Lbp", "Mmp8", "Ccl6", "Casp7", "Itgb1")
filtered_mat <- inflamtarget.mat[rownames(inflamtarget.mat) %in% selected_genes, ]

# 데이터를 전치 (samples를 rows로, genes를 columns로)
transposed_mat <- t(filtered_mat)

# 필터링된 데이터로 히트맵 생성
pheatmap(
  filtered_mat, 
  cluster_rows = TRUE, 
  show_rownames = TRUE, 
  cluster_cols = TRUE, 
  fontsize_row = 8,  # sample 이름 글씨 크기
  fontsize_col = 8,   # gene 이름 글씨 크기
  angle_col=45
)

##KEGG - Juneyoung
#BiocManager::install(c("gage","GO.db","AnnotationDbi","org.Mm.eg.db"))
#BiocManager::install("pathview")
library(gage)
library(pathview)
kg.mmu <- kegg.gsets(species="mmu")
kegg.sigmet.gs <- kg.mmu$kg.sets[kg.mmu$sigmet.idx]
kegg.dise.gs <- kg.mmu$kg.sets[kg.mmu$dise.idx]
# set up go database
go.mmu <- go.gsets(species="mouse")
go.bp.gs <- go.mmu$go.sets[go.mmu$go.subs$BP]
go.mf.gs <- go.mmu$go.sets[go.mmu$go.subs$MF]
go.cc.gs <- go.mmu$go.sets[go.mmu$go.subs$CC]

#https://www.hadriengourle.com/tutorials/rna/
res01.name<- gsub("\\..*","",row.names(res01))
row.names(res01) <- res01.name

res02.name<- gsub("\\..*","",row.names(res02))
row.names(res02) <- res02.name

res03.name<- gsub("\\..*","",row.names(res03))
row.names(res03) <- res03.name

res04.name<- gsub("\\..*","",row.names(res04))
row.names(res04) <- res04.name


res01$symbol <- mapIds(org.Mm.eg.db, keys=row.names(res01), column="SYMBOL", keytype="ENSEMBL", multiVals="first")
res01$entrez <- mapIds(org.Mm.eg.db, keys=row.names(res01), column="ENTREZID", keytype="ENSEMBL", multiVals="first")
res01$name <- mapIds(org.Mm.eg.db, keys=row.names(res01), column="GENENAME",keytype="ENSEMBL",  multiVals="first")

res02$symbol <- mapIds(org.Mm.eg.db, keys=row.names(res02), column="SYMBOL", keytype="ENSEMBL", multiVals="first")
res02$entrez <- mapIds(org.Mm.eg.db, keys=row.names(res02), column="ENTREZID", keytype="ENSEMBL", multiVals="first")
res02$name <- mapIds(org.Mm.eg.db, keys=row.names(res02), column="GENENAME",keytype="ENSEMBL",  multiVals="first")

res03$symbol <- mapIds(org.Mm.eg.db, keys=row.names(res03), column="SYMBOL", keytype="ENSEMBL", multiVals="first")
res03$entrez <- mapIds(org.Mm.eg.db, keys=row.names(res03), column="ENTREZID", keytype="ENSEMBL", multiVals="first")
res03$name <- mapIds(org.Mm.eg.db, keys=row.names(res03), column="GENENAME",keytype="ENSEMBL",  multiVals="first")

res04$symbol <- mapIds(org.Mm.eg.db, keys=row.names(res04), column="SYMBOL", keytype="ENSEMBL", multiVals="first")
res04$entrez <- mapIds(org.Mm.eg.db, keys=row.names(res04), column="ENTREZID", keytype="ENSEMBL", multiVals="first")
res04$name <- mapIds(org.Mm.eg.db, keys=row.names(res04), column="GENENAME",keytype="ENSEMBL",  multiVals="first")

res01.fixna = na.omit(res01)
res02.fixna = na.omit(res02)
res03.fixna = na.omit(res03)
res04.fixna = na.omit(res04)

# grab the log fold changes for everything
res01.fixna.fc <- res01.fixna$log2FoldChange #WT_VEH vs KO_VEH  
res02.fixna.fc <- res02.fixna$log2FoldChange #WT_AKK vs KO_AKK
res03.fixna.fc <- res03.fixna$log2FoldChange #WT_VEH vs WT_AKK
res04.fixna.fc <- res04.fixna$log2FoldChange #KO_VEH vs KO_AKK

names(res01.fixna.fc) <- res01.fixna$entrez #WT_VEH vs KO_VEH
names(res02.fixna.fc) <- res02.fixna$entrez #WT_AKK vs KO_AKK
names(res03.fixna.fc) <- res03.fixna$entrez #WT_VEH vs WT_AKK
names(res04.fixna.fc) <- res04.fixna$entrez #KO_VEH vs KO_AKK

# Run enrichment analysis on all log fc - res01 - 이렇게 할경우 greater가 모두 q.value 1로 나옴.
#all : gsets = c(sigmet.idx,dise.idx)
###
transform.sig.file <- function(x){
  ta <- cbind(as.data.frame(x$greater), data.frame("type" = "greater"))
  tb <- cbind(as.data.frame(x$less), data.frame("type" = "less"))
  result <- rbind(ta, tb)
  return(result)
}

fc.kegg.sigmet.p.res01 <- gage(res01.fixna.fc, gsets = kegg.sigmet.gs, compare="unpaired", use.fold=T, rank.test=T, same.dir = TRUE) #sigaling and metabolism
fc.kegg.dise.p.res01   <- gage(res01.fixna.fc, gsets = kegg.dise.gs, compare="unpaired", use.fold=T, rank.test=T, same.dir = TRUE) # disease
fc.go.bp.p.res01       <- gage(res01.fixna.fc, gsets = go.bp.gs, compare="unpaired", use.fold=T, rank.test=T) #biological process
fc.go.mf.p.res01       <- gage(res01.fixna.fc, gsets = go.mf.gs, compare="unpaired", use.fold=T, rank.test=T) #molecular function
fc.go.cc.p.res01       <- gage(res01.fixna.fc, gsets = go.cc.gs, compare="unpaired", use.fold=T, rank.test=T) #cellular component

fc.kegg.sigmet.p.res02 <- gage(res02.fixna.fc, gsets = kegg.sigmet.gs, compare="unpaired", use.fold=T, rank.test=T)
fc.kegg.dise.p.res02   <- gage(res02.fixna.fc, gsets = kegg.dise.gs, compare="unpaired", use.fold=T, rank.test=T)
fc.go.bp.p.res02       <- gage(res02.fixna.fc, gsets = go.bp.gs, compare="unpaired", use.fold=T, rank.test=T)
fc.go.mf.p.res02       <- gage(res02.fixna.fc, gsets = go.mf.gs, compare="unpaired", use.fold=T, rank.test=T)
fc.go.cc.p.res02       <- gage(res02.fixna.fc, gsets = go.cc.gs, compare="unpaired", use.fold=T, rank.test=T)

fc.kegg.sigmet.p.res03 <- gage(res03.fixna.fc, gsets = kegg.sigmet.gs, compare="unpaired", use.fold=T, rank.test=T)
fc.kegg.dise.p.res03   <- gage(res03.fixna.fc, gsets = kegg.dise.gs, compare="unpaired", use.fold=T, rank.test=T)
fc.go.bp.p.res03       <- gage(res03.fixna.fc, gsets = go.bp.gs, compare="unpaired", use.fold=T, rank.test=T)
fc.go.mf.p.res03       <- gage(res03.fixna.fc, gsets = go.mf.gs, compare="unpaired", use.fold=T, rank.test=T)
fc.go.cc.p.res03       <- gage(res03.fixna.fc, gsets = go.cc.gs, compare="unpaired", use.fold=T, rank.test=T)

fc.kegg.sigmet.p.res04 <- gage(res04.fixna.fc, gsets = kegg.sigmet.gs, compare="unpaired", use.fold=T, rank.test=T)
fc.kegg.dise.p.res04   <- gage(res04.fixna.fc, gsets = kegg.dise.gs, compare="unpaired", use.fold=T, rank.test=T)
fc.go.bp.p.res04       <- gage(res04.fixna.fc, gsets = go.bp.gs, compare="unpaired", use.fold=T, rank.test=T)
fc.go.mf.p.res04       <- gage(res04.fixna.fc, gsets = go.mf.gs, compare="unpaired", use.fold=T, rank.test=T)
fc.go.cc.p.res04       <- gage(res04.fixna.fc, gsets = go.cc.gs, compare="unpaired", use.fold=T, rank.test=T)

write.csv(transform.sig.file(fc.kegg.sigmet.p.res01), file="Pathway/fc.kegg.sigmet.p.res01.csv") #save all pathway
write.csv(transform.sig.file(fc.kegg.dise.p.res01  ), file="Pathway/fc.kegg.dise.p.res01.csv")
write.csv(transform.sig.file(fc.go.bp.p.res01      ), file="Pathway/fc.go.bp.p.res01.csv")
write.csv(transform.sig.file(fc.go.mf.p.res01      ), file="Pathway/fc.go.mf.p.res01.csv")
write.csv(transform.sig.file(fc.go.cc.p.res01      ), file="Pathway/fc.go.cc.p.res01.csv")

write.csv(transform.sig.file(fc.kegg.sigmet.p.res02), file="Pathway/fc.kegg.sigmet.p.res02.csv")
write.csv(transform.sig.file(fc.kegg.dise.p.res02  ), file="Pathway/fc.kegg.dise.p.res02.csv")
write.csv(transform.sig.file(fc.go.bp.p.res02      ), file="Pathway/fc.go.bp.p.res02.csv")
write.csv(transform.sig.file(fc.go.mf.p.res02      ), file="Pathway/fc.go.mf.p.res02.csv")
write.csv(transform.sig.file(fc.go.cc.p.res02      ), file="Pathway/fc.go.cc.p.res02.csv")

write.csv(transform.sig.file(fc.kegg.sigmet.p.res03), file="Pathway/fc.kegg.sigmet.p.res03.csv")
write.csv(transform.sig.file(fc.kegg.dise.p.res03  ), file="Pathway/fc.kegg.dise.p.res03.csv")
write.csv(transform.sig.file(fc.go.bp.p.res03      ), file="Pathway/fc.go.bp.p.res03.csv")
write.csv(transform.sig.file(fc.go.mf.p.res03      ), file="Pathway/fc.go.mf.p.res03.csv")
write.csv(transform.sig.file(fc.go.cc.p.res03      ), file="Pathway/fc.go.cc.p.res03.csv")

write.csv(transform.sig.file(fc.kegg.sigmet.p.res04), file="Pathway/fc.kegg.sigmet.p.res04.csv")
write.csv(transform.sig.file(fc.kegg.dise.p.res04  ), file="Pathway/fc.kegg.dise.p.res04.csv")
write.csv(transform.sig.file(fc.go.bp.p.res04      ), file="Pathway/fc.go.bp.p.res04.csv")
write.csv(transform.sig.file(fc.go.mf.p.res04      ), file="Pathway/fc.go.mf.p.res04.csv")
write.csv(transform.sig.file(fc.go.cc.p.res04      ), file="Pathway/fc.go.cc.p.res04.csv")

#WT_VEH vs KO_VEH
#WT_AKK vs KO_AKK
#WT_VEH vs WT_AKK
#KO_VEH vs KO_AKK

#volcanoplot kegg
cut_lfckegg <- 1
cut_qvaluekegg <- 0.1
par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)

#res01 - kegg sigmet
volcano.fc.kegg.sigmet.p.res01<- read.csv("Pathway_volcano/fc.kegg.sigmet.p.res01.csv", header = T, stringsAsFactors = F )
cut_lfckegg <- 1
cut_qvaluekegg <- 0.1
par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)
volcano.fc.kegg.sigmet.p.res01=na.omit(volcano.fc.kegg.sigmet.p.res01)
volcano.kegg.sigmet.res01 <- as.data.frame(volcano.fc.kegg.sigmet.p.res01)

with(volcano.kegg.sigmet.res01, plot(stat.mean, -log10(q.val), pch=20, main="Volcano plot", col='grey', 
                                     cex=1.0, xlab=bquote(stat.mean), ylab=bquote(~-log[10]~Q~value)))
with(subset(volcano.kegg.sigmet.res01, q.val<cut_qvaluekegg & stat.mean>cut_lfckegg), points(stat.mean, -log10(q.val),
                                                                                             pch=20, col='firebrick1', cex=1.5))
with(subset(volcano.kegg.sigmet.res01, q.val<cut_qvaluekegg & stat.mean<(-cut_lfckegg)), points(stat.mean, -log10(q.val), 
                                                                                                pch=20, col='dodgerblue', cex=1.5))
abline(v=0, col='black', lty=4, lwd=1.0)
abline(v=-cut_lfckegg, col='black', lty=4, lwd=2.0)
abline(v=cut_lfckegg, col='black', lty=4, lwd=2.0)
abline(h=-log10(max(volcano.kegg.sigmet.res01$q.val[volcano.kegg.sigmet.res01$q.val<cut_qvaluekegg], na.rm=TRUE)), col='black', lty=4, lwd=2.0)

#x.extrem.kegg.sigmet.res01 <- na.omit(volcano.kegg.sigmet.res01[volcano.kegg.sigmet.res01$stat.mean < (-2) & -log10(volcano.kegg.sigmet.res01$q.val) > 6, ])
#with(x.extrem.kegg.sigmet.res01, text(stat.mean, -log10(q.val), labels=X, pos=4, cex=0.7) )
#x.extrem2.kegg.sigmet.res01 <- na.omit(volcano.kegg.sigmet.res01[volcano.kegg.sigmet.res01$stat.mean > 2 & -log10(volcano.kegg.sigmet.res01$q.val) > 4, ])
#with(x.extrem2.kegg.sigmet.res01, text(stat.mean, -log10(q.val), labels=X, pos=c(2,3), cex=0.7) )


#res01-kegg disease
volcano.fc.kegg.dise.p.res01<- read.csv("Pathway_volcano/fc.kegg.dise.p.res01.csv", header = T, stringsAsFactors = F )
volcano.fc.kegg.dise.p.res01=na.omit(volcano.fc.kegg.dise.p.res01)
volcano.kegg.dise.res01 <- as.data.frame(volcano.fc.kegg.dise.p.res01)

with(volcano.kegg.dise.res01, plot(stat.mean, -log10(q.val), pch=20, main="Volcano plot", col='grey', 
                                   cex=1.0, xlab=bquote(stat.mean), ylab=bquote(~-log[10]~Q~value)))
with(subset(volcano.kegg.dise.res01, q.val<cut_qvaluekegg & stat.mean>cut_lfckegg), points(stat.mean, -log10(q.val),
                                                                                           pch=20, col='firebrick1', cex=1.5))
with(subset(volcano.kegg.dise.res01, q.val<cut_qvaluekegg & stat.mean<(-cut_lfckegg)), points(stat.mean, -log10(q.val), 
                                                                                              pch=20, col='dodgerblue', cex=1.5))
abline(v=0, col='black', lty=4, lwd=1.0)
abline(v=-cut_lfckegg, col='black', lty=4, lwd=2.0)
abline(v=cut_lfckegg, col='black', lty=4, lwd=2.0)
abline(h=-log10(max(volcano.kegg.dise.res01$q.val[volcano.kegg.dise.res01$q.val<cut_qvaluekegg], na.rm=TRUE)), col='black', lty=4, lwd=2.0)

#x.extrem.kegg.disease.res01 <- na.omit(volcano.kegg.dise.res01[volcano.kegg.dise.res01$stat.mean < (-2) & -log10(volcano.kegg.dise.res01$q.val) > 1, ])
#with(x.extrem.kegg.disease.res01, text(stat.mean, -log10(q.val), labels=X, pos=4, cex=0.7) )
#x.extrem2.kegg.disease.res01 <- na.omit(volcano.kegg.dise.res01[volcano.kegg.dise.res01$stat.mean > 2 & -log10(volcano.kegg.dise.res01$q.val) > 4, ])
#with(x.extrem2.kegg.disease.res01, text(stat.mean, -log10(q.val), labels=X, pos=2, cex=0.7) )


#res02 - kegg sigmet
volcano.fc.kegg.sigmet.p.res02<- read.csv("Pathway_volcano/fc.kegg.sigmet.p.res02.csv", header = T, stringsAsFactors = F )
volcano.fc.kegg.sigmet.p.res02=na.omit(volcano.fc.kegg.sigmet.p.res02)
volcano.kegg.sigmet.res02 <- as.data.frame(volcano.fc.kegg.sigmet.p.res02)

with(volcano.kegg.sigmet.res02, plot(stat.mean, -log10(q.val), pch=20, main="Volcano plot", col='grey', 
                                     cex=1.0, xlab=bquote(stat.mean), ylab=bquote(~-log[10]~Q~value)))
with(subset(volcano.kegg.sigmet.res02, q.val<cut_qvaluekegg & stat.mean>cut_lfckegg), points(stat.mean, -log10(q.val),
                                                                                             pch=20, col='firebrick1', cex=1.5))
with(subset(volcano.kegg.sigmet.res02, q.val<cut_qvaluekegg & stat.mean<(-cut_lfckegg)), points(stat.mean, -log10(q.val), 
                                                                                                pch=20, col='dodgerblue', cex=1.5))
abline(v=0, col='black', lty=4, lwd=1.0)
abline(v=-cut_lfckegg, col='black', lty=4, lwd=2.0)
abline(v=cut_lfckegg, col='black', lty=4, lwd=2.0)
abline(h=-log10(max(volcano.kegg.sigmet.res02$q.val[volcano.kegg.sigmet.res02$q.val<cut_qvaluekegg], na.rm=TRUE)), col='black', lty=4, lwd=2.0)

#x.extrem.kegg.sigmet.res02 <- na.omit(volcano.kegg.sigmet.res02[volcano.kegg.sigmet.res02$stat.mean < (-2) & -log10(volcano.kegg.sigmet.res02$q.val) > 1, ])
#with(x.extrem.kegg.sigmet.res02, text(stat.mean, -log10(q.val), labels=X, pos=4, cex=0.7) )
#x.extrem2.kegg.sigmet.res02 <- na.omit(volcano.kegg.sigmet.res02[volcano.kegg.sigmet.res02$stat.mean > 2 & -log10(volcano.kegg.sigmet.res02$q.val) > 1, ])
#with(x.extrem2.kegg.sigmet.res02, text(stat.mean, -log10(q.val), labels=X, pos=c(2,3), cex=0.7) )

#res02 - kegg disease
volcano.fc.kegg.dise.p.res02<- read.csv("Pathway_volcano/fc.kegg.dise.p.res02.csv", header = T, stringsAsFactors = F )
volcano.fc.kegg.dise.p.res02=na.omit(volcano.fc.kegg.dise.p.res02)
volcano.kegg.dise.res02 <- as.data.frame(volcano.fc.kegg.dise.p.res02)

with(volcano.kegg.dise.res02, plot(stat.mean, -log10(q.val), pch=20, main="Volcano plot", col='grey', 
                                   cex=1.0, xlab=bquote(stat.mean), ylab=bquote(~-log[10]~Q~value)))
with(subset(volcano.kegg.dise.res02, q.val<cut_qvaluekegg & stat.mean>cut_lfckegg), points(stat.mean, -log10(q.val),
                                                                                           pch=20, col='firebrick1', cex=1.5))
with(subset(volcano.kegg.dise.res02, q.val<cut_qvaluekegg & stat.mean<(-cut_lfckegg)), points(stat.mean, -log10(q.val), 
                                                                                              pch=20, col='dodgerblue', cex=1.5))
abline(v=0, col='black', lty=4, lwd=1.0)
abline(v=-cut_lfckegg, col='black', lty=4, lwd=2.0)
abline(v=cut_lfckegg, col='black', lty=4, lwd=2.0)
abline(h=-log10(max(volcano.kegg.dise.res02$q.val[volcano.kegg.dise.res02$q.val<cut_qvaluekegg], na.rm=TRUE)), col='black', lty=4, lwd=2.0)

#res03 - kegg sigmet
volcano.fc.kegg.sigmet.p.res03<- read.csv("Pathway/fc.kegg.sigmet.p.res03.csv", header = T, stringsAsFactors = F )
volcano.fc.kegg.sigmet.p.res03=na.omit(volcano.fc.kegg.sigmet.p.res03)
volcano.kegg.sigmet.res03 <- as.data.frame(volcano.fc.kegg.sigmet.p.res03)

with(volcano.kegg.sigmet.res03, plot(stat.mean, -log10(q.val), pch=20, main="Volcano plot", col='grey', 
                                     cex=1.0, xlab=bquote(stat.mean), ylab=bquote(~-log[10]~Q~value)))
with(subset(volcano.kegg.sigmet.res03, q.val<cut_qvaluekegg & stat.mean>cut_lfckegg), points(stat.mean, -log10(q.val),
                                                                                             pch=20, col='firebrick1', cex=1.5))
with(subset(volcano.kegg.sigmet.res03, q.val<cut_qvaluekegg & stat.mean<(-cut_lfckegg)), points(stat.mean, -log10(q.val), 
                                                                                                pch=20, col='dodgerblue', cex=1.5))
abline(v=0, col='black', lty=4, lwd=1.0)
abline(v=-cut_lfckegg, col='black', lty=4, lwd=2.0)
abline(v=cut_lfckegg, col='black', lty=4, lwd=2.0)
abline(h=-log10(max(volcano.kegg.sigmet.res03$q.val[volcano.kegg.sigmet.res03$q.val<cut_qvaluekegg], na.rm=TRUE)), col='black', lty=4, lwd=2.0)

#x.extrem.kegg.sigmet.res03 <- na.omit(volcano.kegg.sigmet.res03[volcano.kegg.sigmet.res03$stat.mean < (-2) & -log10(volcano.kegg.sigmet.res03$q.val) > 5, ])
#with(x.extrem.kegg.sigmet.res03, text(stat.mean, -log10(q.val), labels=X, pos=4, cex=0.7) )
#x.extrem2.kegg.sigmet.res03 <- na.omit(volcano.kegg.sigmet.res03[volcano.kegg.sigmet.res03$stat.mean > 2 & -log10(volcano.kegg.sigmet.res03$q.val) > 5, ])
#with(x.extrem2.kegg.sigmet.res03, text(stat.mean, -log10(q.val), labels=X, pos=2, cex=0.7) )

#res03 - kegg disease
volcano.fc.kegg.dise.p.res03<- read.csv("Pathway/fc.kegg.dise.p.res03.csv", header = T, stringsAsFactors = F )
volcano.fc.kegg.dise.p.res03=na.omit(volcano.fc.kegg.dise.p.res03)
volcano.kegg.dise.res03 <- as.data.frame(volcano.fc.kegg.dise.p.res03)

with(volcano.kegg.dise.res03, plot(stat.mean, -log10(q.val), pch=20, main="Volcano plot", col='grey', 
                                   cex=1.0, xlab=bquote(stat.mean), ylab=bquote(~-log[10]~Q~value)))
with(subset(volcano.kegg.dise.res03, q.val<cut_qvaluekegg & stat.mean>cut_lfckegg), points(stat.mean, -log10(q.val),
                                                                                           pch=20, col='firebrick1', cex=1.5))
with(subset(volcano.kegg.dise.res03, q.val<cut_qvaluekegg & stat.mean<(-cut_lfckegg)), points(stat.mean, -log10(q.val), 
                                                                                              pch=20, col='dodgerblue', cex=1.5))
abline(v=0, col='black', lty=4, lwd=1.0)
abline(v=-cut_lfckegg, col='black', lty=4, lwd=2.0)
abline(v=cut_lfckegg, col='black', lty=4, lwd=2.0)
abline(h=-log10(max(volcano.kegg.dise.res03$q.val[volcano.kegg.dise.res03$q.val<cut_qvaluekegg], na.rm=TRUE)), col='black', lty=4, lwd=2.0)

#res04 - kegg sigmet
cut_lfckegg <- 1
cut_qvaluekegg <- 0.05
#par(mar=c(10,10,10,10), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)
par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)

volcano.fc.kegg.sigmet.p.res04<- read.csv("Pathway/fc.kegg.sigmet.p.res04.csv", header = T, stringsAsFactors = F )
volcano.fc.kegg.sigmet.p.res04=na.omit(volcano.fc.kegg.sigmet.p.res04)
volcano.kegg.sigmet.res04 <- as.data.frame(volcano.fc.kegg.sigmet.p.res04)

with(volcano.kegg.sigmet.res04, plot(stat.mean, -log10(q.val), pch=20, main="Volcano plot", col='grey', 
                                     cex=1.0, xlab=bquote(stat.mean), ylab=bquote(~-log[10]~Q~value)))
with(subset(volcano.kegg.sigmet.res04, q.val<cut_qvaluekegg & stat.mean>cut_lfckegg), points(stat.mean, -log10(q.val),
                                                                                             pch=20, col='#2A9D8F', cex=1.5))
with(subset(volcano.kegg.sigmet.res04, q.val<cut_qvaluekegg & stat.mean<(-cut_lfckegg)), points(stat.mean, -log10(q.val), 
                                                                                                pch=20, col='#DAA520', cex=1.5))
abline(v=0, col='black', lty=4, lwd=1.0)
abline(v=-cut_lfckegg, col='black', lty=4, lwd=2.0)
abline(v=cut_lfckegg, col='black', lty=4, lwd=2.0)
abline(h=-log10(max(volcano.kegg.sigmet.res04$q.val[volcano.kegg.sigmet.res04$q.val<cut_qvaluekegg], na.rm=TRUE)), col='black', lty=4, lwd=2.0)

#x.extrem.kegg.sigmet.res04 <- na.omit(volcano.kegg.sigmet.res04[volcano.kegg.sigmet.res04$stat.mean < (-2) & -log10(volcano.kegg.sigmet.res04$q.val) > 3, ])
#with(x.extrem.kegg.sigmet.res04, text(stat.mean, -log10(q.val), labels=X, pos=c(4,3,4), cex=0.7) )
#x.extrem2.kegg.sigmet.res04 <- na.omit(volcano.kegg.sigmet.res04[volcano.kegg.sigmet.res04$stat.mean > 5 & -log10(volcano.kegg.sigmet.res04$q.val) > 5, ])
#with(x.extrem2.kegg.sigmet.res04, text(stat.mean, -log10(q.val), labels=X, pos=2, cex=0.7) )


#res04 - kegg disease
volcano.fc.kegg.dise.p.res04<- read.csv("Pathway/fc.kegg.dise.p.res04.csv", header = T, stringsAsFactors = F )
volcano.fc.kegg.dise.p.res04=na.omit(volcano.fc.kegg.dise.p.res04)
volcano.kegg.dise.res04 <- as.data.frame(volcano.fc.kegg.dise.p.res04)

with(volcano.kegg.dise.res04, plot(stat.mean, -log10(q.val), pch=20, main="Volcano plot", col='grey', 
                                   cex=1.0, xlab=bquote(stat.mean), ylab=bquote(~-log[10]~Q~value)))
with(subset(volcano.kegg.dise.res04, q.val<cut_qvaluekegg & stat.mean>cut_lfckegg), points(stat.mean, -log10(q.val),
                                                                                           pch=20, col='firebrick1', cex=1.5))
with(subset(volcano.kegg.dise.res04, q.val<cut_qvaluekegg & stat.mean<(-cut_lfckegg)), points(stat.mean, -log10(q.val), 
                                                                                              pch=20, col='dodgerblue', cex=1.5))
abline(v=0, col='black', lty=4, lwd=1.0)
abline(v=-cut_lfckegg, col='black', lty=4, lwd=2.0)
abline(v=cut_lfckegg, col='black', lty=4, lwd=2.0)
abline(h=-log10(max(volcano.kegg.dise.res04$q.val[volcano.kegg.dise.res04$q.val<cut_qvaluekegg], na.rm=TRUE)), col='black', lty=4, lwd=2.0)


#sigGeneSet(fc.kegg.sigmet.p.res01, cutoff=0.1, qpval="q.val", heatmap=T, pdf.size=c(14, 14), cexRow=0.6, cexCol = 0.6)
fc.kegg.sigmet.p.res01.sig<-sigGeneSet(fc.kegg.sigmet.p.res01, cutoff=0.1, qpval="q.val", heatmap=F, outname = "fc.kegg.sigmet.p.res01.sig")
fc.kegg.dise.p.res01.sig<-sigGeneSet(fc.kegg.dise.p.res01, cutoff=0.1, qpval="q.val", heatmap=F, outname = "fc.kegg.dise.p.res01.sig")
fc.go.bp.p.res01.sig<-sigGeneSet(fc.go.bp.p.res01, cutoff=0.1, qpval="q.val", heatmap=F, outname = "fc.go.bp.p.res01.sig")
fc.go.mf.p.res01.sig<-sigGeneSet(fc.go.mf.p.res01, cutoff=0.1, qpval="q.val", heatmap=F, outname = "fc.go.mf.p.res01.sig")
fc.go.cc.p.res01.sig<-sigGeneSet(fc.go.cc.p.res01, cutoff=0.1, qpval="q.val", heatmap=F, outname = "fc.go.cc.p.res01.sig")

fc.kegg.sigmet.p.res02.sig<-sigGeneSet(fc.kegg.sigmet.p.res02, cutoff=0.1, qpval="q.val", heatmap=F, outname = "fc.kegg.sigmet.p.res02.sig")
fc.kegg.dise.p.res02.sig<-sigGeneSet(fc.kegg.dise.p.res02, cutoff=0.1, qpval="q.val", heatmap=F, outname = "fc.kegg.sigmet.p.res02.sig")
fc.go.bp.p.res02.sig<-sigGeneSet(fc.go.bp.p.res02, cutoff=0.1, qpval="q.val", heatmap=F, outname = "fc.go.bp.p.res02.sig")
fc.go.mf.p.res02.sig<-sigGeneSet(fc.go.mf.p.res02, cutoff=0.1, qpval="q.val", heatmap=F, outname = "fc.go.mf.p.res02.sig")
fc.go.cc.p.res02.sig<-sigGeneSet(fc.go.cc.p.res02, cutoff=0.1, qpval="q.val", heatmap=F, outname = "fc.go.cc.p.res02.sig")

fc.kegg.sigmet.p.res03.sig<-sigGeneSet(fc.kegg.sigmet.p.res03, cutoff=0.1, qpval="q.val", heatmap=F, outname = "fc.kegg.sigmet.p.res03.sig")
fc.kegg.dise.p.res03.sig<-sigGeneSet(fc.kegg.dise.p.res03, cutoff=0.1, qpval="q.val", heatmap=F, outname = "fc.kegg.sigmet.p.res03.sig")
fc.go.bp.p.res03.sig<-sigGeneSet(fc.go.bp.p.res03, cutoff=0.1, qpval="q.val", heatmap=F, outname = "fc.go.bp.p.res03.sig")
fc.go.mf.p.res03.sig<-sigGeneSet(fc.go.mf.p.res03, cutoff=0.1, qpval="q.val", heatmap=F, outname = "fc.go.mf.p.res03.sig")
fc.go.cc.p.res03.sig<-sigGeneSet(fc.go.cc.p.res03, cutoff=0.1, qpval="q.val", heatmap=F, outname = "fc.go.cc.p.res03.sig")

fc.kegg.sigmet.p.res04.sig<-sigGeneSet(fc.kegg.sigmet.p.res04, cutoff=0.1, qpval="q.val", heatmap=F, outname = "fc.kegg.sigmet.p.res04.sig")
fc.kegg.dise.p.res04.sig<-sigGeneSet(fc.kegg.dise.p.res04, cutoff=0.1, qpval="q.val", heatmap=F, outname = "fc.kegg.sigmet.p.res04.sig")
fc.go.bp.p.res04.sig<-sigGeneSet(fc.go.bp.p.res04, cutoff=0.1, qpval="q.val", heatmap=F, outname = "fc.go.bp.p.res04.sig")
fc.go.mf.p.res04.sig<-sigGeneSet(fc.go.mf.p.res04, cutoff=0.1, qpval="q.val", heatmap=F, outname = "fc.go.mf.p.res04.sig")
fc.go.cc.p.res04.sig<-sigGeneSet(fc.go.cc.p.res04, cutoff=0.1, qpval="q.val", heatmap=F, outname = "fc.go.cc.p.res04.sig")

# covert the kegg results to data frames
fc.kegg.sigmet.p.up.res01 <- as.data.frame(fc.kegg.sigmet.p.res01$greater)
fc.kegg.dise.p.up.res01 <- as.data.frame(fc.kegg.dise.p.res01$greater)

fc.kegg.sigmet.p.down.res01 <- as.data.frame(fc.kegg.sigmet.p.res01$less)
fc.kegg.dise.p.down.res01 <- as.data.frame(fc.kegg.dise.p.res01$less)

fc.kegg.sigmet.p.up.res02 <- as.data.frame(fc.kegg.sigmet.p.res02$greater)
fc.kegg.dise.p.up.res02 <- as.data.frame(fc.kegg.dise.p.res02$greater)

fc.kegg.sigmet.p.down.res02 <- as.data.frame(fc.kegg.sigmet.p.res02$less)
fc.kegg.dise.p.down.res02 <- as.data.frame(fc.kegg.dise.p.res02$less)

fc.kegg.sigmet.p.up.res03 <- as.data.frame(fc.kegg.sigmet.p.res03$greater)
fc.kegg.dise.p.up.res03 <- as.data.frame(fc.kegg.dise.p.res03$greater)

fc.kegg.sigmet.p.down.res03 <- as.data.frame(fc.kegg.sigmet.p.res03$less)
fc.kegg.dise.p.down.res03 <- as.data.frame(fc.kegg.dise.p.res03$less)

fc.kegg.sigmet.p.up.res04 <- as.data.frame(fc.kegg.sigmet.p.res04$greater)
fc.kegg.dise.p.up.res04 <- as.data.frame(fc.kegg.dise.p.res04$greater)

fc.kegg.sigmet.p.down.res04 <- as.data.frame(fc.kegg.sigmet.p.res04$less)
fc.kegg.dise.p.down.res04 <- as.data.frame(fc.kegg.dise.p.res04$less)

#write.csv(fc.kegg.sigmet.p.up.res01, file   = "Pathway_updown/1.fc.kegg.sigmet.p.up.res01.csv")
#write.csv(fc.kegg.dise.p.up.res01, file     = "Pathway_updown/1.fc.kegg.dise.p.up.res01.csv")
#write.csv(fc.kegg.sigmet.p.down.res01, file = "Pathway_updown/1.fc.kegg.sigmet.p.down.res01.csv")
#write.csv(fc.kegg.dise.p.down.res01, file   = "Pathway_updown/1.fc.kegg.dise.p.down.res01.csv")
#write.csv(fc.kegg.sigmet.p.up.res02, file   = "Pathway_updown/2.fc.kegg.sigmet.p.up.res02.csv")
#write.csv(fc.kegg.dise.p.up.res02, file     = "Pathway_updown/2.fc.kegg.dise.p.up.res02.csv")
#write.csv(fc.kegg.sigmet.p.down.res02, file = "Pathway_updown/2.fc.kegg.sigmet.p.down.res02.csv")
#write.csv(fc.kegg.dise.p.down.res02, file   = "Pathway_updown/2.fc.kegg.dise.p.down.res02.csv")
#write.csv(fc.kegg.sigmet.p.up.res03, file   = "Pathway_updown/3.fc.kegg.sigmet.p.up.res03.csv")
#write.csv(fc.kegg.dise.p.up.res03, file     = "Pathway_updown/3.fc.kegg.dise.p.up.res03.csv")
#write.csv(fc.kegg.sigmet.p.down.res03, file = "Pathway_updown/3.fc.kegg.sigmet.p.down.res03.csv")
#write.csv(fc.kegg.dise.p.down.res03, file   = "Pathway_updown/3.fc.kegg.dise.p.down.res03.csv")
#write.csv(fc.kegg.sigmet.p.up.res04, file   = "Pathway_updown/4.fc.kegg.sigmet.p.up.res04.csv")
#write.csv(fc.kegg.dise.p.up.res04, file     = "Pathway_updown/4.fc.kegg.dise.p.up.res04.csv")
#write.csv(fc.kegg.sigmet.p.down.res04, file = "Pathway_updown/4.fc.kegg.sigmet.p.down.res04.csv")
#write.csv(fc.kegg.dise.p.down.res04, file   = "Pathway_updown/4.fc.kegg.dise.p.down.res04.csv")


# convert the go results to data frames
fc.go.bp.p.up.res01   <- as.data.frame(fc.go.bp.p.res01$greater)
fc.go.mf.p.up.res01   <- as.data.frame(fc.go.mf.p.res01$greater)
fc.go.cc.p.up.res01   <- as.data.frame(fc.go.cc.p.res01$greater)
fc.go.bp.p.down.res01 <- as.data.frame(fc.go.bp.p.res01$less)
fc.go.mf.p.down.res01 <- as.data.frame(fc.go.mf.p.res01$less)
fc.go.cc.p.down.res01 <- as.data.frame(fc.go.cc.p.res01$less)

fc.go.bp.p.up.res02   <- as.data.frame(fc.go.bp.p.res02$greater)
fc.go.mf.p.up.res02   <- as.data.frame(fc.go.mf.p.res02$greater)
fc.go.cc.p.up.res02   <- as.data.frame(fc.go.cc.p.res02$greater)
fc.go.bp.p.down.res02 <- as.data.frame(fc.go.bp.p.res02$less)
fc.go.mf.p.down.res02 <- as.data.frame(fc.go.mf.p.res02$less)
fc.go.cc.p.down.res02 <- as.data.frame(fc.go.cc.p.res02$less)

fc.go.bp.p.up.res03   <- as.data.frame(fc.go.bp.p.res03$greater)
fc.go.mf.p.up.res03   <- as.data.frame(fc.go.mf.p.res03$greater)
fc.go.cc.p.up.res03   <- as.data.frame(fc.go.cc.p.res03$greater)
fc.go.bp.p.down.res03 <- as.data.frame(fc.go.bp.p.res03$less)
fc.go.mf.p.down.res03 <- as.data.frame(fc.go.mf.p.res03$less)
fc.go.cc.p.down.res03 <- as.data.frame(fc.go.cc.p.res03$less)

fc.go.bp.p.up.res04   <- as.data.frame(fc.go.bp.p.res04$greater)
fc.go.mf.p.up.res04   <- as.data.frame(fc.go.mf.p.res04$greater)
fc.go.cc.p.up.res04   <- as.data.frame(fc.go.cc.p.res04$greater)
fc.go.bp.p.down.res04 <- as.data.frame(fc.go.bp.p.res04$less)
fc.go.mf.p.down.res04 <- as.data.frame(fc.go.mf.p.res04$less)
fc.go.cc.p.down.res04 <- as.data.frame(fc.go.cc.p.res04$less)

#write.csv(fc.go.bp.p.up.res01  , file = "Pathway_GO term/1.fc.go.bp.p.up.res01.csv")
#write.csv(fc.go.mf.p.up.res01  , file = "Pathway_GO term/1.fc.go.mf.p.up.res01.csv")
#write.csv(fc.go.cc.p.up.res01  , file = "Pathway_GO term/1.fc.go.cc.p.up.res01.csv")
#write.csv(fc.go.bp.p.down.res01, file = "Pathway_GO term/1.fc.go.bp.p.down.res01.csv")
#write.csv(fc.go.mf.p.down.res01, file = "Pathway_GO term/1.fc.go.mf.p.down.res01.csv")
#write.csv(fc.go.cc.p.down.res01, file = "Pathway_GO term/1.fc.go.cc.p.down.res01.csv")
#
#write.csv(fc.go.bp.p.up.res02  , file = "Pathway_GO term/2.fc.go.bp.p.up.res02.csv")
#write.csv(fc.go.mf.p.up.res02  , file = "Pathway_GO term/2.fc.go.mf.p.up.res02.csv")
#write.csv(fc.go.cc.p.up.res02  , file = "Pathway_GO term/2.fc.go.cc.p.up.res02.csv")
#write.csv(fc.go.bp.p.down.res02, file = "Pathway_GO term/2.fc.go.bp.p.down.res02.csv")
#write.csv(fc.go.mf.p.down.res02, file = "Pathway_GO term/2.fc.go.mf.p.down.res02.csv")
#write.csv(fc.go.cc.p.down.res02, file = "Pathway_GO term/2.fc.go.cc.p.down.res02.csv")
#
#write.csv(fc.go.bp.p.up.res03  , file = "Pathway_GO term/3.fc.go.bp.p.up.res03.csv")
#write.csv(fc.go.mf.p.up.res03  , file = "Pathway_GO term/3.fc.go.mf.p.up.res03.csv")
#write.csv(fc.go.cc.p.up.res03  , file = "Pathway_GO term/3.fc.go.cc.p.up.res03.csv")
#write.csv(fc.go.bp.p.down.res03, file = "Pathway_GO term/3.fc.go.bp.p.down.res03.csv")
#write.csv(fc.go.mf.p.down.res03, file = "Pathway_GO term/3.fc.go.mf.p.down.res03.csv")
#write.csv(fc.go.cc.p.down.res03, file = "Pathway_GO term/3.fc.go.cc.p.down.res03.csv")
#
#write.csv(fc.go.bp.p.up.res04  , file = "Pathway_GO term/4.fc.go.bp.p.up.res04.csv")
#write.csv(fc.go.mf.p.up.res04  , file = "Pathway_GO term/4.fc.go.mf.p.up.res04.csv")
#write.csv(fc.go.cc.p.up.res04  , file = "Pathway_GO term/4.fc.go.cc.p.up.res04.csv")
#write.csv(fc.go.bp.p.down.res04, file = "Pathway_GO term/4.fc.go.bp.p.down.res04.csv")
#write.csv(fc.go.mf.p.down.res04, file = "Pathway_GO term/4.fc.go.mf.p.down.res04.csv")
#write.csv(fc.go.cc.p.down.res04, file = "Pathway_GO term/4.fc.go.cc.p.down.res04.csv")

#stat.mean graph from KEGG
library(ggplot2)
#res01: KO_veh vs KO_col
#res02: WT_veh vs WT_col
#res03: WT_veh vs KO_veh
#res04: WT_col vs KO_col

sm_res01_sigmet <- read.csv("Pathway_updown/1.fc.kegg.sigmet.p.updown.res01.csv", header = T, stringsAsFactors = F)
sm_res01_sigmet$plus <- ifelse(sm_res0_sigmet$stat.mean > 0, "plus", "minus")
ggplot(data = sm_res01_sigmet, aes(x=reorder(KEGG, stat.mean), y=stat.mean, fill=plus))+
  geom_bar(stat = "identity")+
  coord_flip()+
  scale_fill_manual(values=c("dodgerblue", "firebrick1"), guide="none")+
  ylab("Stat.mean")+
  xlab("KEGG pathway")

sm_res02_sigmet <- read.csv("Pathway_KEGG/2.fc.kegg.sigmet.p.updown.res02.csv", header = T, stringsAsFactors = F)
sm_res02_sigmet$plus <- ifelse(sm_res02_sigmet$stat.mean > 0, "plus", "minus")
ggplot(data = sm_res02_sigmet, aes(x=reorder(KEGG, stat.mean), y=stat.mean, fill=plus))+
  geom_bar(stat = "identity")+
  coord_flip()+
  scale_fill_manual(values=c("dodgerblue", "firebrick1"), guide="none")+
  ylab("Stat.mean")+
  xlab("KEGG pathway")

sm_res03_sigmet <- read.csv("Pathway_KEGG/3.fc.kegg.sigmet.p.updown.res03.csv", header = T, stringsAsFactors = F)
sm_res03_sigmet$plus <- ifelse(sm_res03_sigmet$stat.mean > 0, "plus", "minus")
ggplot(data = sm_res03_sigmet, aes(x=reorder(KEGG, stat.mean), y=stat.mean, fill=plus))+
  geom_bar(stat = "identity")+
  coord_flip()+
  scale_fill_manual(values=c("dodgerblue", "firebrick1"), guide="none")+
  ylab("Stat.mean")+
  xlab("KEGG pathway")

sm_res04_sigmet <- read.csv("Pathway_updown/4.fc.kegg.sigmet.p.up.down.res04.csv", header = T, stringsAsFactors = F)
sm_res04_sigmet$plus <- ifelse(sm_res04_sigmet$stat.mean > 0, "plus", "minus")
ggplot(data = sm_res04_sigmet, aes(x=reorder(KEGG, stat.mean), y=stat.mean, fill=plus))+
  geom_bar(stat = "identity")+
  coord_flip()+
  scale_fill_manual(values=c("#daa520", "#2a9d8f"), guide="none")+
  ylab("Stat.mean")+
  xlab("KEGG pathway")

#stat.mean graph from GO - biological process - too much
sm_res01_bp <- read.csv("Pathway_volcano/1.fc.go.bp.p.updown.res01.csv", header = T, stringsAsFactors = F)
sm_res01_bp$plus <- ifelse(sm_res01_bp$stat.mean > 0, "plus", "minus")
ggplot(data = sm_res01_bp, aes(x=reorder(GO, stat.mean), y=stat.mean, fill=plus))+
  geom_bar(stat = "identity")+
  coord_flip()+
  scale_fill_manual(values=c("dodgerblue", "firebrick1"), guide="none")+
  ylab("Stat.mean")+
  xlab("GO term_biological process")



#KEGG pathway mapping
#https://www.r-bloggers.com/2015/12/tutorial-rna-seq-differential-expression-pathway-analysis-with-sailfish-deseq2-gage-and-pathview/
#https://bioinformaticsandme.tistory.com/104


#res01: KO_veh vs KO_col
#res02: WT_veh vs WT_col
#res03: WT_veh vs KO_veh
#res04: WT_col vs KO_col
setwd("D:/Experiment (실험)/11.MDR2 KO/2.Colesevelam/RNA-seq/salmon/Pathway_KEGG q0.05")
library(ggplot2)
#res01
sm_res01_sigmet <- read.csv("1.updown.res01.csv", header = T, stringsAsFactors = F)
sm_res01_sigmet$plus <- ifelse(sm_res01_sigmet$stat.mean > 0, "plus", "minus")
ggplot(data = sm_res01_sigmet, aes(x=reorder(KEGG, stat.mean), y=stat.mean, fill=plus))+
  geom_bar(stat = "identity")+
  coord_flip()+
  scale_fill_manual(values=c("#E7AC02", "#228833"), guide="none")+
  ylab("Stat.mean")+
  xlab("KEGG pathway")

#res02
sm_res02_sigmet <- read.csv("2.updown.res02.csv", header = T, stringsAsFactors = F)
sm_res02_sigmet$plus <- ifelse(sm_res02_sigmet$stat.mean > 0, "plus", "minus")
ggplot(data = sm_res02_sigmet, aes(x=reorder(KEGG, stat.mean), y=stat.mean, fill=plus))+
  geom_bar(stat = "identity")+
  coord_flip()+
  scale_fill_manual(values=c("dodgerblue", "firebrick1"), guide="none")+
  ylab("Stat.mean")+
  xlab("KEGG pathway")

#res03
sm_res03_sigmet <- read.csv("3.updown.res03.csv", header = T, stringsAsFactors = F)
sm_res03_sigmet$plus <- ifelse(sm_res03_sigmet$stat.mean > 0, "plus", "minus")
ggplot(data = sm_res03_sigmet, aes(x=reorder(KEGG, stat.mean), y=stat.mean, fill=plus))+
  geom_bar(stat = "identity")+
  coord_flip()+
  scale_fill_manual(values=c("#228833", "#EE6677"), guide="none")+
  ylab("Stat.mean")+
  xlab("KEGG pathway")

#res04
sm_res04_sigmet <- read.csv("4.updown.res04.csv", header = T, stringsAsFactors = F)
sm_res04_sigmet$plus <- ifelse(sm_res04_sigmet$stat.mean > 0, "plus", "minus")
ggplot(data = sm_res04_sigmet, aes(x=reorder(KEGG, stat.mean), y=stat.mean, fill=plus))+
  geom_bar(stat = "identity")+
  coord_flip()+
  scale_fill_manual(values=c("dodgerblue", "firebrick1"), guide="none")+
  ylab("Stat.mean")+
  xlab("KEGG pathway")
