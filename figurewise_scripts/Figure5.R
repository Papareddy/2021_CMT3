library(tximport)
library(DESeq2)
library(tidyverse)
library(dplyr)
library(ggplot2)
gene.names<-read.table("/groups/nodine/lab/members/Ranjith/mCHG/figure5_rna/AGI_to_name.txt",header = T)
names(gene.names)[1]<-"GeneId"

--- Gene  Quantification and DeSeq -----
files <- list.files("/scratch-cbe/users/ranjith.papareddy/rCMT3_transcriptome_reanal/rSEM_V44/star_rsem", full.names = T, recursive = T, pattern = ".genes.results")

naming<-c(files)
names(files)<-naming
names(files)<- gsub(names(files),pattern = "/scratch-cbe/users/ranjith.papareddy/rCMT3_transcriptome_reanal/rSEM_V44/star_rsem/|.genes.results",replacement = "")
names(files)<- gsub(names(files),pattern = "_trimmed/abundance.tsv",replacement = "")
names(files)
files<-files[c(10:15,1:3)]
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
rsem.WT_rCMT3<-txi.rsem$abundance

sampleNames <- colnames(txi.rsem$counts)
sampleGroup <- factor(c(rep("TGroup",6),rep("CGroup",3)))
sampleTable <- data.frame(sampleName = sampleNames, type = sampleGroup)
rownames(sampleTable) <- colnames(txi.rsem$counts)
txi.rsem$length[txi.rsem$length == 0] <- 1
ddsTxi <- DESeqDataSetFromTximport(txi.rsem, sampleTable, design = ~ type)
dds<-ddsTxi
ddsTxi <- DESeq(ddsTxi)
res <- results(ddsTxi,independentFiltering=F)
res$GeneId<-rownames(res)
res.tpm<-as.data.frame(rsem.WT_rCMT3)%>%rownames_to_column(var="GeneId")%>%inner_join(as.data.frame(res))%>%
  full_join(gene.names)


DMRsdistnacetogenes<-read_tsv("/groups/nodine/lab/members/Ranjith/mCHG/figure5_rna/rCMT3_WT_BC_DMR.distance2Genes.tsv")

Deseq_DMR_table<- DMRsdistnacetogenes%>%full_join(res.tpm)
write_tsv(Deseq_DMR_table[,-c(9,10)],"/groups/nodine/lab/members/Ranjith/mCHG/figure5_rna/table2-Deseq_DMR_rCMT3vsWT.tsv" )

boxplot(DMRsdistnacetogenes[,c(11:13)])

clusters<-read.table("/groups/nodine/lab/members/Ranjith/mCHG/Gene_Clustering/clustering_results_pcg_1tpm_5Cs.txt")
clusters<- as.data.frame(clusters[,c(4,13)])
res.tpm.clust<-clusters%>%inner_join(as.data.frame(res))

write.table(res.tpm.clust, "/groups/nodine/lab/members/Ranjith/mCHG/figure5_rna/deseq_table.tsv",quote = F,
            sep ="\t")



resSig <- res.tpm[(!is.na(res.tpm$padj) & (res.tpm$padj < 0.010000) & (abs(res.tpm$log2FoldChange)> 1)), ]
dim(resSig)

write.table(res.tpm.clust, "/groups/nodine/lab/members/Ranjith/mCHG/figure5_rna/deg_0.01_log1.tsv",quote = F,
            sep ="\t")



with()



resSig.down<-resSig%>%filter(log2FoldChange < -1)
head(resSig.down)
dim(resSig.down)
res.down.clust<-clusters%>%inner_join(as.data.frame(resSig.down))
head(res.down.clust)
view(resSig.down)



res.down.clust.1<-res.down.clust%>%filter(mClust %in% c(1,2))

resSig.up<-resSig%>%filter(log2FoldChange > 1)
resSig.up.clust<-clusters%>%inner_join(as.data.frame(resSig.up))
table(resSig.up.clust$mClust)
dim(resSig.up)
par(mfrow=c(4,4),las=2)
boxplot(log2(resSig.down[2:10]))
boxplot(log2(resSig.up[2:10]))


pheatmap::pheatmap(log2(resSig.up[2:10]+1),scale = "row",cluster_cols = F)

##### deg features


genes.dist.dmr<- read.table("/groups/nodine/lab/members/Ranjith/mCHG/figure5_rna/genes_dist_todmrs.tsv")
colnames(genes.dist.dmr)[4]<-"GeneId"

down.genes.dist.dmr<-genes.dist.dmr%>%inner_join(as.data.frame(resSig.down), by="GeneId")%>%filter(abs(V14)<1500)
dim(down.genes.dist.dmr)
rownames(down.genes.dist.dmr)<-down.genes.dist.dmr$gene_name
vioplot(down.genes.dist.dmr$V11, down.genes.dist.dmr$V12, down.genes.dist.dmr$V13)
write.table(down.genes.dist.dmr[,c(-c(1:3,8:10,26:28,24))],"/groups/nodine/lab/members/Ranjith/mCHG/figure5_rna/down.genes.dist.dmr.tsv", quote = F, sep="\t" )
library(venn)
grid.newpage()
draw.pairwise.venn(area1 = 4733, area2 = 563, cross.area = 24, category = c("Dog People", 
                                                                         "Cat People"))


meantpm<- read.table("/groups/nodine/lab/members/Ranjith/mCHG/figure5_rna/mean_tpm_table.tsv")%>%rownames_to_column(var="GeneId")
meantpm.drg<-meantpm%>%inner_join(down.genes.dist.dmr)
rownames(meantpm.drg)<-meantpm.drg$gene_name
pheatmap::pheatmap(log2(meantpm.drg[,c(2:12)]+1), scale = "row", cluster_cols = F, cellwidth = 20, cellheight = 15)

freq.up<-c(62, 39, 17, 11 )
freq.down<-c(171,129,58,16)
freq.all<-c(7882,7847,5469,14439)

frq.cbind<-cbind(freq.up,freq.down,freq.all)
rownames(frq.cbind)<-c("clust1","clust2","clust3","clust4")

scale.frq <- as.data.frame(scale(frq.cbind , center = FALSE, scale = colSums(frq.cbind)))

par(mfrow=c(3,3),las=2)
scale.frq$log2.down<- log2(scale.frq$freq.down/scale.frq$freq.all)
barplot(as.matrix(scale.frq$log2),beside=T,ylim=c(-2,2),space = 0.15)

####


with(res.tpm, plot(log2FoldChange,-log10(padj), xlim=c(-11,11), pch=20, col="gainsboro"))

title( main="multi",xlab = "-log10(padj)" , ylab="log2 FC")

# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(as.data.frame(res.tpm), padj<.01 & log2FoldChange>=1 ), points(-log10(padj), log2FoldChange, pch=20, col="slateblue", bg= '#1c75bc',cex=1))
with(subset(as.data.frame(res.tpm), padj<.01 & log2FoldChange <=-1 ), points(-log10(padj), log2FoldChange, pch=20, col="lightcoral", bg= 'black',cex=1))
box(bty="o")
abline(h=0,lty=3)



pie(freq.down)


Volcano <-res.tpm %>%
  mutate(padj = replace(padj, padj <= 0.00000001, 0.00000001)) %>%
  mutate(log2FoldChange = replace(log2FoldChange, log2FoldChange <= -4.5, -4.5))%>%
  mutate(log2FoldChange = replace(log2FoldChange, log2FoldChange >= 4.5, 4.5))%>%
  mutate(threshold = ifelse(log2FoldChange >= 1,"A", ifelse(log2FoldChange<=-1 , "B", "C")))

V.plot <- ggplot(Volcano, aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(aes(colour = threshold), alpha=0.6, size=1.5) +
  scale_colour_manual(values = c("A"= "#113559", "C"="#e4e3e2",  "B"= "#591111"))+
  xlim(c(-5, 5)) + ylim(c(0, 8))+theme_classic()+theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))+
  geom_vline(xintercept=c(-1.5,1.5),linetype="dashed", color = "black",size=0.5)


library(ggplot2)

Deseq_DMR_table.test<-as.data.frame(Deseq_DMR_table)%>%filter(padj<.01 & log2FoldChange< -1)%>%filter(abs(Distance2gene)<1500)%>%
 filter(Hypo %in% c("WT"))

names(Deseq_DMR_table.test)[c(15:17)]<-c("rCMT3_L1_R1","rCMT3_L1_R2","rCMT3_L1_R3")
names(Deseq_DMR_table.test)[c(18:20)]<-c("rCMT3_L3_R1","rCMT3_L3_R2","rCMT3_L3_R3")
write_tsv(Deseq_DMR_table.test[-c(9,10)],"/groups/nodine/lab/members/Ranjith/mCHG/figure5_rna/downregulated.hypomethylated.tsv")


meantpm.drg<-meantpm%>%inner_join(as.data.frame(Deseq_DMR_table.test[c(4,30)]))
rownames(meantpm.drg)<-meantpm.drg$gene_name
breaksList = c(-2,seq(-1.5,1,0.075),seq(1.5,2,0.05), 2.5)
library(RColorBrewer)
pheatmap::pheatmap(log2(meantpm.drg[,c(5:12)]+1), scale = "row", cluster_cols = F, cellwidth = 20, cellheight = 15,border_color = "white",
                   color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
                   breaks = breaksList)

boxplot(Deseq_DMR_table.test[,c(11:13)])




