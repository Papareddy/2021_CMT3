ibrary(tximport)
library(DESeq2)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(tm)
library(factoextra)


#### LOAD NODINE DATA
rsemcounts.MD.list<- list.files("/scratch-cbe/users/ranjith.papareddy/rCMT3_transcriptome_reanal/nodine_embryonic/results/rsem/",
                                      pattern = "gene_counts.csv",full.names=T)

rsemcounts.MD.list<-rsemcounts.MD.list[-28]

rsemcounts.MN = as.data.frame(do.call(cbind,lapply(rsemcounts.MD.list , function(x) {
  DF <- read.csv(x)
  colnames(DF)<- c(paste0(basename(x),'.GeneID'),paste0(basename(x),'count'))
  return(DF)})))

colnames(rsemcounts.MN)<-gsub(colnames(rsemcounts.MN),pattern = "_rsem_gene_counts.csv|_|count",replacement = "")
colnames(rsemcounts.MN)[1]<-"GeneId"
rsemcounts.MN<-rsemcounts.MN%>%dplyr::select(-matches(".GeneID"))
factor_table.MD = data.frame('name'=colnames(rsemcounts.MN[-1]),'genotype'=removeNumbers(colnames(rsemcounts.MN[-1])))

#### LOAD RP DATA
rsemcounts.RP.list<- list.files("/scratch-cbe/users/ranjith.papareddy/rCMT3_transcriptome_reanal/results/rsem/",
                                  pattern = "gene_counts.csv",full.names=T)

rsemcounts.RP.list<-rsemcounts.RP.list[-70]


rsemcounts.RP = as.data.frame(do.call(cbind,lapply(rsemcounts.RP.list , function(x) {
  DF <- read.csv(x)
  colnames(DF)<- c(paste0(basename(x),'.GeneID'),paste0(basename(x),'count'))
  return(DF)})))

colnames(rsemcounts.RP)<-gsub(colnames(rsemcounts.RP),pattern = "_rsem_gene_counts.csv|_|count",replacement = "")
colnames(rsemcounts.RP)[1]<-"GeneId"
rsemcounts.RP<-rsemcounts.RP%>%dplyr::select(-matches(".GeneID"))%>%dplyr::select(matches("GeneId|BAMpooled"))
colnames(rsemcounts.RP)<-gsub(colnames(rsemcounts.RP),pattern = ".BAMpooled",replacement = "")
factor_table.Rp = data.frame('name'=colnames(rsemcounts.RP[-1]),'genotype'=c(rep('RP.bc',3),rep('gCMT3',6),rep('rCMT3',6)))



full.counts.withplastids<-rsemcounts.RP%>%inner_join(rsemcounts.MN)
### filter pcg
pc.genes <- scan('/groups/nodine/lab/members/Ranjith/annotations_ranj/gene_types/protein_coding.txt','character')
mc.genes = scan('/groups/nodine/lab/members/Ranjith/annotations_ranj/gene_types/mitochondrial_chloroplast.txt','character')
nm.genes = scan('/groups/nodine/lab/members/Ranjith/annotations_ranj/gene_types/nuclear_mitochondria.txt','character')
ercc.genes = scan('/groups/nodine/lab/members/Ranjith/annotations_ranj/gene_types/ercc_spike_ins.txt','character')
protein.coding.genes <- pc.genes[!pc.genes %in% c(mc.genes,nm.genes)]
protein.coding.genes <- as.data.frame(protein.coding.genes[protein.coding.genes %in% full.counts.withplastids$GeneId])
colnames(protein.coding.genes)[1]<-"GeneId"
full.counts =  as.data.frame(protein.coding.genes%>%inner_join(full.counts.withplastids))
full.counts= full.counts%>%column_to_rownames(var="GeneId")

## VST normalization
factor_table.ful<-rbind(factor_table.Rp,factor_table.MD)
counts = round(full.counts[-1],0)
dm = DESeqDataSetFromMatrix(counts,colData = factor_table.ful, design = ~ genotype)
vst = varianceStabilizingTransformation(dm,blind=T)
v = assay(vst)
vst_for_pca = apply(v,2,function(x){
  new_x = x - 5
  new_x[new_x < 0] = 0
  # new_x = new_x / max(new_x)
  return(new_x)
})
vst_for_pca = round(vst_for_pca,2)



##### ploting
pheatmap(cor(vst_for_pca),border_color = NA)

prcomp.pca<-prcomp(vst_for_pca[,-c(4:9)],scale. = F)
## PC1 89.411105007
#PC2  6.357028988
prcomp.pca<-prcomp.pca$rotation
prcomp.pca<-prcomp.pca[order(rownames(prcomp.pca)),]
prcomp.pca<-as.data.frame(prcomp.pca)
prcomp.pca$tissue<-rownames(prcomp.pca)

par(mfrow=c(4,4),las=2,bty="n")
dev.off()
plot(prcomp.pca[,1]~prcomp.pca[,2],data=prcomp.pca,col=color1,pch=pchs, cex=1.25)
text(prcomp.pca[,1]~prcomp.pca[,2],data=prcomp.pca, labels=rownames(prcomp.pca))

prcomp.pca$color<-c(rep("1ng",6),rep("RP.bc",3),rep("eh",3),rep("et",3),
                    rep("gl",3),rep("lh",3),rep("lt",3), rep("mg",3),rep("pg",3),rep("R5",3),rep("R12",3))

prcomp.pca$forshape<-c(rep(8,9),rep(19,21),rep(8,6))

ggplot(prcomp.pca, aes(PC2,PC1, color=color,group=forshape))+geom_point(size=4,aes(shape=forshape))+
  scale_shape_manual(values=c(17, 16))+
scale_color_manual(values=c('#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','grey25','#b15928'))+
  theme_classic()



color1<-c(rep("#a6cee3",6),rep("#1f78b4",3),rep("#b2df8a",3),rep("#fb9a99",3),
         rep("#e31a1c",3),rep("#fdbf6f",3),rep("#ff7f00",3), rep("#cab2d6",3),rep("#6a3d9a",3),rep("#d9d9d9",3),rep("#636363",3))

pchs<-c(rep(8,9),rep(19,21),rep(8,6))




##### DEGS remove CA and R12C

degs.rsemRP<-rsemcounts.RP
rownames(degs.rsemRP)<-degs.rsemRP$GeneId
degs.rsemRP<-degs.rsemRP[,c(2:4,11:13)]
head(degs.rsemRP)
groups <- factor(c(rep("CGroup",3),rep("TGroup",3)))

data<-round(degs.rsemRP)
min_read <- 5
data <- data[apply(data,1,function(x){max(x)}) > min_read,]
sampleInfo <- data.frame(groups,row.names=colnames(data))
library(DESeq2, quietly=T)
dds <- DESeqDataSetFromMatrix(countData = data, colData = sampleInfo, design = ~ groups)
dds$condition = relevel(dds$groups,"CGroup")
dds <- DESeq(dds)
res <- results(dds)




res <- results(dds,independentFiltering=F)
countdesq=counts(dds,normalize=T)
#write.table(res, file="multi.P4.WT_gene_TE_analysis.txt", sep="\t",quote=F)
resSig <- res[(!is.na(res$padj) & (res$padj < 0.010000) & (abs(res$log2FoldChange)> 2)), ]
#write.table(resSig, file="multi.P4.WT_sigdiff_gene_TE.txt",sep="\t", quote=F)
dim(resSig)
with(resSig, plot(log2FoldChange,-log10(padj)))
abline(v=c(2,-2),col="red")

resSig.tpm<-as.data.frame(resSig)%>%rownames_to_column(var="GeneId")%>%inner_join(rsemcounts.RP)

boxplot(log2(resSig.tpm[,c(8:22)]))
