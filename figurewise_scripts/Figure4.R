library(tidyverse)
library(dplyr)
library(pkgsearch)
library(factoextra)
library(NbClust)
--- Identify K clusters -------
clusters.R12methrate<- read.table("/groups/nodine/lab/members/Ranjith/mCHG/BS_all/pooledTsvs/BC_R12.CHG.bed.allClusters")
clusters.WTmethrate<- read.table("/groups/nodine/lab/members/Ranjith/mCHG/BS_all/pooledTsvs/BC_WT.CHG.bed.allClusters")

clusters.R12_WT.methrate<-as.data.frame(clusters.R12methrate[,c(4:7)])%>%inner_join(clusters.WTmethrate[,c(4:7)],by="V4")
clusters.R12_WT.methrate<-clusters.R12_WT.methrate%>%filter()

median.diff<-as.data.frame(aggregate((V6.x-V6.y)~V5.x,clusters.R12_WT.methrate,median))
par(mfrow=c(3,3),las=2, tcl=-0.35)
mp<-barplot(median.diff[,2]*100,border=NA,names.arg=median.diff[,1],xaxt="n", ylim=c(0,17))
axis(1, at=mp, labels=c(1:4), las=2)


#pc.genes <- read.table('/groups/nodine/lab/members/Ranjith/annotations_ranj/gene_types/protein_coding.txt')
#mc.genes = read.table('/groups/nodine/lab/members/Ranjith/annotations_ranj/gene_types/mitochondrial_chloroplast.txt')
#nm.genes = read.table('/groups/nodine/lab/members/Ranjith/annotations_ranj/gene_types/nuclear_mitochondria.txt')
#ercc.genes = read.table('/groups/nodine/lab/members/Ranjith/annotations_ranj/gene_types/ercc_spike_ins.txt')
#protein.coding.genes <- pc.genes[!pc.genes %in% c(mc.genes,nm.genes)]
#names(protein.coding.genes)<-"V4"


#genes.bed<-read.table("/groups/nodine/lab/members/Ranjith/mCHG/BS_all/pooledTsvs/nonAth_TAIR10_protein_coding.bed.ID")
#genes.bed<-as.data.frame(genes.bed)
#PCG.genes.bed<-as.data.frame(protein.coding.genes)%>%inner_join(genes.bed,by="V4")
#PCG.genes.bed<-PCG.genes.bed[,c(2:4,1,5:7)]
#write.table(PCG.genes.bed,"/groups/nodine/lab/members/Ranjith/mCHG/BS_all/pooledTsvs/nonAth_TAIR10_protein_coding.bed",quote=F)

###3 Cluster
siRNA.TE<- list.files("/groups/nodine/lab/members/Ranjith/mCHG/BS_all/pooledTsvs",pattern = "CHG.Genemapped",
                      full.names=T)
siRNA.TE<-siRNA.TE[c(3,6,7,10,13,14)]
siRNA.TEmap = as.data.frame(do.call(cbind,lapply(siRNA.TE , function(x) {
  DF <- read.table(x)
  colnames(DF)<- c(paste0(basename(x),'.chr'),paste0(basename(x),'.start'),
                   paste0(basename(x),'.end'),paste0(basename(x),'.AGI'),
                   paste0(basename(x),'.V5'),paste0(basename(x),'.Strand'),paste0(basename(x),'.V6'),
                   paste0(basename(x),'.count'),paste0(basename(x),'.methRate'))
  return(DF)})))

names(siRNA.TEmap)<-gsub(x = names(siRNA.TEmap), pattern = ".CHG.Genemapped", replacement = "")
siRNA.TEmap1<-siRNA.TEmap%>%dplyr::select(matches("BC_Col_old.chr|BC_Col_old.start|BC_Col_old.end|BC_Col_old.AGI|BC_Col_old.V5|BC_Col_old.Strand|BC_Col_old.V6|methRate"))
names(siRNA.TEmap1)

siRNA.TEmap1$r12.wt.diff.emb<-siRNA.TEmap1$BC_R12.methRate-siRNA.TEmap1$BC_Col_old.methRate
mclusted.sirna<- Mclust(siRNA.TEmap1$r12.wt.diff.emb,G=4, modelNames="V")

siRNA.TEmap1 <- cbind(siRNA.TEmap1,mclusted.sirna$classification)
names(siRNA.TEmap1)
colnames(siRNA.TEmap1)[15] <-"mClust"
table(siRNA.TEmap1$mClust)
boxplot(siRNA.TEmap1$r12.wt.diff.emb~siRNA.TEmap1$mClust,outline=F)
boxplot((siRNA.TEmap1$BC_Col_old.end-siRNA.TEmap1$BC_Col_old.start)~siRNA.TEmap1$mClust,outline=F)
#siRNA.TEmap1$ntiles<- ntile(siRNA.TEmap1$r12.wt.diff.emb, 5)
#boxplot(siRNA.TEmap1$r12.wt.diff.emb~siRNA.TEmap1$ntiles,outline=F)

#km.res <- kmeans(siRNA.TEmap1$r12.wt.diff.emb, 4, nstart = 25)
#siRNA.TEmap1$clusters<-km.res$cluster
#table(siRNA.TEmap1$clusters)
#boxplot(siRNA.TEmap1$r12.wt.diff.emb~siRNA.TEmap1$clusters,outline=F)
#set.seed(123)

siRNA.TEmap1$TAIR10<-"TAIR10"
siRNA.TEmap1$gene<-"gene"
forgtfs<-siRNA.TEmap1[,c(1,16,17,2:3,5:6,15,4)]
head(forgtfs)
#write.table(forgtfs,"/groups/nodine/lab/members/Ranjith/mCHG/main.allclusters.txt",quote=F)


library(mclust)
BIC.mCHGdiff <- mclustBIC(te.24nt.mclust.mean1,G=seq(2,20),by=2)
#plot(BIC.siRNA,las=2)
#abline(v=8)

siRNA.TEmap1$TAIR10<-"TAIR10"
siRNA.TEmap1$gene<-"gene"
forgtfs<-siRNA.TEmap1[,c(1,16,17,2:3,5:6,18,4)]
