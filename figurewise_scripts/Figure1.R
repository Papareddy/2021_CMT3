"AT5G49160", MET1 
"AT1G57820", VIM1
"AT1G66050", VIM2
"AT5G39550", VIM3


------ 1A For bar chart of MET1, VIM123 Transcript levels ------

### Import TPM Table

tpm = read.table('~/embryo_timecourse_TPM.tsv',header=T,row.names = 1,stringsAsFactors = F)%>%
  rownames_to_column(var="GeneId")%>%filter(GeneId %in% c("AT5G49160","AT1G57820","AT1G66050","AT5G39550"))%>%
  dplyr::select(matches('GeneId|floralbud|pg|gl|eh|lh|et|lt|bc|mg|leaf'))%>%melt()
tpm$variable <- gsub("\\_.*","",tpm$variable)



ggbarplot(tpm, x = "variable", y = "value",add = "mean_sd",fill="#ffaaaa",width = 0.8, size=0.1,facet.by ="GeneId")+
theme(axis.text.x=element_text(angle = 90,hjust = 0.5),text = element_text(size=10),
      axis.ticks = element_line(colour = "black", size = 0.5),
      axis.line.x = element_line(colour = "black", size=0.5),
      panel.border = element_rect(colour = "black", fill="transparent", size=0.5),
      axis.ticks.length=unit(.15, "cm"),strip.background = element_blank(),
      panel.spacing = unit(0, "lines"), strip.text = element_blank())+
  stat_compare_means(ref.group = "leaf", label = "p.signif",method = "t.test")+facet_wrap(~GeneId)

------ 1B  Co-Varing genes ------

### Import scripts from other locations
source('~/nearest.R')
genes = c("AT5G49160","AT1G57820","AT1G66050","AT5G39550")# MET1, VIM1,VIM2 and VIM3
weight = distances(genes, TPM, weighted = TRUE)

namedweights<- as.data.frame(weight)%>%rownames_to_column(var="ensembl_gene_id")
geneorder_byweight<- inner_join(namedweights,araGenes[,c(1,6)])
geneorder_byweight<-geneorder_byweight[order(geneorder_byweight$weight),]


#LOAD DATA TABLES
tpm = read.table('~/embryo_timecourse_TPM.tsv',header=T,row.names = 1,stringsAsFactors = F)
sample_factors = read.table('~/embryo_timecourse_factors.tsv',header = T,stringsAsFactors = F)

#CALCULATE MEANS
mean_values = sapply(
  as.character(unique(sample_factors[,'type'])),
  function(i){
    picked_columns = sample_factors[sample_factors[,'type'] == i, 'name']
    return(round(rowMeans(tpm[,picked_columns]),2))
  }
)

#### chose stages and genes ###
picked_stages = c("fb",'pg','gl','eh','lh','et','lt','bc','mg','lf')
####randomization
#forSamp<-as.data.frame(mean_values[,picked_stages])
#df<- c()
#for (i in 1:1000){
 # dat <-colMeans(sample_n(forSamp, 25))
  #df <- bind_rows(df, dat)
  #}
#randomized<-colMedians(as.matrix(df))

###
picked_genes=c(genes,names(weight[1:25]))
weigted.tpm<-mean_values[picked_genes,picked_stages]
weigted.tpm<-as.data.frame(weigted.tpm)%>%rownames_to_column(var="ensembl_gene_id")

Covarient.Genes <- inner_join(weigted.tpm,araGenes[,c(1,6)],by="ensembl_gene_id")
Covarient.Genes$description<-gsub("\\[.*","",Covarient.Genes$description)
Covarient.Genes<-Covarient.Genes[order(Covarient.Genes$description),]
rownames(Covarient.Genes)<-Covarient.Genes$description
#Covarient.Genes<-Covarient.Genes[c(11,17:19,30,2,1,3:10,12:16,20:29),]
Covarient.Genes<-Covarient.Genes[c(11,17:19,2:10,12:16,1,20:29),]
Covarient.Genes<-Covarient.Genes[,-c(1,12)]
Covarient.Genes[30,]<-c(randomized)

Covarient.Genes<-Covarient.Genes[c(1:4,30,5:30),]

breaksList = c(-2, seq(-1, -0.5, by = 0.05),seq(-0.35, 2, by = 0.05))
color = colorRampPalette(c("black",rev(brewer.pal(n = 11, name ="RdYlBu"))),bias=0.95)(length(breaksList))

color = colorRampPalette(c("black","#222568",rev(brewer.pal(n = 11, name ="RdYlBu")),"#73001a","#520013"),bias=0.95)(50)
pheatmap::pheatmap(Covarient.Genes,scale = "row",cluster_cols = F,cluster_rows = F, border_color = "white",
                   cellwidth = 10,cellheight = 10,gaps_row = c(4,5,30),
                    fontsize_row = 4,color = color)


dummt=c(5,100)
basemean<-as.matrix( rbind(cbind(rowMeans(Covarient.Genes),rowMeans(Covarient.Genes)),dummt))
basemean<-basemean[c(1:4,30,5:29,31),]
breaksList = c(0,seq(3, 6, by = 0.05),seq(6.5, 8, by = 0.05))
color = colorRampPalette(c("white","#542788",rev(brewer.pal(n = 11, name ="Spectral"))[-c(3,6)]),bias=1)(length(breaksList))


#color = colorRampPalette(c("black",rev(brewer.pal(n = 11, name ="Spectral"))),bias=0.75)(10)
pheatmap::pheatmap(log2(basemean),show_rownames = T,cellwidth = 10,cellheight = 10,
                   fontsize_row = 4,gaps_row = c(4,5,30),
                   cluster_cols = F,cluster_rows = F, border_color = "grey25", color = color)
                   


------ 1C DMC  BIT MAP
DMS<- read.table("/Volumes/nodine/lab/members/Ranjith/mCHG/figure1/CGDMRs/Developmental_DMRs/allC.DMS.sorted.bed")%>%filter(V5<=0.01)
DMS<- separate(DMS,col="V4",into = c("Context", "strant"),sep="\\;")#ignore strand information

ggplot() + geom_logo(DMS$Context,method='p') +theme_classic()+theme(axis.ticks.length=unit(.15, "cm"))+
  geom_hline(yintercept=seq(0,1,0.1), lty=2, color="grey50",lwd=1)+
  scale_y_continuous(breaks=seq(0,1,0.1))+theme(legend.title = element_blank())



------ 1D DMR annotation ------

### Load Pear script fro DMR annotation interscection


CG<- read.table("/Volumes/nodine/lab/members/Ranjith/mCHG/figure1/CGDMRs/Developmental_DMRs/annotation/sig.CG.DMRs.coord.bed.Anno",stringsAsFactors=F,header=T)
CG.S1      <- sum(CG[,4])
CG    <- cbind(CG[,9],CG[,7],CG[,5],CG[,6],CG[,8],CG[,10],CG[,11:15],CG[,17],CG[,20],CG[,16],CG[,19],CG[,18],CG[,21])
colnames(CG) <- c("2kb-US","5pUTR","CDS","Intron","3pUTR","2kb-DS","aslncRNA","lncRNA","miRNA","pri-miRNA","ncRNA","snoRNA","tRNA","pseudo","TE gene","TE","Intergenic")
CG <- CG%>%select( -matches('gene|pri-miRNA|aslncRNA|ncRNA|miRNA'))
CG.D <- colSums(CG)
CG.D.rel <- CG.D/CG.S1
CG.D.rel<- as.data.frame(CG.D.rel)
CG.freq <- as.data.frame(scale(CG.D.rel , center = FALSE, scale = colSums(CG.D.rel)))


CHG<- read.table("/Volumes/nodine/lab/members/Ranjith/mCHG/figure1/CGDMRs/Developmental_DMRs/annotation/sig.CHG.DMRs.coord.bed.Anno",stringsAsFactors=F,header=T)
CHG.S1      <- sum(CHG[,4])
CHG    <- cbind(CHG[,9],CHG[,7],CHG[,5],CHG[,6],CHG[,8],CHG[,10],CHG[,11:15],CHG[,17],CHG[,20],CHG[,16],CHG[,19],CHG[,18],CHG[,21])
colnames(CHG) <- c("2kb-US","5pUTR","CDS","Intron","3pUTR","2kb-DS","aslncRNA","lncRNA","miRNA","pri-miRNA","ncRNA","snoRNA","tRNA","pseudo","TE gene","TE","Intergenic")
CHG <- CHG%>%select( -matches('gene|pri-miRNA|aslncRNA|ncRNA|miRNA'))
CHG.D <- colSums(CHG)
CHG.D.rel <- CHG.D/CHG.S1
CHG.D.rel<- as.data.frame(CHG.D.rel)
CHG.freq <- as.data.frame(scale(CHG.D.rel , center = FALSE, scale = colSums(CHG.D.rel)))

symmetric.features<-cbind(CG.freq,-CHG.freq)
symmetric.features<- symmetric.features[order(-symmetric.features[,1]),]
Tree<-barplot(symmetric.features$CHG.D.rel,ylim=c(-0.41,0.41),las=2, border = "NA",col="#673b8e",yaxt="n",xaxt="n")
axis(1,at=Tree,labels = row.names(symmetric.features),las=2)
barplot(symmetric.features$CG.D.rel,add=T, names.arg = row.names(symmetric.features),col = "#a70f3d",border = NA,xaxt="n",las=2)


------ 1E and 1F DMR annotation ------

overlaps<-read.table("~/Desktop/Manustrips/mCHG/figure1/CHGin_Overlaps.txt",header = T,stringsAsFactors = T)
vioplot(overlaps$methylation_level_MET1-overlaps$methylation_level_WT,
      overlaps$methylation_level_SJ.VIM-overlaps$methylation_level_SJ.WT,
        overlaps$methylation_level_CMT3-overlaps$methylation_level_WT,
        main="overlaped",las=2,outline=T,notch=T,names=c("MET1","VIM123","CMT3"),
      col="lightcoral",lwd=1,boxcol="grey50",whiskcolr="grey50",outcol="grey40",staplecol="grey50")

nonoverlaps<-read.table("~/Desktop/Manustrips/mCHG/figure1/CHGinNonOverlaps.txt",header = T,stringsAsFactors = T)
vioplot(nonoverlaps$methylation_level_MET1-nonoverlaps$methylation_level_WT,
        nonoverlaps$methylation_level_SJ.VIM-nonoverlaps$methylation_level_SJ.WT,
        nonoverlaps$methylation_level_CMT3-nonoverlaps$methylation_level_WT,
        main="Non_overlaped",las=2,outline=T,notch=T,names=c("MET1","VIM123","CMT3"),
        col="#9088d4",lwd=1,boxcol="grey50",whiskcolr="grey50",outcol="grey40",staplecol="grey50")



------ 1G and 1H  DMRs ------

col=c("#8091c1","#7766cc","#f785b2","#e0b45c","#f7a59a","#9f8a76","#9ad485","#2b8977")
CG.DMRs<- read.delim("/Volumes/nodine/lab/members/Ranjith/mCHG/figure1/CGDMRs/Developmental_DMRs/500sim.100bp.CG_Main_Dev_series_rms_results_collapsed.tsv",header=T)
          %>%filter(number_of_dms>=8)
dim(CG.DMRs)
boxplot(CG.DMRs[,c(7:14)],outline=F,col=col,medcol="white",border="grey15",lty=1,lwd=1.5,medlwd=3,range=1.5)


        

CHG.DMRs<- read.delim("/Volumes/nodine/lab/members/Ranjith/mCHG/figure1/CGDMRs/Developmental_DMRs/100sim.100bp.CHG_Main_Dev_series_rms_results_collapsed.tsv",header=T)
          %>%filter(number_of_dms>=4)
dim(CHG.DMRs)
boxplot(CHG.DMRs[,c(7:14)],outline=F,col=col,medcol="white",border="grey15",lty=1,lwd=1.5,medlwd=3,range=1.5)



------ S 1A GO Analysis ------
S1A ### GO Analysis
int.genes <- names(weight[1:25])
int.genes <- factor(as.integer(all.genes %in% int.genes))
names(int.genes) = all.genes
GOdata <- new("topGOdata", ontology='BP', allGenes = int.genes, annot = annFUN.gene2GO, gene2GO = geneID2GO)
results <- runTest(GOdata, algorithm = "elim", statistic = "fisher") ##running the GO Test

goEnrichment <- GenTable(object = GOdata, elimFisher = results,topNodes=15)
goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
write.xlsx(goEnrichment,"~/Desktop/Manustrips/mCHG/figure1/GoTerms_withhCGmMTas.xlsx")
goEnrichment$Term <- paste(goEnrichment$Term,goEnrichment$GO.ID,  sep="-")

goEnrichment$log10<- -log10(as.numeric(goEnrichment$elimFisher))


------ S 1B Methylation rate and methyltransferases TPM correlation ------
CGmean<- as.data.frame(colMeans(CG.DMRs[,c(7:11,13,14)], na.rm = T, dims = 1))
CGmean<-CGmean*100
cormet<-as.data.frame(t(Covarient.Genes[1,c(1,2,4,5,8,9,10)])) ## cormat means correlation matrix
CGm_MET<-cbind(cormet,CGmean)
names(CGm_MET)<-c("met1","mCG")
CGm_MET<-CGm_MET%>%rownames_to_column(var = "stage")

CGm_MET<-ggscatter(CGm_MET, x = "mCG", y = "met1",add = "reg.line", color = "stage",size =3.5,add.params = list(color = "black", fill = "lightgray"))+
  stat_cor(method = "pearson")+xlim(30, 60)+geom_hline(yintercept=1, lty=2, color="grey50")

HGmean<- as.data.frame(colMeans(CHG.DMRs[,c(7:11,13,14)], na.rm = T, dims = 1))
CHGmean<-CHGmean*100
cormet<-as.data.frame(t(Covarient.Genes[16,c(1,2,4,5,8,9,10)]))

CHGm_MET<-cbind(cormet,CHGmean)
names(CHGm_MET)<-c("CMT3","mCHG")
CHGm_MET<-CHGm_MET%>%rownames_to_column(var = "stage")

CHGm_MET<-ggscatter(CHGm_MET, x = "mCHG", y = "CMT3",add = "reg.line", color = "stage",fill="stage",size =3.5,add.params = list(color = "black", fill = "lightgray"))+
  stat_cor(method = "pearson")+scale_y_continuous(breaks=seq(0,60,30))+geom_hline(yintercept=1, lty=2, color="grey50")



