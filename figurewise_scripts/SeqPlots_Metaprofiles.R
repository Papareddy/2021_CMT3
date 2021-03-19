library(seqplots)

library("BSgenome.Athaliana.TAIR.TAIR9")
library(Biostrings)
library(IRanges)
library(GenomeInfoDb)
library(rtracklayer)
bw <- list.files("~/",pattern = ".bw", full.names=T)

RdDMfeatures<- "/groups/nodine/lab/members/Ranjith/RdDM_TE.bed.gtf"
CMT2features<- "/groups/nodine/lab/members/Ranjith/CMT2_TE.bed.gtf"
genes<-"/groups/nodine/lab/members/Ranjith/nonAth_TAIR10_protein_coding.gtf"

Genes.Clusters <-list.files("/groups/nodine/lab/members/Ranjith/mCHG/Gene_Clustering/",pattern = ".gtf",full.names=T)
genes.meta <- getPlotSetArray(bw, Genes.Clusters, 'TAIR9',xmin = 1500L, xmax = 1500L,xanchored = 4000L, type = "af",bin = 100L)
plot(genes.meta,error.estimates = T,legend=T)
