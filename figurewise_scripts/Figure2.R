
-------- 2c, 2D, 2G and 2H ---------------

CHG<- list.files("/Volumes/nodine/lab/members/Ranjith/mCHG/P4effect/",pattern = "CHG.TEmapped",full.names=T)
CHG = as.data.frame(do.call(cbind,lapply(CHG , function(x) {                 
  DF <- read.table(x, stringsAsFactors = T)
  DF<-DF[,c(4,6,7)]
  colnames(DF)<- c(paste0(basename(x),'.GeneID'),paste0(basename(x),'.count'),paste0(basename(x),'.meth'))
  return(DF)})))
head(CHG) 

CHG1<- data.frame(CHG[,-1],row.names = CHG[,1])
CHG1<-CHG1%>%dplyr::select(-matches('GeneID'))%>%rownames_to_column(var="V4")%>%
  inner_join(TE.clasif.clusterId[,c(4,6)])%>%
  filter(Clusters!="Zero")%>%filter_at(vars(ends_with(".count")), all_vars(. > 4))%>%dplyr::select(-matches('count'))

CHG1[,-c(1,10)]<- CHG1[,-c(1,10)]*100

CHG.A<-CHG1%>%filter(Clusters %in% c("7","1","2"))%>%dplyr::select(-matches('Clusters'))
CHG.B<-CHG1%>%filter(Clusters %in% c("5","6","8","3"))%>%dplyr::select(-matches('Clusters'))


par(mfrow=c(4,4),las=2)
boxplot(CHG.A$FB.P4.4x.CHG-CHG.A$FB.WT.4x.CHG,
        CHG.A$Heart.P4.4x.CHG-CHG.A$Heart.WT.4x.CHG,
        CHG.A$BC.P4.4x.CHG-CHG.A$BC.WT.4x.CHG,
        CHG.A$Leaf.P4.4x.CHG-CHG.A$RP.Leaf.WT.4x.CHG,
        outline=F,main="Euchrom_WTvsP4_CHG",col=col,medcol="white",border="grey15",lty=1,lwd=1.5,medlwd=3,range=1.5)

boxplot(CHG.B$FB.P4.4x.CHG-CHG.B$FB.WT.4x.CHG,
        CHG.B$Heart.P4.4x.CHG-CHG.B$Heart.WT.4x.CHG,
        CHG.B$BC.P4.4x.CHG-CHG.B$BC.WT.4x.CHG,
        CHG.B$Leaf.P4.4x.CHG-CHG.B$RP.Leaf.WT.4x.CHG,
        outline=F,main="Hetchrom_WTvsP4_CHG",col=col,medcol="white",border="grey15",lty=1,lwd=1.5,medlwd=3,range=1.5)
        
        
        
#### same for CG TEmapedd files in the same location


