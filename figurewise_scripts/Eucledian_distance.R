euc.dist=function(A, B, W=rep(1,length(A))){
  # Calculate a Euclidean distance between array A and array B.
  # Optional W is an array of weights, where larger weights have stronger effects
  return(sqrt(sum((((A - B)*W) ^ 2))))
}

zscore = function(v){
  (v - mean(v))/sd(v)
}

TPM = read.table("~/Desktop/Manustrips/mCHG/figure1/covariance/embryo_mean_tpm.tsv",stringsAsFactors = F,header = T,row.names = 1)
TPM = TPM[rowSums(TPM > 1)>=1,] # Subset for genes that are expressed > 1TPM in >= 1 stage of embryogenesis
zscores = round(t(apply(TPM,1,zscore)),2)
TPM = round(TPM,2)

distances = function(gene, dataset=zscores, weighted=FALSE, include_input=FALSE){
  # Returns a sorted array of genes by their distance from the input gene(s)
  gene_subset = gene[gene %in% rownames(dataset)] # Subset the list of input genes for those that exist in the dataset
  weights = rep(1,length(gene_subset)) # Default weight of 1 for every point
  
  if(length(gene_subset)>1){ # More than 1 gene was provided. Calculate the centroid
    search_values = colMeans(dataset[gene_subset,])
    if(weighted){ # Calculate the inverse of the standard deviation for each dimension in the input data
      weights = (apply(dataset[gene_subset,], 2, sd))^-1
    }
  }else{
    if(length(gene_subset)==0){ # No genes were found in the dataset
      return()
    }
    # Exactly 1 gene makes up the search values
    search_values = dataset[gene_subset,]
  }
  
  dists = apply(dataset, 1, function(gene_values){
    return(euc.dist(as.numeric(search_values),as.numeric(gene_values),as.numeric(weights)))
  })
  
  if(!include_input){
    dists = dists[!names(dists) %in% gene]
  }
  return(sort(dists))
}



