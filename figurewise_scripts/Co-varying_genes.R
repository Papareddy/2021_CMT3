source('Eucledian_distance.R')

genes = c('AT4G19020','AT5G14620')
no_weight = distances(genes, TPM)
weight = distances(genes, TPM, weighted = TRUE)
z_no_weight = distances(genes, zscores)
z_weight = distances(genes, zscores, weighted = T)
