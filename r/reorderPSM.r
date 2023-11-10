#!/usr/bin/env Rscript
args = commandArgs( trailingOnly = TRUE )

library( rhdf5 )

fileName = args[1]
psm <- h5read( fileName, "posteriorSimilarityMatrix")

n = nrow( psm )
distance = as.dist( 640 - psm )

hc = hclust( distance, method = "complete", members = NULL)
psm_hc = psm

psm_hc[ 1:n, ] = psm_hc[ hc$order, ]
psm_hc[ , 1:n] = psm_hc[ , hc$order]

saveFileName = paste( fileName, "_reoder.txt", sep = "")
write.table( psm_hc, file = saveFileName, row.names = FALSE, col.names = FALSE)
