#! /usr/bin/Rscript

library(pheatmap)
args <-  commandArgs(trailingOnly = T)
a <- read.table(args[1], header = T, row.names =1)
a <- 1-a # Plot 1-PPAU
#head(a)
a <- apply(a, 1, function(x) x=x-x[1])
a <- t(a)
na_count <- apply(a, 1, function(x) sum(is.na(x)))
a <- a[na_count < 3,]
print(paste(sum(na_count > 3),"record(s) were removed"))

pheatmap(a, color = colorRampPalette(c("blue","white","red"))(11),
         breaks = seq(from = -1, to=1, length.out = 12),
         border_color = NA,
         show_rownames = F, cluster_cols = F,
         filename = args[2], width = 4, height = 8)
