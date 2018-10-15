#! /usr/bin/Rscript

                                        # This script is to plot bar plot of GO enrichment

library(ggplot2)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
go <- read.delim(args[1], header = TRUE,
                 colClasses = c("character"))

go2 <- go[,c(2,4,6,8)]

go2 <- data.frame(GOterms=paste(go2[,1]," (",go2[,2],")",sep=""),
                  Pvalue=go2[,3],
                  Aspect=go2[,4], stringsAsFactors = FALSE)

go2$Pvalue <- as.numeric(go2$Pvalue)
go2$Pvalue <- -log10(go2$Pvalue)
#head(go2)
go2$Aspect <- factor(go2$Aspect, levels = c("P", "F", "C"))
#head(go2)
#str(go2)
go2 <- arrange(go2, Aspect, Pvalue)

go2$GOterms <- as.character(go2$GOterms)
#print(unique(go2$GOterm))
go2$GOterms <- factor(go2$GOterms, levels = unique(go2$GOterms))
p <- ggplot(go2, aes(x=GOterms, y=Pvalue))
p <- p+geom_bar(stat="identity", width = 0.75, fill="#2ca25f",
                color="black", size=0.4, alpha = 0.85)
p <- p+facet_grid(Aspect~., scales = "free", space="free")
p <- p+ylab("-log10(Pvalue)")
p <- p+xlab("GO terms")
p <- p+geom_hline(yintercept=1.3, linetype = "dashed",
                  alpha = 0.8, size = 0.4)
p <- p+coord_flip()

p <- p+theme_classic()
p <- p+theme(axis.text.y = element_text(face = "bold"))
#pdf("tmp.pdf")
#p
#dev.off()
ggsave(filename = args[2], width = as.numeric(args[3]), height = as.numeric(args[4]))
