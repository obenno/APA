#! /usr/bin/Rscript

library(ggplot2)
library(gridExtra)

args <- commandArgs(trailingOnly = T)
a <- read.table(args[1], header=T)
a$Diff_PPAU <- 0-a$Diff_PPAU

a$CountLabel <- rep(0, nrow(a)) # make CountLabel
a$CountLabel[a$Diff_PPAU >=0.2 & a$DE_Status=="Up"] <- "lengthening_up"
a$CountLabel[a$Diff_PPAU >=0.2 & a$DE_Status=="Down"] <- "lengthening_down"
#a$CountLabel[a$Diff_PPAU <=-0.2 & a$DE_Status=="Up"] <- "shortening_up"
#a$CountLabel[a$Diff_PPAU <=-0.2 & a$DE_Status=="Down"] <- "shortening_down"

mini_right <- ggplotGrob(ggplot(a, aes(CountLabel))+geom_bar()+theme(axis.text = element_blank(), axis.title = element_blank()))
a$Comparison <- factor(a$Comparison, levels = c("1hvs0h","2hvs0h","4hvs0h","24hvs0h","48hvs0h"))
a$DE_Status <- factor(a$DE_Status, levels = c("Down", "Up", "Not_significant"))
p <- ggplot(a, aes(Diff_PPAU, log2FC, color=DE_Status))
p <- p+geom_point(alpha=0.6)+facet_grid(Tissue~Comparison)
p <- p+geom_vline(xintercept = c(-0.2, 0.2), linetype="dashed", size=0.5, alpha=0.5)
p <- p+xlim(c(-1,1))
p <- p+xlab(expression("Time Point -"~Delta~"PPAU"))+ylab("Log2 Fold Change")
p <- p+scale_color_manual(values=c("#0571b0","#ca0020","grey"))
p <- p+theme_bw()+theme(legend.position = "bottom")
p <- p+annotation_custom(
           grob = mini_right,
           xmin = -0.75,
           xmax = 1,
           ymin = -2.5,
           ymax = 7.5)
ggsave(p, filename = args[2], width=12, height=8)
