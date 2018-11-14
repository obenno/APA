#! /usr/bin/Rscript

library(ggplot2)
library(gridExtra)
library(dplyr)

args <- commandArgs(trailingOnly = T)
a <- read.table(args[1], header=T)
a$Diff_PPAU <- 0-a$Diff_PPAU
nrow(a)
a <- a[a$Diff_PPAU >=0.2 | a$Diff_PPAU <=-0.2,] # Only plot |Diff_PPAU|>=0.2
nrow(a)

a$CountLabel <- rep(0, nrow(a)) # make CountLabel
a$CountLabel[a$Diff_PPAU > 0 & a$DE_Status=="Up"] <- "lengthening_up"
a$CountLabel[a$Diff_PPAU > 0 & a$DE_Status=="Down"] <- "lengthening_down"
a$CountLabel[a$Diff_PPAU < 0 & a$DE_Status=="Up"] <- "shortening_up"
a$CountLabel[a$Diff_PPAU < 0 & a$DE_Status=="Down"] <- "shortening_down"

a$Comparison <- factor(a$Comparison, levels = c("1hvs0h","2hvs0h","4hvs0h","24hvs0h","48hvs0h"))
a$DE_Status <- factor(a$DE_Status, levels = c("Down", "Up", "Not_significant"))

make_miniPlot <- function(d){
    mini <- ggplot(d, aes(x=condition, y=counts))
    mini <- mini+geom_bar(stat = "identity", width = 0.8, fill=c("#0571b0","#ca0020"))
    mini <- mini+geom_text(aes(label = paste(counts, "\n", condition)),
                           vjust = 0.9, size=2, lineheight = 0.8)
    mini <- mini+theme_bw()
    mini <- mini+theme(axis.text = element_blank(), axis.title = element_blank(),
                       axis.ticks = element_blank(), panel.grid = element_blank(),
                       panel.border = element_blank(),
                       plot.margin = margin(0,0,0,0))
    mini <- ggplotGrob(mini)
    return(mini)
}


for (i in unique(a$Tissue)){
    for (k in unique(a$Comparison)){
        b <- data.frame()
        b <- a[a$Tissue == i & a$Comparison == k,]
        c <- data.frame()
        c <- as.data.frame(b %>% count(CountLabel))
        
        d_right <- data.frame(condition=c("up", "down"),
                              counts=c(c[c$CountLabel=="lengthening_up",2], c[c$CountLabel=="lengthening_down",2]))
        #d_right$condition <- factor(d_right$condition, levels = c("up", "down"))
        #print(d_right)
        #str(d_right)
        d_left <- data.frame(condition=c("up", "down"),
                             counts=c(c[c$CountLabel=="shortening_up",2], c[c$CountLabel=="shortening_down",2]))
        #d_left$condition <- factor(d_left$condition, levels = c("up", "down"))
        p <- ggplot(b, aes(Diff_PPAU, log2FC, color=DE_Status))
        p <- p+geom_point(alpha=0.6)
                                        #p <- p+facet_grid(Tissue~Comparison)
        p <- p+geom_vline(xintercept = c(0), linetype="dashed", size=0.5, alpha=0.5)
        p <- p+xlim(c(-1,1))
        p <- p+xlab(expression("Time Point -"~Delta~"PPAU"))+ylab("Log2 Fold Change")
        p <- p+scale_color_manual(values=c("#0571b0","#ca0020","grey"))
        p <- p+theme_bw()+theme(legend.position = "bottom")
        #print(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range[2])
                                        # make mini plot object
        mini_left <- make_miniPlot(d_left)
        mini_right <- make_miniPlot(d_right)
        mini_y_min <- ggplot_build(p)$layout$panel_scales_y[[1]]$range$range[2]-(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range[2]-ggplot_build(p)$layout$panel_scales_y[[1]]$range$range[1])/3
        mini_y_max <- ggplot_build(p)$layout$panel_scales_y[[1]]$range$range[2]
        
        p <- p + annotation_custom(
                     grob = mini_right,
                     xmin = 0.6,
                     xmax = 1,
                     ymin = mini_y_min,
                     ymax = mini_y_max)
        p <- p + annotation_custom(
                     grob = mini_left,
                     xmin = -1,
                     xmax = -0.6,
                     ymin = mini_y_min,
                     ymax = mini_y_max)
        
        outName <- paste(i,k,"pdf", sep = ".")
        Pvalue <- fisher.test(cbind(d_left$counts, d_right$counts))$p.value
        Pvalue <- signif(Pvalue, digits = 2)
        p <- p + annotate("text",
                          label = paste("P-value =", Pvalue, "\n", "Fisher's Exact Test"),
                          x = 0.8, y = ggplot_build(p)$layout$panel_scales_y[[1]]$range$range[1],
                          size = 2)
        ggsave(plot=p, filename =outName , width=4, height=4)
    }
}
