#! /usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))

# Make Options ##

args <- commandArgs(trailingOnly = T)


option.list <- list(
    make_option(c("-r", "--ref"), type="character", default=NULL,
                help="label of gene tracks, could be a comma separated list [%default]"),
    make_option(c("-g", "--gff"), type="character", default=NULL,
                help="Input gene model (GFF File), or TxDB file is also supported [%default]"),
    make_option(c("-b", "--bam"), type="character", default=NULL,
                help="Input mapping bam file, could be a comma separated list [%default]"),
    make_option(c("-m", "--mode"), type="character", default="exp",
                help="Select illustration mode, exp or sashimi [%default]"),
    make_option(c("-c", "--chr"), type="character", default=NULL,
                help="Chromosome"),
    make_option(c("-s", "--start"), type="numeric", default=NULL,
                help="start position"),
    make_option(c("-e", "--end"), type="numeric", default=NULL,
                help="end position"),
    make_option(c("-o", "--out"), type="character", default=NULL,
                help="output pdf file"),
    make_option(c("--width"), type="numeric", default=7,
                help="output file width [%default]"),
    make_option(c("--height"), type="numeric", default=5,
                help="output file height [%default]"),
    make_option(c("--label"), type="character", default=NULL,
                help="label of samples, could be a comma separated list [%default]"),
    make_option(c("--log"), action="store_true", default=FALSE,
                help="Whether log-transform coverage"),
    make_option(c("--size"), type="character", default=NULL,
                help="Whether manually input track size"),
    make_option(c("--peptide"), type="character", default=NULL,
                help="provide peptide position file for visualization"),
    make_option(c("--PAS"), type="character", default=NULL,
                help="provide PAS bed file for visualization"),
    make_option(c("--PSI"), type="character", default=NULL,
                help="provide PSI data for visualization"),
    make_option(c("--Trans"), type="character", default=NULL,
                help="provide transcript ID to extract PSI data"),
    make_option(c("--Tissue"), type="character", default=NULL,
                help="provide tissue to extract PSI data"),
    make_option(c("--color"), type="character", default=NULL,
                help="set color of gene tracks, could be a comma separated list [%default]")
)
# It's sad that alignment track doesn't suppport wig or bw file, only data track
# support, but data track doesn't support coverage style plot 

desc <- "Plot genome and RNA-seq mapping track"
parser <- OptionParser(option_list=option.list,
                      description = desc,
                      usage="usage: %prog -g gff.input -b mapping.bam -c chromosome_name -s start_position -e end_position -o out.pdf\n")
opt <- parse_args(parser, args=args, positional_arguments=T)

#opt
#length(opt$args)

if(length(opt$options) == 0)
    stop(print_help(parser))


suppressPackageStartupMessages(library(Gviz))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggplotify))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(tidyverse))
options(ucscChromosomeNames=FALSE) # Soybean chromosome name does not stick to ucsc nameing convention

## Modify scheme
customScheme <- getScheme()
customScheme$GeneRegionTrack$fill <- "#8282d2"
customScheme$GeneRegionTrack$col <- NULL
customScheme$GeneRegionTrack$background.title <- "transparent"
customScheme$AnnotationTrack$background.title <- "transparent"
customScheme$AlignmentsTrack$background.title <- "transparent"
addScheme(customScheme, "customScheme")
options(Gviz.scheme="customScheme")

tracks <- list()

# Make GRanges object of Gene Model Annotation
if(!is.null(opt$options$gff)){
    gtf_File <- unlist(strsplit(opt$options$gff, ","))
}else{
    stop("Please provide at least one gff/gtf file.")
}
## First set colors for multiple GeneRegionTrack
if(!is.null(opt$options$color)){
    geneColor <- unlist(strsplit(opt$option$color, ","))
    if(length(geneColor)!=length(gtf_File)){
        stop("Color number is not equal to gtf number")
    }
}else{
    geneColor <- rep("#8282d2", length(gtf_File))
}


if(!is.null(opt$options$ref)){
    geneLabels <- unlist(strsplit(opt$options$ref, ","))
}else{
    geneLabels <- rep(NULL, length(gtf_File))
}
geneTracks <- list()
for(i in 1:length(gtf_File)){
    if(str_detect(gtf_File[i], "\\.sqlite$")){
        gene_TxDb <- loadDb(gtf_File[i])
    }else if(str_detect(gtf_File[i], "\\.(gtf|gff3)$")){
        gene_TxDb <- makeTxDbFromGFF(gtf_File[i], format="auto")
    }else{
        stop("Please provide GTF/GFF file as input.")
    }
    ## Extract intron for assembled dataset
    if(i==length(gtf_File)){
            introns <- unlist(intronsByTranscript(gene_TxDb))
    }        
    gTrack <- GeneRegionTrack(gene_TxDb, name = geneLabels[i],
                              fill = geneColor[i], col.line = NULL,
                              transcriptAnnotation = "mRNA",
                              showTitle = FALSE, alpha.title=0,
                              col.border.title="transparent",
                              col = NULL)
    geneTracks[[i]] <- gTrack
    
}
message("Finished load GTF...")



##loci <- paste(opt$options$chr, ":",
##              opt$options$start-1000, "-",
##              opt$options$end+1000, sep="")
##
##which <- GRanges(loci)
##number_of_transcripts <- sum(countOverlaps(transcripts(gmx_TxDb), which))
##gmx_gr <- import(con = opt$options$gff)
##gmx_db <- makeTxDbFromGRanges(gmx_gr)
genomeTrack <- GenomeAxisTrack() # Position Ruler Track
##geneTrack <- GeneRegionTrack(gmx_TxDb, name = "Gene Model",
##                             fill = "#8282d2", col.line = NULL,
##                             transcriptAnnotation = "transcript",
##                             col = NULL)
#Fst.data <- read.table(opt$options$fst, header=T)
#Fst.data <- cbind(Fst.data[,1:3], Fst.data[,5])
#Fst.gr <- makeGRangesFromDataFrame(Fst.data, seqnames.field = "CHROM", start.field = "BIN_START", end.field = "BIN_END", keep.extra.columns = T)
#Fst.track <- DataTrack(range = Fst.gr, type = "p", name = "Fst", ylim=c(0,0.8))
if(!is.null(opt$options$bam)){
    inputFile <- unlist(strsplit(opt$options$bam, ","))
    ## Make expTracks

    if(!is.null(opt$options$label)){
        sampleLabels <- unlist(strsplit(opt$options$label, ","))
    }else{
        sampleLabels <- rep(NULL, length(inputFile))
    }
    
    expTracks <- list()
    if(opt$options$mode == "exp"){
        exp_type <- "coverage"
    }else if(opt$options$mode == "sashimi"){
        exp_type <- c("coverage", "sashimi")
    }else{
        print("mode arg must be either exp or sashimi")
    }
    
    for(i in 1:length(inputFile)){
        newTrack <- AlignmentsTrack(range = inputFile[i],
                                    name = sampleLabels[i],
                                    isPaired = TRUE)
        expTracks <- c(expTracks, newTrack)
    }
    tracks <- c(tracks, expTracks)
}
message("Finished adding alignment tracks...")

if(!is.null(opt$options$PAS)){
    if(is.null(opt$options$size)){
        sizes <- c(sizes, 0.2)
    }
    PAS <- read.table(opt$options$PAS, header = F)
    colnames(PAS) <- c("chr","start",
                       "end", "id", "col5", "strand")
    PAS_Track <- AnnotationTrack(start=PAS$start, width=15,
                                 chromosome=PAS$chr,
                                 strand=PAS$strand, group=PAS$id,
                                 name = "PAS",
                                 col.line = NULL, col = NULL,
                                 fill = "#7fbf7b",
                                 stacking="dense")
    tracks <- c(tracks, genomeTrack, PAS_Track, geneTracks)
    message("Added PAS track...")
}else{
    tracks <- c(tracks, genomeTrack, geneTracks)
}
### Disable auto size

##auto_size <- function(number_of_tracks, stackHeight, number_of_transcripts, ref){
##    ## Since the auto sizing funtion in Gviz is not friendly and functional
##    ## I write my own autosizing function with tuning sizes parameter
##    track_size <- c(rep(1, number_of_tracks), 0.8)
##    gTrack_size <- 0.24*stackHeight*number_of_transcripts+(1-stackHeight)*0.24
##    if(!is.null(ref)){
##        ref_trans_num <- sum(countOverlaps(transcripts(ref_TxDb), which))
##        ratio1 <- ref_trans_num/(ref_trans_num+number_of_transcripts)
##        ratio2 <- number_of_transcripts/(ref_trans_num+number_of_transcripts)
##        gTrack_size <- gTrack_size*2
##        track_size <- c(track_size, 2*gTrack_size*ratio1, 2*gTrack_size*ratio2)
##    }else{
##        track_size <- c(track_size, gTrack_size)
##    }
##    return(track_size)
##}

if(opt$options$log){
    transform_function <- function(x){log(x+1, 2)}
}else{
    transform_function <- function(x){x}
}

if(!is.null(opt$options$size)){
    sizes <- unlist(strsplit(opt$option$size, ","))
    sizes <- as.numeric(sizes)
}else{
    sizes <- c(rep(0.8, length(expTracks)), 0.6, rep(1, length(geneTracks)))
}

if(!is.null(opt$options$peptide)){
    if(is.null(opt$options$size)){
        sizes <- c(sizes, 0.8)
    }
    peptide <- read.table(opt$options$peptide, header=T)
    peptide$feature <- as.character(peptide$feature)
    pepTrack <- GeneRegionTrack(peptide, name = "Peptides",
                                col.line = NULL, col = NULL,
                                fill = "#fdb863",
                                transcriptAnnotation="symbol")
    tracks <- c(tracks, pepTrack)
}


if(!is.null(opt$options$bam)){
    message("Track sizes are ", sizes)
    f <- function() {
        plotTracks(tracks,
                   chromosome = opt$options$chr,
                   from= as.numeric(opt$options$start),
                   to= as.numeric(opt$options$end),
                   type = exp_type,
                   ##background.title = "white",
                   ##col.axis = "lightgrey",
                   ##col.title = "black",
                   coverageHeight = 0.01, # default 0.1
                   minCoverageHeight = 0, # default 50
                   sashimiHeight = 0.01, # default 0.1
                   lwd.sashimiMax = 2,  # line width of sashimi, default 10, too wide
                   minSashimiHeight = 0, # default 50
                   sashimiScore = 5, # default 1
                   stackHeight = 0.5,
                   background.title = "transparent",
                   col.axis= "black",
                   col.title="black",
                   fontcolor="black",
                   fontface.title=2,
                   fontsize.title=10,
                   rotation.title=90,
                   sashimiFilter=introns,
                   ##sashimiNumbers = TRUE,
                   sizes = sizes,
                   transformation = transform_function)
    }
    p1 <- as.ggplot(f)
    
}else{
    f <- function() {
        plotTracks(tracks,
                   chromosome = opt$options$chr,
                   from= as.numeric(opt$options$start),
                   to= as.numeric(opt$options$end),
                   sizes = sizes)
    }
    p1 <- as.ggplot(f)
}
p1 <-  p1 + theme(plot.margin= margin(0,0,0,0, "mm"))
message("Finishe major plot...")

if(!is.null(opt$options$PSI)){
    PSI <- read.table(file=opt$options$PSI, header = T)
    PSI <- PSI %>%
        filter(Trans==opt$options$Trans,
               Tissue==opt$options$Tissue) %>%
        mutate(Sample=factor(Sample, levels=c("0h","1h","2h","4h","24h","48h")))
    p2 <- ggplot(PSI, aes(x=Sample, y=PSI, ymin=down, ymax=up))
    p2 <- p2 + geom_errorbar(width=0.25) + geom_point(shape=21, size=4, fill="white")
    p2 <- p2 + scale_x_discrete(limits=rev(levels(PSI$Sample))) + coord_flip()
    p2 <- p2 + theme_classic() + xlab("") + theme(axis.text.y=element_blank())

    right_up <- sum(sizes[1:6])
    right_down <- sum(sizes)-sum(sizes[1:6])

    ## Align by x axis
    plots <- align_plots(p1, p2, align = "h", axis="bt")
    right_p <- plot_grid(plots[[2]], NULL, ncol = 1, rel_heights = c(right_up,right_down))
    p <- plot_grid(plots[[1]], right_p, nrow = 1, rel_widths = c(5,1))
}else{
    p <- p1
}
p <- p + theme(plot.margin=unit(c(0,0,0,0),"mm"))
ggsave(p, file=opt$options$out,
       width=opt$options$width, height=opt$options$height)
