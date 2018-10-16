#! /usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))

# Make Options ##

args <- commandArgs(trailingOnly = T)


option.list <- list(
    make_option(c("-g", "--gff"), type="character", default=NULL,
                help="Input gene model (GFF File) [%default]"),
    make_option(c("-w", "--wig"), type="character", default=NULL,
                help="Input mapping wig/bigwig file  [%default]"),
    make_option(c("-c", "--chr"), type="character", default=NULL,
                help="Chromosome"),
    make_option(c("-s", "--start"), type="numeric", default=NULL,
                help="start position"),
    make_option(c("-e", "--end"), type="numeric", default=NULL,
                help="end position"),
    make_option(c("-o", "--out"), type="character", default=NULL,
                help="output pdf file")
)


desc <- "Plot genome and RNA-seq mapping track"
parser <- OptionParser(option_list=option.list,
                      description = desc,
                      usage="usage: %prog -g gff.input -w mapping.bw/wig -c chromosome_name -s start_position -e end_position -o out.pdf\n")
opt <- parse_args(parser, args=args, positional_arguments=T)

opt
length(opt$args)

if(length(opt$options) == 0)
    stop(print_help(parser))


suppressPackageStartupMessages(library(Gviz))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(GenomicFeatures))
options(ucscChromosomeNames=FALSE) # Soybean chromosome name does not stick to ucsc nameing convention



# Make GRanges object of Gene Model Annotation


gmx_TxDb <- makeTxDbFromGFF(opt$options$gff)

genomeTrack <- GenomeAxisTrack() # Position Ruler Track
geneTrack <- GeneRegionTrack(gmx_TxDb, name = "Genes")
#Fst.data <- read.table(opt$options$fst, header=T)
#Fst.data <- cbind(Fst.data[,1:3], Fst.data[,5])
#Fst.gr <- makeGRangesFromDataFrame(Fst.data, seqnames.field = "CHROM", start.field = "BIN_START", end.field = "BIN_END", keep.extra.columns = T)
#Fst.track <- DataTrack(range = Fst.gr, type = "p", name = "Fst", ylim=c(0,0.8))

expTrack <- DataTrack(range = opt$options$wig,
                      type = "h",
                      chromosome = opt$options$chr,
                      name = "Exp")


pdf(opt$options$out, height=3)
plotTracks(list(genomeTrack, geneTrack, expTrack),
           chromosome = opt$options$chr,
           from= as.numeric(opt$options$start),
           to= as.numeric(opt$options$end))
dev.off()
