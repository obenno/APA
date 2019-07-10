#! /usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))

# Make Options ##

args <- commandArgs(trailingOnly = T)


option.list <- list(
    make_option(c("-r", "--ref"), type="character", default=NULL,
                help="Whether include reference gene model (GFF File) [%default]"),
    make_option(c("-g", "--gff"), type="character", default=NULL,
                help="Input gene model (GFF File) [%default]"),
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
                help="output file width"),
    make_option(c("--height"), type="numeric", default=5,
                help="output file height"),
    make_option(c("--label"), type="character", default=NULL,
                help="label of samples, could be a comma separated list [%default]")
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
options(ucscChromosomeNames=FALSE) # Soybean chromosome name does not stick to ucsc nameing convention

# Make GRanges object of Gene Model Annotation
if(!is.null(opt$options$ref)){
    ref_TxDb <- makeTxDbFromGFF(opt$option$ref, format="auto")
    refTrack <- GeneRegionTrack(ref_TxDb, name = "Ref Gene Model",
                                fill = "#8282d2", col.line = NULL,
                                transcriptAnnotation = "mRNA",
                                col = NULL)
}

gmx_TxDb <- makeTxDbFromGFF(opt$options$gff, format="auto")

loci <- paste(opt$options$chr, ":",
              opt$options$start-1000, "-",
              opt$options$end+1000, sep="")

which <- GRanges(loci)
number_of_transcripts <- sum(countOverlaps(transcripts(gmx_TxDb), which))
#gmx_gr <- import(con = opt$options$gff)
#gmx_db <- makeTxDbFromGRanges(gmx_gr)
genomeTrack <- GenomeAxisTrack() # Position Ruler Track
geneTrack <- GeneRegionTrack(gmx_TxDb, name = "Gene Model",
                             fill = "#8282d2", col.line = NULL,
                             transcriptAnnotation = "transcript",
                             col = NULL)
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
}

auto_size <- function(number_of_tracks, stackHeight, number_of_transcripts){
    # Since the auto sizing funtion in Gviz is not friendly and functional
    # I write my own autosizing function with tuning sizes parameter
    track_size <- c(rep(1, number_of_tracks), 0.6)
    gTrack_size <- 0.3*stackHeight*number_of_transcripts+(1-stackHeight)*0.3
    track_size <- c(track_size, gTrack_size)
    return(track_size)
}

pdf(opt$options$out,
    height=opt$options$height,
    width=opt$options$width)
tracks <- list()
if(!is.null(opt$options$bam)){
    if(!is.null(opt$options$ref)){
        tracks <- c(tracks, expTracks, genomeTrack, refTrack, geneTrack)
    }else{
        tracks <- c(tracks, expTracks,genomeTrack, geneTrack)
    }
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
               sashimiScore = 1, # default 1
               stackHeight = 0.5,
               ##sashimiNumbers = TRUE,
               sizes = auto_size(length(expTracks), 0.5, number_of_transcripts))
               ##sizes=c(rep(1, length(expTracks)), 0.6, 0.3))
}else{
    if(!is.null(opt$options$ref)){
        tracks <- c(tracks, genomeTrack, refTrack, geneTrack)
    }else{
        tracks <- c(tracks, genomeTrack, geneTrack)
    }
    plotTracks(tracks,
               chromosome = opt$options$chr,
               from= as.numeric(opt$options$start),
               to= as.numeric(opt$options$end),
               sizes = c(0.3,0.7))
}    
dev.off()
