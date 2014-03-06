options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE) #take commandline args
if(length(args) < 2) {
  args <- c("--help")
}
Dat_file <- as.character(unlist(args[1])) #sample 1 file
print(Dat_file) 
Dat_name <- as.character(unlist(args[1])) #sample 1 file

results_graph <- unlist(paste("circos_",args[2],".pdf",sep='')) #graph file
print(results_graph) 
library(methylKit)
library("graphics")
ideoDMC <- function(methylDiff.obj, chrom.length, difference = 25, 
    qvalue = 0.01, circos = FALSE, title = "test", hyper.col = "magenta", 
    hypo.col = "green") {
    require(methylKit)
    require(GenomicRanges)
    require(ggbio)

    # chrom.length
    myIdeo <- GRanges(seqnames = names(chrom.length), ranges = IRanges(start = 1, 
        width = chrom.length))
    seqlevels(myIdeo) = names(chrom.length)
    seqlengths(myIdeo) = (chrom.length)


    hypo = get.methylDiff(methylDiff.obj, difference = difference, qvalue = qvalue, 
        type = "hypo")
    hyper = get.methylDiff(methylDiff.obj, difference = difference, qvalue = qvalue, 
        type = "hyper")

    g.per = as(hyper, "GRanges")
    seqlevels(g.per, force=TRUE) = seqlevels(myIdeo)
    seqlengths(g.per)=(chrom.length)

    g.po = as(hypo, "GRanges")
    seqlevels(g.po, force=TRUE) = seqlevels(myIdeo)
    seqlengths(g.po)=(chrom.length)

    values(g.po)$id = "hypo"
    values(g.per)$id = "hyper"

    if (circos) {

        p <- ggplot() + layout_circle(myIdeo, geom = "ideo", fill = "gray70", 
            radius = 39, trackWidth = 2)


        p <- p + layout_circle(c(g.po, g.per), geom = "point", 
                 size = 1, aes(x = midpoint, 
            y = meth.diff, color = id), radius = 25, trackWidth = 30) +              
            scale_colour_manual(values = c(hyper.col, hypo.col))
        p + layout_circle(myIdeo, geom = "text", aes(label = seqnames), 
            vjust = 0, radius = 55, trackWidth = 7) + opts(title = title)

    } else {

        p <- ggplot() + layout_karyogram(myIdeo, cytoband = FALSE)
        p + layout_karyogram(c(g.po, g.per), geom = "point", size = 1, 
         aes(x = midpoint, 
            y = meth.diff, color = id)) + scale_colour_manual(values = c(hyper.col, 
            hypo.col)) + opts(title = title)

    }
}
library(BSgenome)
load(Dat_file)
gene.obj=read.transcript.features("gmax_bed.bed")
seqlengths(gene.obj)<-c(55915595,51656713,47781076,49243852,41936504,50722821,44683157,46995532,46843750,50969635,39172790,40113140,44408971,49711204,50939160,37397385,41906774,62308140,50589441,46773167)
chr.len = seqlengths(gene.obj)
pdf(file=results_graph)  
ideoDMC(myDiff, chrom.length = chr.len, difference = 25, qvalue = 0.01, 
    circos = TRUE, title = Dat_name, hyper.col = "magenta", hypo.col = "green")
dev.off()
