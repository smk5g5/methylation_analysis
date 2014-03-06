options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE) #take commandline args
if(length(args) < 2) {
  args <- c("--help")
}
my_output_directory <- "/home/skhan/data/soybeandata/5x/"
sample_1 <- as.character(unlist(args[1])) #sample 1 file
sample_2 <- as.character(unlist(args[2])) #sample 2 file
sample_name_1 <- as.character(unlist(args[3])) # sample 1 name
sample_name_2 <- as.character(unlist(args[4])) #sample 2 name
con_text <- as.character(unlist(args[6])) #context
results_graph <- unlist(paste(my_output_directory,args[5],".pdf",sep='')) #graph file
results_data <- unlist(paste(my_output_directory,args[5],".RData",sep='')) 
results_hyper <- unlist(paste(my_output_directory,'hyper_dmr_',args[5],".txt",sep=''))
results_hypo <- unlist(paste(my_output_directory,'hypo_dmr_',args[5],".txt",sep=''))	
other_info <- unlist(paste(my_output_directory,'Stats_',args[5],".txt",sep=''))	
all_diff <- unlist(paste(my_output_directory,'all_dmr_',args[5],".txt",sep=''))
library(methylKit)
library("graphics")

annotate_df <- function(methdiff,ann,output) {
all_mat=getMembers(ann)
mat_df <- as.data.frame(all_mat)
mat_df$id <-methdiff$id
dframe <- getData(methdiff)
m <-merge(dframe,mat_df,by="id")
write.table(m,file=output,sep="\t", col.names = T, row.names = F,quote=FALSE)
}


file.list <- list(sample_1, sample_2)
print(file.list) 
clist <- read(file.list, sample.id = list(sample_name_1, sample_name_2), assembly = "gmax", treatment = c(0,1), context = con_text)
filtered.clist <- filterByCoverage(clist, lo.count = 5,lo.perc = NULL, hi.count = NULL, hi.perc = 99.9)
# tiles=tileMethylCounts(filtered.clist,win.size=50,step.size=50)
newMeth <- unite(filtered.clist,destrand=FALSE)
gene.obj=read.transcript.features("gmax_bed.bed")
pdf(file=results_graph)
getCorrelation(newMeth, plot = T)
# clusterSamples(newMeth, dist = "correlation", method = "ward",plot = TRUE)
# PCASamples(newMeth, screeplot = TRUE)
# PCASamples(newMeth)
myDiff <- calculateDiffMeth(newMeth,num.cores = 6)
myDiff25p.hyper <- get.methylDiff(myDiff, difference = 25,qvalue = 0.01, type = "hyper")
myDiff25p.hypo <- get.methylDiff(myDiff, difference = 25,qvalue = 0.01, type = "hypo")
myDiff25p <- get.methylDiff(myDiff, difference = 25,qvalue = 0.01)
diffMethPerChr(myDiff, plot = TRUE, qvalue.cutoff = 0.01, newMeth.cutoff = 25) ##per chromosome plot
annotate.WithGenicParts(myDiff25p,gene.obj)
diffAnn=annotate.WithGenicParts(myDiff25p,gene.obj)
diff_hyper=annotate.WithGenicParts(myDiff25p.hyper,gene.obj)
diff_hypo=annotate.WithGenicParts(myDiff25p.hypo,gene.obj)
plotTargetAnnotation(diffAnn,precedence=TRUE,main="differential methylation annotation")
annotate_df(myDiff25p,diffAnn,all_diff)
annotate_df(myDiff25p.hyper,diff_hyper,results_hyper)
annotate_df(myDiff25p.hypo,diff_hypo,results_hypo)
# head(getAssociationWithTSS(diffAnn))
ds <- getTargetAnnotationStats(diffAnn,percentage=TRUE,precedence=TRUE)
write(ds, other_info)
# getFeatsWithTargetsStats(diffAnn,percentage=TRUE)
# all_mat=getMembers(diffAnn)
# hyper_mat=getMembers(diff_hyper)
# hypo_mat=getMembers(diff_hypo)
# write.table(hyper_mat,file=results_hyper,sep="\t", col.names = T, row.names = F)
# write.table(hypo_mat,file=results_hypo,sep="\t", col.names = T, row.names = F)
# write.table(all_mat,file=all_diff,sep="\t", col.names = T, row.names = F)
save(myDiff, file = results_data)
circos_graph <- unlist(paste("circos_",args[5],".pdf",sep='')) #graph file
# system("Rscript ideoDMC.r results_data args[5]")
ideoDMC <- function(methylDiff.obj, chrom.length, difference = 25, 
    qvalue = 0.01, circos = FALSE, title = "test", hyper.col = "magenta", 
    hypo.col = "green") {
    require(methylKit)
    require(GenomicRanges)
    require(ggbio)
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
# load(Dat_file)
# gene.obj=read.transcript.features("gmax_bed.bed")
seqlengths(gene.obj)<-c(55915595,51656713,47781076,49243852,41936504,50722821,44683157,46995532,46843750,50969635,39172790,40113140,44408971,49711204,50939160,37397385,41906774,62308140,50589441,46773167)
chr.len = seqlengths(gene.obj)
ideoDMC(myDiff, chrom.length = chr.len, difference = 25, qvalue = 0.01, 
circos = TRUE, title = args[5], hyper.col = "magenta", hypo.col = "green")
dev.off()