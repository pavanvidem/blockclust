#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

suppressPackageStartupMessages(library(dendextend))

if(args[1]=="clust"){
    mydata<-read.table(args[2], sep=" ", header=FALSE, row.names=1) #"hclust_input.mtx"
    annot<-read.table(args[3], sep="\t", header=FALSE, row.names=1) #"blockgroup_annotations.txt"
    d <- as.dendrogram(hclust(dist(mydata, method="euclidean"), method="average"))
    pdf(args[4], w=20, h=7)
    par(mar = c(3,3,1,1))
    d <- set(d, "labels_cex", 0.001)
    cols <- c("lightskyblue", "red1", "tomato", "gray10", "blue", "cyan", "darkmagenta", "orange", "yellow", "grey", "red4", "pink",  "darkgreen",  "navajowhite2", "palegreen", "darkviolet", "chartreuse", "plum2", "sienna1", "orangered1")
    rna_type <- factor(annot$V2)
    col_rna_type <- cols[rna_type]
    plot(d)
    colored_bars(colors=col_rna_type, dend=d, rowLabels="rna_class", cex.rowLabels=1.5, y_shift=-0.1, y_scale=0.2)
    legend("topright", legend=levels(rna_type), fill=cols)
    dev.off()
}else if(args[1]=="hist"){
    mydat<-read.table(args[2], sep="\t", header=TRUE)
    par(xpd=TRUE)
    library('ggplot2') #suppressPackageStartupMessages(
#    pdf(args[3], width=15, height=10)
    cols <- c("lightskyblue", "red1", "tomato", "gray10", "blue", "cyan", "darkmagenta", "orange", "yellow", "grey", "red4", "pink",  "darkgreen",  "navajowhite2", "palegreen", "darkviolet", "chartreuse", "plum2", "sienna1", "orangered1") #,
    print(length(unique(mydat$RNA_type)))
    p<-qplot(Clusters, data=mydat, geom="bar", fill=factor(RNA_type))+theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_fill_manual(values=cols)
    ggsave(args[3], p, width=15, height=10)
}else{
    write("Undefined plotting mode!", stderr())
}

