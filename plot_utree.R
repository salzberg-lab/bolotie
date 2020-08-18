#! /usr/bin/env Rscript

# loads the tree and cluster
# creates a color scheme for clusters for subsequent use
# save the plot

# apeTree.R treefile.nwk clusterfile.tsv palette.tsv output.png paths

args = commandArgs(trailingOnly=TRUE)

suppressPackageStartupMessages(library(ape,quietly=TRUE,warn.conflicts=FALSE,verbose=FALSE))
suppressPackageStartupMessages(library(plyr,quietly=TRUE,warn.conflicts=FALSE,verbose=FALSE))

tree_fname <- args[1]
cluster_fname <-args[2]
palette_fname <-args[3]
out_fname <- args[4]
paths_fname <- args[5]

# load palette
cur_palette_df = read.table(palette_fname,header=FALSE,col.names=c("clu_id","cur_palette"),sep="\t",comment.char = "-")
if(min(cur_palette_df$clu_id)==0){
    cur_palette_df$clu_id <- cur_palette_df$clu_id+1
}
cur_palette_df$cur_palette <- as.character(cur_palette_df$cur_palette)

# load clusters and assign colors
clusters = read.table(cluster_fname,header=FALSE,col.names=c("SequenceName","ClusterNumber"))
if(min(clusters$ClusterNumber)==0){
    clusters$ClusterNumber <- clusters$ClusterNumber+1
}

clusters_palette = merge(clusters,cur_palette_df,by.x='ClusterNumber',by.y='clu_id')[c("SequenceName","cur_palette")]

# load tree
tree <- read.tree(file = tree_fname)
tree_tips_df = data.frame(tree$tip.label)
colnames(tree_tips_df) <- c("SequenceName")
tipcols = join(tree_tips_df,clusters_palette,by='SequenceName',type="left")["cur_palette"]
tipcols = as.vector(tipcols$cur_palette)

# load paths
paths <- read.table(paths_fname,header=FALSE,col.names=c("seqid","ll","br"),sep="\t")
paths <- paths[paths$seqid %in% tree$tip.label,]

# find edges were both nodes belong to the same cluster and set the color
edgecols <- rep('black', nrow(tree$edge))
ec=0
for(row in 1:nrow(tree$edge)) {
    ec=ec+1
    edge2 = tree$edge[row,2]
    
    if(edge2<=length(tipcols)){
        edgecols[ec]=tipcols[edge2]
    }
}

png(out_fname,height=1200,width=1200)
plot(tree,"u",
     use.edge.length = TRUE,
     show.tip.label = FALSE,
     tip.color=tipcols,
     edge.color=edgecols,
     edge.width=5)

seqids<-paths$seqid
ii<-sapply(seqids,grep,tree$tip.label)
for(i in ii){
    tiplabels(text="",tip=i, pch=21, col="black", bg='red', cex=4)
}
dev.off()
