source("00_importation.R")

all.clean.info <- readRDS("../../RDS_files/Gene_analysis/all.clean.info.RDS")
all.clean.100k <- readRDS("../../RDS_files/Gene_analysis/all.clean.100k.RDS")

## ----2. Normalize read counts by gene length and make species-specific list of dataframes------------

gene.counts.list<-list()
ngl.gene.counts.list<-list()

for (i in sort(unique(all.clean.info[!all.clean.info$Taxa == "", "Taxa"]))){
  ## Make table of samples by genes for the taxon
  print(i)
  gene.count<-as.data.frame(all.clean.100k[all.clean.info$Taxa == i, ])
  taxa.info<-all.clean.info[all.clean.info$Taxa == i, ]
  gene.count.s<-gene.count[colSums(gene.count) > 0]
  gene.count.s<-gene.count.s[rowSums(gene.count) > 0, ]
  gene.count.s.ngl<-gene.count.s ## Make a cp of the dataframe to fill in with ngl
  ## normalize read counts by gene length -> ngl
  for (n in rownames(gene.count.s))
  {
    gL<-taxa.info[taxa.info$Gene == n, "Length"]
    gene.count.s.ngl[n, ]<-apply(gene.count.s[n, ], 2, function(x) (x*150)/gL) ## normalize nReads to gene length. 
  }
  gene.counts.list[[i]]<-gene.count.s ## place the raw gene count into the list 
  ngl.gene.counts.list[[i]]<-gene.count.s.ngl
  rm(gene.count.s)
  rm(gene.count.s.ngl)
}

## remove df's that don't have enough genes (any)
gene.counts.list<-gene.counts.list[sapply(gene.counts.list, function(x) any(nrow(x) > 0))]
saveRDS(object = gene.counts.list, file = "../../RDS_files/Gene_analysis/gene.counts.list.RDS")

ngl.gene.counts.list<-ngl.gene.counts.list[sapply(ngl.gene.counts.list, function(x) any(nrow(x) > 0))]
saveRDS(object = ngl.gene.counts.list, file = "../../RDS_files/Gene_analysis/ngl.gene.counts.list.RDS")

print("02_normalization.R has run")
