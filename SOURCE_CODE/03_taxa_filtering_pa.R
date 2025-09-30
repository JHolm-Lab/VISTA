source("00_importation.R")

#### NEW VERSION BASED ON GTDB PROTEIN COUNT

gene.counts.list <- readRDS("../../RDS_files/Gene_analysis/gene.counts.list.RDS")
ngl.gene.counts.list <- readRDS("../../RDS_files/Gene_analysis/ngl.gene.counts.list.RDS")
gen.sizes <- read_xlsx("../../SOURCE_DATA/gen.sizes2.xlsx")
ngl.gene.counts.list.mgSs<-list()
mgss.taxa<-vector()
other.taxa<-vector()


for (i in gen.sizes$Taxon){
  
  print(i)
  l<-gen.sizes[gen.sizes$Taxon == i, "protein_count"]$protein_count
  print(l)
  if (l == 0) { l<-1000}
  print(l)

  samples.w.nGenes_0.5X<-apply(ngl.gene.counts.list[[i]], 2, function(x) sum (x > 0.5))
  print("samples.w.nGenes_0.5X ok ")
  test<-as.data.frame(ngl.gene.counts.list[[i]][, names(ngl.gene.counts.list[[i]]) %in% names(samples.w.nGenes_0.5X[samples.w.nGenes_0.5X >= (l*.75)])])
  print("test ok")
  
  if (ncol(test) > 19){
    ngl.gene.counts.list.mgSs[[i]]<-test
    mgss.taxa<-append(mgss.taxa, i)
  } else {
    other.taxa<-append(other.taxa, i)
  }
}

ngl.gene.counts.list.mgSs.pa<-lapply(ngl.gene.counts.list.mgSs, function(x) replace(x, x > 0, 1))
ngl.gene.counts.list.mgSs.pa.clean<-lapply(ngl.gene.counts.list.mgSs.pa, function(x) x[which(rowSums(x) > 0), ])

saveRDS(object = ngl.gene.counts.list.mgSs.pa.clean, file = "../../RDS_files/Gene_analysis/ngl.gene.counts.list.mgSs.pa.clean1.RDS")

ngl.gene.counts.list.mgSs.pa.clean <- readRDS("../../RDS_files/Gene_analysis/ngl.gene.counts.list.mgSs.pa.clean1.RDS")

mgss.candidates.df <- data.frame()

for (i in names(ngl.gene.counts.list.mgSs.pa.clean)){
  protein.count <- gen.sizes[gen.sizes$Taxon == i, "protein_count"]$protein_count
  if (protein.count == 0) {protein.count <- 1000}
  print(paste(i, protein.count))
  a <- colSums(ngl.gene.counts.list.mgSs.pa.clean[[i]])
  mgss.candidates.df <- rbind(mgss.candidates.df, cbind("90"=sum(a > 0.9*protein.count) / dim(ngl.gene.counts.list.mgSs.pa.clean[[i]])[2] * 100,
                                                        "80"=sum(a > 0.8*protein.count) / dim(ngl.gene.counts.list.mgSs.pa.clean[[i]])[2] * 100,
                                                        "70"=sum(a > 0.7*protein.count) / dim(ngl.gene.counts.list.mgSs.pa.clean[[i]])[2] * 100,
                                                        "60"=sum(a > 0.6*protein.count) / dim(ngl.gene.counts.list.mgSs.pa.clean[[i]])[2] * 100,
                                                        "50"=sum(a > 0.5*protein.count) / dim(ngl.gene.counts.list.mgSs.pa.clean[[i]])[2] * 100)
  )
}


mgss.candidates.df$species <- names(ngl.gene.counts.list.mgSs.pa.clean)

species <- mgss.candidates.df[mgss.candidates.df$`90` >= 80,]$species

ngl.gene.counts.list.mgSs.pa.clean <- ngl.gene.counts.list.mgSs.pa.clean[names(ngl.gene.counts.list.mgSs.pa.clean) %in% species]
saveRDS(object = ngl.gene.counts.list.mgSs.pa.clean, file = "../../RDS_files/Gene_analysis/ngl.gene.counts.list.mgSs.pa.clean.RDS")

# write_csv(mgss.candidates.df, "CSV_files/potential.mgss.candidates.csv")

ngl.gene.counts.list.colSums<-lapply(ngl.gene.counts.list.mgSs.pa.clean, as.data.frame(colSums))
# Make table of samples by Taxa where values are nGenes (using ngl Reads)
library(tidyverse)
ngl.nGenes<-ngl.gene.counts.list.colSums %>% map( ~ .x %>% rownames_to_column('rn')) %>% reduce(full_join, by = 'rn')
ngl.nGenes<-ngl.gene.counts.list.colSums %>% map( ~ .x %>% rownames_to_column('rn')) %>% reduce(full_join, by = 'rn') %>% rename_at(2:ncol(ngl.nGenes), ~ names(ngl.gene.counts.list.colSums))
ngl.nGenes[is.na(ngl.nGenes)]<-0
rownames(ngl.nGenes)<-ngl.nGenes$rn
ngl.nGenes$rn<-NULL
saveRDS(ngl.nGenes, "../../RDS_files/Gene_analysis/ngl.nGenes.RDS")


## ----SANITY CHECK SUMMARY DATA-----------------------------------------------------------------------
## Make a table Taxon nSamples Summary-nGenes/Sample

ngl.gene.counts.list.mgSs.pa.clean <- readRDS("../../RDS_files/Gene_analysis/ngl.gene.counts.list.mgSs.pa.clean.RDS")
nSamples.summary<-plyr::ldply(ngl.gene.counts.list.mgSs.pa.clean, function(x) ncol(x), .id = "Taxon")
colnames(nSamples.summary)<-c("Taxon", "nSamples")

gene.summary<-plyr::ldply(ngl.gene.counts.list.mgSs.pa.clean, function(x) round(summary(colSums(x))), .id = "Taxon")
ngl.gene.counts.list.pa.clean.summary<-merge(nSamples.summary, gene.summary, all=TRUE)

write.table(ngl.gene.counts.list.pa.clean.summary, "../../RDS_files/Gene_analysis/mgSs_Taxon_Summary_nSamples_nGenes.txt", quote=F, row.names = F, sep="\t")


print("03_taxa_filtering_pa.R has run")


#gene.counts.list <- readRDS("RDS_files/gene.counts.list.RDS")
#ngl.gene.counts.list <- readRDS("RDS_files/ngl.gene.counts.list.RDS")
#gen.sizes.sum <- readRDS("RDS_files/gen.sizes.sum.RDS")
#
### ----3.1 Filter taxa for samples w/ 0.75 and 0.5X coverage cutoff to determine if mgSs should be performed.----
#
### Bacillus_cereus is NULL in both gene.counts.list and ngl.gene.counts.list... rm.
#
#gene.counts.list$Bacillus_cereus<-NULL
#ngl.gene.counts.list$Bacillus_cereus<-NULL
#
#saveRDS(gene.counts.list, "RDS_files/gene.counts.list.RDS")
#saveRDS(ngl.gene.counts.list, "RDS_files/ngl.gene.counts.list.RDS")
#
#ngl.gene.counts.list.mgSs<-list()
#mgss.taxa<-vector()
#other.taxa<-vector()
#for (i in names(ngl.gene.counts.list)){
#  l<-gen.sizes.sum[gen.sizes.sum$taxonomy == i, "median"]
#  if (i == "Acidaminococcus_intestini"){ l<-2167 }
#  if (i == "Acinetobacter_baumannii"){ l<-3715 }
#  if (i == "BVAB1"){ l<-1421 }
#  if (i == "BVAB3"){ l<-1474.5 }
#  if (i == "Megasphaera_genomosp."){l<-1521}
#  if (i == "Prevotella_tannerae"){l<-2170}
#  if (i == "Dialister_microaerophilus"){l<-1208}
#  if (i == "Clostridium_difficile"){l<-3829}
#  if (length(l) == 0) { l<-1000}  #gene.counts.list[[i]
#  #sub<-apply(ngl.gene.counts.list[[i]], 2, function(x) sum (x > 0)) ## For each column, get the total number of rows with a non-zero value. 
#  samples.w.nGenes_0.5X<-apply(ngl.gene.counts.list[[i]], 2, function(x) sum (x > 0.5)) ## nGenes in samples w/ngl gene counts > 0.5X
#  test<-as.data.frame(ngl.gene.counts.list[[i]][, names(ngl.gene.counts.list[[i]]) %in% names(samples.w.nGenes_0.5X[samples.w.nGenes_0.5X >= (l*.75)])])
#  #if (is.null(ncol(test)) == FALSE){
#  #samples.keep<-names(nGenes_0.5X[nGenes_0.5X > 0.75*l])
#  if (ncol(test) > 19){ ## If enough samples (> 19) w/enough high-coverage (>0.5X) genes ( > 0.75*medGenome), then do mgSs.
#    ngl.gene.counts.list.mgSs[[i]]<-test
#    mgss.taxa<-append(mgss.taxa, i)
#  } else {
#    other.taxa<-append(other.taxa, i)
#  }
#}
#
### ----3.2 Convert to Presence Absence + rm rows (genes) summing to zero--------------------------------
#
#ngl.gene.counts.list.mgSs.pa<-lapply(ngl.gene.counts.list.mgSs, function(x) replace(x, x > 0, 1))
#ngl.gene.counts.list.mgSs.pa.clean<-lapply(ngl.gene.counts.list.mgSs.pa, function(x) x[which(rowSums(x) > 0), ])
#
#saveRDS(object = ngl.gene.counts.list.mgSs.pa.clean, file = "RDS_files/ngl.gene.counts.list.mgSs.pa.clean.RDS")
#
#ngl.gene.counts.list.colSums<-lapply(ngl.gene.counts.list.mgSs.pa.clean, as.data.frame(colSums))
#
## Make table of samples by Taxa where values are nGenes (using ngl Reads)
#library(tidyverse)
#
#ngl.nGenes<-ngl.gene.counts.list.colSums %>% map( ~ .x %>% rownames_to_column('rn')) %>% reduce(full_join, by = 'rn')
#ngl.nGenes<-ngl.gene.counts.list.colSums %>% map( ~ .x %>% rownames_to_column('rn')) %>% reduce(full_join, by = 'rn') %>% rename_at(2:ncol(ngl.nGenes), ~ names(ngl.gene.counts.list.colSums))
#
#ngl.nGenes[is.na(ngl.nGenes)]<-0
#rownames(ngl.nGenes)<-ngl.nGenes$rn
#ngl.nGenes$rn<-NULL
#saveRDS(ngl.nGenes, "RDS_files/ngl.nGenes.RDS")
#
### ----SANITY CHECK SUMMARY DATA-----------------------------------------------------------------------
### Make a table Taxon nSamples Summary-nGenes/Sample
#
### All taxa should have > 20 samples. All samples should have > 75% of the estimated nGenes in genome
#nSamples.summary<-plyr::ldply(ngl.gene.counts.list.mgSs.pa.clean, function(x) ncol(x), .id = "Taxon")
#colnames(nSamples.summary)<-c("Taxon", "nSamples")
#
#gene.summary<-plyr::ldply(ngl.gene.counts.list.mgSs.pa.clean, function(x) round(summary(colSums(x))), .id = "Taxon")
#ngl.gene.counts.list.pa.clean.summary<-merge(nSamples.summary, gene.summary, all=TRUE)
#
#write.table(ngl.gene.counts.list.pa.clean.summary, "CSV_files/mgSs_Taxon_Summary_nSamples_nGenes.txt", quote=F, row.names = F, sep="\t")
#
#
### Result - all taxa have enough samples, but some samples have few genes (like G.vag samples).
### This is because the 75% of the genome filtering step occurred prior to normalization by gene length.
### Following this step, it may have become apparent that using these rules, some samples do not have enough genes.
### However, it isn't a clear cut step to simply remove samples below a certain cutoff because in some taxa this could remove samples unnecessarily.
### Thus, proceed with dist. matrix and clustering.
#
#print("03_taxa_filtering_pa.R has run")
