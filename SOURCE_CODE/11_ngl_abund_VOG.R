source("00_importation.R")

ngl.gene.counts.list <- readRDS("../../RDS_files/VOG_analysis_updated/vog.count.list.RDS")

## ----6. Produce relative abundance table from ngl counts -----------

ngl.gene.counts.list.colSums<-lapply(ngl.gene.counts.list, function(x) as.data.frame(colSums(x)))

library(tidyverse)

ngl.abund<-ngl.gene.counts.list.colSums %>% map( ~ .x %>% rownames_to_column('rn')) %>% reduce(full_join, by = 'rn')
ngl.abund<-ngl.gene.counts.list.colSums %>% map( ~ .x %>% rownames_to_column('rn')) %>% reduce(full_join, by = 'rn') %>% rename_at(2:ncol(ngl.abund), ~ names(ngl.gene.counts.list.colSums))

ngl.abund[is.na(ngl.abund)]<-0
rownames(ngl.abund)<-ngl.abund$rn
ngl.abund$rn<-NULL

saveRDS(ngl.abund, "../../RDS_files/VOG_analysis_updated/ngl.abund.vog.RDS")

ngl.gene.counts.list.pa<-lapply(ngl.gene.counts.list, function(x) replace(x, x >= 0.5, 1))
ngl.gene.counts.list.pa<-lapply(ngl.gene.counts.list.pa, function(x) replace(x, x < 0.5, 0))
ngl.gene.counts.list.colSums<-lapply(ngl.gene.counts.list.pa, as.data.frame(colSums))
ngl.gene.counts.list.colSums<-lapply(ngl.gene.counts.list.colSums, function(x) x["sampleID"]<-rownames(x))

# Make table of samples by Taxa where values are nGenes (using ngl Reads)

library(tidyverse)

ngl.nGenes<-ngl.gene.counts.list.colSums %>% map( ~ .x %>% as.data.frame() %>% rownames_to_column('rn')) %>% reduce(full_join, by = 'rn')
ngl.nGenes<-ngl.gene.counts.list.colSums %>% map( ~ .x %>% as.data.frame() %>% rownames_to_column('rn')) %>% reduce(full_join, by = 'rn') %>% rename_at(2:ncol(ngl.nGenes), ~ names(ngl.gene.counts.list.colSums))

ngl.nGenes[is.na(ngl.nGenes)]<-0
rownames(ngl.nGenes)<-ngl.nGenes$rn
ngl.nGenes.m<-reshape2::melt(ngl.nGenes, id.vars="rn", variable.name="taxon", value.name = "Coverage")
names(ngl.nGenes.m)[1]<-"sampleID"
ngl.nGenes$rn<-NULL

saveRDS(ngl.nGenes.m, "../../RDS_files/VOG_analysis_updated/ngl.nGenes.m.vog.RDS")

print("11_ngl_abund_VOG.R has run")
