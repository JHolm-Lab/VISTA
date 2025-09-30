source("00_importation.R")

ngl.abund <- readRDS("../../RDS_files/VOG_analysis_updated/ngl.abund.vog.RDS")
ngl.nGenes.m <- readRDS("../../RDS_files/VOG_analysis_updated/ngl.nGenes.m.vog.RDS")
#mgss.candidates <- read_csv("potential.mgss.candidates.csv")


#species <- mgss.candidates[mgss.candidates$`90` >= 80,]$species

## ----7. transfer mgss numbers------------------------------------------------------------------------
## Make a df of samples x taxa and values are the mgSs for that sample.

if (exists("test.df") == TRUE){
 rm(test.df)
} else {
 print("no test.df")
}

samples.clusters <- readRDS("../../RDS_files/VOG_analysis_updated/samples.clusters.updated.vog.RDS")
#samples.clusters <- list()

#for (i in species){
#	samples.clusters[[i]] <- samples.clusters.list[[i]]
#}
#saveRDS(samples.clusters, "RDS_files/samples.clusters.RDS")
#samples.clusters <- readRDS("RDS_files/samples.clusters.RDS")

print("Data OK")

#samples.clusters <- readRDS("RDS_files/samples.clusters2.RDS")
test <- samples.clusters
test<-lapply(test, function(x) as.data.frame(x))

for (i in names(test)){
  names(test[[i]])<-c("sampleID", i)
}

for (i in names(test)){
  if (length(table(samples.clusters[[i]][, 2])) > 1) {
    if (exists("test.df") == FALSE) {
      test.df <- test[[i]]
    } else {
      test.df <- merge(test.df, test[[i]], all = TRUE)
    }
  }
}

samples.clusters.df<-test.df ## A df of samples x taxa and values are the mgSs for that sample.

## Make the resulting dataframe numeric and replace NA w/0.. req'd prior to melting 
rownames(samples.clusters.df)<-samples.clusters.df$sampleID
samples.clusters.df$sampleID<-NULL
test<-sapply(samples.clusters.df, function(x) as.numeric(x))
samples.clusters.num.df<-as.data.frame(test)
samples.clusters.num.df$sampleID<-test.df$sampleID
samples.clusters.num.df[is.na(samples.clusters.num.df)]<-0

## Make melted dataframe of sample ID to taxon to ss# 
samples.clusters.num.df.m<-reshape2::melt(samples.clusters.num.df, id.vars=c("sampleID"), variable.name="taxon", value.name="ss")
samples.clusters.num.df.m$MG_ss<-paste(samples.clusters.num.df.m$taxon, samples.clusters.num.df.m$ss, sep="_")

samples.clusters.num.df.m.d<-reshape2::dcast(samples.clusters.num.df.m, sampleID~MG_ss, value.var = "MG_ss", length)
n<-samples.clusters.num.df.m.d$sampleID

ngl.abund$sampleID<-rownames(ngl.abund)
ngl.abund.m<-reshape2::melt(ngl.abund, id.vars="sampleID", variable.name="taxon", value.name="abundance")
ngl.abund$sampleID<-NULL
ngl.abund.m$sam_tax<-paste(ngl.abund.m$sampleID, ngl.abund.m$taxon, sep="_")
samples.clusters.num.df.m$sam_tax<-paste(samples.clusters.num.df.m$sampleID, samples.clusters.num.df.m$taxon, sep="_")

ngl.abund.clusters<-merge(ngl.abund.m, samples.clusters.num.df.m, all=TRUE)
ngl.abund.clusters[is.na(ngl.abund.clusters)]<-0
ngl.abund.clusters$MG_ss<-ifelse(ngl.abund.clusters$MG_ss == "0", as.character(ngl.abund.clusters$taxon), as.character(ngl.abund.clusters$MG_ss))

ngl.abund.clusters<-merge(ngl.abund.clusters, ngl.nGenes.m, all.x=TRUE)

if (any(is.na(ngl.abund.clusters$nGenes))) {
  ngl.abund.clusters$nGenes[is.na(ngl.abund.clusters$nGenes)] <- 0
}

## Remove base taxa names that have been made into ss
list.of.ss.taxa<-as.vector(unique(ngl.abund.clusters[ngl.abund.clusters$ss == 1, 2]))
test<-ngl.abund.clusters
ngl.abund.clusters$MG_ss<-ifelse(ngl.abund.clusters$taxon %in% list.of.ss.taxa, yes = paste(ngl.abund.clusters$taxon, ngl.abund.clusters$ss, sep="_"),no=as.character(ngl.abund.clusters$MG_ss))
ngl.abund.clusters.cast<-reshape2::dcast(ngl.abund.clusters, sampleID~MG_ss, value.var = "abundance")
ngl.abund.clusters.cast[is.na(ngl.abund.clusters.cast)]<-0
rownames(ngl.abund.clusters.cast)<-ngl.abund.clusters.cast$sampleID
ngl.abund.clusters.cast$sampleID<-NULL
ngl.abund.clusters.cast<-ngl.abund.clusters.cast[, colSums(ngl.abund.clusters.cast) > 0]
ngl.abund.clusters.cast<-ngl.abund.clusters.cast[rowSums(ngl.abund.clusters.cast) > 0, ]

saveRDS(ngl.abund.clusters.cast, "../../RDS_files/VOG_analysis_updated/ngl.abund.clusters.cast.vog.RDS")
saveRDS(samples.clusters.num.df, "../../RDS_files/VOG_analysis_updated/samples.clusters.num.df.vog.RDS")

print("12_transfer_mgss_number_VOG.R has run")
