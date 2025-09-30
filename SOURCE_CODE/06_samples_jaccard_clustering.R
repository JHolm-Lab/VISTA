source("00_importation.R")

vog.list.mgss.pa <- readRDS("../../RDS_files/VOG_analysis_updated/vog.list.mgss.pa.RDS")

#### Perform Jaccard clustering of Samples

samples.dist<-lapply(vog.list.mgss.pa, function(x) vegdist(t(x), method="jaccard", binary = TRUE))
saveRDS(samples.dist, "../../RDS_files/VOG_analysis_updated/samples.dist.vog.RDS")

samples.hc<-lapply(samples.dist, function(x) hclust(x, method="ward.D"))
saveRDS(samples.hc, "../../RDS_files/VOG_analysis_updated/samples.hc.vog.RDS")

print("06_samples_jaccard_clustering.R has run")
