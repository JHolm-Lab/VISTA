source("00_importation.R")

vog.list.mgss.pa <- readRDS("../../RDS_files/VOG_analysis_updated/vog.list.mgss.pa.RDS")

###  Perform Jaccard clustering of VOGs

vog.dist<-lapply(vog.list.mgss.pa, function(x) vegdist(x, method="jaccard", binary = TRUE))
saveRDS(vog.dist, "../../RDS_files/VOG_analysis_updated/vog.dist.RDS")

vog.hc<-lapply(vog.dist, function(x) hclust(x, method="ward.D"))
saveRDS(vog.hc, "../../RDS_files/VOG_analysis_updated/vog.hc.RDS")

print("05_vog_jaccard_clustering.R has run")
