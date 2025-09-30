source("00_importation.R")

ngl.abund.clusters.cast <- readRDS("../../RDS_files/VOG_analysis_updated/ngl.abund.clusters.cast.vog.RDS")
#potential_mgss <- read_csv("CSV_files/potential.mgss.candidates.csv")

print("read ok")

# ----8. Make mgCSTs & official mgCST tables----------------------------------------------------------

## Remove taxa in which most hits go to a few genes, and these genes have bestblast matches (nr/nt) to human
ngl.abund.clusters.cast$Burkholderia_mallei<-NULL
ngl.abund.clusters.cast$Coprobacillus_sp.<-NULL
ngl.abund.clusters.cast$Enterococcus_faecium<-NULL
ngl.abund.clusters.cast$Beggiatoa_sp.<-NULL
ngl.abund.clusters.cast$Candidatus_Pelagibacter<-NULL

print("OK")
## Keep taxa with at leat 80% of the samples containing at leat 90% of the genome size
#not_considered_species  <- potential_mgss[potential_mgss$`90` < 80,]$species
#print("OK2")
#for (i in not_considered_species){
#	print(i)
#	ngl.abund.clusters.cast$i <- NULL
#}

ngl.abund.clusters.cast<-ngl.abund.clusters.cast[ ,colSums(ngl.abund.clusters.cast) > 0]
ngl.abund.clusters.cast<-ngl.abund.clusters.cast[rowSums(ngl.abund.clusters.cast) > 0, ]

#ngl.abund.clusters.cast<-ngl.abund.clusters.cast[colSums(ngl.abund.clusters.cast) > 0, ]

## Remove those samples that were dominated by the above taxa (mgCST27)
#ngl.abund.clusters.cast<-ngl.abund.clusters.cast[!rownames(ngl.abund.clusters.cast) %in% mgCST27, ]
relabund<-ngl.abund.clusters.cast/rowSums(ngl.abund.clusters.cast) ### Currently does not include reads from non-mgss species. How to add? 
mgCST.dist<-philentropy::JSD(as.matrix(relabund))

saveRDS(mgCST.dist, "../../RDS_files/VOG_analysis_updated/mgCST.dist.vog.RDS")
saveRDS(ngl.abund.clusters.cast, "../../RDS_files/VOG_analysis_updated/ngl.abund.clusters.cast.vog.RDS")

rownames(mgCST.dist)<-rownames(relabund)
colnames(mgCST.dist)<-rownames(relabund)
mgCST.hclust <- hclust(as.dist(mgCST.dist), method="ward.D")

saveRDS(mgCST.hclust, "../../RDS_files/VOG_analysis_updated/mgCST.hclust.vog.RDS")


print("13_mgCSTs_vog.R has run")

