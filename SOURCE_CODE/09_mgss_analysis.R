source("00_importation.R")

samples.hc<-readRDS("../../RDS_files/VOG_analysis_updated/samples.hc.vog.RDS")
genes.hc<-readRDS("../../RDS_files/VOG_analysis_updated/vog.hc.RDS")
ngl.gene.counts.list.mgSs.pa.clean<-readRDS("../../RDS_files/VOG_analysis_updated/vog.list.mgss.pa.RDS")
ngl.gene.counts.list <- readRDS("../../RDS_files/VOG_analysis_updated/vog.count.list.RDS")
samples.dist <- readRDS("../../RDS_files/VOG_analysis_updated/samples.dist.vog.RDS")
genes.dist <- readRDS("../../RDS_files/VOG_analysis_updated/vog.dist.RDS")

ngl.nGenes <- readRDS("../../RDS_files/VOG_analysis_updated/ngl.nVOG.RDS")

print("Data importation ok")

# ## ----5. mgSs analysis: Define clusters of subset of samples for mgss and genes to explain the mgss----

samples.clusters <- readRDS("../../RDS_files/VOG_analysis_updated/samples.clusters.updated.vog.RDS")
genes.clusters <- readRDS("../../RDS_files/VOG_analysis_updated/vog.clusters.RDS")


ngl.gene.counts.list.colSums<-lapply(ngl.gene.counts.list, function(x) as.data.frame(colSums(x)))
ngl.abund<-ngl.gene.counts.list.colSums %>% map( ~ .x %>% rownames_to_column('rn')) %>% reduce(full_join, by = 'rn')
ngl.abund<-ngl.gene.counts.list.colSums %>% map( ~ .x %>% rownames_to_column('rn')) %>% reduce(full_join, by = 'rn') %>% rename_at(2:ncol(ngl.abund), ~ names(ngl.gene.counts.list.colSums))
ngl.abund[is.na(ngl.abund)]<-0
rownames(ngl.abund)<-ngl.abund$rn
ngl.abund$rn<-NULL

saveRDS(ngl.abund, "../../RDS_files/VOG_analysis_updated/ngl.abund.vog.RDS")

mgSs.stats<-list()
library(ggplot2)

#for (i in names(samples.clusters)[!(names(samples.clusters) %in% not_to_cluster)]){
for (i in names(samples.clusters)){
  l<-as.data.frame(ngl.abund[i])
  l$sampleID<-rownames(ngl.abund)
  names(l)<-c("species_coverage", "sampleID")
  mgSs.stats[[i]]<-merge(l, samples.clusters[[i]], all.y=TRUE)
  l<-ngl.nGenes[i]
  l$sampleID<-rownames(l)
  names(l)<-c("No_VOG", "sampleID")
  mgSs.stats[[i]]<-merge(l, mgSs.stats[[i]], all.y=TRUE)
  l<-as.data.frame(rowSums(ngl.abund))
  l$sampleID<-rownames(ngl.abund)
  names(l)<-c("sample_coverage", "sampleID")
  mgSs.stats[[i]]<-merge(l, mgSs.stats[[i]], all.y=TRUE)
  mgSs.stats[[i]]$sample_cluster<-factor(mgSs.stats[[i]]$sample_cluster, levels=c(1:max(as.numeric(mgSs.stats[[i]]$sample_cluster))))
  mgSs.stats[[i]]$species_abund<-mgSs.stats[[i]]$species_coverage/mgSs.stats[[i]]$sample_coverage
}
saveRDS(mgSs.stats, "../../RDS_files/VOG_analysis_updated/vog.mgSs.stats.RDS")

### Summarize(coverage)  =>  Summarize(species_coverage) or sample_coverage ???
t<-lapply(mgSs.stats, function(x) Summarize(species_coverage ~ as.factor(sample_cluster), data=x, digits=0))
mgSs.stats.df<-do.call(rbind, t)
# t<-lapply(mgSs.stats, function(x) x %>% group_by(sample_cluster) %>% summarize())
write.csv(mgSs.stats.df, "../../RDS_files/VOG_analysis_updated/vog.mgSs.coverage.stats.csv", row.names = T, quote=F)

### y=log10(coverage)  =>  y=log10(species_coverage)

for (i in names(mgSs.stats)) {
print(ggplot(mgSs.stats[[i]], aes(x=No_VOG, y=log10(species_coverage), color=sample_cluster))+geom_point()+scale_color_brewer(palette="Set3")+theme_classic()+xlab("No. VOG")+ylab("log10(Coverage)")+theme(legend.position = "none")+ggtitle(i))
print(ggsave(paste("../../FIGURES/VOG_analysis_updated/mgSs_coverage/mgss_coverage_pdf/", paste(i, "subspecies_coverage_by_NoVOG.pdf", sep="_"), sep=""), width=4, height=3))
print(ggplot(mgSs.stats[[i]], aes(x=as.factor(sample_cluster), y=log10(species_coverage),fill=as.factor(sample_cluster), group=as.factor(sample_cluster),color=as.factor(sample_cluster)))+geom_violin()+scale_fill_brewer(palette="Set3")+scale_color_brewer(palette="Set3")+theme_classic()+theme(legend.position = "none")+ggtitle(i)+xlab("mgSs")+ylab("log10(Coverage)")+geom_jitter(size=0.05, color="black", width=0.2))
print(ggsave(paste("../../FIGURES/VOG_analysis_updated/mgSs_coverage/mgss_coverage_pdf/", paste(i, "subspecies_coverage_boxplot.pdf", sep="_"), sep=""), width=4, height=3))
}

#for (i in names(samples.clusters)[!(names(samples.clusters) %in% not_to_cluster)]){
for (i in names(samples.clusters)){
  v<-mgSs.stats[[i]]
  rownames(v)<-v$sampleID
  v.1<-v[names(ngl.gene.counts.list.mgSs.pa.clean[[i]]),]
  v.2<-v.1[order.dendrogram(as.dendrogram(samples.hc[[i]])),]
  pdf(paste("../../FIGURES/VOG_analysis_updated/mgSs_coverage/mgss_coverage_pdf/",paste(i, "mgCST_coverage_Barplot.pdf", sep=""), sep=""),height=0.287)
  par(lwd=0.0026, mai=c(0,0,0,0), omi=c(0,0,0,0))
  barplot(log10(v.2$species_coverage), cex.axis = 0.4, space = c(0,0), axes = F, col="black")
  dev.off()
}

for (i in names(mgSs.stats)) {
  ggplot(mgSs.stats[[i]], aes(x=No_VOG, y=log10(species_coverage), color=sample_cluster))+geom_point()+scale_color_brewer(palette="Set3") + 
    theme_classic() + 
    xlab("No. VOG") + 
    ylab("log10(Coverage)") + 
    theme(legend.position = "none") + 
    ggtitle(i)
  ggsave(paste("../../FIGURES/VOG_analysis_updated/mgSs_coverage/mgss_coverage_png/", i, "_subspecies_coverage_by_NoVOG.png", sep = ""), width=4, height=3)
  print(i)
  
  ggplot(mgSs.stats[[i]], aes(x=as.factor(sample_cluster), y=log10(species_coverage),fill=as.factor(sample_cluster), group=as.factor(sample_cluster),color=as.factor(sample_cluster))) + 
    geom_violin()+scale_fill_brewer(palette="Set3") + 
    scale_color_brewer(palette="Set3") + 
    theme_classic() + 
    theme(legend.position = "none") + 
    ggtitle(i) + 
    xlab("mgSs") + 
    ylab("log10(Coverage)") +
    geom_jitter(size=0.05, color="black", width=0.2)
  
  ggsave(paste("../../FIGURES/VOG_analysis_updated/mgSs_coverage/mgss_coverage_png/", i, "_subspecies_coverage_boxplot.png", sep = ""), width=4, height=3)
  
  print(i)
}



## ----SANITY CHECK sequencing coverage for each mgSs--------------------------------------------------

sample.cluster.stats<-list()
library(ggplot2)

#for (i in names(samples.clusters)[!(names(samples.clusters) %in% not_to_cluster)]){
for (i in names(samples.clusters)){
  #l<-as.data.frame(rowSums(ngl.abund))
  l<-as.data.frame(ngl.abund[i])
  l$sampleID<-rownames(ngl.abund)
  names(l)<-c("species_coverage", "sampleID")
  sample.cluster.stats[[i]]<-merge(l, samples.clusters[[i]], all.y=TRUE)
  l<-ngl.nGenes[i]
  l$sampleID<-rownames(l)
  names(l)<-c("No_VOG", "sampleID")
  sample.cluster.stats[[i]]<-merge(l, sample.cluster.stats[[i]], all.y=TRUE)
  l<-as.data.frame(rowSums(ngl.abund))
  l$sampleID<-rownames(ngl.abund)
  names(l)<-c("sample_coverage", "sampleID")
  sample.cluster.stats[[i]]<-merge(l, sample.cluster.stats[[i]], all.y=TRUE)
  sample.cluster.stats[[i]]$sample_cluster<-factor(sample.cluster.stats[[i]]$sample_cluster, levels=c(1:max(as.numeric(sample.cluster.stats[[i]]$sample_cluster))))
  sample.cluster.stats[[i]]$species_abund<-sample.cluster.stats[[i]]$species_coverage/sample.cluster.stats[[i]]$sample_coverage
}

saveRDS(sample.cluster.stats, "../../RDS_files/VOG_analysis_updated/sample.cluster.stats.vog.RDS")

library(FSA)
t<-lapply(sample.cluster.stats, function(x) Summarize(species_coverage ~ as.factor(sample_cluster), data=x, digits=0))
sample.cluster.stats.df<-do.call(rbind, t)
#t<-lapply(sample.cluster.stats, function(x) x %>% group_by(sample_cluster) %>% summarize())
write.csv(sample.cluster.stats.df, "../../RDS_files/VOG_analysis_updated/vog.mgSs.coverage.stats.csv", row.names = T, quote=F)


print("09_mgss_analysis.R has run")

