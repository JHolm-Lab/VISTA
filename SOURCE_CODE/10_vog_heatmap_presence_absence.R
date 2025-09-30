source("00_importation.R")

samples.clusters <- readRDS("../../RDS_files/VOG_analysis_updated/samples.clusters.updated.vog.RDS")
ngl.gene.counts.list.mgSs.pa.clean <- readRDS("../../RDS_files/VOG_analysis_updated/vog.list.mgss.pa.RDS")
samples.hc <- readRDS("../../RDS_files/VOG_analysis_updated/samples.hc.vog.RDS")
genes.clusters <- readRDS("../../RDS_files/VOG_analysis_updated/vog.clusters.RDS")
genes.hc <- readRDS("../../RDS_files/VOG_analysis_updated/vog.hc.RDS")
samples.dist <- readRDS("../../RDS_files/VOG_analysis_updated/samples.dist.vog.RDS")

for (i in names(ngl.gene.counts.list.mgSs.pa.clean)){

  x<-ngl.gene.counts.list.mgSs.pa.clean[[i]]
  
  f = gplots:::heatmap.2
  ## dend row margins
  body(f)[[79]][["lwd"]]<-0.0001
  body(f)[[79]][["cex"]]<-0.2
  body(f)[[79]][[2]][[2]]<-7
  body(f)[[79]][[2]][[3]]<-7
  body(f)[[80]][[3]][[2]][[3]][[2]][[4]]<-TRUE
  
  ## dend col margins
  body(f)[[81]][["lwd"]]<-0.0001
  body(f)[[81]][["cex"]]<-0.2
  body(f)[[82]][[3]][[2]][[3]][[2]][[3]]<-TRUE

  
  main<-paste(gsub("_", " ", i), paste("\nNumber of Samples=", paste(ncol(x)), paste(", Number of VOG=", paste(nrow(x), sep=""), sep=""), sep=""), sep="")
  cols<-c(RColorBrewer::brewer.pal(12,"Set3"), RColorBrewer::brewer.pal(8,"Set1"))
  colfunc <- colorRampPalette(c("antiquewhite", "mediumblue"))

  png(paste("../../FIGURES/VOG_analysis_updated/vog_heatmap_presence_absence/", i,"heatmap_presence_absence.png", sep = "_"),width=3, height=3, bg="white", units="in", res=1200)
  par(cex.main=1)
  f(as.matrix(x), 
    Colv = set(as.dendrogram(samples.hc[[i]]), "branches_lwd", 0.5),
    Rowv = set(as.dendrogram(genes.hc[[i]]), "branches_lwd", 0.5),
    col=colfunc(2), 
    main = main,
    RowSideColors = cols[as.factor(as.numeric(genes.clusters[[i]]$vog_cluster))], 
    ColSideColors = cols[as.factor(as.numeric(samples.clusters[[i]]$sample_cluster))],
    labRow = "", 
    labCol = "",
    srtCol = 0,
    key = FALSE, 
    trace="none",
    margins = c(2, 0.1),
    cexCol = 0.7,
    offsetCol=-0.5,
    adjCol = c(NA,1)
  )

  legend(xpd = TRUE, x=1.1, y=1.55,
         legend = as.factor(sort(unique(as.numeric(samples.clusters[[i]]$sample_cluster)))),
         col = cols[as.factor(sort(unique(as.numeric(samples.clusters[[i]]$sample_cluster))))],
         lty = 1, lwd = 1, cex=.15, bty="n", title = "mgSs", ncol=2)
  
  legend(xpd = TRUE, x=-0.2, y=1.1,
         legend = as.factor(sort(unique(as.numeric(genes.clusters[[i]]$vog_cluster)))),
         col = cols[as.factor(sort(unique(as.numeric(genes.clusters[[i]]$vog_cluster))))],
         lty = 1, lwd = 1, cex=.15, bty="n", title = "VOG Clusters", ncol = 2)
  
  dev.off()
}


print("10_vog_heatmap_presence_absence.R has run")

