source("00_importation.R")

ngl.abund.clusters.cast <- readRDS("../../RDS_files/VOG_analysis_updated/ngl.abund.clusters.cast.vog.RDS")
mgCST.hclust <- readRDS("../../RDS_files/VOG_analysis_updated/mgCST.hclust.vog.RDS")
mgCST.dist <- readRDS("../../RDS_files/VOG_analysis_updated/mgCST.dist.vog.RDS")

#potential_mgss <- read_csv("CSV_files/potential.mgss.candidates.csv")

#species <- potential_mgss[potential_mgss$`90` >= 80,]$species

relabund<-ngl.abund.clusters.cast/rowSums(ngl.abund.clusters.cast)

print(sum(is.na(relabund)))
rownames(mgCST.dist)<-rownames(relabund)
colnames(mgCST.dist)<-rownames(relabund)

mgCSTs.sort.df <- data.frame(dtc = numeric(0),
                             meanRelabund = numeric(0),
                             domTaxa = character(0),
                             mgCST = numeric(0),
                             minClusterSize = numeric(0),
                             deepSplit = numeric(0)
)

mgCSTs.samples.df <- data.frame(dtc = numeric(0),
                                sampleID = character(0),
                                domTaxa = character(0),
                                relabund = numeric(0),
                                mgCST = numeric(0),
                                minClusterSize = numeric(0),
                                deepSplit = numeric(0)
)

mgss <- read.csv("../../RDS_files/VOG_analysis_updated/mgss.clustering.parameters.vog.csv")
species <- mgss$X

#for.sort<-c("Lactobacillus_crispatus", "Lactobacillus_gasseri", "Lactobacillus_iners", "Lactobacillus_jensenii", "Lactobacillus_paragasseri", "Lactobacillus_mulieris","UBA629_sp005465875", "Gardnerella_vaginalis", "Gardnerella_swidsinkii", "Gardnerella_piotii", "Gardnerella_leopoldii", "Streptococcus", "Prevotella", "Prevotella_bivia", "MultiGenera","Bifidobacterium_breve","Limosilactobacillus_oris","Bifidobacterium_dentium","Trichomonas_vaginalis")
for.sort<-c("Lactobacillus_crispatus", "Lactobacillus_gasseri", "Lactobacillus_paragasseri", "Lactobacillus_iners", "Lactobacillus_jensenii", "Lactobacillus_mulieris", 
            "UBA629_sp005465875", 
            "Gardnerella_vaginalis", "Gardnerella_swidsinkii", "Gardnerella_piotii", "Gardnerella_leopoldi", 
            "Limosilactobacillus_oris", "Bifidobacterium_breve","Bifidobacterium_dentium", "Prevotella_bivia", "MultiGenera", "Trichomonas_vaginalis")

#[35] "Bifidobacterium_breve_2"   "Bifidobacterium_dentium"  
#[37] "Prevotella_bivia" 
#for (k in species){
#      if (!(k %in% for.sort)){
#	      for.sort <- append(for.sort, k)
#      }
#    }

#print(sum(is.na(mgCST.hclust)))

for (i in c(10:50)){
  for (j in c(0:4)){
    print(paste(i,j))
    
    dtc <- dynamicTreeCut::cutreeDynamic(mgCST.hclust, distM = mgCST.dist, method="hybrid", minClusterSize = i, deepSplit = j)
    print(table(dtc))
    dtc.df<-data.frame(dtc = dtc, sampleID = mgCST.hclust$labels)
    

    # Create the initial DataFrame with cluster assignments and sample IDs
    # dtc.df <- data.frame(dtc = dtc, sampleID = mgCST.hclust$labels)
    # Calculate the maximum relative abundance and the corresponding taxon
    max_relabund <- apply(relabund, 1, max)
    domTaxa <- colnames(relabund)[max.col(relabund, ties.method = "first")]
    # Identify where domTaxa is "MultiGenera"
    is_multigenera <- domTaxa == "MultiGenera"
    # For rows where domTaxa is "MultiGenera", find the second most abundant taxon
    second_max_idx <- apply(relabund, 1, function(x) order(x, decreasing = TRUE)[2])
    second_domTaxa <- colnames(relabund)[second_max_idx]
    # Replace domTaxa with the second most abundant taxon where it is "MultiGenera"
    domTaxa[is_multigenera] <- second_domTaxa[is_multigenera]
    max_relabund[is_multigenera] <- apply(relabund, 1, function(x) sort(x, decreasing = TRUE)[2])[is_multigenera]
    # Combine everything into the final DataFrame
    mgCSTs.samples <- cbind(dtc.df, domTaxa = domTaxa, relabund = max_relabund)
    
    # mgCSTs.samples<-cbind(dtc.df, domTaxa=colnames(relabund)[max.col(relabund, ties.method="first")], relabund=apply(relabund, 1, max))
    mgCSTs<-as.data.frame(mgCSTs.samples %>% group_by(dtc) %>% dplyr::summarise(meanRelabund=mean(relabund)))
    mgCSTs$domTaxa<-as.data.frame(mgCSTs.samples %>% count(dtc, domTaxa) %>% group_by(dtc) %>% slice(which.max(n)))[,2]
    mgCSTs$domTaxa<-as.character(mgCSTs$domTaxa)
    mgCSTs$minClusterSize <- i
    mgCSTs$deepSplit <- j
    
    mgCSTs$domTaxa<-as.character(mgCSTs$domTaxa)
    print(mgCSTs)
    
    v<-vector()
    
    for (l in for.sort){
      v<-append(v, sort(as.vector(mgCSTs[grepl(pattern = l, x = mgCSTs$domTaxa), "domTaxa"])))
    }
    print(paste("v",sum(is.na(v))))
    print(v)
    
    #v <- unique(v)
    #print(v)
    # mgCSTs.sort<-mgCSTs[pmatch(v, mgCSTs$domTaxa), ]
    
    mgCSTs.sort<-mgCSTs[pmatch(v, mgCSTs$domTaxa), ]
    
    print(mgCSTs.sort)
    mgCSTs.sort$mgCST <- c(1:length(table(dtc)))
    print("you are here")
    mgCSTs.sort$minClusterSize <- i
    mgCSTs.sort$deepSplit <- j
    print(mgCSTs.sort) 
    mgCSTs.samples<-merge(mgCSTs.samples, mgCSTs.sort[,c("dtc", "mgCST")], all.x=TRUE)
    mgCSTs.samples$minClusterSize <- i
    mgCSTs.samples$deepSplit <- j
    
    mgCSTs.samples.df <- rbind(mgCSTs.samples.df, mgCSTs.samples)
    mgCSTs.sort.df <- rbind(mgCSTs.sort.df, mgCSTs.sort)
    
  }
}



write_csv(mgCSTs.samples.df, "../../RDS_files/VOG_analysis_updated/vog.mgCSTs.samples.df.updated.csv")
write_csv(mgCSTs.sort.df, "../../RDS_files/VOG_analysis_updated/vog.mgCSTs.sort.df.updated.csv")

print("14_mgCST_clustering_vog.R has run")
