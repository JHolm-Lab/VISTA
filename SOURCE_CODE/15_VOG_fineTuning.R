source("00_importation.R")
require(randomForestSRC)

ngl.gene.counts.list.mgSs.pa.clean<-readRDS("../../RDS_files/VOG_analysis_updated/vog.list.mgss.pa.RDS")
samples.clusters <- readRDS("../../RDS_files/VOG_analysis_updated/samples.clusters.updated.vog.RDS")

# all.clean.100k<-readRDS("../../RDS_files/Gene_analysis/all.clean.100k.RDS")
# nReads<-colSums(all.clean.100k)
# saveRDS(nReads, "../../RDS_files/VOG_analysis_updated/nReads.RDS")

nReads <- readRDS("../../RDS_files/VOG_analysis_updated/nReads.RDS")

tuning<-list()

species <- c("Gardnerella_vaginalis", "Gardnerella_vaginalis_C", "Gardnerella_swidsinkii", "Gardnerella_piotii",
	     "Lactobacillus_crispatus", "Lactobacillus_gasseri", "Lactobacillus_jensenii", "Lactobacillus_iners")
#,"Trichomonas_vaginalis","UBA629_sp005465875", "Bifidobacterium_breve", "Lactobacillus_paragasseri")

for (i in species){
  print(i)
  the.data<-as.data.frame(t(ngl.gene.counts.list.mgSs.pa.clean[[i]]))
  print(i)
  if(i %in% c("Trichomonas_vaginalis", "Gardnerella_vaginalis", "Prevotella_bivia", "Megasphaera_genomosp.", "Lactobacillus_iners", "Lactobacillus_crispatus", "Atopobium_vaginae")){
    the.data<-as.data.frame(t(ngl.gene.counts.list.mgSs.pa.clean[[i]][,colnames(ngl.gene.counts.list.mgSs.pa.clean[[i]]) %in% names(nReads[nReads > 1000000])]))
  }
  print(i)
  sample.clusters1<-samples.clusters[[i]]
  print(i)
  the.data<-the.data[, colSums(the.data) > min(table(sample.clusters1[["sample_cluster"]]))]
  print(i)
  the.data[["mgss"]]<-as.vector(sample.clusters1[match(rownames(the.data), sample.clusters1$sampleID), "sample_cluster"])
  print(i)
  the.data<-the.data[!the.data$mgss %in% "untyped", ]
  print(i)
  the.data[["mgss"]]<-factor(the.data[["mgss"]])
  print(i)
  the.data<-the.data[!is.na(the.data$mgss), colSums(the.data[, !names(the.data) %in% "mgss"]) > 0]
  print(i)
  tuning<-tune.rfsrc(mgss~., the.data,
                            mtryStart = ncol(the.data) / 2,
                            nodesizeTry = c(1:9, seq(10, 100, by = 5)), ntreeTry = ncol(the.data)*0.1,
                            sampsize = function(x){min(x * .632, max(150, x ^ (3/4)))}, 
                            nsplit = 1, stepFactor = 1.25, improve = 1e-3, strikeout = 3, maxIter = 25,
                            trace = FALSE, doBest = TRUE)

  saveRDS(tuning, paste("../../RDS_files/VOG_analysis_updated/fineTuning_RDS/",i, "_fineTuning.vog.RDS", sep=""))
}

print("15_VOG_fine_tuning.R has run")
