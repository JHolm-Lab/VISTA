source("00_importation.R")

##### Date
today <- strsplit(date(), " ")
month <- today[[1]][2]
day <- today[[1]][3]
year <- today[[1]][5]
today2 <- paste(day,month,year, sep="")

##### ---1. Data importation
ngl.gene.counts.list.mgSs.pa.clean<-readRDS("../../RDS_files/VOG_analysis_updated/vog.list.mgss.pa.RDS")
samples.clusters <- readRDS("../../RDS_files/VOG_analysis_updated/samples.clusters.updated.vog.RDS")
nReads <- readRDS("../../RDS_files/VOG_analysis_updated/nReads.RDS")
ngl.abund.clusters.cast<-readRDS("../../RDS_files/VOG_analysis_updated/ngl.abund.clusters.cast.vog.RDS")
mgcsts.samples.df <- read_csv("../../RDS_files/VOG_analysis_updated/vog.mgCSTs.samples.df.updated.csv", show_col_types = FALSE)

deepsplit <- 4
minclustersize <- 10
original.mgCST <- mgcsts.samples.df[mgcsts.samples.df$deepSplit == 4 & mgcsts.samples.df$minClusterSize == 10, ]


##### ---2. Build table w/tuning parameters for major taxa --- see 01_fineTuning.R

## -2.1 List of species with tuning parameters - Limosilobactocillus_oris not considered because not in ngl.gene.counts.list.mgSs.pa.clean
#species <- c("Gardnerella_vaginalis", "Gardnerella_vaginalis_C", "Gardnerella_swidsinkii", "Gardnerella_piotii",
#             "Lactobacillus_crispatus", "Lactobacillus_gasseri", "Lactobacillus_jensenii", "Lactobacillus_iners", "Lactobacillus_mulieris",
#             "UBA629_sp005465875", "Bifidobacterium_breve", "Lactobacillus_paragasseri")

species <- c("Gardnerella_vaginalis", "Gardnerella_vaginalis_C", "Gardnerella_swidsinkii", "Gardnerella_piotii",
             "Lactobacillus_crispatus", "Lactobacillus_gasseri", "Lactobacillus_jensenii", "Lactobacillus_iners")

df<- list()
for (i in species) {
  df2 <- data.frame((t(as.data.frame(readRDS(paste("../../RDS_files/VOG_analysis_updated/fineTuning_RDS/", i, "_fineTuning.vog.RDS", sep=""))$optimal))))
  df2$Taxa <- i
  df[[i]] <- df2
}

optimal.tbl <- do.call(rbind, df)
write.csv(optimal.tbl, "../../RDS_files/VOG_analysis_updated/vog_optimal_tbl.csv")


##### ---3. Classifier for Mgss

## -3.1 Get species that are not clustered
one_cluster_species <- c()
for (i in names(samples.clusters)) {
  if (length(table(samples.clusters[[i]][["sample_cluster"]])) == 1){
    one_cluster_species <- c(one_cluster_species, i)
  }
}

classifiers<-list()
predicted<-list()
set.seed(123)
mgss.cv<-data.frame()

## -3.2 Perform 10 fold cross validation for each mgSs.
for(i in names(samples.clusters)[!names(samples.clusters) %in% one_cluster_species]){
  print(i)

  the.data<-as.data.frame(t(ngl.gene.counts.list.mgSs.pa.clean[[i]]))
  if(i %in% c("Gardnerella_vaginalis", "Prevotella_bivia", "Megasphaera_genomosp.", "Lactobacillus_iners", "Lactobacillus_crispatus", "Atopobium_vaginae")){
    the.data<-as.data.frame(t(ngl.gene.counts.list.mgSs.pa.clean[[i]][,colnames(ngl.gene.counts.list.mgSs.pa.clean[[i]]) %in% names(nReads[nReads > 5.5e5])]))
  }
  sample.clusters1<-samples.clusters[[i]]
  the.data<-the.data[, colSums(the.data) > min(table(sample.clusters1[["sample_cluster"]]))]
  the.data[["mgss"]]<-as.vector(sample.clusters1[match(rownames(the.data), sample.clusters1$sampleID), "sample_cluster"])
  the.data[["mgss"]]<-factor(the.data[["mgss"]])
  the.data<-the.data[!is.na(the.data$mgss), colSums(the.data[, !names(the.data) %in% "mgss"]) > 0]
  the.node<-as.numeric(optimal.tbl[optimal.tbl$Taxa %in% i, "nodesize"])
  the.mtry<-as.numeric(optimal.tbl[optimal.tbl$Taxa %in% i, "mtry"])
  #Randomly shuffle the data
  the.data<-the.data[sample(nrow(the.data)),]
  #Create 10 equally size folds
  folds <- cut(seq(1,nrow(the.data)),breaks=10,labels=FALSE)
  for(n in 1:10){
    print("Enter 2nd For Loop ok")
    #Segement your data by fold using the which() function
    testIndexes <- which(folds==n, arr.ind=TRUE)
    testing <- the.data[testIndexes, ]
    training <- the.data[-testIndexes, ]
    #Use the test and train data partitions however you desire..
    if(length(the.mtry) > 0){
      classifiers[[i]]<- rfsrc.fast(mgss~., training, ntree = ncol(the.data)*0.5, forest = TRUE, nodes=the.node, mtry=the.mtry)
    }else{
      classifiers[[i]]<- rfsrc.fast(mgss~., training, ntree = ncol(the.data)*0.5, forest = TRUE)
    }
    predicted[[i]]<-predict.rfsrc(classifiers[[i]], testing)
    mgss.cv<-rbind(mgss.cv, cbind(Taxon=i, CV=n, err=predicted[[i]]$err.rate[nrow(predicted[[i]]$err.rate)]))
  }
}

mgss.cv$Taxon<-gsub("_", " ", mgss.cv$Taxon)


saveRDS(classifiers, "../../RDS_files/VOG_analysis_updated/vog.classifier.RDS")
saveRDS(predicted, "../../RDS_files/VOG_analysis_updated/vog.predicted.RDS")
write_csv(mgss.cv, "../../RDS_files/VOG_analysis_updated/vog.mgss.cv.csv")

mgss.cv <- read_csv("../../RDS_files/VOG_analysis_updated/vog.mgss.cv.csv", show_col_types = FALSE)
classifiers <- readRDS("../../RDS_files/VOG_analysis_updated/vog.classifier.RDS")
predicted <- readRDS("../../RDS_files/VOG_analysis_updated/vog.predicted.RDS")

## -3.3 Out of Bag Error Rate representation
ggplot(mgss.cv, aes(y=reorder(Taxon, as.numeric(err)), x=as.numeric(err)))+
  geom_boxplot(notch=FALSE, outlier.shape = NA, lwd=0.1, fill="gray")+
  geom_jitter(width=0.001, size=0.01)+
  xlab("Out of Bag Error Rate")+
  theme_bw()+
  ylab("")+
  theme(text=element_text(size=5))

ggsave("../../FIGURES/VOG_analysis_updated/vog_mgSs_oob_err.pdf", height=7, width=6)

## -3.4 Full classifier for use
mgss.classifier<-list()

for(i in names(samples.clusters)[!names(samples.clusters) %in% one_cluster_species]){
  the.data<-as.data.frame(t(ngl.gene.counts.list.mgSs.pa.clean[[i]]))
  if(i %in% c("Gardnerella_vaginalis", "Lactobacillus_iners", "Atopobium_vaginae")){
    the.data<-as.data.frame(t(ngl.gene.counts.list.mgSs.pa.clean[[i]][,colnames(ngl.gene.counts.list.mgSs.pa.clean[[i]]) %in% names(nReads[nReads > 5.5e5])]))
  }
  sample.clusters1<-samples.clusters[[i]]
  the.data<-the.data[, colSums(the.data) > min(table(sample.clusters1[["sample_cluster"]]))]
  the.data[["mgss"]]<-as.vector(sample.clusters1[match(rownames(the.data), sample.clusters1$sampleID), "sample_cluster"])
  the.data<-the.data[!the.data$mgss %in% "untyped", ]
  the.data[["mgss"]]<-factor(the.data[["mgss"]])
  the.mtry<-as.numeric(optimal.tbl[optimal.tbl$Taxa %in% i, "mtry"])
  the.node<-as.numeric(optimal.tbl[optimal.tbl$Taxa %in% i, "nodesize"])
  the.data<-the.data[!is.na(the.data$mgss), colSums(the.data[, !names(the.data) %in% "mgss"]) > 0]
  if(length(the.mtry) > 0){
    mgss.classifier[[i]] <- rfsrc.fast(mgss~., the.data, forest = TRUE, nodes=the.node, mtry=the.mtry)
  }else{
    mgss.classifier[[i]] <- rfsrc.fast(mgss~., the.data, forest = TRUE)
  }
}

saveRDS(mgss.classifier, "../../RDS_files/VOG_analysis_updated/vog.mgss_classifier.RDS")

print("16_mgss_VOG_classifier.R has run")
