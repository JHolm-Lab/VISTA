source("00_importation.R")

##### Date
today <- strsplit(date(), " ")
month <- today[[1]][2]
day <- today[[1]][3]
year <- today[[1]][5]
today2 <- paste(day,month,year, sep="")


## Moveme function used in the next section
moveme <- function (invec, movecommand) {
  movecommand <- lapply(strsplit(strsplit(movecommand, ";")[[1]],
                                 ",|\\s+"), function(x) x[x != ""])
  movelist <- lapply(movecommand, function(x) {
    Where <- x[which(x %in% c("before", "after", "first",
                              "last")):length(x)]
    ToMove <- setdiff(x, Where)
    list(ToMove, Where)
  })
  myVec <- invec
  for (i in seq_along(movelist)) {
    temp <- setdiff(myVec, movelist[[i]][[1]])
    A <- movelist[[i]][[2]][1]
    if (A %in% c("before", "after")) {
      ba <- movelist[[i]][[2]][2]
      if (A == "before") {
        after <- match(ba, temp) - 1
      }
      else if (A == "after") {
        after <- match(ba, temp)
      }
    }
    else if (A == "first") {
      after <- 0
    }
    else if (A == "last") {
      after <- length(myVec)
    }
    myVec <- append(temp, values = movelist[[i]][[1]], after = after)
  }
  myVec
}


ngl.gene.counts.list.mgSs.pa.clean<-readRDS("../../RDS_files/VOG_analysis_updated/vog.list.mgss.pa.RDS")
samples.clusters <- readRDS("../../RDS_files/VOG_analysis_updated/samples.clusters.updated.vog.RDS")
nReads <- readRDS("../../RDS_files/nReads.RDS")
ngl.abund.clusters.cast<-readRDS("../../RDS_files/VOG_analysis_updated/ngl.abund.clusters.cast.vog.RDS")
mgcsts.samples.df <- read_csv("../../RDS_files/VOG_analysis_updated/vog.mgCSTs.samples.df.updated.csv", show_col_types = FALSE)

deepsplit <- 4
minclustersize <- 10
original.mgCST <- mgcsts.samples.df[mgcsts.samples.df$deepSplit == 4 & mgcsts.samples.df$minClusterSize == 10, ]

relabund<-ngl.abund.clusters.cast/rowSums(ngl.abund.clusters.cast)
relabund$mgCST<-original.mgCST[match(rownames(relabund), original.mgCST$sampleID), "mgCST"]$mgCST
relabund$mgCST<-factor(relabund$mgCST)

folds <- cut(seq(1,nrow(relabund)),breaks=10,labels=FALSE)

# -1 Centroid based classifier
for(n in 1:10){
 # Segment your data by fold using the which() function
 testIndexes <- which(folds==n, arr.ind=TRUE)
 # testing has to be read counts
 testing <- ngl.abund.clusters.cast[testIndexes, ]
 # training has to be relative abundances
 training <- relabund[-testIndexes, ]
 training$SID<-rownames(training)
 the.data.m<-reshape2::melt(training, id.vars = c("SID", "mgCST"), value.name = "relabund", variable.name = "Taxon")
 centroids<-reshape2::dcast(the.data.m, mgCST~Taxon, value.var = "relabund", fun.aggregate = mean)
 centroids$mgCST<-paste("mgCST ", centroids$mgCST, sep="")
 write.csv(centroids, paste("../../RDS_files/VOG_analysis_updated/CENTROIDS/vog_mgCST_centroids_cv", n, ".csv", sep=""), row.names = F, quote=F)
 testing$mgCST<-NULL
 testing$read_count<-rowSums(testing)
 testing$SID<-rownames(testing)
 testing<-testing[moveme(names(testing), "read_count first")]
 testing<-testing[moveme(names(testing), "SID first")]
 write.csv(testing, paste("../../RDS_files/VOG_analysis_updated/CENTROIDS/vog_mgCST_testing_cv", n, ".csv", sep=""), row.names = F, quote=F)
}

print("Centroid based classifier done")

## -2 Run Python code from Valencia
# cd ../../RDS_files/VOG_analysis_updated/CENTROIDS
# terminal : for i in {1..10}; do python vog_mgCST.py -ref vog_mgCST_centroids_cv"$i".csv -i vog_mgCST_testing_cv"$i".csv -o cv_"$i"; done

n_mgCST <- 25
file.list <- list.files(path = "../../RDS_files/VOG_analysis_updated/CENTROIDS", pattern='^cv_*', full.names = T)
df.list <- lapply(file.list, function(x) read.csv(file=x))
maps<-do.call("rbind", df.list)
cv<-maps[,c("SID", "read_count", "mgCST", "score")]
cv$mgCST_pred<-gsub(pattern="mgCST ", replacement = "", cv$mgCST)
cv$mgCST_obs<-original.mgCST[match(cv$SID, original.mgCST$sampleID), "mgCST"]$mgCST
kappa_result <- cohen.kappa(x = cv[,5:6]) ## report the weighted kappa.

cv.df<-reshape2::dcast(cv, mgCST_obs~mgCST_pred, value.var = "mgCST_obs", fun.aggregate = sum)
cv.err.centroids<-as.data.frame(cbind(mgCST=1:n_mgCST, accuracy=apply(cv.df[,2:ncol(cv.df)], 1, function(x) max(x)/sum(x))))
print(mean(cv.err.centroids$accuracy))
# [1] 0.8677283


## -3 mgCST confusion plot relabund centroids
relabund<-cv.df
relabund<-relabund[order(relabund$mgCST_obs), ]
relabund<-relabund[,c("mgCST_obs", as.character(1:n_mgCST))]
relabund[,2:ncol(relabund)]<-relabund[,2:ncol(relabund)]/rowSums(relabund[,2:ncol(relabund)])
cv.df.m<-reshape2::melt(relabund, id.vars="mgCST_obs", variable.name="predicted", value.name="n")
cv.df.m$predicted<-factor(cv.df.m$predicted, levels = as.character(1:n_mgCST))
ggplot(data =  cv.df.m, mapping = aes(y = as.factor(mgCST_obs), x = as.factor(predicted))) +
  xlab("Predicted mgCST")+
  ylab("Observed mgCST")+
  geom_tile(aes(fill = n), colour = "white") +
  scale_fill_gradient(low = "white", high = "black")+
  theme_classic()
ggsave("../../FIGURES/VOG_analysis_updated/vog_mgCST_Confusion_plot_relabund_centroids.pdf", height=5, width=7)

## -4 Write out centroids
relabund<-ngl.abund.clusters.cast/rowSums(ngl.abund.clusters.cast)
relabund$SID<-rownames(relabund)
relabund$mgCST<-original.mgCST[match(rownames(relabund), original.mgCST$sampleID), "mgCST"]$mgCST
relabund$mgCST<-factor(relabund$mgCST)
the.data.m<-reshape2::melt(relabund, id.vars = c("SID", "mgCST"), value.name = "relabund", variable.name = "Taxon")
centroids<-reshape2::dcast(the.data.m, mgCST~Taxon, value.var = "relabund", fun.aggregate = mean)
centroids$mgCST<-paste("mgCST ", centroids$mgCST, sep="")
write.csv(centroids, paste("../../RDS_files/VOG_analysis_updated/CENTROIDS/vog_mgCST_centroids_", today2, ".csv", sep=""), row.names = F, quote=F)


print("17_mgCSTs_VOG_classifier.R has run")
