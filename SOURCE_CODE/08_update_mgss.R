## Libraries
library(fpc)
library(ggplot2)
library(tidyverse)



## Import data
original_clustering_results <- read.csv("../mgss_randomization/mgss.clustering.parameters.vog.csv")
randomized_clustering_results_5percent <- read.csv("../mgss_randomization/mgss_randomization_clustering_resultss_all_species.csv")
randomized_clustering_results_10percent <- read.csv("../mgss_randomization/mgss_randomization_clustering_resultss_all_species_10percent.csv")
randomized.results.5percent <- readRDS("../mgss_randomization/mgss_randomization_clustering_results_all_species.RDS")
randomized.results.10percent <- readRDS("../mgss_randomization/mgss_randomization_clustering_resultss_all_species_10percent.RDS")


samples.hc <- readRDS("../../RDS_files/VOG_analysis/samples.hc.vog.RDS")
samples.clusters <- readRDS("../../RDS_files/VOG_analysis/samples.clusters.vog.RDS")
samples.dist <- readRDS("../../RDS_files/VOG_analysis/samples.dist.vog.RDS")


randomized.results <- randomized.results.10percent
## Filter data : taxa where original results is not greater than all randomized results
names(original_clustering_results)[names(original_clustering_results) == "X"] <- "Taxa"
taxa.to.redo.clustering <- c()
for (i in names(randomized.results)){
  randomized.results.df <- as.data.frame(randomized.results[[i]])
  max.avg.silwidth.random <- max(randomized.results.df$AverageSilhouetteWidth, na.rm = TRUE)
  original.results.df <- original_clustering_results[original_clustering_results$Taxa == i,]
  avg.silwidth.original <- max(original.results.df$AverageSilhouetteWidth, na.rm = TRUE)
  print(paste(i, max.avg.silwidth.random, avg.silwidth.original))
  if (avg.silwidth.original < max.avg.silwidth.random){
    taxa.to.redo.clustering <- c(taxa.to.redo.clustering, i)
  }
}


## Filter original_clustering_results by keeping only taxa.to.redo.clustering with original k>1 
df <- original_clustering_results[original_clustering_results$Taxa %in% taxa.to.redo.clustering, ]
df <- df[df$Number_of_cluster > 1, ]

samples.clusters.to.update <- samples.clusters[names(samples.clusters) %in% df$Taxa]
updated_clustering_results <- original_clustering_results

k.cluster <- list()
for (i in df$Taxa){
  print(i)
  # Initial setup
  data <- df[df$Taxa == i,]
  samples.hc.taxa <- samples.hc[[i]]
  samples.dist.taxa <- samples.dist[[i]]
  
  k.original <- data$Number_of_cluster
  avg.silwidth.original <- data$AverageSilhouetteWidth
  
  max.avg.silwidth.random <- max(randomized.results[[i]]$AverageSilhouetteWidth, na.rm = TRUE)
  
  # Initialize variables for the while loop
  k.current <- k.original
  avg.silwidth.new <- -2  # Start with a value smaller than any possible silhouette width
  
  # Perform clustering and check silhouette width
  while (avg.silwidth.new < max.avg.silwidth.random && k.current > 1) {
    # Perform clustering
    dtc.taxa <- cutree(samples.hc.taxa, k = k.current)
    
    # Compute clustering stats
    clustering_stats <- cluster.stats(as.matrix(samples.dist.taxa), dtc.taxa)
    avg.silwidth.new <- clustering_stats$avg.silwidth
    
    
    # Check if we need to redo clustering
    if (avg.silwidth.new < max.avg.silwidth.random){
      print(paste("Redo clustering with k =", k.current))
      k.current <- k.current - 1
      
    } else if (avg.silwidth.new > max.avg.silwidth.random && min(table(dtc.taxa)) >= 10) {
      print("Clustering ok")
      k.cluster[[i]] <- k.current
      samples.clusters.to.update[[i]] <- as.data.frame(cbind(sampleID=samples.hc.taxa$labels, sample_cluster = dtc.taxa))
      clustering_stats <- cluster.stats(as.matrix(samples.dist.taxa), dtc.taxa)
      # Modify the values for that row
      row_index <- which(updated_clustering_results$Taxa == i)
      updated_clustering_results[row_index, "minClusterSize"] <- NA
      updated_clustering_results[row_index, "deepSplit"] <- NA
      updated_clustering_results[row_index, "AverageSilhouetteWidth"] <- clustering_stats$avg.silwidth
      updated_clustering_results[row_index, "Number_of_cluster"] <- k.current
      updated_clustering_results[row_index, "Number_of_unassigned"] <- NA
      
    } else if (avg.silwidth.new > max.avg.silwidth.random && min(table(dtc.taxa)) < 10) {
      k.cluster[[i]] <- k.current - 1
      dtc.taxa <- cutree(samples.hc.taxa, k = k.current - 1)
      samples.clusters.to.update[[i]] <- as.data.frame(cbind(sampleID=samples.hc.taxa$labels, sample_cluster = dtc.taxa))
      clustering_stats <- cluster.stats(as.matrix(samples.dist.taxa), dtc.taxa)
      # Modify the values for that row
      row_index <- which(updated_clustering_results$Taxa == i)
      updated_clustering_results[row_index, "minClusterSize"] <- NA
      updated_clustering_results[row_index, "deepSplit"] <- NA
      updated_clustering_results[row_index, "AverageSilhouetteWidth"] <- clustering_stats$avg.silwidth
      updated_clustering_results[row_index, "Number_of_cluster"] <- k.current - 1
      updated_clustering_results[row_index, "Number_of_unassigned"] <- NA
    }
  }
  
  # If loop exits and k.current has become too small
  if (k.current <= 1) {
    print("Minimum k reached. Final clustering done.")
    k.cluster[[i]] <- 1
    samples.clusters.to.update[[i]] <- as.data.frame(cbind(sampleID=samples.hc.taxa$labels, sample_cluster = 1))
    row_index <- which(updated_clustering_results$Taxa == i)
    updated_clustering_results[row_index, "minClusterSize"] <- NA
    updated_clustering_results[row_index, "deepSplit"] <- NA
    updated_clustering_results[row_index, "AverageSilhouetteWidth"] <- 0
    updated_clustering_results[row_index, "Number_of_cluster"] <- 1
    updated_clustering_results[row_index, "Number_of_unassigned"] <- NA
  }
  
}

## k.cluster.df contains Taxa and k to use
k.cluster.df <- as.data.frame(do.call(rbind, k.cluster))
k.cluster.df$Taxa <- rownames(k.cluster.df)
rownames(k.cluster.df) <- NULL
names(k.cluster.df) <- c("k","Taxa")

## Update samples.clusters
for (i in names(samples.clusters.to.update)){
  samples.clusters[[i]] <- samples.clusters.to.update[[i]]
}

# saveRDS(samples.clusters, "../../RDS_files/VOG_analysis_updated/samples.clusters.updated.vog.RDS")

# write.csv(updated_clustering_results, "../mgss_randomization/updated_clustering_results_w_while.csv", row.names = FALSE)


## Updates visualization
updated.results.df <- read.csv("../mgss_randomization/updated_clustering_results_w_while.csv")
updated.results.df <- updated.results.df[,c("Taxa", "AverageSilhouetteWidth", "Number_of_cluster")]
updated.results.df$Type <- "Observed"

randomized.results.df <- read.csv("../mgss_randomization/mgss_randomization_clustering_resultss_all_species.csv")
names(randomized.results.df)[names(randomized.results.df) == "species"] <- "Taxa"
randomized.results.df <- randomized.results.df[,c("Taxa", "AverageSilhouetteWidth", "Number_of_cluster")]
randomized.results.df$Type <- "Randomized"

# Filter data considering only mgSs (k > 1)
original.mgss <- updated.results.df[updated.results.df$Number_of_cluster > 1,]
randomized.mgss <- randomized.results.df[randomized.results.df$Taxa %in% original.mgss$Taxa,]

# combine orginial and randomized data
all.data <- rbind(original.mgss, randomized.mgss)

# Fix Taxa names and data type
all.data$Taxa <- gsub("_", " ", all.data$Taxa)
all.data$Taxa <- gsub("(^[A-Za-z])[a-z]+\\s([a-z]+)", "\\1. \\2", all.data$Taxa)

# Create a lookup table for "Observed" rows
observed_data <- all.data %>%
  filter(Type == "Observed") %>%
  select(Taxa, Number_of_cluster)

# Left join the observed data with the original data to get the Number_of_cluster
all.data <- all.data %>%
  left_join(observed_data, by = "Taxa", suffix = c("", "_observed")) %>%
  mutate(Taxa_cluster = paste(Taxa, ", k = ", Number_of_cluster_observed))

# Set plot characteristics
random.pal <- c("Observed" = "#D73027", "Randomized" = "lightblue")
cust.theme.2 <- theme_bw() + theme(text = element_text(size = 12),
                                   axis.text.x = element_text(size = 8, 
                                                              angle = 90, 
                                                              hjust = 1, 
                                                              face = "italic"),
                                   aspect.ratio = 5/11)

# Plot all results
ggplot(all.data, aes(x = Taxa_cluster, y = AverageSilhouetteWidth, fill = Type)) +
  geom_point(data = filter(all.data, Type == "Randomized"), 
             size = 2,        # Adjust point size as needed
             shape = 21,      # Shape 21 allows for border color
             stroke = 0.5,    # Width of the border
             color = "black") +  # Border color
  geom_point(data = filter(all.data, Type == "Observed"), 
             size = 2,        # Adjust point size as needed
             shape = 21,      # Shape 21 allows for border color
             stroke = 0.5,    # Width of the border
             color = "black") +  # Border color
  scale_fill_manual(name = "Data type", values = random.pal) +
  labs(x = "", y = "Avg. silhouette width") +
  cust.theme.2






# uba <- all.data[all.data$Taxa == "UBA629 sp005465875",]
# a <- uba[uba$Type == "Observed",]$AverageSilhouetteWidth
# b <- max(uba[uba$Type == "Randomized",]$AverageSilhouetteWidth)
# ratio <- 100-100*b/a
# # [1] 0.5462126


