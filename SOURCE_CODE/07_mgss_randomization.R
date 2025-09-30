library(fpc)
library(vegan)
library(ggplot2)
library(tidyverse)

# setwd("../mgss_randomization/")

# 1. ----- Species to consider for randomization comparison
# species <- c("Lactobacillus_crispatus", "Lactobacillus_iners", "Lactobacillus_jensenii", "Lactobacillus_gasseri", "Lactobacillus_mulieris",
#              "UBA629_sp005465875","Gardnerella_vaginalis", "Gardnerella_vaginalis_A","Gardnerella_vaginalis_C","Gardnerella_vaginalis_D","Gardnerella_vaginalis_E",
#              "Gardnerella_vaginalis_F","Gardnerella_vaginalis_H","Gardnerella_leopoldii","Gardnerella_swidsinkii", "Gardnerella_piotii", 
#              "Prevotella_amnii", "Bifidobacterium_breve","Sneathia_vaginalis", "Aerococcus_christensenii")

mgss.pa <- readRDS("../../RDS_files/VOG_analysis_updated/vog.list.mgss.pa.RDS")

mgss.pa <- mgss.pa[names(mgss.pa) %in% species]

# 2. ----- Creation of 10 randomized presence absence table for each species

percent.to.change <- 10

pa.random.10.times <- list()

for (i in names(mgss.pa)){

  ## Select data range for randomization
  data <- mgss.pa[[i]]
  data<- as.matrix(data)
  # data <- as.matrix(data[,-1])
  # rownames(data) <- data[,1]

  ## Set data dimensions
  dimensions <- dim(data)
  print(paste("Loaded matrix with ",dimensions[1]," rows and ",dimensions[2]," columns",sep=''))

  ## Convert matrix to vector
  data.vector <- c(data)
  num.elements <- length(data.vector)
  num.elements.to.change <- num.elements*(percent.to.change/100)
  min.value <- min(data.vector)
  max.value <- max(data.vector)

  print(paste("Total elements in matrix: ",num.elements,sep=''))
  print(paste("Number of elements in matrix to change: ",num.elements.to.change,sep=''))
  print(paste("Min value in matrix: ",min.value,sep = ''))
  print(paste("Max value in matrix: ",max.value,sep = ''))

  ## Sample without replacement
  index.reordered <- sample(1:num.elements,replace=FALSE)
  
  # Set the number of times to iterate
  num_iterations <- 10
  
  # Set min and max values for the random numbers
  min_value <- 0
  max_value <- 1
  
  sp_pa <- list()
  # Iterate X times
  for (iteration in 1:num_iterations) {
    # Iterate over the indices to change
    for (j in 1:num.elements.to.change) {
      index.to.change <- index.reordered[j]
      rand.value <- round(runif(1, min.value, max.value), digits = 3)
      data.vector[index.to.change] <- rand.value
      print(paste("Round: ", iteration, "; Index to change: ", index.to.change, "; Random value to add: ", rand.value, sep = ''))
    }

    # Convert vector into matrix using dimensions of original matrix
    data.final <- matrix(data.vector, ncol = dimensions[2])
    data.final[data.final > 0.5] <- 1
    data.final[data.final <= 0.5] <- 0
    colnames(data.final) <- colnames(data)
    rownames(data.final) <- rownames(data)
  
    # Add extra column to matrix with the row names
    to.write <- data.frame("Name" = rownames(data.final), data.final)
    
    sp_pa[[as.character(iteration)]] <- data.final
  }
  pa.random.10.times[[i]] <- sp_pa
}

saveRDS(pa.random.10.times, "../mgss_randomization/randomized_pa_table_all_species_10percent.RDS")


# 2 ----- Mgss clustering protocol tested on each randomized table for each species

not_to_cluster <- c("Alterileibacterium","Anaerococcus", "Bacteroides","Campylobacter","Corynebacterium", "Gardnerella","Gulosibacter", "Lactobacillus","Limosilactobacillus","MultiGenera","Porphyromonas", "Prevotella", "Streptococcus")

species <- names(mgss.pa)[!(names(mgss.pa) %in% not_to_cluster)]

clustering.results <- list()

for (sp in species) {
  
  print(sp)
  vog.list.mgss.pa <- pa.random.10.times[[sp]]
  
  #### Perform Jaccard clustering of Samples
  samples.dist<-lapply(vog.list.mgss.pa, function(x) vegdist(t(x), method="jaccard", binary = TRUE))
  samples.hc<-lapply(samples.dist, function(x) hclust(x, method="ward.D"))
  
  ### 1. Clustering with all combination of parameters
  samples.params <- list()
  samples.metrics <- list()
  samples.dtc <- list()
  
  for (i in sort(names(samples.hc))){
    
    print(paste(i, which(sort(names(samples.hc))==i)))
    
    # Empty DataFrame to store the parameters of the clustering for each subspecies
    samples.params.df <- data.frame(minClusterSize = numeric(0),
                                    deepSplit = numeric(0),
                                    AvgSilhouetteWidth = numeric(0)
    )
    
    dist <- as.matrix(samples.dist[[i]])
    samples.dtc[[i]] <- list()
    samples.metrics[[i]] <- list()
    
    # Clustering for each subspecies with variation of parameters : minClusterSize and deepSplit
    for (j in c(10:20)) {  # each value for minClusterSize (adjust the range as needed)
      samples.dtc[[i]][[j]] <- list()
      samples.metrics[[i]][[j]] <- list()
      
      for(k in c(0:4)){   # each value for deepSplit
        print(paste(i,j,k))

        # Clustering
        try <- tryCatch({
          cut_tree <- dynamicTreeCut::cutreeDynamic(samples.hc[[i]],
                                                    distM = dist,
                                                    method = "hybrid",
                                                    minClusterSize = j,
                                                    deepSplit = k)
          samples.dtc[[i]][[j]][[k+1]] <- data.frame(dtc = cut_tree, sampleID = samples.hc[[i]]$labels)
          
          # Metrics of the clustering
          stats <- cluster.stats(dist, cut_tree)
          samples.metrics[[i]][[j]][[k+1]] <- stats
          
          # Insertion of the metrics for each combination of parameters into a row of the dataframe
          row <- data.frame(
            minClusterSize = j,
            deepSplit = k,
            AverageSilhouetteWidth = stats$avg.silwidth
          )
          
          # Merging the row to the dataframe
          samples.params.df <- rbind(samples.params.df, row)
          
        }, error = function(e) {
          print(e)
        })
      }
    }
    samples.params[[i]] <- samples.params.df
  }

  ### 2. Filtering follwing this protocol :       a/ Keep the rows where Average Sil Width is maximum
  #                                               b/ Keep the maximum value of minClusterSize
  #                                               c/ Keep the min value of deepSplit
  
  params <- list()
  
  ## a / Keep the rows where the AvgSilWidth is maximum
  for (i in names(samples.params)){
    params[[i]] <- samples.params[[i]][samples.params[[i]]$AverageSilhouetteWidth == max(samples.params[[i]]$AverageSilhouetteWidth), ]
  }
  
  ## b / Keep the max value for minClusterSize
  params_min_minClusterSize <- list()
  for (i in names(params)){
    params_min_minClusterSize[[i]] <- params[[i]][params[[i]]$minClusterSize == max(params[[i]]$minClusterSize),]
  }
  
  ## c / Keep the min value for deepSplit
  best_combination_per_species <- list()
  for (i in names(params_min_minClusterSize)){
    best_combination_per_species[[i]] <- params_min_minClusterSize[[i]][params_min_minClusterSize[[i]]$deepSplit == min(params_min_minClusterSize[[i]]$deepSplit),]
  }
  

  ### 3. Look at if 10 is a possible minclustersize for each species
  
  # perform a clustering with minsize of 10 and each values of deepSplit. If a cluster is created with a size < 10 despite of the condition "minClusterSize = 10", it is corresponding to "unassigned" samples (this is allowed by DynamicTreeCut)
  # We consider those unassigned as a cluster if they are closed from each others and if they are > 10.
  # If not, we consider that the specie can't be clustered.
  
  dtc.10 <- list()
  stats.10 <- list()
  for (i in names(samples.hc)) {
    
    dtc.10[[i]] <- list()
    stats.10[[i]] <- list()
    
    for (deep.split in 0:4) {
      
      dtc.10[[i]][[deep.split+1]] <- list()
      stats.10[[i]][[deep.split+1]] <- list()
      dtc1 <- dynamicTreeCut::cutreeDynamic(samples.hc[[i]],
                                            distM = as.matrix(samples.dist[[i]]),
                                            method="hybrid",
                                            minClusterSize = 10,
                                            deepSplit = deep.split
      )
      dtc.10[[i]][[deep.split+1]] <- dtc1
      stats.10[[i]][[deep.split+1]] <- cluster.stats(as.matrix(samples.dist[[i]]), dtc1)
    }
  }

  # Look at all the cluster's size
  
  min.size <- list()
  for (i in names(stats.10)){
    min.size[[i]] <- c()
    for (k in 1:5) {
      min.size[[i]] <- c(min.size[[i]], stats.10[[i]][[k]]$min.cluster.size)
    }  
  }

  ### 4. Function Clustering
  
  # take in inputs a species and a minclustersize
  # return 3 elements : 1) The cluster vector processed with the optimal deepsplit : we keep here the maximum value of deepsplit because we increase a lot our minClusterSize
  #                     2) The value of deepsplit used
  #                     3) the cluster.stats object from the library fpc
  #                     4) the vector with all values of avg.silwidth for each deepsplit (just in case)
  
  cluster_params <- function(species, minsize, deepsplit) {
    
    stats.results <- c()
    hc <- samples.hc[[species]]
    dist <- samples.dist[[species]]
    
    for (d in c(0,1,2,3,4)) {
      cut_tree <- dynamicTreeCut::cutreeDynamic(hc, distM = as.matrix(dist), method = "hybrid", minClusterSize = minsize, deepSplit = d)
      stats <- cluster.stats(as.matrix(dist), cut_tree)
      stats.results <- c(stats.results, stats$avg.silwidth)
    }
    
    ds <- max(table(stats.results)) -1
    dtc <- dynamicTreeCut::cutreeDynamic(hc, distM = as.matrix(dist), method = "hybrid", minClusterSize = minsize, deepSplit = ds)
    stats <- cluster.stats(as.matrix(dist), dtc)
    
    result <- list("dtc"=dtc, "deepsplit" = ds, "cluster.stats" = stats, "avg.silwidth.all" = stats.results)
    
    return(result)
  }
  
  ### 5. Clutering : if 10 is actually a valid minClusterSize, cluster with find optimal parameters. If not, don't cluster
  
  samples.clusters <- list()
  samples.clusters.metrics <- list()
  mgss.clustering.parameters <- list()
  
  not_to_cluster <- c("Alterileibacterium","Anaerococcus", "Bacteroides","Campylobacter","Corynebacterium", "Gardnerella","Gulosibacter", "Lactobacillus","Limosilactobacillus","MultiGenera","Porphyromonas", "Prevotella", "Streptococcus")
  
  for (i in names(min.size)) {
    
    # If 10 is a valid  minimum cluster size for the specie
    if (min(min.size[[i]]) >= 10 && !(i %in% not_to_cluster)) {
      #if (min(min.size[[i]]) >= 10) {
      # Initial clustering with the optimal parameters determined in "2. Filtering"
      dtc <- dynamicTreeCut::cutreeDynamic(samples.hc[[i]],
                                           distM = as.matrix(samples.dist[[i]]),
                                           method="hybrid",
                                           minClusterSize = best_combination_per_species[[i]]$minClusterSize,
                                           deepSplit = best_combination_per_species[[i]]$deepSplit
      )
      
      # Transform "unassigned" samples to a cluster into one cluster (unassigned are labelled 0 by default)
      if (0 %in% dtc) { dtc <- dtc + 1 }
      
      
      # if the initial clustering return more than 10 clusters
      if (length(unique(dtc)) > 10) {
        k <- 100
        d <- best_combination_per_species[[i]]$deepSplit - 1
        
        while (k > 10 && d > 0){
          dtc <- dynamicTreeCut::cutreeDynamic(samples.hc[[i]],
                                               distM = as.matrix(samples.dist[[i]]),
                                               method="hybrid",
                                               minClusterSize = best_combination_per_species[[i]]$minClusterSize,
                                               deepSplit = d
          )
          k <- length(unique(dtc))
          if (k > 10){
            d <- d - 1
          }
          print(paste(i, d, k))
        }
        
        # # utilisation of the function "cluster_params" from part "4. Function Clustering" with input of minsize = 10% of the number of samples for the specie
        # a <- round(length(samples.hc[[i]]$labels)*0.1)
        # dtc2 <- cluster_params(i, a)
        
        # creation of 3 lists of dataframe : samples.clusters ; samples.clusters.metrics ; mgss.clustering.parameters
        samples.clusters[[i]]<-as.data.frame(cbind(sampleID=samples.hc[[i]]$labels, sample_cluster = dtc))
        samples.clusters.metrics[[i]] <- cluster.stats(as.matrix(samples.dist[[i]]), dtc)
        
        mgss.clustering.parameters[[i]] <- as.data.frame(cbind(minClusterSize = best_combination_per_species[[i]]$minClusterSize,
                                                               deepSplit = d,
                                                               AverageSilhouetteWidth = cluster.stats(as.matrix(samples.dist[[i]]), dtc)$avg.silwidth,
                                                               Number_of_cluster = length(unique(dtc)),
                                                               Number_of_unassigned = ifelse(is.na(table(samples.clusters[[i]][['sample_cluster']])["0"]) == FALSE, table(samples.clusters[[i]][['sample_cluster']])["0"], 0),
                                                               Number_of_sample = dim(samples.clusters[[i]])[1]
        )
        )
        
        # if the initial clustering return less than 10 clusters, we keep the initial clustering parameters & results
      } else {
        
        samples.clusters[[i]]<-as.data.frame(cbind(sampleID=samples.hc[[i]]$labels, sample_cluster = dtc))
        samples.clusters.metrics[[i]] <- cluster.stats(as.matrix(samples.dist[[i]]), dtc)
        
        mgss.clustering.parameters[[i]] <- as.data.frame(cbind(minClusterSize = best_combination_per_species[[i]]$minClusterSize,
                                                               deepSplit = best_combination_per_species[[i]]$deepSplit,
                                                               AverageSilhouetteWidth = best_combination_per_species[[i]]$AverageSilhouetteWidth,
                                                               Number_of_cluster = length(unique(dtc)),
                                                               Number_of_unassigned = ifelse(is.na(table(samples.clusters[[i]][['sample_cluster']])["0"]) == FALSE, table(samples.clusters[[i]][['sample_cluster']])["0"], 0),
                                                               Number_of_sample = dim(samples.clusters[[i]])[1]
        )
        )
      }
      
      # If 10 is not a valid minimum cluster size : don't cluster
    } else {
      
      no.dtc <- numeric(length(samples.hc[[i]]$labels))
      no.dtc[] <- 1
      samples.clusters[[i]]<-as.data.frame(cbind(sampleID=samples.hc[[i]]$labels, sample_cluster = no.dtc))
      samples.clusters.metrics[[i]] <- 0
      
      mgss.clustering.parameters[[i]] <- as.data.frame(cbind(minClusterSize = 0,
                                                             deepSplit = 0,
                                                             AverageSilhouetteWidth = 0,
                                                             Number_of_cluster = 1,
                                                             Number_of_unassigned = 0,
                                                             Number_of_sample = dim(samples.clusters[[i]])[1]
      )
      )
    }
  }
  
  # Store clustering stats for each species in a list
  clustering.results[[sp]] <- do.call(rbind, mgss.clustering.parameters)

}

paste("mgss_randomization_clustering_resultss_all_species_", percent.to.change,"percent.csv", sep="")

df <- do.call(rbind, clustering.results)
df$species <- rownames(df)
df$species <- sub("\\.\\d{1,2}$", "", df$species)
rownames(df) <- NULL
# write.csv(df, paste("../mgss_randomization/mgss_randomization_clustering_resultss_all_species_", percent.to.change,"percent.csv", sep=""), row.names = F)
# saveRDS(clustering.results, paste("../mgss_randomization/mgss_randomization_clustering_resultss_all_species_", percent.to.change,"percent.RDS", sep=""))
# 
# 





# Function to read and process clustering results
process_clustering_results <- function(file_path, type_label, original_taxa) {
  results <- read.csv(file_path) %>%
    rename(Taxa = species) %>%
    mutate(Type = type_label) %>%
    filter(Taxa %in% original_taxa) %>%
    select(AverageSilhouetteWidth, Taxa, Number_of_cluster, Type)
  return(results)
}

# Read in original clustering results
original_clustering_results <- read.csv("../mgss_randomization/mgss.clustering.parameters.vog.csv") %>%
  mutate(Taxa = X, Type = "Observed") %>%
  select(-X) %>%
  filter(Number_of_cluster > 1) %>%
  select(AverageSilhouetteWidth, Taxa, Number_of_cluster, Type)

# Process randomized results
mgss_randomization_clustering_results_5percent <- process_clustering_results(
  "../mgss_randomization/mgss_randomization_clustering_resultss_all_species.csv", "Randomized", original_clustering_results$Taxa
)

mgss_randomization_clustering_results_10percent <- process_clustering_results(
  "../mgss_randomization/mgss_randomization_clustering_resultss_all_species_10percent.csv", "Randomized", original_clustering_results$Taxa
)


# Final data preparation for 5% and 10%
all.data.5percent <- bind_rows(mgss_randomization_clustering_results_5percent, original_clustering_results) %>%
  left_join(select(original_clustering_results, Taxa, Number_of_cluster), by = "Taxa", suffix = c("", "_observed")) %>%
  mutate(Taxa_cluster = paste(Taxa, "k=", Number_of_cluster_observed))

all.data.10percent <- bind_rows(mgss_randomization_clustering_results_10percent, original_clustering_results) %>%
  left_join(select(original_clustering_results, Taxa, Number_of_cluster), by = "Taxa", suffix = c("", "_observed")) %>%
  mutate(Taxa_cluster = paste(Taxa, "k=", Number_of_cluster_observed))

# Function to plot the results
plot_results <- function(data, title) {
  random.pal <- c("Observed" = "#D73027", "Randomized" = "lightblue")
  cust.theme.2 <- theme_bw() + theme(text = element_text(size = 12),
                                     axis.text.x = element_text(size = 8, angle = 90, hjust = 1, face = "italic"),
                                     aspect.ratio = 5/11)
  
  ggplot(data, aes(x = Taxa_cluster, y = AverageSilhouetteWidth, fill = Type)) +
    geom_point(data = filter(data, Type == "Randomized"), 
               size = 2, shape = 21, stroke = 0.5, color = "black") +
    geom_point(data = filter(data, Type == "Observed"), 
               size = 2, shape = 21, stroke = 0.5, color = "black") +
    scale_fill_manual(name = "Data type", values = random.pal) +
    labs(x = "", y = "Avg. silhouette width", title = title) +
    cust.theme.2
}

# Plot the results for 5% and 10%
plot_results(all.data.5percent, "Results for 5% Randomization")
plot_results(all.data.10percent, "Results for 10% Randomization")





