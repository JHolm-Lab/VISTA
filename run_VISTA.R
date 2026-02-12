#!/usr/bin/env Rscript

# Print the current date and time
current_time <- Sys.time()
print(paste("Start time:", current_time, sep = " "))

## IMPORT PACKAGES
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!require(pheatmap, quietly = TRUE, warn.conflicts = FALSE)) install.packages("pheatmap", quiet = TRUE)
if (!require(dplyr, quietly = TRUE)) install.packages("dplyr", quiet = TRUE)
if (!require(data.table, quietly = TRUE)) install.packages("data.table", quiet = TRUE)
if (!require(R.utils, quietly = TRUE)) install.packages("R.utils", quiet = TRUE)
if (!require(randomForestSRC, quietly = TRUE)) install.packages("randomForestSRC", quiet = TRUE)
# This script requires approximately 20 GB of RAM to run efficiently.

#######################################################     IMPORT PACKAGES     ###################################################################################################
library(randomForestSRC)
library(pheatmap)
library(dplyr)
library(parallel)
library(data.table)

args = commandArgs(trailingOnly=TRUE)
wd <- getwd()

#######################################################     TEST ARGUMENTS     ####################################################################################################
# Check if paths exist
if (!file.exists(args[1]) || !dir.exists(args[2])) {
  stop("Full paths to the following files/directories must be supplied and must exist:\n
       (1) VIRGO2_Compiled.summary.NR.txt (your VIRGO output file)\n
       (2) VISTA-master directory", call. = FALSE)
}

#######################################################     CAPTURE THE DATE FOR THE OUTPUT     ###################################################################################
today <- strsplit(date(), " ")
month <- today[[1]][2]
if (today[[1]][3] %in% "") {
  day <- today[[1]][4]
  year <- today[[1]][6]
} else {
  day <- today[[1]][3]
  year <- today[[1]][5]
}
today2 <- paste(day, month, year, sep="")
cat(paste(today2, "\n"))

cat(paste("Current wd:", wd, "\n"))
cat(paste("Input File:", args[1], "\n"))
cat(paste("VISTA-master directory:", args[2], "\n"))

#######################################################     READ THE MGSS.CLASSIFIER RDS FILE     #################################################################################
cat("Reading in mgSs classification models\n")
mgss.classifiers <- readRDS(paste(args[2], "/VISTA_data/volume/VISTA_files/vog.mgss_classifier.RDS", sep=""))
if (!exists("mgss.classifiers")) {
  print("Path to VISTA-master not correct")
  stop()
}

## READ IN REFERENCE CENTROIDS
cat("Reading in VISTA mgCST reference centroids\n")
reference_centroids <- as.data.frame(fread(paste(args[2], "/VISTA_data/volume/VISTA_files/vog_mgCST_centroids_25Oct2024.csv", sep="")))
rownames(reference_centroids) <- reference_centroids$mgCST
reference_centroids$mgCST <- NULL

mgCST.centroids <- as.data.frame(reference_centroids[, 2:ncol(reference_centroids)])
rownames(mgCST.centroids) <- reference_centroids$mgCST

#######################################################     READ THE VIRGO2 COMPILED OUTPUT FILE     ##############################################################################
file_path <- file.path(args[1])

# Check for gzipped version if plain text doesn't exist
if (!file.exists(file_path)) {
  gz_path <- paste0(file_path, ".gz")
  if (file.exists(gz_path)) {
    file_path <- gz_path
  } else {
    stop("Neither plain text nor gzipped version of the file was found at: ", file_path)
  }
}
counts.genes <- fread(file_path, fill=TRUE)
names(counts.genes)[1] <- "Gene"

## Added sample info to print for user upfront 20May2025
cat(paste0("From ", ncol(counts.genes)-1, " samples, ", sum(colSums(counts.genes[, -1]) > 1), " had > 1 read.\nRead distribution:\n"))
print(summary(colSums(counts.genes[, -1])))

#######################################################     READ THE VIRGO2 ANNOTATION FILES     ##################################################################################
# Read in VIRGO2 Gene length table
file_path <- file.path(args[2], "/VISTA_data/volume/virgo_data/0.VIRGO2.geneLength.txt")

# Check for gzipped version if plain text doesn't exist
if (!file.exists(file_path)) {
  gz_path <- paste0(file_path, ".gz")
  if (file.exists(gz_path)) {
    file_path <- gz_path
  } else {
    stop("Neither plain text nor gzipped version of the file was found at: ", file_path)
  }
}
gene.length <- fread(file_path)
names(gene.length) <- c("Gene", "Length")

# Read in VIRGO2 Taxon table
file_path <- file.path(args[2], "/VISTA_data/volume/virgo_data/1.VIRGO2.taxon.txt")

# Check for gzipped version if plain text doesn't exist
if (!file.exists(file_path)) {
  gz_path <- paste0(file_path, ".gz")
  if (file.exists(gz_path)) {
    file_path <- gz_path
  } else {
    stop("Neither plain text nor gzipped version of the file was found at: ", file_path)
  }
}
taxon.tbl <- fread(file_path)
names(taxon.tbl) <- c("Cluster", "Gene", "Taxa", "Cat")
taxon.tbl <- taxon.tbl[Taxa != ""]

# Read in VIRGO2 VOG table
file_path <- file.path(args[2], "/VISTA_data/volume/virgo_data/VIRGO2_VOGkey.txt")

# Check for gzipped version if plain text doesn't exist
if (!file.exists(file_path)) {
  gz_path <- paste0(file_path, ".gz")
  if (file.exists(gz_path)) {
    file_path <- gz_path
  } else {
    stop("Neither plain text nor gzipped version of the file was found at: ", file_path)
  }
}
vog.tbl <- fread(file_path)
names(vog.tbl) <- c("Gene", "Taxa", "VOG")

cat("Data and libraries importation : OK\n")

#######################################################     MERGE ANNOTATION FILES     ############################################################################################
genes <- merge(gene.length, taxon.tbl[, .(Gene, Taxa)], by='Gene', all.x=FALSE, all.y=FALSE)
genes <- merge(genes, vog.tbl[, .(Gene, VOG)], by='Gene', all.x=FALSE, all.y=FALSE)
genes <- merge(genes, counts.genes, by='Gene', all.x=FALSE, all.y=FALSE)
genes <- na.omit(genes)

#######################################################     NORMALIZATION OF COUNTS BY 150/geneLength     #########################################################################
exclude_cols <- c("Gene", "Length", "Taxa", "VOG")      # Specify the columns to exclude
a <- which(!colnames(genes) %in% exclude_cols)[1]       # Find the index of the first column that is not in the exclude list

genes.ngl <- copy(genes)
genes.ngl[, (a:ncol(genes.ngl)) := lapply(.SD, function(x) x * 150 / Length), .SDcols=a:ncol(genes.ngl)]
fwrite(genes.ngl, paste(wd, "/norm_counts_genes_", today2, ".csv", sep=""))

#######################################################     CREATE GENE PRESENCE/ABSENCE TABLE     ################################################################################
genes.ngl.pa <- copy(genes)
genes.ngl.pa[, (a:ncol(genes.ngl.pa)) := lapply(.SD, function(x) ifelse((x * 150 / Length) >= 0.5, 1, 0)), .SDcols=a:ncol(genes.ngl.pa)]

rm(genes)
#######################################################     CREATE SAMPLES BY TAXA TABLE WITH NORMALIZED COUNTS (INDEPENDANT OF GENE/VOG)     ####################################
cat("Normalizing counts by gene length ...\n")
counts.mgss <- genes.ngl[, !names(genes.ngl) %in% c("Gene", "Length", "Cat", "VOG"), with=FALSE]
counts.mgss <- counts.mgss[Taxa != ""]
counts.mgss <- counts.mgss[, lapply(.SD, sum), by=Taxa]
old.names <- colnames(counts.mgss)
counts.mgss <- data.table::transpose(counts.mgss, keep.names="Sample")
names(counts.mgss) <- c(as.character(counts.mgss[1]))
names(counts.mgss)[1] <- "Sample"
counts.mgss <- counts.mgss[-1]
counts.mgss[is.na(counts.mgss)] <- 0
counts.mgss[, (2:ncol(counts.mgss)) := lapply(.SD, as.numeric), .SDcols = 2:ncol(counts.mgss)]

counts.mgss <- as.data.frame(counts.mgss)
rownames(counts.mgss) <- counts.mgss$Sample
## Added in printing out taxa table with NGL counts, 20May2025
fwrite(counts.mgss, paste(wd, "/norm_counts_taxa_", today2, ".csv", sep=""), row.names = FALSE)
counts.mgss$Sample <- NULL

#######################################################     GROUP GENES TABLES BY VOGS TO RUN THE CLASSIFIER     ##################################################################
concat_text <- function(x) {
  paste(x, collapse=", ")
}

vog.ngl <- genes.ngl[, c(.(Gene = concat_text(Gene)),lapply(.SD, sum)), by = .(VOG, Taxa), .SDcols = where(is.numeric)]
vog.ngl.pa <- genes.ngl.pa[, c(.(Gene = concat_text(Gene)),lapply(.SD, sum)), by = .(VOG, Taxa), .SDcols = where(is.numeric)]
# Convert values > 0 to 1 from the 5th column onward -- added 20Mar2025, JBH. (Without this, the PA table contained integers > 1 representing the number of genes in a VOG present in the sampel)
vog.ngl.pa[, (5:ncol(vog.ngl.pa)) := lapply(.SD, function(x) ifelse(x > 0, 1, 0)), .SDcols = 5:ncol(vog.ngl.pa)] 

cat("Classifying mgSs ...")

#######################################################     CLASSIFY MGSS     #####################################################################################################
counts.mgss.ngl <- copy(counts.mgss)
rm(genes.ngl)
rm(genes.ngl.pa)

run_classifier <- function(taxon) {
  tryCatch({
    print(taxon)
    
    a <- which(!colnames(vog.ngl.pa) %in% exclude_cols)[1]
    
    table <- as.data.frame(t(vog.ngl.pa[vog.ngl.pa[["Taxa"]] %in% taxon, a:ncol(vog.ngl.pa)]))
    names(table) <- vog.ngl.pa[vog.ngl.pa[["Taxa"]] %in% taxon, ]$VOG
    
    ngl.sum <- rowSums(as.data.frame(t(vog.ngl[vog.ngl[["Taxa"]] %in% taxon, a:ncol(vog.ngl)])))
    
    names(counts.mgss.ngl)[names(counts.mgss.ngl) %in% taxon] <- paste(taxon, 0, sep="_")
    samples.for.0 <- rownames(table)[rowSums(table) < 500]
    
    taxon_0 <- paste(taxon, 0, sep="_")
    
    # Initialize a local data.table to store results
    local_counts <- counts.mgss.ngl
    local_counts <- as.data.frame(local_counts[[taxon_0]])
    rownames(local_counts) <- rownames(counts.mgss.ngl)
    names(local_counts) <- taxon_0
    
    if (length(samples.for.0) < nrow(table)) {
      for (gene in mgss.classifiers[[taxon]]$xvar.names) {
        if (is.null(table[[gene]])) {
          table[[gene]] <- 0
        }
      }
      
      table <- table[, mgss.classifiers[[taxon]]$xvar.names]
      predicted.mgss <- predict.rfsrc(mgss.classifiers[[taxon]], table)$predicted
      ss <- as.data.frame(cbind(mgss = colnames(predicted.mgss)[apply(predicted.mgss, 1, which.max)]), row.names = rownames(table))
      
      
      for (sample in rownames(table)) {
        new <- ifelse(sample %in% samples.for.0, paste(taxon, 0, sep="_"), paste(taxon, ss[sample, "mgss"], sep="_"))
        if (!new %in% paste(taxon, 0, sep="_")) {
          local_counts[sample, new] <- as.numeric(ngl.sum[sample])
          local_counts[sample, paste(taxon, 0, sep="_")] <- 0
        }
      }
    }
    
    for (sample in samples.for.0) {
      local_counts[sample, paste(taxon, 0, sep="_")] <- as.numeric(ngl.sum[sample])
    }
    
    # Return the local counts data.table
    return(local_counts)
    
  }, error = function(e) {
    message(sprintf("Error processing taxon %s: %s", taxon, e$message))
    return(NULL)  # Return NULL in case of error
  })
}

# Set the number of cores for parallel processing
# num_cores <- detectCores() - 1  # Use one less than the total number of cores
# Process taxa in parallel and collect results
taxa <- names(counts.mgss)[names(counts.mgss) %in% names(mgss.classifiers)]
#results.1 <- mclapply(taxa, run_classifier, mc.cores = num_cores) ## mclapply seemed to slow grid processing -- altered with below loop
library(parallel)

cat("Processing", length(taxa), " taxa. \nCombining results")
results.1 <- list()
#for (chunk in taxa_chunks) {
for (the.taxon in taxa){
  cat("Processing taxon:", paste(the.taxon, collapse = ", "), "\n")
  time_taken <- system.time({
    chunk_results <- run_classifier(the.taxon)
  })
  print(time_taken)
  results.1 <- c(results.1, chunk_results)
}

cat("Completed processing", length(results.1), " taxa. \nCombining results")
counts.mgss.ngl.1 <- do.call(cbind, results.1)

#######################################################     TRANSFER ALL NON MGSS TAXON COUNTS TO NGL     ##################################################################
#### How do we include the species that HAVE classifiers, but not enough genes present to
non_mgss_taxa <- names(counts.mgss)[!names(counts.mgss) %in% names(mgss.classifiers)] 
results.2 <- lapply(non_mgss_taxa, function(taxon) {
  tryCatch({
    # Initialize a local data.table to store results
    local_counts <- counts.mgss
    local_counts <- as.data.frame(local_counts[[taxon]])
    rownames(local_counts) <- rownames(counts.mgss.ngl)
    names(local_counts) <- taxon
    
    return(local_counts)
    
  }, error = function(e) {
    message(sprintf("Error processing taxon %s: %s", taxon, e$message))
    return(NULL)
  })
})


counts.mgss.ngl.2 <- do.call(cbind, results.2)

counts.mgss.ngl <- cbind(counts.mgss.ngl.1, counts.mgss.ngl.2)
rownames(counts.mgss.ngl) <- rownames(counts.mgss)

counts.mgss.ngl[is.na(counts.mgss.ngl)] <- 0
fwrite(counts.mgss.ngl, paste(wd, "/norm_counts_mgSs_mgCST_", today2, ".csv", sep=""), row.names = TRUE)

print("norm_counts_mgSs_mgCST_ done")

#######################################################     CENTROID CLASSIFIER     #####################################################################################################
# Defining function to determine yue-clayton theta
## A vectorized version of yue-distance significantly increases the speed. 
yue_distance <- function(row, median) {
  # Convert to numeric vectors if not already
  row <- as.numeric(row)
  median <- as.numeric(median)
  
  # Vectorized calculations
  product <- sum(row * median)
  diff_sq <- sum((row - median)^2)
  
  # Compute Yue distance
  yue_med_dist <- product / (product + diff_sq)
  return(yue_med_dist)
}


## MAKE RELABUND TABLE FRESH
relabund<-counts.mgss.ngl/rowSums(counts.mgss.ngl)

## REFORMAT RELABUND TO INCLUDE ALL EXPECTED COLUMN NAMES (xvar.names)
cat("Reformatting relative abundance table to include all expected column names\n")
relabund<-relabund[, names(relabund) %in% names(reference_centroids)]
for(taxa in names(reference_centroids)){
  if(is.null(relabund[[taxa]])){
    relabund[[taxa]]<-0
  }
}
print("Dimension of relabund:")
relabund<-relabund[,names(reference_centroids)]
n<-ncol(relabund)
print(dim(relabund))

## FOR EACH MGCST, MEASURE THE SIMILARITY OF EACH SAMPLE TO EACH MGCST CENTROID USING YUE + CLAYTON THETA
## A vectorized version of run_yue_distance significantly increases the speed. 
run_yue_distance <- function(i) {
  tryCatch({
    cat("Running run_yue_distance\n")
    mgCST <- paste("mgCST", i, sep = " ")
    
    # Extract reference vector once
    ref_vec <- as.numeric(reference_centroids[mgCST, ])
    
    # Convert relabund to matrix if not already
    relabund_matrix <- as.matrix(relabund[, 1:n])
    
    # Vectorized apply using vapply for speed and type safety
    result <- vapply(
      X = seq_len(nrow(relabund_matrix)),
      FUN = function(j) yue_distance(relabund_matrix[j, ], ref_vec),
      FUN.VALUE = numeric(1)
    )
    
    relabund[[mgCST]] <- result
    return(result)
  }, error = function(e) {
    message(sprintf("Error processing taxon %s: %s", mgCST, e$message))
    return(NULL)
  })
}


i <- 1:25

# Use mclapply for parallel processing
#results.yue.ditance <- mclapply(i, run_yue_distance, mc.cores = num_cores)
## BOTTLENECKS OF PARALLELIZATION HERE:
# 1. Memory Usage: Each parallel process may duplicate large objects like relabund or reference_centroids, increasing memory demand.
# 2. Overhead: For small or fast tasks, the overhead of managing parallel processes may outweigh the speedup.
# 3. SLURM Configuration: If your SLURM job doesn't request enough CPUs (e.g., via --cpus-per-task), mclapply may not be able to use the specified number of cores.
# saveRDS(results.yue.ditance, paste(wd, "/results.yue.distance.RDS", sep=""))
## Given this and updates to vectorized functions, there may be no need for multiple cores.
cat("Assigning VISTA mgCSTs\n")
results.yue.ditance <- lapply(i, run_yue_distance)

# Combine the list of numeric vectors into a single data frame
relabund.yue.distance <- as.data.frame(t(do.call(rbind, lapply(results.yue.ditance, function(x) as.data.frame(t(x))))))

# Rename the columns to reflect the mgCST names
rownames(relabund.yue.distance) <- rownames(relabund)
colnames(relabund.yue.distance) <- paste("mgCST", 1:length(results.yue.ditance), sep = " ")

m<-n+1
relabund <- cbind(relabund, relabund.yue.distance)
relabund[is.na(relabund)] <- 0
relabund[["mgCST"]]<-colnames(relabund[,m:which(colnames(relabund) %in% "mgCST 25")])[apply(relabund[,m:which(colnames(relabund) %in% "mgCST 25")],1,which.max)]

write.csv(relabund, paste(wd, "/relabund_w_mgCSTs_", today2, ".csv", sep=""), row.names = TRUE, quote=F)
write.csv(cbind(mgCST=relabund["mgCST"], max_YC_theta=apply(relabund[(ncol(relabund)-1):(ncol(relabund)-25)], 1, max)), paste(wd, "/mgCSTs_", today2, ".csv", sep=""), row.names = TRUE, quote=F)


## PLOT HEATMAP
cat("Making heatmap\n")
mgCST<-as.data.frame(rbind(c("1", "#FE0308"), c("2", "#F54C5E"), c("3", "#F07084"), c("4", "#EC94A5"),c("5", "#F0BCCC"),c("6", "#F6D3DA"),
                           c("7", "#86C61A"), c("8", "#B4DB29"), 
                           c("9", "#F68A11"), c("10", "#FF981C"),c("11", "#FFA435"),
                           c("12", "#FAE727"),c("13", "#FBEA3F"), c("14", "#FBED58"),
                           c("15", "#E1C775"),
                           c("16", "#589682"),c("17", "#6BA290"),
                           c("18", "#2C31A0"),c("19", "#3C44A8"),c("20", "#444DAC"),c("21", "#676EBC"),c("22", "#6B7EC0"),c("23", "#829CCD"),
                           c("24", "#C7FFC7"), c("25", "#8c8c8c"), c("", "white"), c("NA", "white")))

names(mgCST)<-c("mgCST", "color")
colfunc <- colorRampPalette(c("khaki", "limegreen", "darkslategray1", "mediumblue", "magenta", "red"))
relabund.mgCST<-relabund[,1:n]
relabund.mgCST<-relabund.mgCST[,order(colSums(relabund.mgCST), decreasing = TRUE)]
relabund.mgCST$mgCST<-gsub("mgCST ", "", relabund$mgCST)
relabund.mgCST<-relabund.mgCST[order(as.numeric(relabund.mgCST[["mgCST"]])),]
relabund.mgCST[["color"]]<-mgCST[match(relabund.mgCST[["mgCST"]], mgCST$mgCST), "color"]
names(relabund.mgCST)<-gsub("_", " ", names(relabund.mgCST))

# Prepare matrix and annotations
mat <- t(as.matrix(relabund.mgCST[, 1:50]))
annotation_col <- data.frame(mgCST = relabund.mgCST[["mgCST"]])
rownames(annotation_col) <- rownames(relabund.mgCST)
annotation_colors <- list(mgCST = setNames(mgCST$color[mgCST$mgCST %in% annotation_col$mgCST], mgCST$mgCST[mgCST$mgCST %in% annotation_col$mgCST]))

# Save to PDF
pdf(paste0(wd, "/mgCST_heatmap_", today2, ".pdf"), width = 7, height = 10)
pheatmap(
  mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colfunc(100),
  main = paste("mgCST Heatmap\nnSamples =", nrow(relabund)),
  fontsize_row = 7,
  fontsize_col = 1,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  legend = TRUE,
  border_color = NA
)
dev.off()

print(paste("mgCST_heatmap_", today2, ".csv has been saved", sep=""))


current_time <- Sys.time()
print(paste("End time:", current_time, sep = " "))

