#!/usr/bin/env Rscript

# Print the current date and time
current_time <- Sys.time()
print(paste("Start time:", current_time, sep = " "))

## IMPORT PACKAGES
if (!require(randomForestSRC, quietly = TRUE)) install.packages("randomForestSRC", repos = "https://cran.us.r-project.org")
if (!require(gplots, quietly = TRUE, warn.conflicts = FALSE)) install.packages("gplots", quiet = TRUE)
if (!require(dplyr, quietly = TRUE)) install.packages("dplyr", quiet = TRUE)
if (!require(data.table, quietly = TRUE)) install.packages("data.table", quiet = TRUE)

### Added parallel processing (mc.cores) + created functions to speed up processing time
##### Requires submitting script with a 4th argument, which is the number of cores/threads



#######################################################     IMPORT PACKAGES     ###################################################################################################

library(randomForestSRC)
library(gplots)
library(dplyr)
library(parallel)
library(data.table)

args = commandArgs(trailingOnly=TRUE)
wd <- getwd()

#######################################################     TEST ARGUMENTS     ####################################################################################################
if (length(args) < 4) {
  stop("Full paths to the following files/directories must be supplied in this order:\n(1) VIRGO2_Compiled.summary.NR.txt\n(2) VIRGO2-master directory\n(3) mgCST-classifier-master directory\n(4) Number of cores", call.=FALSE)
}

num_cores <- as.integer(args[4])
if (is.na(num_cores) || num_cores < 1) {
  stop("Number of cores must be a positive integer", call.=FALSE)
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
print(today2)

#######################################################     READ THE MGSS.CLASSIFIER RDS FILE     #################################################################################

mgss.classifiers <- readRDS(paste(args[3], "/vog.mgss_classifier.RDS", sep=""))
if (!exists("mgss.classifiers")) {
  print("Path to mgCST-classifier-master not correct")
  stop()
}


#######################################################     READ THE VIRGO2 COMPILED OUTPUT FILE     ##############################################################################

counts.genes <- fread(args[1])
names(counts.genes)[1] <- "Gene"

#######################################################     READ THE VIRGO2 ANNOTATION FILES     ##################################################################################

gene.length <- fread(paste(args[2], "/0.VIRGO2.geneLength.txt.gz", sep=""))
names(gene.length) <- c("Gene", "Length")

taxon.tbl <- fread(paste(args[2], "/1.VIRGO2.taxon.txt.gz", sep=""))
names(taxon.tbl) <- c("Cluster", "Gene", "Taxa", "Cat")
taxon.tbl <- taxon.tbl[Taxa != ""]

vog.tbl <- fread(paste(args[2], "/VIRGO2_VOGkey.txt", sep=""))
names(vog.tbl) <- c("Gene", "Taxa", "VOG")

print("Data and libraries importation : OK")

#######################################################     MERGE ANNOTATION FILES     ############################################################################################

genes <- merge(gene.length, taxon.tbl[, .(Gene, Taxa)], by='Gene', all.x=FALSE, all.y=FALSE)
genes <- merge(genes, vog.tbl[, .(Gene, VOG)], by='Gene', all.x=FALSE, all.y=FALSE)
genes <- merge(genes, counts.genes, by='Gene', all.x=FALSE, all.y=FALSE)
genes <- na.omit(genes)

print("Dimension of the final count table:")
print(dim(genes))

#######################################################     NORMALIZATION OF COUNTS BY 150/geneLength     #########################################################################

exclude_cols <- c("Gene", "Length", "Taxa", "VOG")      # Specify the columns to exclude
a <- which(!colnames(genes) %in% exclude_cols)[1]       # Find the index of the first column that is not in the exclude list

genes.ngl <- copy(genes)
genes.ngl[, (a:ncol(genes.ngl)) := lapply(.SD, function(x) x * 150 / Length), .SDcols=a:ncol(genes.ngl)]
# fwrite(genes.ngl, paste(wd, "/norm_counts_genes_", today2, ".csv", sep=""))

#######################################################     CREATE GENE PRESENCE/ABSENCE TABLE     ################################################################################

genes.ngl.pa <- copy(genes)
genes.ngl.pa[, (a:ncol(genes.ngl.pa)) := lapply(.SD, function(x) ifelse(x * 150 / Length >= 0.5, 1, 0)), .SDcols=a:ncol(genes.ngl.pa)]


#######################################################     CREATE SAMPLES BY TAXA TABLE WITH NORMALIZED COUNTS (INDEPENDANT OF GENE/VOG)     ####################################

print("Process counts.mgss ...")
counts.mgss <- genes.ngl[, !names(genes.ngl) %in% c("Gene", "Length", "Cat", "VOG"), with=FALSE]
counts.mgss <- counts.mgss[Taxa != ""]
counts.mgss <- counts.mgss[, lapply(.SD, sum), by=Taxa]
counts.mgss <- transpose(counts.mgss, keep.names="Sample")
old.names <- colnames(counts.mgss)
setnames(counts.mgss, old=old.names, new=as.character(counts.mgss[1]))
setnames(counts.mgss, old = names(counts.mgss)[1], new = "Sample")
counts.mgss <- counts.mgss[-1]
counts.mgss[is.na(counts.mgss)] <- 0
counts.mgss[, (2:ncol(counts.mgss)) := lapply(.SD, as.numeric), .SDcols = 2:ncol(counts.mgss)]

counts.mgss <- as.data.frame(counts.mgss)
rownames(counts.mgss) <- counts.mgss$Sample
counts.mgss$Sample <- NULL

#######################################################     GROUP GENES TABLES BY VOGS TO RUN THE CLASSIFIER     ##################################################################

concat_text <- function(x) {
  paste(x, collapse=", ")
}


vog.ngl <- genes.ngl[, c(.(Gene = concat_text(Gene)),lapply(.SD, sum)), by = .(VOG, Taxa), .SDcols = where(is.numeric)]
vog.ngl.pa <- genes.ngl.pa[, c(.(Gene = concat_text(Gene)),lapply(.SD, sum)), by = .(VOG, Taxa), .SDcols = where(is.numeric)]


print("Classify MGSS ...")

#######################################################     CLASSIFY MGSS     #####################################################################################################

counts.mgss.ngl <- copy(counts.mgss)
genes.ngl <- vog.ngl
genes.ngl.pa <- vog.ngl.pa

# taxon <- "Gardnerella_vaginalis_A"

run_classifier <- function(taxon) {
  tryCatch({
    print(taxon)
    
    a <- which(!colnames(genes.ngl.pa) %in% exclude_cols)[1]
    
    table <- as.data.frame(t(genes.ngl.pa[genes.ngl.pa[["Taxa"]] %in% taxon, a:ncol(genes.ngl.pa)]))
    names(table) <- genes.ngl.pa[genes.ngl.pa[["Taxa"]] %in% taxon, ]$VOG
    
    ngl.sum <- rowSums(as.data.frame(t(genes.ngl[genes.ngl[["Taxa"]] %in% taxon, a:ncol(genes.ngl)])))
    
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
results.1 <- mclapply(taxa, run_classifier, mc.cores = num_cores)
counts.mgss.ngl.1 <- do.call(cbind, results.1)


#######################################################     TRANSFER ALL NON MGSS TAXON COUNTS TO NGL     ##################################################################

non_mgss_taxa <- names(counts.mgss)[!names(counts.mgss) %in% names(mgss.classifiers)]

results.2 <- mclapply(non_mgss_taxa,
                      
                      function(taxon) {
                        tryCatch({
                          # Initialize a local data.table to store results
                          local_counts <- counts.mgss
                          local_counts <- as.data.frame(local_counts[[taxon]])
                          rownames(local_counts) <- rownames(counts.mgss.ngl)
                          names(local_counts) <- taxon
                          
                          a <- which(!colnames(genes.ngl) %in% exclude_cols)[1]
                          ngl.sum <- rowSums(transpose(genes.ngl[Taxa == taxon, .SD, .SDcols=a:ncol(genes.ngl)]))
                          
                          a <- which(!colnames(genes.ngl.pa) %in% exclude_cols)[1]
                          # table <- transpose(genes.ngl.pa[Taxa == taxon, .SD, .SDcols=a:ncol(genes.ngl.pa)])
                          table<-as.data.frame(t(genes.ngl.pa[genes.ngl.pa[["Taxa"]] %in% taxon, a:ncol(genes.ngl.pa)]))
                          
                          for (sample in rownames(table)) {
                            local_counts[sample, taxon] <- as.numeric(ngl.sum[sample])
                            # counts.mgss.ngl[sample, (taxon) := as.numeric(ngl.sum[sample]), on = .(sample)]
                          }
                          
                          # Return the local counts data.table
                          return(local_counts)
                          
                        }, error = function(e) {
                          message(sprintf("Error processing taxon %s: %s", taxon, e$message))
                          return(NULL)  # Return NULL in case of error
                        })
                      },
                      
                      mc.cores=num_cores)

counts.mgss.ngl.2 <- do.call(cbind, results.2)

counts.mgss.ngl <- cbind(counts.mgss.ngl.1, counts.mgss.ngl.2)
rownames(counts.mgss.ngl) <- rownames(counts.mgss.ngl.1)

counts.mgss.ngl[is.na(counts.mgss.ngl)] <- 0
fwrite(counts.mgss.ngl, paste(wd, "/norm_counts_mgSs_mgCST_", today2, ".csv", sep=""), row.names = TRUE)

print("norm_counts_mgSs_mgCST_ done")



#######################################################     CENTROID CLASSIFIER     #####################################################################################################


# Defining function to determine yue-clayton theta
yue_distance<-function(row, median){
  #creating a counting variable to index the median list
  taxon_count = 1
  #creating lists to iteratively store output
  median_times_obs <- vector()
  median_minus_obs_sq <- vector()
  #looping through the row and calculating product and difference squared between row data and median data
  for (taxon_abund in row) {
    #calculate p * q
    median_times_obs<-append(median_times_obs, as.numeric(median[taxon_count])*taxon_abund)
    #calculate p-q squared
    median_minus_obs_sq<-append(median_minus_obs_sq, as.numeric((median[taxon_count]-taxon_abund)**2))
    taxon_count <- taxon_count+1
  }
  #calculate sum p* q
  product = sum(median_times_obs)
  #calculate sum p-q squared
  diff_sq = sum(median_minus_obs_sq)
  #calculate yue_med_dist
  yue_med_dist = product / (diff_sq + product)
  #return the value of yue distance
  return(yue_med_dist)
}

## READ IN REFERENCE CENTROIDS
reference_centroids <- as.data.frame(fread(paste(args[3], "/vog_mgCST_centroids_24Sep2024.csv", sep="")))
rownames(reference_centroids) <- reference_centroids$mgCST
reference_centroids$mgCST <- NULL

mgCST.centroids <- as.data.frame(reference_centroids[, 2:ncol(reference_centroids)])
rownames(mgCST.centroids) <- reference_centroids$mgCST

## MAKE RELABUND TABLE FRESH
relabund<-counts.mgss.ngl/rowSums(counts.mgss.ngl)

## REFORMAT RELABUND TO INCLUDE ALL EXPECTED COLUMN NAMES (xvar.names)
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
# 
# 
# for(i in 1:26){
#   print(i)
#   mgCST<-paste("mgCST", i, sep=" ")
#   relabund[[mgCST]]<-apply(relabund[,1:n], 1, function(x) yue_distance(x, reference_centroids[mgCST,]))
# }


## FOR EACH MGCST, MEASURE THE SIMILARITY OF EACH SAMPLE TO EACH MGCST CENTROID USING YUE + CLAYTON THETA

run_yue_distance <- function(i){
  tryCatch({
    mgCST<-paste("mgCST", i, sep=" ")
    relabund[[mgCST]] <- apply(relabund[,1:n], 1, function(x) yue_distance(x, reference_centroids[mgCST,]))
    return(relabund[[mgCST]])
  }, error = function(e) {
    message(sprintf("Error processing taxon %s: %s", mgCST, e$message))
    return(NULL)  # Return NULL in case of error
  })
}

# Set the number of cores for parallel processing
# num_cores <- detectCores() - 1  # Use one less than the total number of cores

i <- 1:26
results.yue.ditance <- mclapply(i, run_yue_distance, mc.cores = num_cores)
# saveRDS(results.yue.ditance, paste(wd, "/results.yue.distance.RDS", sep=""))

# Combine the list of numeric vectors into a single data frame
relabund.yue.distance <- as.data.frame(t(do.call(rbind, lapply(results.yue.ditance, function(x) as.data.frame(t(x))))))

# Rename the columns to reflect the mgCST names
rownames(relabund.yue.distance) <- rownames(relabund)
colnames(relabund.yue.distance) <- paste("mgCST", 1:length(results.yue.ditance), sep = " ")

# relabund.yue.distance <- do.call(cbind, results.yue.ditance)

# write.csv(relabund.yue.distance, paste(wd, "/relabund_test", today2, ".csv", sep=""), row.names = TRUE, quote=F)

m<-n+1
relabund <- cbind(relabund, relabund.yue.distance)
relabund[is.na(relabund)] <- 0
relabund[["mgCST"]]<-colnames(relabund[,m:which(colnames(relabund) %in% "mgCST 26")])[apply(relabund[,m:which(colnames(relabund) %in% "mgCST 26")],1,which.max)]

write.csv(relabund, paste(wd, "/relabund_w_mgCSTs_", today2, ".csv", sep=""), row.names = TRUE, quote=F)
write.csv(relabund["mgCST"], paste(wd, "/mgCSTs_", today2, ".csv", sep=""), row.names = TRUE, quote=F)


## PLOT HEATMAP
mgCST<-as.data.frame(rbind(c("1", "#FE0308"), c("2", "#F54C5E"), c("3", "#F07084"), c("4", "#EC94A5"),c("5", "#F0BCCC"),c("6", "#F6D3DA"),
                           c("7", "#86C61A"), c("8", "#B4DB29"),
                           c("9", "#FF7200"), c("10", "#F68A11"), c("11", "#F8A40E"),
                           c("12", "#FAE50D"),c("13", "#F3F40E"),c("14", "#F3F45F"), 
                           c("15", "#448A73"),c("16", "#BCD6CD"),
                           c("17", "#221886"),c("18", "#3E3792"), c("19", "#5D579E"),c("20", "#7C76AC"),c("21", "#9A98BF"),
                           c("22", "#c997cc"),c("23", "#0d4018"),c("24", "#a16060"),c("25", "#c1ffc1"),c("26", "#8c8c8c"), c("", "white"), c("NA", "white")))

names(mgCST)<-c("mgCST", "color")
colfunc <- colorRampPalette(c("khaki", "limegreen", "darkslategray1", "mediumblue", "magenta", "red"))
relabund.mgCST<-relabund[,1:n]
relabund.mgCST<-relabund.mgCST[,order(colSums(relabund.mgCST), decreasing = TRUE)]
relabund.mgCST$mgCST<-gsub("mgCST ", "", relabund$mgCST)
relabund.mgCST<-relabund.mgCST[order(as.numeric(relabund.mgCST[["mgCST"]])),]
relabund.mgCST[["color"]]<-mgCST[match(relabund.mgCST[["mgCST"]], mgCST$mgCST), "color"]
names(relabund.mgCST)<-gsub("_", " ", names(relabund.mgCST))

pdf(paste(wd, "/mgCST_heatmap_", today2, ".pdf", sep=""), width=7, height=10)
gplots::heatmap.2(t(as.matrix(relabund.mgCST[,1:50])), Colv = FALSE, Rowv = FALSE, col=colfunc(100), keysize= 1.0, densadj=0, density.info='none', key = TRUE, key.ylab=NA, key.title=NA, key.ytickfun=NA, key.xlab="Relative Abundance", trace="none", cexRow = 0.7, cexCol = 0.1, adjRow = c(1, NA),offsetRow = -38,  main = paste("mgCST Heatmap\nnSamples=", paste(nrow(relabund))), title(main = paste("mgCST Heatmap\nnSamples=", paste(nrow(relabund))), line = -2), ColSideColors = as.vector(relabund.mgCST[["color"]]), lhei = c(1,7), dendrogram = "none")
dev.off()

print(paste("mgCST_heatmap_", today2, ".csv has been saved", sep=""))


current_time <- Sys.time()
print(paste("End time:", current_time, sep = " "))

