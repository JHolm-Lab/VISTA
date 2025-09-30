
library(tidyverse)
library(pheatmap)
library(data.table)
library(grid)
library(gridExtra)

# Import data----
mgCST.vog.table <- readRDS("../../RDS_files/VOG_analysis_updated/mgCST.vog.table.RDS")
samples_w_mgCST <- read.csv("../../SOURCE_DATA/samples_w_mgCSTs.csv")
VOG <- read.table("../../../VIRGO2_sharing/VIRGO2_VOGkey.txt", sep="\t", header=TRUE)
CAZy <- read.csv("../../../VIRGO2_sharing/8.VIRGO2.CAZy.txt", sep="\t")
eggNOG <- read.csv("../../../VIRGO2_sharing/3.VIRGO2.eggNog.txt", sep="\t")
AMR <- read.csv("../../../VIRGO2_sharing/9.VIRGO2.AMR.txt", sep="\t")

## CAZy with VOG annotation----
CAZy_w_vog <- merge(CAZy, VOG, by = "Gene", all.x = TRUE)
CAZy_w_vog <- CAZy_w_vog[, c("VOG", "CAZy")] %>%
  mutate(CAZy = strsplit(as.character(CAZy), ",")) %>%
  unnest(CAZy) %>%
  distinct()

## COG category with VOG annotation----
eggNOG_w_vog <- eggNOG %>%
  select(Gene, COG_category) %>%
  mutate(COG_category = strsplit(as.character(COG_category), "")) %>%  # Split COG_category with multiple letter
  unnest(COG_category) %>%
  left_join(VOG[,c("Gene", "VOG")], by="Gene") %>%
  select(VOG, COG_category) %>%
  distinct()

## AMR with VOG annotation----
AMR_w_vog <- AMR %>%
  select(Gene, GeneSymbol) %>%
  left_join(VOG[,c("Gene", "VOG")], by="Gene") %>%
  select(VOG, GeneSymbol) %>%
  distinct()

## Dominant taxa w mgCSTs numbers----
mgCST.numbers <- list()
for (i in names(mgCST.vog.table)){
  mgCST.numbers[[i]] <- as.integer(strsplit(i, "_")[[1]][2])
}


# Functions----
generate_annotation_x_mgCST <- function(annotation_type, mgCSTs){
  
  if (annotation_type == "CAZy"){
    annot_w_VOG <- CAZy_w_vog
    col.to.consider <- "CAZy"  # Use string representation
  } else if (annotation_type == "COG"){
    annot_w_VOG <- eggNOG_w_vog
    col.to.consider <- "COG_category"  # Use string representation
  } else if (annotation_type == "AMR"){
    annot_w_VOG <- AMR_w_vog
    col.to.consider <- "GeneSymbol"  # Use string representation
  }
  # Extract VOG presence/absence for considered taxa
  mgCST_data <- mgCST.vog.table[[mgCSTs]] %>% rownames_to_column("VOG")
  mgCST.id <- strsplit(mgCSTs, "_")[[1]][2]
  # Extract samples belonging to mgCSTs dominated by considered taxa
  samples <- samples_w_mgCST[samples_w_mgCST$mgCST %in% mgCST.numbers[[mgCSTs]], ]
  samples <- samples[, c("sampleID", "mgCST")]
  
  # Make annotation x SampleID table
  mgCST_data_w_annot <- merge(mgCST_data, annot_w_VOG, by = "VOG", all.x = TRUE)
  mgCST_data_w_annot <- mgCST_data_w_annot[, -which(colnames(mgCST_data_w_annot) == "VOG")] %>%
    aggregate(as.formula(paste(". ~", col.to.consider)), data = ., sum) %>%
    mutate(across(-all_of(col.to.consider), ~ ifelse(. > 0, 1, 0)))
  
  # Make CAZy x mgCST table
  df <- t(mgCST_data_w_annot)
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  colnames(df) <- df[1,]
  df <- df[-1,]
  df$sampleID <- rownames(df)
  rownames(df) <- NULL
  df <- merge(df, samples, by = "sampleID")
  df <- df[, -which(colnames(df) == "sampleID")]
  df[, -which(colnames(df) == "mgCST")] <- lapply(df[, -which(colnames(df) == "mgCST")], as.numeric)
  df <- aggregate(. ~ mgCST, data = df, sum)
  mgCSTs <- df$mgCST
  df <- t(df[, -1])
  colnames(df) <- mgCSTs
  
  # Number of samples in an mgCST
  count.samples <- as.data.frame(table(samples$mgCST)) %>%
    rename(mgCST = Var1, count = Freq)
  
  # Get proportion of a CAZy in a mgCST
  proportions_annot <- as.data.frame(df)
  colnames(proportions_annot) <- as.character(colnames(proportions_annot))
  for (i in seq_along(proportions_annot)) {
    proportions_annot[[i]] <- proportions_annot[[i]] / count.samples$count[i]
  }
  
  return(list(mgCST_data_w_annot, df, proportions_annot))
}

generate_table <- function(annotation_type){
  annot_w_mgCST <- list()
  for (mgCST in names(mgCST.vog.table)) {
    temp <- generate_annotation_x_mgCST(annotation_type, mgCST)
    df <- as.data.frame(temp[[3]])
    df$rownames <- rownames(df)
    annot_w_mgCST[[mgCST]] <- df
  }
  
  merged_annot <- annot_w_mgCST[[1]]
  for (i in 2:length(annot_w_mgCST)) {
    merged_annot <- merge(merged_annot, annot_w_mgCST[[i]], by = "rownames", all = TRUE)
  }
  
  merged_annot[is.na(merged_annot)] <- 0
  rownames(merged_annot) <- merged_annot$rownames
  merged_annot$rownames <- NULL
  colnames(merged_annot) <- as.character(colnames(merged_annot))
  
  return(merged_annot)
}

generate_heatmap <- function(data, size, title){
  
  mgCST<-as.data.frame(rbind(c("1", "#FE0308"), c("2", "#F54C5E"), c("3", "#F07084"), c("4", "#EC94A5"),c("5", "#F0BCCC"),c("6", "#F6D3DA"),
                             c("7", "#86C61A"), c("8", "#B4DB29"), 
                             c("9", "#F68A11"), c("10", "#FF981C"),c("11", "#FFA435"),
                             c("12", "#FAE727"),c("13", "#FBEA3F"), c("14", "#FBED58"),
                             c("15", "#E1C775"),
                             c("16", "#589682"),c("17", "#6BA290"),
                             c("18", "#2C31A0"),c("19", "#3C44A8"),c("20", "#444DAC"),c("21", "#676EBC"),c("22", "#6B7EC0"),c("23", "#829CCD"),
                             c("24", "#C7FFC7"), c("25", "#8c8c8c"))) %>% rename(mgCST = V1, color = V2)
  
  # columns annotation (mgCST color)
  annotation_colors_mgCST <- setNames(mgCST$color, mgCST$mgCST)
  col.annotation <- data.frame(mgCST = as.character(c(1:25)))
  # italic_row_labels<-rownames(data)
  # if(grepl(", ", italic_row_labels[1])){
  #   name1<-str_split_fixed(string = rownames(data), ", ", 2)[,1]
  #   name2<-str_split_fixed(string = rownames(data), ", ", 2)[,2]
  #   italic_row_labels <- sapply(1:length(name1), function(i) {
  #     bquote(italic(.(name1[i])) ~ "-" ~ .(name2[i]))
  #   })
  # }
  
  number.of.samples <- as.data.frame(table(samples_w_mgCST[,c("mgCST")]))
  names(number.of.samples) <- c("mgCST", "Number_of_samples")
  col_labels <- number.of.samples$Number_of_samples
  
  # colfunc <- colorRampPalette(c("khaki", "limegreen", "darkslategray1", "mediumblue", "magenta", "red"))
  pheatmap(data,
           cluster_rows = F,
           cluster_cols = F,
           show_colnames = TRUE,
           show_rownames = TRUE,
           annotation_col = col.annotation,
           annotation_colors = list(mgCST = annotation_colors_mgCST),
           fontsize_row = size,
           border_color = "black",
           lwd = 0.5,
           angle_col = 0,
           main = title,
           labels_col = col_labels,
           fontsize_col = 8
           # labels_row = do.call(expression, italic_row_labels)
  )
}

# Visualization----
data.CAZy <- generate_table("CAZy")
CAZy_filtered <- data.CAZy[apply(data.CAZy, 1, function(row) any(row > 0.5)), ]

#data.AMR <- generate_table("AMR")
#AMR_filtered <- data.AMR[apply(data.AMR, 1, function(row) any(row > 0.5)), ]

data.eggNOG <- generate_table("COG")
eggNOG_filtered <- data.eggNOG[apply(data.eggNOG, 1, function(row) any(row > 0.5)), ]

p1 <- generate_heatmap(CAZy_filtered, 5, "CAZy")
p2 <- generate_heatmap(eggNOG_filtered, 5, "COG Category")

grid.arrange(p1$gtable, p2$gtable, ncol=2)



#' old code non optmized
#' generate_CAZy_x_mgCST <- function(taxa){
#'   
#'   #' This function creates a CAZy (rows) x mgCST (columns) table for a given dominant taxa 
#'   #' defined in the mgCSTs. It also calculates the proportions of each CAZy category relative 
#'   #' to the number of samples in each mgCST.
#'   #'
#'   #' @param taxa character string representing the dominant taxa for which the CAZy x mgCST 
#'   #'             table is to be generated.
#'   #'
#'   #' @return list containing two elements:
#'   #'   \item{df}{A data frame representing the number of CAZy categories found in each mgCST.}
#'   #'   \item{proportions_CAZy}{A data frame representing the proportions of each CAZy category relative to the number of samples in each mgCST.}
#'   #'
#'   #' @export
#'   
#'   # Extract VOG presence/absence for considered taxa
#'   taxa_data_w_CAZy <- mgCST.vog.pa[[taxa]]
#'   
#'   # Extract samples belonging to mgCSTs dominated by considered taxa
#'   samples <- samples_w_mgCST[samples_w_mgCST$mgCST %in% mgCST.numbers[[taxa]], ]
#'   samples <- samples[, c("sampleID", "mgCST")]
#'   
#'   # Make CAZy x SampleID table
#'   taxa_data_w_CAZy <- merge(taxa_data_w_CAZy, CAZy_w_vog, by = "VOG")
#'   taxa_data_w_CAZy <- taxa_data_w_CAZy[, -which(colnames(taxa_data_w_CAZy) == "VOG")] %>%
#'     aggregate(. ~ CAZy, data = ., sum) %>%
#'     mutate(across(-CAZy, ~ ifelse(. > 0, 1, 0)))
#'   
#'   # Make CAZy x mgCST table
#'   df <- t(taxa_data_w_CAZy)
#'   df <- as.data.frame(df, stringsAsFactors = FALSE)
#'   colnames(df) <- df[1,]
#'   df <- df[-1,]
#'   df$sampleID <- rownames(df)
#'   rownames(df) <- NULL
#'   df <- merge(df, samples, by = "sampleID")
#'   df <- df[, -which(colnames(df) == "sampleID")]
#'   df[, -which(colnames(df) == "mgCST")] <- lapply(df[, -which(colnames(df) == "mgCST")], as.numeric)
#'   df <- aggregate(. ~ mgCST, data = df, sum)
#'   mgCSTs <- df$mgCST
#'   df <- t(df[, -1])
#'   colnames(df) <- mgCSTs
#'   
#'   # Number of sample in an mgCST
#'   count.samples <- as.data.frame(table(samples$mgCST)) %>%
#'     rename(mgCST = Var1, count = Freq)
#'   
#'   # Get proportion of a CAZy in a mgCST
#'   proportions_CAZy <- as.data.frame(df)
#'   colnames(proportions_CAZy) <- as.character(colnames(proportions_CAZy))
#'   for (i in seq_along(proportions_CAZy)) {
#'     proportions_CAZy[[i]] <- proportions_CAZy[[i]] / count.samples$count[i]
#'   }
#'   
#'   return(list(taxa_data_w_CAZy, df, proportions_CAZy))
#' }
#' 
#' generate_AMR_x_mgCST <- function(taxa) {
#'   
#'   #' This function creates a COG_category (rows) x mgCST (columns) table for a given dominant taxa 
#'   #' defined in the mgCSTs. It also calculates the proportions of each COG category relative 
#'   #' to the number of samples in each mgCST.
#'   #'
#'   #' @param taxa character string representing the dominant taxa for which the COG x mgCST 
#'   #'             table is to be generated.
#'   #'
#'   #' @return list containing two elements:
#'   #'   \item{df}{A data frame representing the number of COG categories found in each mgCST.}
#'   #'   \item{proportions_COG}{A data frame representing the proportions of each COG category relative to the number of samples in each mgCST.}
#'   #'
#'   #' @export
#'   
#'   # Extract VOG presence/absence for considered taxa
#'   taxa_data_w_AMR <- mgCST.vog.pa[[taxa]]
#'   
#'   # Extract samples belonging to mgCSTs dominated by considered taxa
#'   samples <- samples_w_mgCST[samples_w_mgCST$mgCST %in% mgCST.numbers[[taxa]], ]
#'   samples <- samples[, c("sampleID", "mgCST")]
#'   
#'   # Make COG x SampleID table
#'   taxa_data_w_AMR <- merge(taxa_data_w_AMR, AMR_w_vog, by = "VOG")
#'   taxa_data_w_AMR <- taxa_data_w_AMR[, -which(colnames(taxa_data_w_AMR) == "VOG")] %>%
#'     aggregate(. ~ GeneSymbol, data = ., sum) %>%
#'     mutate(across(-GeneSymbol, ~ ifelse(. > 0, 1, 0)))
#'   
#'   # Make COG x mgCST table
#'   df <- t(taxa_data_w_AMR)
#'   df <- as.data.frame(df, stringsAsFactors = FALSE)
#'   colnames(df) <- df[1,]
#'   df <- df[-1,]
#'   df$sampleID <- rownames(df)
#'   rownames(df) <- NULL
#'   df <- merge(df, samples, by = "sampleID")
#'   df <- df[, -which(colnames(df) == "sampleID")]
#'   df[, -which(colnames(df) == "mgCST")] <- lapply(df[, -which(colnames(df) == "mgCST")], as.numeric)
#'   df <- aggregate(. ~ mgCST, data = df, sum)
#'   mgCSTs <- df$mgCST
#'   df <- t(df[, -1])
#'   colnames(df) <- mgCSTs
#'   
#'   # Number of sample in an mgCST
#'   count.samples <- as.data.frame(table(samples$mgCST)) %>%
#'     rename(mgCST = Var1, count = Freq)
#'   
#'   # Get proportion of a CAZy in a mgCST
#'   proportions_AMR <- as.data.frame(df)
#'   colnames(proportions_AMR) <- as.character(colnames(proportions_AMR))
#'   for (i in seq_along(proportions_AMR)) {
#'     proportions_AMR[[i]] <- proportions_AMR[[i]] / count.samples$count[i]
#'   }
#'   
#'   return(list(taxa_data_w_AMR, df, proportions_AMR))
#'   
#' }
#' 
#' generate_COG_x_mgCST <- function(taxa) {
#'   
#'   #' This function creates a COG_category (rows) x mgCST (columns) table for a given dominant taxa 
#'   #' defined in the mgCSTs. It also calculates the proportions of each COG category relative 
#'   #' to the number of samples in each mgCST.
#'   #'
#'   #' @param taxa character string representing the dominant taxa for which the COG x mgCST 
#'   #'             table is to be generated.
#'   #'
#'   #' @return list containing two elements:
#'   #'   \item{df}{A data frame representing the number of COG categories found in each mgCST.}
#'   #'   \item{proportions_COG}{A data frame representing the proportions of each COG category relative to the number of samples in each mgCST.}
#'   #'
#'   #' @export
#'   
#'   # Extract VOG presence/absence for considered taxa
#'   taxa_data_w_COG <- mgCST.vog.pa[[taxa]]
#'   
#'   # Extract samples belonging to mgCSTs dominated by considered taxa
#'   samples <- samples_w_mgCST[samples_w_mgCST$mgCST %in% mgCST.numbers[[taxa]], ]
#'   samples <- samples[, c("sampleID", "mgCST")]
#'   
#'   # Filter COG annotation file
#'   COG <- eggNOG_w_vog %>%
#'     filter(VOG %in% taxa_data_w_COG$VOG) %>%
#'     distinct()
#'   
#'   # Make COG x SampleID table
#'   taxa_data_w_COG <- merge(taxa_data_w_COG, COG, by = "VOG")
#'   taxa_data_w_COG <- taxa_data_w_COG[, -which(colnames(taxa_data_w_COG) == "VOG")] %>%
#'     aggregate(. ~ COG_category, data = ., sum) %>%
#'     mutate(across(-COG_category, ~ ifelse(. > 0, 1, 0)))
#'   
#'   # Make COG x mgCST table
#'   df <- t(taxa_data_w_COG)
#'   df <- as.data.frame(df, stringsAsFactors = FALSE)
#'   colnames(df) <- df[1,]
#'   df <- df[-1,]
#'   df$sampleID <- rownames(df)
#'   rownames(df) <- NULL
#'   df <- merge(df, samples, by = "sampleID")
#'   df <- df[, -which(colnames(df) == "sampleID")]
#'   df[, -which(colnames(df) == "mgCST")] <- lapply(df[, -which(colnames(df) == "mgCST")], as.numeric)
#'   df <- aggregate(. ~ mgCST, data = df, sum)
#'   mgCSTs <- df$mgCST
#'   df <- t(df[, -1])
#'   colnames(df) <- mgCSTs
#'   
#'   # Number of sample in an mgCST
#'   count.samples <- as.data.frame(table(samples$mgCST)) %>%
#'     rename(mgCST = Var1, count = Freq)
#'   
#'   # Get proportion of a CAZy in a mgCST
#'   proportions_COG <- as.data.frame(df)
#'   colnames(proportions_COG) <- as.character(colnames(proportions_COG))
#'   for (i in seq_along(proportions_COG)) {
#'     proportions_COG[[i]] <- proportions_COG[[i]] / count.samples$count[i]
#'   }
#'   
#'   return(list(taxa_data_w_COG, df, proportions_COG))
#'   
#' }
#' 
#' # Visualization----
#' 
#' ## Make CAZy x mgCST table----
#' CAZy_w_mgCST <- list()
#' for (taxon in names(mgCST.vog.pa)) {
#'   temp <- generate_CAZy_x_mgCST(taxon)
#'   df <- as.data.frame(temp[[3]])
#'   df$rownames <- rownames(df)
#'   CAZy_w_mgCST[[taxon]] <- df
#' }
#' 
#' merged_CAZy <- CAZy_w_mgCST[[1]]
#' for (i in 2:length(CAZy_w_mgCST)) {
#'   merged_CAZy <- merge(merged_CAZy, CAZy_w_mgCST[[i]], by = "rownames", all = TRUE)
#' }
#' 
#' merged_CAZy[is.na(merged_CAZy)] <- 0
#' rownames(merged_CAZy) <- merged_CAZy$rownames
#' merged_CAZy$rownames <- NULL
#' colnames(merged_CAZy) <- as.character(colnames(merged_CAZy))
#' 
#' ## Make AMR x mgCST table----
#' AMR_w_mgCST <- list()
#' for (taxon in names(mgCST.vog.pa)) {
#'   temp <- generate_AMR_x_mgCST(taxon)
#'   df <- as.data.frame(temp[[3]])
#'   df$rownames <- rownames(df)
#'   AMR_w_mgCST[[taxon]] <- df
#' }
#' 
#' merged_AMR <- AMR_w_mgCST[[1]]
#' for (i in 2:length(AMR_w_mgCST)) {
#'   merged_AMR <- merge(merged_AMR, AMR_w_mgCST[[i]], by = "rownames", all = TRUE)
#' }
#' 
#' merged_AMR[is.na(merged_AMR)] <- 0
#' rownames(merged_AMR) <- merged_AMR$rownames
#' merged_AMR$rownames <- NULL
#' colnames(merged_AMR) <- as.character(colnames(merged_AMR))
#' 
#' ## Make COG x mgCST table----
#' COG_w_mgCST <- list()
#' for (taxon in names(mgCST.vog.pa)) {
#'   temp <- generate_COG_x_mgCST(taxon)
#'   df <- as.data.frame(temp[[3]])
#'   df$rownames <- rownames(df)
#'   COG_w_mgCST[[taxon]] <- df
#' }
#' 
#' merged_COG <- COG_w_mgCST[[1]]
#' for (i in 2:length(COG_w_mgCST)) {
#'   merged_COG <- merge(merged_COG, COG_w_mgCST[[i]], by = "rownames", all = TRUE)
#' }
#' 
#' merged_COG[is.na(merged_COG)] <- 0
#' rownames(merged_COG) <- merged_COG$rownames
#' merged_COG$rownames <- NULL  # Remove the rownames column if not needed
#' colnames(merged_COG) <- as.character(colnames(merged_COG))
#' 
#' ## Heatmaps----
#' 
#' # Filter data : CAZy categories found in at least 90% of any mgCST’s samples
#' CAZy_filtered <- merged_CAZy[apply(merged_CAZy, 1, function(row) any(row > 0.9)), ]
#' generate_heatmap(CAZy_filtered, 5, "CAZy")
#' 
#' AMR_filtered <- merged_AMR[apply(merged_AMR, 1, function(row) any(row > 0.9)), ]
#' generate_heatmap(AMR_filtered, 5, "AMR")
#' 
#' COG_filtered <- merged_COG[apply(merged_COG, 1, function(row) any(row > 0.9)), ]
#' generate_heatmap(COG_filtered, 5, "COG")
#' 
#' 
