require(reshape2)
require(data.table)
require(dplyr)
require(stringr)
library(tidyverse)
library(pheatmap)
library(gridExtra)

# VIRGO2 annotation files----
VOG <- read.table("../../../VIRGO2_sharing/VIRGO2_VOGkey.txt", sep="\t", header=TRUE) %>%
  filter(Taxonomy != "unknown")

cazy <- read.table("../../../VIRGO2_sharing/8.VIRGO2.CAZy.txt", sep="\t", header=TRUE)
cazy <- merge(cazy, VOG, all.x=TRUE)
cazy <- cazy %>%
  separate_rows(CAZy, sep = ",")

amr <- read.table("../../../VIRGO2_sharing/9.VIRGO2.AMR.txt", sep="\t", header=TRUE, fill = TRUE)
amr <- merge(amr, VOG, all.x=TRUE)
amr <- amr %>%
  separate_rows(GeneSymbol, sep = ",")

cog <- read.table("../../../VIRGO2_sharing/3.VIRGO2.eggNog.txt", sep="\t", header=TRUE, fill = TRUE)
cog <- merge(cog, VOG, all.x=TRUE)
cog <- cog %>%
  # mutate(COG_category = strsplit(as.character(COG_category), "")) %>%  # Split COG_category with multiple letter
  # unnest(COG_category) %>%
  # filter(COG_category %in% c("W", "U")) %>%
  select(Gene, VOG, COG_category, Taxonomy) %>%   # Select columns you want
  distinct() 

eggNOG <- read.csv("../../../VIRGO2_sharing/3.VIRGO2.eggNog.txt", sep = "\t")
eggNOG <- merge(eggNOG, VOG, all.x = TRUE)

result.list <- list()
cogs_to_keep <- c("COG0674", "COG1014", "COG0370", "COG1918", "COG1974", "COG0514", "COG0845", "COG3467", "CPG0778")
for (cog in cogs_to_keep){
  result.list[[cog]] <- eggNOG[grep(cog, eggNOG$eggNOG_OGs, ignore.case = T), ]
}

eggNOG.filtered <- do.call(rbind, result.list)
pattern <- "COG\\d+"
eggNOG.filtered <- eggNOG.filtered %>%
  mutate(cog_matches = str_extract_all(eggNOG_OGs, pattern)) %>%
  mutate(cog_matches = sapply(cog_matches, function(x) paste(unique(x), collapse = ", "))) %>%
  mutate(cog_matches = strsplit(cog_matches, ", ")) %>%
  unnest(cog_matches) %>%
  mutate(cog_matches_label = case_when(
    cog_matches == "COG0674" ~ "pfoR",
    cog_matches == "COG1014" ~ "pfoR",
    cog_matches == "COG0370" ~ "FeoAB",
    cog_matches == "COG1918" ~ "FeoAB",
    cog_matches == "COG0845" ~ "BME",
    cog_matches == "COG3467" ~ "NIM",
    cog_matches == "COG0778" ~ "NIM",
    cog_matches == "COG1974" ~ "recA",
    cog_matches == "COG0514" ~ "recA",
    
    TRUE ~ cog_matches)) %>%
  filter(cog_matches %in% cogs_to_keep)
  # filter(cog_matches %in% c("pfoR", "FeoAB", "BME", "NIM", "recA"))


# Make table----
#VOG_SampleID_mgCST_Taxonomy

ngl.abund.clusters.cast.vog <- readRDS("../../RDS_files/VOG_analysis_updated/ngl.abund.clusters.cast.vog.RDS")
ngl.abund.clusters.cast.vog<-ngl.abund.clusters.cast.vog/rowSums(ngl.abund.clusters.cast.vog)
ngl.abund.clusters.cast.vog$sampleID<-rownames(ngl.abund.clusters.cast.vog)
ngl.abund.clusters.cast.vog.m<-reshape2::melt(ngl.abund.clusters.cast.vog, id.vars="sampleID", variable.name="mgSs", value.name="proportion")
ngl.abund.clusters.cast.vog.m<-ngl.abund.clusters.cast.vog.m[ngl.abund.clusters.cast.vog.m$proportion > 0, ]
ngl.abund.clusters.cast.vog.m$species<-ifelse(str_count(ngl.abund.clusters.cast.vog.m$mgSs, "_") > 2, 
                                              paste(str_split_fixed(ngl.abund.clusters.cast.vog.m$mgSs, "_", n=4)[,1], 
                                                    str_split_fixed(ngl.abund.clusters.cast.vog.m$mgSs, "_", n=4)[,2], 
                                                    str_split_fixed(ngl.abund.clusters.cast.vog.m$mgSs, "_", n=4)[,3], sep="_"), 
                                              ifelse(str_count(ngl.abund.clusters.cast.vog.m$mgSs, "_") > 1,
                                                     paste(str_split_fixed(ngl.abund.clusters.cast.vog.m$mgSs, "_", n=4)[,1], 
                                                           str_split_fixed(ngl.abund.clusters.cast.vog.m$mgSs, "_", n=4)[,2], sep="_"),
                                                     as.character(ngl.abund.clusters.cast.vog.m$mgSs)))

samples_w_mgCST <- read.csv("../../SOURCE_DATA/samples_w_mgCSTs.csv")[,c(1,4)]
samples_w_mgCST<-merge(samples_w_mgCST, ngl.abund.clusters.cast.vog.m, all=TRUE)

mgCST.vog.table <- readRDS("../../RDS_files/VOG_analysis_updated/mgCST.vog.table.RDS")
filtered_list <- lapply(mgCST.vog.table, function(df) {
  melted_df <- reshape2::melt(cbind(VOG=rownames(df), df), id.vars = "VOG", variable.name = "sampleID", value.name = "pa")
  melted_df <- subset(melted_df, pa != 0)
  melted_df$mgCST <- samples_w_mgCST$mgCST[match(melted_df$sampleID, samples_w_mgCST$sampleID)]
  melted_df$Taxonomy <- VOG$Taxonomy[match(melted_df$VOG, VOG$VOG)]
  return(melted_df)
})

mgCST.vog.pa.df<-as.data.frame(rbindlist(filtered_list))
total.samples<-samples_w_mgCST %>% group_by(mgCST) %>% summarise(tSamples=length(unique(sampleID)))

# Functions----
generate_table <- function(annotation_type, mgcsts) {
  
  # Initialization : annotation file to use
  if (annotation_type == "CAZy") {
    annot_w_VOG <- cazy
    col.to.consider <- "CAZy"
  } else if (annotation_type == "COG_category") {
    annot_w_VOG <- cog
    col.to.consider <- "COG_category"
  } else if (annotation_type == "GeneSymbol") {
    annot_w_VOG <- amr
    col.to.consider <- "GeneSymbol"
  } else if (annotation_type == "eggNOG") {
    annot_w_VOG <- eggNOG.filtered
    #col.to.consider <- "cog_matches"
    col.to.consider <- "cog_matches_label"
  }
  
  # Filter mgCST.vog.pa.df for matching VOGs with annotation file
  annotation.df <- mgCST.vog.pa.df[mgCST.vog.pa.df$VOG %in% unique(annot_w_VOG$VOG), ]
  annotation.df <- annotation.df %>%
    filter(mgCST %in% mgcsts)
  
  # Process data
  annotation.df$col_to_consider <- annot_w_VOG[[col.to.consider]][match(annotation.df$VOG, annot_w_VOG$VOG)]
  to.plot.annotation <- annotation.df[!is.na(annotation.df$col_to_consider) & !annotation.df$col_to_consider %in% "-" & !annotation.df$Taxonomy %in% "MultiGenera", ] %>%
    group_by(mgCST, Taxonomy, col_to_consider) %>%
    summarise(nSamples = length(unique(sampleID)), .groups = 'drop')
  
  # Merge with total samples data and calculate proportion
  to.plot.annotation <- merge(to.plot.annotation, total.samples, all.x = TRUE)
  to.plot.annotation$prop <- to.plot.annotation$nSamples / to.plot.annotation$tSamples
  
  # Reshape data for annotation
  data.annotation <- reshape2::dcast(to.plot.annotation, Taxonomy + col_to_consider ~ mgCST, value.var = "prop", fill = 0)
  data.annotation$Tax_annotation <- paste(data.annotation$Taxonomy, data.annotation$col_to_consider, sep = ", ")
  rownames(data.annotation) <- data.annotation$Tax_annotation
  data.annotation <- data.annotation[order(data.annotation$col_to_consider, data.annotation$Taxonomy), ]
  data.annotation$Taxonomy <- NULL
  data.annotation$col_to_consider <- NULL
  data.annotation$Tax_annotation <- NULL
  
  # # Filter rows based on percentage
  # annotation_filtered <- data.annotation[apply(data.annotation, 1, function(row) any(row > percentage)), ]
  # rownames(annotation_filtered) <- gsub("_", " ", rownames(annotation_filtered))
  
  return(data.annotation)
}

generate_heatmap <- function(data, title){
  
  mgcsts <- colnames(data)
  
  mgCST<-as.data.frame(rbind(c("1", "#FE0308"), c("2", "#F54C5E"), c("3", "#F07084"), c("4", "#EC94A5"),c("5", "#F0BCCC"),c("6", "#F6D3DA"),
                             c("7", "#86C61A"), c("8", "#B4DB29"), c("9", "#0F411A"), 
                             c("10", "#F68A11"), c("11", "#FF981C"),c("12", "#FFA435"),
                             c("13", "#FAE727"),c("14", "#FBEA3F"), c("15", "#FBED58"),
                             c("16", "#E1C775"),
                             c("17", "#589682"),c("18", "#6BA290"),
                             c("19", "#444DAC"),c("20", "#6B7EC0"),c("21", "#829CCD"),c("22", "#3C44A8"),c("23", "#676EBC"),c("24", "#2C31A0"),
                             c("25", "#C7FFC7"))) %>% rename(mgCST = V1, color = V2) %>%
    filter(mgCST %in% mgcsts)
  
  # columns annotation (mgCST color)
  annotation_colors_mgCST <- setNames(mgCST$color, mgCST$mgCST)
  # col.annotation <- data.frame(mgCST = as.character(c(1:25)))
  col.annotation <- data.frame(mgCST = mgcsts)
  rownames(col.annotation) <- mgcsts
  
  italic_row_labels<-rownames(data)
  if(str_count(italic_row_labels[1], ", ") > 1){
    if(grepl(", ", italic_row_labels[1])){
      name1<-gsub("_", " ", str_split_fixed(string = rownames(data), ", ", 3)[,1])
      name2<-str_split_fixed(string = rownames(data), ", ", 3)[,2]
      name3<-str_split_fixed(string = rownames(data), ", ", 3)[,3]
      italic_row_labels <- sapply(1:length(name1), function(i) {
        bquote(italic(.(name1[i])) ~ "-" ~ .(name2[i]) ~ "-" ~ .(name3[i]))
      })
    }
  }else{
    if(grepl(", ", italic_row_labels[1])){
      name1<-gsub("_", " ", str_split_fixed(string = rownames(data), ", ", 2)[,1])
      name2<-str_split_fixed(string = rownames(data), ", ", 2)[,2]
      italic_row_labels <- sapply(1:length(name1), function(i) {
        bquote(italic(.(name1[i])) ~ "-" ~ .(name2[i]))
      })
    }
  }
  
  
  nSamples <- total.samples %>% filter(mgCST %in% mgcsts)
  
  pheatmap(data,
           cluster_rows = F,
           cluster_cols = F,
           show_colnames = TRUE,
           show_rownames = TRUE,
           annotation_legend = F,
           fontsize = 6,
           annotation_col = col.annotation,
           annotation_colors = list(mgCST = annotation_colors_mgCST),
           fontsize_row = 7,
           border_color = "black",
           lwd = 0.5,
           angle_col = 0,
           main = title,
           labels_col = nSamples$tSamples,
           fontsize_col = 6, 
           gaps_row = c(0, 3, 10, 21),
           cellwidth = 14,
           labels_row = do.call(expression, italic_row_labels)
  )
}

# Plot data----

# Use functions to create table and plot heatmap: select annotation data and list of mgCST to consider
cazy.to.plot <- generate_table("CAZy", c(17:24))
eggNOG.to.plot <- generate_table("eggNOG", c(1:25))

# Filtering
percentage.cazy <- 0.8
cazy.to.plot.filtered <- cazy.to.plot[apply(cazy.to.plot, 1, function(row) any(row > percentage.cazy)), ]

percentage.eggNOG <- 0.9
eggNOG.to.plot.filtered <- eggNOG.to.plot[apply(eggNOG.to.plot, 1, function(row) any(row > percentage.eggNOG)), ]

p1 <- generate_heatmap(cazy.to.plot.filtered,
                       paste("CAZy - ",percentage.cazy*100,"%", sep=""))

p2 <- generate_heatmap(eggNOG.to.plot.filtered, 
                       paste("eggNOG - ",percentage.eggNOG*100,"%", sep=""))

grid.arrange(p1$gtable, p2$gtable, ncol = 2)


# cog.to.plot <- generate_table("COG_category", c(17:24))
# cog.to.plot.filtered <- cog.to.plot[apply(cog.to.plot, 1, function(row) any(row > 0.8)), ]
# generate_heatmap(cog.to.plot.filtered, "COG")
# 
# amr.to.plot <- generate_table("GeneSymbol", c(17:24))
# amr.to.plot.filtered <- amr.to.plot[apply(amr.to.plot, 1, function(row) any(row > 0)), ]
# generate_heatmap(amr.to.plot.filtered, "AMR")


# Iners mgSs analysis----
mgss <- readRDS("../../RDS_files/VOG_analysis_updated/vog.list.mgss.pa.RDS") 
samples.clusters <- readRDS("../../RDS_files/VOG_analysis_updated/samples.clusters.updated.vog.RDS")

data.iners <-mgss[["Lactobacillus_iners"]]
iners.clusters <- samples.clusters[['Lactobacillus_iners']]

iners.1.samples <- iners.clusters[iners.clusters$sample_cluster == 1,]$sampleID
iners.2.samples <- iners.clusters[iners.clusters$sample_cluster == 2,]$sampleID

iners.mgss1 <- data.iners %>%
  select(all_of(iners.1.samples))

iners.mgss2 <- data.iners %>%
  select(all_of(iners.2.samples))

# Get proportion of VOG in samples
vog.count.iners1 <- data.frame(count = rowSums(iners.mgss1))
vog.count.iners1$proportion <- vog.count.iners1$count / length(iners.1.samples)      

vog.count.iners2 <- data.frame(count = rowSums(iners.mgss2))
vog.count.iners2$proportion <- vog.count.iners2$count / length(iners.2.samples)      

# Add VOG as a column in each DataFrame
vog.count.iners1$VOG <- rownames(vog.count.iners1)
vog.count.iners2$VOG <- rownames(vog.count.iners2)

# Merge the two DataFrames by VOG
vog.count.combined <- merge(vog.count.iners1[, c("VOG", "count", "proportion")], vog.count.iners2[, c("VOG", "count", "proportion")], 
                            by = "VOG", suffixes = c("_1", "_2"))

# Get delta of proportion : mgss1 - mgss2
vog.count.combined$delta <- vog.count.combined$proportion_1 - vog.count.combined$proportion_2
vog.count.combined <- vog.count.combined %>%
  arrange(delta)

## CAZy----
vog.cazy.iners <- vog.count.combined %>%
  left_join(cazy[,c("VOG","CAZy")], by="VOG") %>%
  filter(!(is.na(CAZy))) %>%
  group_by(VOG, delta) %>%
  summarise(CAZy = paste(unique(CAZy), collapse = ", ")) %>%
  distinct()

# Reshape the data to have proportion_1 and proportion_2 as separate columns
vog.cazy.heatmap <- vog.cazy.iners %>%
  select(VOG, CAZy, delta)%>%
  pivot_wider(names_from = CAZy, values_from = delta) %>%
  mutate(across(everything(), ~ replace_na(.x, 0)))

# Plot
vog.cazy.matrix <- as.matrix(vog.cazy.heatmap[,-1]) 
rownames(vog.cazy.matrix) <- vog.cazy.heatmap$VOG
vog.cazy.matrix_transposed <- t(vog.cazy.matrix)
color_gradient <- colorRampPalette(c("blue", "white"))(100)
pheatmap(vog.cazy.matrix_transposed, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         main = "Heatmap of Delta Values (VOG vs CAZy)", 
         color = color_gradient)
    
## eggNOG----   
pattern <- "COG\\d+"
vog.eggNOG.iners <- as.data.frame(eggNOG) %>%
  filter(VOG %in% vog.count.combined$VOG) %>%
  mutate(cog_matches = str_extract_all(eggNOG_OGs, pattern)) %>%
  mutate(cog_matches = sapply(cog_matches, function(x) paste(unique(x), collapse = ", "))) %>%
  mutate(cog_matches = strsplit(cog_matches, ", ")) %>%
  unnest(cog_matches) %>%
  distinct() %>%
  select(VOG, cog_matches, COG_category) %>%
  left_join(vog.count.combined, by="VOG") %>%
  filter(abs(delta) > 0.5)

vog.eggNOG.heatmap <- vog.eggNOG.iners %>%
  select(VOG, cog_matches, delta) %>%
  pivot_wider(
    names_from = cog_matches, 
    values_from = delta,
    values_fn = mean
  ) %>%
  mutate(across(everything(), ~ replace_na(.x, 0)))

vog.eggNOG.matrix <- as.matrix(vog.eggNOG.heatmap[,-1])
rownames(vog.eggNOG.matrix) <- vog.eggNOG.heatmap$VOG

pheatmap(t(vog.eggNOG.matrix), 
         cluster_rows = F,
         cluster_cols = F,
         main = "Iners1 - Iners2",
         border_color = "black",
         lwd = 0.5, 
         color = color_gradient)


## KEGG----
kegg <- as.data.frame(read.table("../../../VIRGO2_sharing/7.VIRGO2.kegg.txt.gz", sep = "\t", header = TRUE))
kegg <- merge(kegg, VOG, by="Gene", all.x = TRUE)

vog.kegg.iners <- vog.count.combined %>%
  left_join(kegg[,c("VOG", "KEGG")], by="VOG") %>%
  filter(!(is.na(KEGG))) %>%
  mutate(KEGG = strsplit(KEGG, ",")) %>%
  unnest(KEGG) %>%
  distinct()

vog.kegg.iners.filtered <- vog.kegg.iners %>%
  filter(abs(delta) > 0.01)

vog.kegg.heatmap <- vog.kegg.iners.filtered %>%
  select(VOG, KEGG, delta) %>%
  pivot_wider(
    names_from = KEGG, 
    values_from = delta,
    # values_fn = mean
  ) %>%
  mutate(across(everything(), ~ replace_na(.x, 0)))

vog.kegg.matrix <- as.matrix(vog.kegg.heatmap[,-1])
rownames(vog.kegg.matrix) <- vog.kegg.heatmap$VOG

pheatmap(t(vog.kegg.matrix), 
         cluster_rows = F,
         cluster_cols = F,
         main = "Iners1 - Iners2",
         border_color = "black",
         lwd = 0.5, 
         color = color_gradient)
  

kegg.file <- kegg %>%
  filter(VOG %in% vog.kegg.iners[abs(vog.kegg.iners$delta) > 0.01,]$VOG) %>%
  select(VOG, KEGG) %>%
  mutate(KEGG = strsplit(KEGG, ",")) %>%
  unnest(KEGG) %>%
  distinct()
write.table(kegg.file, "../../RDS_files/VOG_analysis_updated/iners.mgss.kegg.txt", sep="\t", row.names = F, col.names = F, quote = F)
  



# # Filter the data for iner1
# iner1.file <- vog.kegg.iners.filtered %>%
#   filter(delta > 0.01) %>%
#   select(VOG, KEGG)
# # Filter the data for iner2
# iner2.file <- vog.kegg.iners.filtered %>%
#   filter(delta < -0.01) %>%
#   select(VOG, KEGG)
# # Specify the output file path
# output_file <- "~/Desktop/iners_combined.txt"
# cat("# iners 1\n", file = output_file)
# write.table(iner1.file, output_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
# cat("# iners 2\n", file = output_file, append = TRUE)
# write.table(iner2.file, output_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)


# find most diffrent vog between mgss
# ex cog present in mgss1 or 2

# keep cazy and eggNOG
# keep cog label
# mgcst 17 to 24
# 
# 
# diagram for iners:
#   identfy which cog are unique betweeen mgss from a same major domTaxa


# Old code not working
## Iners accross all mgCSTs and COG_categories
# cog <- read.table("../../../VIRGO2_sharing/3.VIRGO2.eggNog.txt", sep="\t", header=TRUE, fill = TRUE)
# cog <- merge(cog, VOG, all.x=TRUE)
# cog <- cog %>%
#   mutate(COG_category = strsplit(as.character(COG_category), "")) %>%
#   unnest(COG_category) %>%
#   # filter(COG_category %in% c("N","W","U","X","I","Q")) %>%
#   select(Gene, VOG, COG_category, Taxonomy) %>%
#   distinct() 
# 
# data.iners<-mgCST.vog.pa.df[mgCST.vog.pa.df$VOG %in% unique(cog$VOG), ] %>%
#   filter(Taxonomy == "Lactobacillus_iners")
# data.iners$COG <- cog$COG_category[match(data.iners$VOG, cog$VOG)]
# data.iners<-data.iners[!is.na(data.iners$COG) & !data.iners$COG %in% "-" & !data.iners$Taxonomy %in% "MultiGenera", ] %>% group_by(mgCST, Taxonomy, COG) %>% summarise(nSamples=length(unique(sampleID)))
# data.iners<-merge(data.iners, total.samples, all.x=TRUE)
# data.iners$prop<-data.iners$nSamples/data.iners$tSamples
# 
# data.iners.to.plot<-reshape2::dcast(data.iners, Taxonomy+COG~mgCST, value.var = "prop", fill = 0)
# data.iners.to.plot$Tax_COG<-paste(data.iners.to.plot$Taxonomy, data.iners.to.plot$COG, sep=", ")
# rownames(data.iners.to.plot)<-data.iners.to.plot$Tax_COG
# data.iners.to.plot<-data.iners.to.plot[order(data.iners.to.plot$Taxonomy, data.iners.to.plot$COG), ]
# data.iners.to.plot$Taxonomy<-NULL
# data.iners.to.plot$COG<-NULL
# data.iners.to.plot$Tax_COG<-NULL
# 
# iners_filtered <- data.iners.to.plot[apply(data.iners.to.plot, 1, function(row) any(row > 0.5)), ] # Change filtering value here
# 
# rownames(iners_filtered)<-gsub("_", " ", rownames(iners_filtered))
# col.annotation.iners <- data.frame(mgCST = as.character(c(1:25)))
# rownames(col.annotation.iners) <- col.annotation.iners$mgCST
# annotation_colors_iners <- setNames(mgCST$color, mgCST$mgCST)
# 
# pheatmap(iners_filtered,
#          cluster_rows = F,
#          cluster_cols = F,
#          show_colnames = TRUE,
#          show_rownames = TRUE,
#          annotation_col = col.annotation.iners,
#          annotation_colors = list(mgCST = annotation_colors_iners),
#          fontsize_row = 5,
#          border_color = "black",
#          lwd = 0.5,
#          angle_col = 0,
#          labels_col =total.samples$tSamples,
#          main="Iners x COG",
#          fontsize_col = 5)
# 
