## Libraries----
require(reshape2)
require(data.table)
require(dplyr)
require(stringr)
library(tidyverse)
library(pheatmap)
library(gridExtra)
library(grid)
library(ggh4x)
library(glue)


path.to.VIRGO2.sharing <- "~/bin/VIRGO2_sharing"


## Data importation----

VOG <- read.table(glue("{path.to.VIRGO2.sharing}/VIRGO2_VOGkey.txt"), sep="\t", header=TRUE) %>%
  filter(Taxonomy != "unknown")

gp <- read.table(glue("{path.to.VIRGO2.sharing}/6.VIRGO2.geneProduct.txt"), sep = "\t", header = TRUE, quote = "")
gp <- merge(gp, VOG, all=TRUE)
gp.filtered <- gp[grepl("lactate dehydrogenase", gp$GeneProduct, ignore.case = TRUE), ]

vog.clust <- readRDS("../../RDS_files/VOG_analysis_updated/vog.clusters.RDS")
vog.clust <- as.data.frame(rbindlist(vog.clust, use.names = T, fill = T, idcol = T))
names(vog.clust) <- c("Taxonomy", "VOG", "Cluster")

vog.clust.gp <- merge(vog.clust, gp, all=TRUE) 


cazy <- read.table(glue("{path.to.VIRGO2.sharing}/8.VIRGO2.CAZy.txt"), sep="\t", header=TRUE)
cazy <- merge(cazy, VOG, all.x=TRUE)
cazy <- cazy %>%
  separate_rows(CAZy, sep = ",") %>% 
  filter()


amr <- read.table(glue("{path.to.VIRGO2.sharing}/9.VIRGO2.AMR.txt"), sep="\t", header=TRUE, fill = TRUE)
amr <- merge(amr, VOG, all.x=TRUE)
amr <- amr %>%
  separate_rows(GeneSymbol, sep = ",")


ec <- read.table(glue("{path.to.VIRGO2.sharing}/5.VIRGO2.EC.txt"), sep="\t", header=TRUE, fill = TRUE)
ec <- merge(ec, VOG, all.x=TRUE)
ec <- ec %>%
  separate_rows(EC, sep = ",") %>% 
  select(Gene, VOG, EC, Taxonomy) %>%   # Select columns you want
  distinct() 
#ec.ba<- as.data.frame(cbind(EC=c("4.1.1.18","2.5.1.6","4.1.1.50","3.5.3.1","4.1.1.17","3.5.3.11","3.5.3.12","3.5.1.53","2.1.3.6","2.3.1.57", "1.5.3.13","4.1.1.19","2.5.1.16","1.5.3.16","2.3.1.57","1.5.3.13","2.5.1.22","1.7.2.3","4.3.99.4","1.1.99.1", "1.21.4.4", "4.1.1.25"), 
#            Description=c("Lysine decarboxylase", "Methionine adenosyltransferase","S-adenosylmethionine decarboxylase","Arginase", "Ornithine decarboxylase", "Agmatinase", "Agmatine deiminase", "N-Carbamoylputrescine amidohydrolase","Putrescine transcarbamylase", "Spermidine acetyltransferase","Acetylpolyamine oxidase", "Arginine decarboxylase", "Spermidine synthase", "Spermine oxidase", "Spermine acetyltransferase","Acetylpolyamine oxidase","Spermine synthase","Trimethylamine N-oxide reductase", "Choline trimethylamine-lyase","Choline dehydrogenase", "Betaine reductase", "Tyrosine decarboxylase")))
ec.ba<- as.data.frame(cbind(EC=c("4.1.1.18","3.5.3.11","4.1.1.19","2.5.1.16","1.5.3.16","2.5.1.22","1.7.2.3","4.3.99.4","4.1.1.25"), 
                            Description=c("Lysine decarboxylase", "Agmatinase", "Arginine decarboxylase","Spermidine synthase","Spermine oxidase","Spermine synthase", "Trimethylamine N-oxide reductase","Choline trimethylamine-lyase","Tyrosine decarboxylase"),
                            Product=c("Cadaverine","Putrescine","Agmatine","Spermidine","Spermidine","Spermine","Trimethylamine","Trimethylamine","Tyramine")))

ec.filtered<-merge(ec, ec.ba, all.y=TRUE)
ec.filtered$EC<-paste0(ec.filtered$Product, ", EC:", ec.filtered$EC)
#ec.filtered$EC<-factor(ec.filtered$EC, levels=c("Lysine decarboxylase, EC:4.1.1.18", "Agmatinase, EC:3.5.3.11","Arginine decarboxylase, EC:4.1.1.19","Spermidine synthase, EC:2.5.1.16","Spermine oxidase, EC:1.5.3.16","Spermine synthase, EC:2.5.1.22","Trimethylamine N-oxide reductase, EC:1.7.2.3","Choline trimethylamine-lyase, EC:4.3.99.4","Tyrosine decarboxylase, EC:4.1.1.25"))
#ec.filtered$EC<-paste0(ec.filtered$Description, ", EC:", ec.filtered$EC)
#ec.filtered$EC<-factor(ec.filtered$EC, levels=c("Lysine decarboxylase, EC:4.1.1.18", "Agmatinase, EC:3.5.3.11","Arginine decarboxylase, EC:4.1.1.19","Spermidine synthase, EC:2.5.1.16","Spermine oxidase, EC:1.5.3.16","Spermine synthase, EC:2.5.1.22","Trimethylamine N-oxide reductase, EC:1.7.2.3","Choline trimethylamine-lyase, EC:4.3.99.4","Tyrosine decarboxylase, EC:4.1.1.25"))

cog <- read.table(glue("{path.to.VIRGO2.sharing}/3.VIRGO2.eggNog.txt"), sep="\t", header=TRUE, fill = TRUE)
cog <- merge(cog, VOG, all.x=TRUE)
cog <- cog %>%
  # mutate(COG_category = strsplit(as.character(COG_category), "")) %>%  # Split COG_category with multiple letter
  # unnest(COG_category) %>%
  # filter(COG_category %in% c("W", "U")) %>%
  select(Gene, VOG, COG_category, Taxonomy) %>%   # Select columns you want
  distinct() 

eggNOG <- read.csv(glue("{path.to.VIRGO2.sharing}/3.VIRGO2.eggNog.txt"), sep = "\t")
eggNOG <- merge(eggNOG, VOG, all.x = TRUE)

result.list <- list()
#cogs_to_keep <- c("COG0674", "COG1014", "COG0370", "COG1918", "COG1974", "COG0514", "COG0845", "COG3467", "CPG0778")
cogs_to_keep <- c("COG3467")
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
    #cog_matches == "COG0674" ~ "porA",
    #cog_matches == "COG1013" ~ "porB",
    #cog_matches == "COG1014" ~ "porC",
    #cog_matches == "COG1144" ~ "porD",
    #cog_matches == "COG1918" ~ "feoA",
    #cog_matches == "COG0370" ~ "feoB",
    #cog_matches == "COG0845" ~ "acrA",
    #cog_matches == "COG0841" ~ "acrB",
    cog_matches == "COG3467" ~ "nim",
    #cog_matches == "COG0778" ~ "nfnB",
    #cog_matches == "COG1974" ~ "recA",
    #cog_matches == "COG0514" ~ "recQ",
    
    TRUE ~ cog_matches)) %>%
  filter(cog_matches %in% cogs_to_keep)
# filter(cog_matches %in% c("porA", "porB", "porC", "porD", "recA"))



ngl.abund.clusters.cast.vog <- readRDS("../../RDS_files/VOG_analysis_updated/ngl.abund.clusters.cast.vog.RDS")
ngl.abund.clusters.cast.vog<-ngl.abund.clusters.cast.vog/rowSums(ngl.abund.clusters.cast.vog)
ngl.abund.clusters.cast.vog$sampleID<-rownames(ngl.abund.clusters.cast.vog)
ngl.abund.clusters.cast.vog.m<-reshape2::melt(ngl.abund.clusters.cast.vog, id.vars="sampleID", variable.name="mgSs", value.name="proportion")
#ngl.abund.clusters.cast.vog.m<-ngl.abund.clusters.cast.vog.m[ngl.abund.clusters.cast.vog.m$proportion > 0, ]
ngl.abund.clusters.cast.vog.m$species<-ifelse(str_count(ngl.abund.clusters.cast.vog.m$mgSs, "_") > 2, 
                                              paste(str_split_fixed(ngl.abund.clusters.cast.vog.m$mgSs, "_", n=4)[,1], 
                                                    str_split_fixed(ngl.abund.clusters.cast.vog.m$mgSs, "_", n=4)[,2], 
                                                    str_split_fixed(ngl.abund.clusters.cast.vog.m$mgSs, "_", n=4)[,3], sep="_"), 
                                              ifelse(str_count(ngl.abund.clusters.cast.vog.m$mgSs, "_") > 1,
                                                     paste(str_split_fixed(ngl.abund.clusters.cast.vog.m$mgSs, "_", n=4)[,1], 
                                                           str_split_fixed(ngl.abund.clusters.cast.vog.m$mgSs, "_", n=4)[,2], sep="_"),
                                                     as.character(ngl.abund.clusters.cast.vog.m$mgSs)))

samples_w_mgCST <- read.csv("../../SOURCE_DATA/samples_w_mgCSTs.csv")[,c(1,4)]
total.samples<-samples_w_mgCST %>% group_by(mgCST) %>% summarise(tSamples=length(unique(sampleID)))
samples_w_mgCST<-merge(samples_w_mgCST, ngl.abund.clusters.cast.vog.m, all=TRUE)
samples_w_mgCST$mgCST <- factor(samples_w_mgCST$mgCST, levels=1:25)

mgCST.vog.table <- readRDS("../../RDS_files/VOG_analysis_updated/mgCST.vog.table.RDS")
filtered_list <- lapply(mgCST.vog.table, function(df) {
  melted_df <- reshape2::melt(cbind(VOG=rownames(df), df), id.vars = "VOG", variable.name = "sampleID", value.name = "pa")
  melted_df <- subset(melted_df, pa != 0)
  melted_df$mgCST <- samples_w_mgCST$mgCST[match(melted_df$sampleID, samples_w_mgCST$sampleID)]
  melted_df$Taxonomy <- VOG$Taxonomy[match(melted_df$VOG, VOG$VOG)]
  return(melted_df)
})

mgCST.vog.pa.df<-as.data.frame(rbindlist(filtered_list))
mgCST.vog.pa.df$mgCST <- samples_w_mgCST$mgCST[match(mgCST.vog.pa.df$sampleID, samples_w_mgCST$sampleID)]

## Function----
generate_table <- function(annotation_type, mgcsts, taxon) {
  
  # Initialization : annotation file to use
  if (annotation_type == "CAZy") {
    annot_w_VOG <- cazy
    col.to.consider <- "CAZy"
  } else if (annotation_type == "COG_category") {
    annot_w_VOG <- cog
    col.to.consider <- "COG_category"
  } else if (annotation_type == "amr") {
    annot_w_VOG <- amr
    col.to.consider <- "GeneSymbol"
  } else if (annotation_type == "eggNOG") {
    annot_w_VOG <- eggNOG
    col.to.consider <- "Preferred_name"
    #col.to.consider <- "cog_matches_label"
  } else if (annotation_type == "ec") {
    annot_w_VOG <- ec.filtered
    col.to.consider <- "EC"
  } else if (annotation_type == "gp") {
    annot_w_VOG <- gp
    col.to.consider <- "GeneProduct"
  } else if (annotation_type == "vog.clust") {
    annot_w_VOG <- vog.clust
    col.to.consider <- "Cluster"
  } 
  
  
  # Filter mgCST.vog.pa.df for matching VOGs with annotation file
  annotation.df <- mgCST.vog.pa.df[mgCST.vog.pa.df$VOG %in% unique(annot_w_VOG$VOG), ]
  annotation.df <- annotation.df %>%
    filter(mgCST %in% mgcsts)
  
  # Process data
  annotation.df$col_to_consider <- annot_w_VOG[[col.to.consider]][match(annotation.df$VOG, annot_w_VOG$VOG)]
  if(exists("taxon")){
    to.plot.annotation <- annotation.df[!is.na(annotation.df$col_to_consider) & !annotation.df$col_to_consider %in% "-" & annotation.df$Taxonomy %in% taxon, ] %>%
      group_by(mgCST, Taxonomy, col_to_consider) %>%
      summarise(nSamples = length(unique(sampleID)), .groups = 'drop')
  }else{ 
    to.plot.annotation <- annotation.df[!is.na(annotation.df$col_to_consider) & !annotation.df$col_to_consider %in% "-" & !annotation.df$Taxonomy %in% "MultiGenera", ] %>%
      group_by(mgCST, Taxonomy, col_to_consider) %>%
      summarise(nSamples = length(unique(sampleID)), .groups = 'drop')
  }
  
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

## Data processing----

# vog.to.plot <- generate_table("vog.clust", c(1:25), c("Prevotella_timonensis", "Lactobacillus_crispatus", "Lactobacillus_gasseri", "Lactobacillus_paragasseri", "Lactobacillus_iners", "Lactobacillus_jensenii", "Lactobacillus_mulieris", "UBA629_sp005465875", as.vector(unique(VOG$Taxonomy[grepl("Gardnerella", VOG$Taxonomy)])), "Bifidobacterium_breve", as.vector(unique(VOG$Taxonomy[grepl("Streptococcus", VOG$Taxonomy)]))))


taxa.to.grep <- "Gardnerella"

vog.to.plot.gv <- generate_table("vog.clust", c(1:25), as.vector(unique(VOG$Taxonomy[grepl(taxa.to.grep, VOG$Taxonomy)])))

vog.to.plot.filtered <- vog.to.plot.gv

vog.to.plot.filtered$tax_group <- rownames(vog.to.plot.filtered)
vog.to.plot.m <- reshape2::melt(vog.to.plot.filtered, id.vars="tax_group", variable.name="mgCST", value.name="Prevalence")
vog.to.plot.m$Taxonomy <- str_split_fixed(vog.to.plot.m$tax_group, ", ", 2)[,1]
vog.to.plot.m$VOG_cluster <- str_split_fixed(vog.to.plot.m$tax_group, ", ", 2)[,2]
vog.to.plot.m$VOG_cluster <- factor(vog.to.plot.m$VOG_cluster, levels=10:1)
vog.to.plot.m$tax.plot <- gsub("Gardnerella_", "G.", vog.to.plot.m$Taxonomy)
vog.to.plot.m$tax.plot <- gsub("Lactobacillus_", "L. ", vog.to.plot.m$Taxonomy)
vog.to.plot.m$tax.plot <- gsub("Prevotella_", "P. ", vog.to.plot.m$Taxonomy)
vog.to.plot.m$tax.plot <- gsub("_", " ", vog.to.plot.m$tax.plot)

to.plot <- reshape2::dcast(vog.to.plot.m, tax_group~mgCST, value.var = "Prevalence")
rownames(to.plot) <- to.plot$tax_group
to.plot$tax_group <- NULL
to.plot <- to.plot[,order(colSums(to.plot))]
rownames(to.plot) <- gsub("Gardnerella_", "G. ", rownames(to.plot))
rownames(to.plot) <- gsub("Lactobacillus_", "L. ", rownames(to.plot))
rownames(to.plot) <- gsub("Prevotella_", "P. ", rownames(to.plot))

gard.order <- c("7", "12", "1", "8", "15", "6", "3", "14", "13", "5", "9", "4", "2", "10", "11", "18", "23", "22", "21", "19", "20", "16", "17", "24", "25")
to.plot <- to.plot[,gard.order]

## Plot----
pheatmap(to.plot, 
         color = colorRampPalette(c("antiquewhite", "limegreen", "darkslategray1", "mediumblue", "magenta", "red"))(100), 
         clustering_method = "ward.D2", 
         cluster_cols = F, 
         cluster_rows = T,
         fontsize = 5, 
         angle_col = 0, 
         fontsize_col = 6, 
         border_color = NA,
         main = glue("The prevalence of {taxa.to.grep} VOGs in mgCSTs", sep=" "))



