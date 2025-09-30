## Load libraries
library(tidyverse)
library(vegan)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)


# Import data----
ngl.abund.clusters.cast <- readRDS("../../RDS_files/VOG_analysis_updated/ngl.abund.clusters.cast.vog.RDS")
mgCST.hclust <- readRDS("../../RDS_files/VOG_analysis_updated/mgCST.hclust.vog.RDS")
mgCST.dist <- readRDS("../../RDS_files/VOG_analysis_updated/mgCST.dist.vog.RDS")

mgCSTs.sort <- read.csv("../../SOURCE_DATA/mgCSTs.csv")
mgCSTs.samples <- read.csv("../../SOURCE_DATA/samples_w_mgCSTs.csv")

mgss.colors <- as.data.frame(read.csv("../../SOURCE_DATA/color.csv"))
names(mgss.colors) <- c("domTaxa", "color")

mgCST.col <- mgCSTs.sort[,c("mgCST", "color")]
names(mgCST.col)<-c("mgCST", "color")

samples_w_mgcst <- read.csv("../../SOURCE_DATA/samples_w_mgCSTs.csv")

participants <- read.csv("../../SOURCE_DATA/VIRGO_participants.csv") %>%
  rename(sampleID = Sample) %>%
  select(sampleID, Geography, CST)

# Calculate relative abundance----
## df == relabund table (samples, rows x taxa, columns)
relabund<-ngl.abund.clusters.cast/rowSums(ngl.abund.clusters.cast)
df <- relabund

samples_w_mgcst <- mgCSTs.samples
# reorder df by sampleID (increasing mgCST)
sample.order <- samples_w_mgcst %>%
  arrange(mgCST) %>%
  pull(sampleID)

df <- df[match(sample.order, rownames(df)),]

# reorder by dominant taxon
taxa.order <- mgCSTs.sort$domTaxa[1:24]
taxa.not.in.order <- setdiff(names(df), taxa.order)

df.DomTaxa <- df %>%
  select(all_of(taxa.order))

df.Other <- df %>%
  select(all_of(taxa.not.in.order))

df.to.plot <- merge(df.DomTaxa, df.Other, by = "row.names", all.x = TRUE) %>%
  column_to_rownames("Row.names")

# group non dominant taxon under "Other" label by summing relabund
n.to.plot <- length(unique(taxa.order))
df.to.plot$Other<-rowSums(df.to.plot[,(n.to.plot+1):ncol(df.to.plot)])
df.to.plot<-df.to.plot[, c(1:n.to.plot, which(names(df.to.plot) %in% "Other"))]
df.to.plot <- df.to.plot[match(sample.order, rownames(df.to.plot)),]

# Shannon diversity annotation
shannon.data <- data.frame("Shannon_diversity" = diversity(df.to.plot, index = "shannon"))

# annotation colors for mgCST
mgCST<-as.data.frame(rbind(c("1", "#FE0308"), c("2", "#F54C5E"), c("3", "#F07084"), c("4", "#EC94A5"),c("5", "#F0BCCC"),c("6", "#F6D3DA"),
                           c("7", "#86C61A"), c("8", "#B4DB29"), 
                           c("9", "#F68A11"), c("10", "#FF981C"),c("11", "#FFA435"),
                           c("12", "#FAE727"),c("13", "#FBEA3F"), c("14", "#FBED58"),
                           c("15", "#E1C775"),
                           c("16", "#589682"),c("17", "#6BA290"),
                           c("18", "#2C31A0"),c("19", "#3C44A8"),c("20", "#444DAC"),c("21", "#676EBC"),c("22", "#6B7EC0"),c("23", "#829CCD"),
                           c("24", "#C7FFC7"), c("25", "#8c8c8c"))) %>% rename(mgCST = V1, color = V2)

annotation_colors_mgCST <- setNames(mgCST$color, mgCST$mgCST)

# Combine Shannon index data with the existing annotation
col.annotation <- samples_w_mgcst %>%
  select(sampleID, mgCST) %>%
  left_join(participants, by="sampleID") %>%
  column_to_rownames("sampleID")%>%
  mutate(mgCST = as.factor(mgCST)) %>%
  rownames_to_column("sampleID") %>%
  left_join(shannon.data %>% rownames_to_column("sampleID"), by = "sampleID") %>%
  column_to_rownames("sampleID") %>%
  select(mgCST, CST, Geography, Shannon_diversity)

# Create columns annotations for heatmap
shannon_colors <- colorRampPalette(c("white", "black"))(4)

CST.color <- data.frame(CST = c('I-A', 'III-A', 'IV-B','III-B','I-B','II','IV-A','IV-C0','IV-C1','IV-C2','IV-C3','IV-C4','V'), 
                        CST_color = c('#fe0308', '#ff7200', '#221886', '#f8a40e', '#f6d3da', '#86c61f','#448a73','#989898','#ef53a7','#a7ddcf','#98c999','#7f0b7c','#fae50d')) %>%
  arrange(CST)
annotation_colors_CST <- setNames(CST.color$CST_color, CST.color$CST)

geography.color <- data.frame(Geography = c("Africa","Asia","Europe","NorthAmer","Oceania"),Geography_color = c("#EE4266", "#000000", "#FFD23F", "#599cd9", "#0EAD69"))
annotation_colors_Geo <- setNames(geography.color$Geography_color, geography.color$Geography)

# Fix taxa names
taxa.names.to.plot <- data.frame(Taxa = colnames(df.to.plot))
taxa.names.to.plot$Taxa <- gsub("(^[A-Za-z])[a-z]+_([a-z]+)(_[A-Za-z]*)?_?([0-9]*)", "\\1. \\2 \\3 \\4", taxa.names.to.plot$Taxa)
taxa.names.to.plot$Taxa <- gsub("\\s+", " ", taxa.names.to.plot$Taxa)
taxa.names.to.plot$Taxa <- gsub("_", "", taxa.names.to.plot$Taxa)
taxa.names.to.plot$Taxa <- gsub("UBA629sp005465875", "Ca. L. vaginae ", taxa.names.to.plot$Taxa)
names(df.to.plot) <- taxa.names.to.plot$Taxa

italic_row_labels <- lapply(names(df.to.plot), function(label) {
  if(grepl("^f ", label)){
    return(label)
  }
  if(grepl("^o ", label)){
    return(label)
  }
  if(grepl("^v ", label)){
    return(label)
  }
  if(grepl("^c ", label)){
    return(label)
  }
  if(grepl("^p ", label)){
    return(label)
  }
  if(grepl("Other", label)){
    return(label)
  }
  else{
    return(bquote(italic(.(label))))
  }
})

colfunc <- colorRampPalette(c("khaki", "limegreen", "darkslategray1", "mediumblue", "magenta", "red"))
png("mgCST_heatmap.png", width = 3, height = 3, units = "in", res = 600)
pheatmap(t(as.matrix(df.to.plot)),
         color = colfunc(50),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_colnames = F,
         fontsize = 2,
         fontsize_col = 5,
         legend = F,
         annotation_row = NA,
         annotation_col = col.annotation,
         annotation_colors = list(mgCST = annotation_colors_mgCST,
                                  Shannon_diversity = shannon_colors,
                                  CST = annotation_colors_CST,
                                  Geography = annotation_colors_Geo),
         labels_row = as.expression(italic_row_labels)
)
dev.off()

