source("00_importation.R")

## VOG count list and VOG presence absence table

gene.counts.list <- readRDS("../../RDS_files/Gene_analysis/gene.counts.list.RDS")
# ngl.gene.counts.list.mgSs.pa.clean <- readRDS("../../RDS_files/Gene_analysis/ngl.gene.counts.list.mgSs.pa.clean.RDS")
all.clean.100k <- readRDS("../../RDS_files/Gene_analysis/all.clean.100k.RDS")
vog <- read.delim("../../SOURCE_DATA/VIRGO2_VOGkey.txt", row.names = 1)
vog$Taxonomy <- NULL

### VOG based mgss list
vog.count.list <- list()
for (i in names(gene.counts.list)){
 print(i)
 df <- gene.counts.list[[i]]
 df.vog <- merge(df, vog, by = 'row.names', all.x = TRUE)
 df.vog$Row.names <- NULL
 df.vog <- df.vog %>% group_by(VOG) %>% summarise_all(sum) %>% as.data.frame()
 df.vog <- df.vog %>% filter(!is.na(VOG))
 rownames(df.vog) <- df.vog$VOG
 df.vog$VOG <- NULL
#  df.vog <- df.vog %>% mutate_all(~ ifelse(. > 0, 1, .))
 vog.count.list[[i]] <- df.vog
}
saveRDS(vog.count.list, "../../RDS_files/VOG_analysis_updated/vog.count.list.RDS")
# saveRDS(vog.count.list, "../../RDS_files/VOG_analysis_updated/vog.list.mgss.pa.RDS")


print("genes to VOG has run")

