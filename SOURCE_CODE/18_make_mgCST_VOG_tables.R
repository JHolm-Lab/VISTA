# csv files are created from python script located in RDS_files/VOG_analysis_updated/COG_CAZy_analysis/

# 1-6
crispatus <- read.csv("../../RDS_files/VOG_analysis_updated/COG_CAZy_analysis/crispatus.pa.vog.csv", check.names = F)
# 7-8
gasseri <- read.csv("../../RDS_files/VOG_analysis_updated/COG_CAZy_analysis/gasseri.pa.vog.csv", check.names = F)
# 9
paragasseri <- read.csv("../../RDS_files/VOG_analysis_updated/COG_CAZy_analysis/paragasseri.pa.vog.csv", check.names = F)
# 10-12
iners <- read.csv("../../RDS_files/VOG_analysis_updated/COG_CAZy_analysis/iners.pa.vog.csv", check.names = F)
# 13-15
jensenii <- read.csv("../../RDS_files/VOG_analysis_updated/COG_CAZy_analysis/jensenii.pa.vog.csv", check.names = F)
# 16
mulieris <- read.csv("../../RDS_files/VOG_analysis_updated/COG_CAZy_analysis/mulieris.pa.vog.csv", check.names = F)
# 17-18 UBA629
uba629 <- read.csv("../../RDS_files/VOG_analysis_updated/COG_CAZy_analysis/uba629.pa.vog.csv", check.names = F)
# 19-24 Gardnerella
gardnerella <- read.csv("../../RDS_files/VOG_analysis_updated/COG_CAZy_analysis/gardnerella.pa.vog.csv", check.names = F)
# 25
breve <- read.csv("../../RDS_files/VOG_analysis_updated/COG_CAZy_analysis/breve.pa.vog.csv", check.names = F)


mgCST.vog.pa <- list()

mgCST.vog.pa[["Lactobacillus_crispatus"]] <- crispatus
mgCST.vog.pa[["Lactobacillus_gasseri"]] <- gasseri
mgCST.vog.pa[["Lactobacillus_paragasseri"]] <- paragasseri
mgCST.vog.pa[["Lactobacillus_iners"]] <- iners
mgCST.vog.pa[["Lactobacillus_jensenii"]] <- jensenii
mgCST.vog.pa[["Lactobacillis_mulieris"]] <- mulieris
mgCST.vog.pa[["UBA629"]] <- uba629
mgCST.vog.pa[["Gardnerella"]] <- gardnerella
mgCST.vog.pa[["Bifidobacterium_breve"]] <- breve

saveRDS(mgCST.vog.pa, "../../RDS_files/VOG_analysis_updated/COG_CAZy_analysis/mgCST.vog.pa.RDS")




mgCST <- read.csv("../../SOURCE_DATA/mgCSTs.csv")[,c(2,3)]
samples_w_mgCST <- read.csv("../../SOURCE_DATA/samples_w_mgCSTs.csv")[,c(1,4)] %>%
  left_join(mgCST, by = "mgCST") %>%
  arrange(mgCST)

mgCST.vog.pa <- readRDS("../../RDS_files/VOG_analysis_updated/COG_CAZy_analysis/mgCST.vog.pa.RDS")
names(mgCST.vog.pa)[names(mgCST.vog.pa) == "Lactobacillis_mulieris"] <- "Lactobacillus_mulieris"

mgCST.vog.list <- list()
for (i in unique(samples_w_mgCST$mgCST)) {
  # Get 'domTaxa' value for 'mgCST'
  domTaxa <- samples_w_mgCST %>%
    filter(mgCST == i) %>%
    select(domTaxa) %>%
    distinct() %>%
    slice(1) %>%
    pull()
  
  # Find a matching key in 'mgCST.vog.pa' where the key is part of 'domTaxa'
  matching_key <- names(mgCST.vog.pa)[sapply(names(mgCST.vog.pa), function(x) grepl(x, domTaxa))]
  print(paste(i, domTaxa, matching_key))
  
  if (length(matching_key) > 0) {
    corresponding_data <- mgCST.vog.pa[[matching_key]] %>% 
      column_to_rownames("VOG") %>% 
      select(samples_w_mgCST[samples_w_mgCST$mgCST == i,]$sampleID)
    filtered_data <- corresponding_data[rowSums(corresponding_data) > 0, ]
    
    mgCST.label <- paste("mgCST_", i, sep="")
    mgCST.vog.list[[mgCST.label]] <- filtered_data
  } 
}

for (i in names(mgCST.vog.list)){
  mgCST.a <- mgCST.vog.list[[i]]
  print(sum(colSums(mgCST.a)[colSums(mgCST.a) == 0]))
}


saveRDS(mgCST.vog.list, "../../RDS_files/VOG_analysis_updated/mgCST.vog.table.RDS")

