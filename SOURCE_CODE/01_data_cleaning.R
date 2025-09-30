source("00_importation.R")

## ----1. Data cleaning ##

gzipped_file <- "../../SOURCE_DATA/VIRGO2_072823.summary.NR.txt.gz"

all<-data.table::fread(cmd = sprintf("gzip -dc %s", gzipped_file), header=TRUE, fill=TRUE, stringsAsFactors = FALSE, data.table = FALSE, showProgress = TRUE)

gen.sizes.sum <- readRDS("../../RDS_files/Gene_analysis/gen.sizes.sum.RDS")

all[is.na(all)] <- 0

print("Dim all :")
print(dim(all))

#[1] "Dim all :"
#[1] 1631749    2569

## Duplicate Gene cleanup
dups <- all[duplicated(all$Gene), "Gene"]
#rowSums(all[all$Gene %in% dups, 5:ncol(all)])
rownames(all) <- all$Gene

r <- rowSums(all[, 5:ncol(all)])
c <- colSums(all[, 5:ncol(all)])

print("Row summary :")
print(summary(r))

#[1] "Row summary :"
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#       1      120      590    22862     3043 85033096 

print("Col summary :")
print(summary(c))

#[1] "Col summary :"
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#     8548   3154136   8101208  14543985  18829998 277478534

## After fixing dups: remove all rows sums of 0
rm.genes <- r[r == 0]
all.clean <- all[!all$Gene %in% names(rm.genes), ]
print("Dim all.clean :")
print(dim(all.clean))

#[1] "Dim all.clean :"
#[1] 1631749    2569

rm(all)

## Remove all samples w/nReads < 100k
s100k <- c[c > 100000]
all.clean.100k  <- all.clean[, names(all.clean) %in% names(s100k)]
all.clean.info <- all.clean[, 1:4]
saveRDS(all.clean.100k, "../../RDS_files/Gene_analysis/all.clean.100k.RDS")
saveRDS(all.clean.info, "../../RDS_files/Gene_analysis/all.clean.info.RDS")

print("01_data_cleaning.R has run")
