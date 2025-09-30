# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

## ----Package install---------------------------------------------------------------------------------

packages_to_install <- c("philentropy",
			 "dynamicTreeCut",
			 "plyr",
			 "dplyr",
			 "tidyverse",
			 "gplots",
			 "ggplot2",
			 "reshape2",
			 "RColorBrewer",
			 "vegan",
			 "FSA",
			 "readxl",
			 "wesanderson",
			 "fastcluster",
			 "cluster",
			 "fpc",
			 "dendextend",
			 "randomForestSRC",
                         "vimp", 
			 "psych")

for (package in packages_to_install) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package, dependencies = TRUE)
  }
}

library(philentropy)
library(dynamicTreeCut)
library(plyr)
library(dplyr)
library(tidyverse)
library(gplots)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(vegan)
library(FSA)
require(readxl)
library(wesanderson)
library(fastcluster)
library(cluster)
library(fpc)
library(dendextend)
library(randomForestSRC)
library(vimp)
library(psych)

print("00_importation.R has run")
