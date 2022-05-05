# Startup -----------------------------------------------------------------

pkg_list <- c("data.table","survey","purrr","loo")

for (pkg in pkg_list) {
  if (!(pkg %in% installed.packages())) {
    install.packages(pkg)
  }
}


# install.packages("frugalSim")
# causl_path <- file.path("C:/Users/Stacey/Desktop/PhD/R/Causl - Simulate Data/Package/causl_0.1.6.tar.gz")
# install.packages(causl_path, repos = NULL, type = "source")
# devtools::install_github("rje42/causl")

library(causl)
library(data.table)
library(survey)
library(mvtnorm)
library(ggplot2)
library(loo)
library(coda)
library(gridExtra)
library(coda)

# File locations ----------------------------------------------------------

folder_dir<-"D:/xlin/R_code/WIP-SMI/Scenarios/Baseline/GPC"
data_dir <- file.path(folder_dir,"/data")
results_dir <- file.path(folder_dir,"/results")
