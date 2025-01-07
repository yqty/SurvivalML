init.config <- function() {
  ## Set CRAN and Bioconductor repositories
  options("repos" = c(CRAN = "https://mirrors.westlake.edu.cn/CRAN/"))
  options(BioC_mirror = "https://mirrors.westlake.edu.cn/bioconductor")

  message("[++] It will cost a while for the first time.")

  ## Load pacman and devtools if not installed
  if (!requireNamespace("pacman", quietly = TRUE)) {
    utils::install.packages("pacman")
  }

  if (!requireNamespace("devtools", quietly = TRUE)) {
    utils::install.packages("devtools")
  }

  ## Load packages from GitHub with error handling
  install_github_package <- function(pkg_name, repo) {
    if (!requireNamespace(pkg_name, quietly = TRUE)) {
      message(paste("[++] Installing", pkg_name, "from GitHub"))
      tryCatch(
        devtools::install_github(repo),
        error = function(e) {
          message(paste("[++] ERROR: R package [", pkg_name, "] from GitHub can NOT be installed."))
        }
      )
    }
  }

  ## Install necessary GitHub packages
  install_github_package("ComplexHeatmap", "jokergoo/ComplexHeatmap")
  suppressPackageStartupMessages(library(ComplexHeatmap))
  install_github_package("summaryBox", "deepanshu88/summaryBox")
  suppressPackageStartupMessages(library(summaryBox))
  install_github_package("bsplus", "ijlyttle/bsplus")
  suppressPackageStartupMessages(library(bsplus))
  install_github_package("Github-Yilei/ggcor", NULL)
  suppressPackageStartupMessages(library(ggcor))
  install_github_package("YuLab-SMU/clusterProfiler", NULL)
  suppressPackageStartupMessages(library(clusterProfiler))
  install_github_package("binderh/CoxBoost", NULL)
  suppressPackageStartupMessages(library(CoxBoost))

  ## Load the packages quietly
  load_package <- function(pkg_name) {
    if (!requireNamespace(pkg_name, quietly = TRUE)) {
      message(paste("[++] ERROR: R package [", pkg_name, "] is not available."))
    } else {
      library(pkg_name, character.only = TRUE, quietly = TRUE)
    }
  }

  ## List of CRAN and Bioconductor packages
  pkgs.to.check <- c(
    "shiny", "shinyjs", "bslib", "waiter", "openxlsx", "plotly", "ggplot2", "stringr",
    "slickR", "shinyvalidate", "markdown", "bs4Dash", "DT", "shinyWidgets", "readr",
    "GetoptLong", "gplots", "shinyFeedback", "patchwork", "cowplot", "ggsci", "tidyverse",
    "RColorBrewer", "ggpubr", "BiocParallel", "shinypop", "remotes", "circlize", "ggridges",
    "RobustRankAggreg", "ggpubr", "org.Hs.eg.db", "GSVA", "survminer", "pROC", "WGCNA",
    "randomForestSRC", "glmnet", "gbm", "survivalsvm", "plsRcox", "RSpectra", "rARPACK",
    "mixOmics", "survcomp", "superpc", "timeROC", "rms", "dcurves", "snowfall"
  )

  ## Check and install packages from CRAN or Bioconductor
  for (pkg.name in pkgs.to.check) {
    if (!requireNamespace(pkg.name, quietly = TRUE)) {
      message(paste("[++] Installing", pkg.name))
      tryCatch(
        pacman::p_install(pkg.name, character.only = TRUE),
        error = function(e) {
          message(paste("[++] ERROR: R package [", pkg.name, "] can NOT be installed."))
        }
      )
    }
    load_package(pkg.name)
  }

  message("[++] All required packages are loaded successfully.")
}
