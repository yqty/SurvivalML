# ---------------------------------------------------------------------------
# load_packages.R – SurvivalML dependency loader  (rev‑f 2025‑05‑14)
#   * fix: GitHub install loop (names() bug)                                     
# ---------------------------------------------------------------------------

## --0.  GitHub‑only packages -------------------------------------------------
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")

gh_pkgs <- c(
  CoxBoost   = "binderh/CoxBoost",
  shinypop   = "dreamRs/shinypop",
  summaryBox = "deepanshu88/summaryBox",
  ggcor      = "Github-Yilei/ggcor"
)

for (pkg in names(gh_pkgs)) {
  if (!requireNamespace(pkg, quietly = TRUE))
    remotes::install_github(gh_pkgs[[pkg]], upgrade = "never")
}

## --1.  Mirror settings ------------------------------------------------------
cran.mirror <- "https://mirrors.westlake.edu.cn/CRAN/"
options(repos = c(CRAN = cran.mirror))                       # 只改 CRAN 行
options(BioC_mirror = "https://mirrors.westlake.edu.cn/bioconductor")
options(install.packages.compile.from.source = "always")     # 无 binary 时回退源码

msg <- function(...) message("[load_packages] ", ...)
ensure_ns <- function(pkg) requireNamespace(pkg, quietly = TRUE)

## helper: binary → source fallback ------------------------------------------
install_bioc <- function(pkgs) {
  if (!ensure_ns("BiocManager")) install.packages("BiocManager")
  for (p in pkgs) {
    if (ensure_ns(p)) next
    # --- try binary first ----------------------------------------------------
    tryCatch({
      BiocManager::install(p, update = FALSE, ask = FALSE, type = "binary")
    }, error = function(e) {
      msg("Binary install failed for ", p, ": ", e$message)
    })
    # --- if still missing → compile from source -----------------------------
    if (!ensure_ns(p)) {
      msg("Compiling ", p, " from source …")
      BiocManager::install(p, update = FALSE, ask = FALSE, type = "source")
    }
  }
}

## --2.  BiocParallel (需 binary) --------------------------------------------
if (!ensure_ns("BiocParallel")) install_bioc("BiocParallel")

## --3.  Data & heavy Bioconductor packages ----------------------------------
bioc_pkgs <- c("GenomeInfoDbData","GenomeInfoDb","GO.db",
               "org.Hs.eg.db","clusterProfiler")
install_bioc(bioc_pkgs)

## --4.  Routine CRAN / Bioc packages ---------------------------------------
cran_pkgs <- c(
  "shiny","shinyjs","bslib","waiter","openxlsx","plotly","ggplot2",
  "stringr","slickR","shinyvalidate","markdown","bs4Dash","DT",
  "shinyWidgets","readr","GetoptLong","gplots","shinyFeedback",
  "patchwork","cowplot","ggsci","tidyverse","RColorBrewer","ggpubr",
  "circlize","ggridges","RobustRankAggreg","GSVA","survminer","pROC",
  "WGCNA","randomForestSRC","glmnet","gbm","survivalsvm","plsRcox",
  "RSpectra","rARPACK","mixOmics","survcomp","superpc","timeROC",
  "rms","dcurves","snowfall","bsplus","ComplexHeatmap"
)
install_bioc(cran_pkgs)

## --5.  Silent load ---------------------------------------------------------
all_pkgs <- c(names(gh_pkgs),"BiocParallel",bioc_pkgs,cran_pkgs)

invisible(lapply(all_pkgs, function(p)
  suppressPackageStartupMessages(library(p, character.only = TRUE))))
msg("All packages are ready – happy Shiny! 🎉")
