# -----------------------------------------------------------------------------
# load_packages.R – Optimised dependency loader for SurvivalML Shiny app
# 2025‑05‑13  (rev‑c – fix BiocParallel binary; add CoxBoost & summaryBox from GitHub)
# -----------------------------------------------------------------------------
# 1. Mirror settings -----------------------------------------------------------
options("repos" = c(CRAN = "https://mirrors.westlake.edu.cn/CRAN/"))
options(BioC_mirror = "https://mirrors.westlake.edu.cn/bioconductor")

msg <- function(...) message("[load_packages] ", ...)
msg("Starting dependency check – this may take a while on first run …")

# 2. Helper utilities ----------------------------------------------------------
ensure_ns <- function(pkg) requireNamespace(pkg, quietly = TRUE)

install_cran_bioc <- function(pkg) {
  if (!ensure_ns(pkg)) {
    msg("Installing ", pkg, " from CRAN / Bioconductor …")
    if (!ensure_ns("BiocManager")) install.packages("BiocManager")
    tryCatch(
      BiocManager::install(pkg, update = FALSE, ask = FALSE, type = "binary"),
      error = function(e) {
        msg("Bioc install failed – fallback to CRAN for ", pkg)
        install.packages(pkg)
      }
    )
  }
}

install_github_pkg <- function(pkg, repo) {
  if (!ensure_ns(pkg)) {
    msg("Installing ", pkg, " from GitHub (", repo, ") …")
    if (!ensure_ns("remotes")) install.packages("remotes")
    remotes::install_github(repo, upgrade = "never")
  }
}

load_pkg <- function(pkg) {
  suppressPackageStartupMessages(
    library(pkg, character.only = TRUE, quietly = TRUE)
  )
}

# 3. Pre‑install BiocParallel binary ------------------------------------------
if (!ensure_ns("BiocParallel")) {
  msg("Installing BiocParallel binary …")
  if (!ensure_ns("BiocManager")) install.packages("BiocManager")
  tryCatch(
    BiocManager::install("BiocParallel", update = FALSE, ask = FALSE, type = "binary"),
    error = function(e) msg("BiocParallel binary install failed: ", e$message)
  )
}

# 4. GitHub packages -----------------------------------------------------------
install_github_pkg("ggcor", "Github-Yilei/ggcor")
install_github_pkg("shinypop", "dreamRs/shinypop")   # CRAN/Bioc 无此包
install_github_pkg("CoxBoost", "binderh/CoxBoost")
install_github_pkg("summaryBox", "deepanshu88/summaryBox")

# 5. Bioconductor heavyweights & core data ------------------------------------
bioc_pkgs <- c("GenomeInfoDbData", "GenomeInfoDb", "GO.db", "clusterProfiler")
lapply(bioc_pkgs, install_cran_bioc)

# 6. CRAN / Bioc routine packages --------------------------------------------
cran_pkgs <- c(
  "shiny", "shinyjs", "bslib", "waiter", "openxlsx", "plotly", "ggplot2",
  "stringr", "slickR", "shinyvalidate", "markdown", "bs4Dash", "DT",
  "shinyWidgets", "readr", "GetoptLong", "gplots", "shinyFeedback",
  "patchwork", "cowplot", "ggsci", "tidyverse", "RColorBrewer", "ggpubr",
  "remotes", "circlize", "ggridges", "RobustRankAggreg", "org.Hs.eg.db",
  "GSVA", "survminer", "pROC", "WGCNA", "randomForestSRC", "glmnet", "gbm",
  "survivalsvm", "plsRcox", "RSpectra", "rARPACK", "mixOmics", "survcomp",
  "superpc", "timeROC", "rms", "dcurves", "snowfall", "bsplus",
  "ComplexHeatmap"
)

lapply(cran_pkgs, install_cran_bioc)

# 7. Load everything silently --------------------------------------------------
all_pkgs <- c("ggcor", "shinypop", "CoxBoost", "summaryBox", "BiocParallel", bioc_pkgs, cran_pkgs)

invisible(lapply(all_pkgs, load_pkg))

msg("All packages are ready – happy Shiny! 🎉")
