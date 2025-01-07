rm(list = ls())
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
source('lzq_alter_landscape.R')

## 二 score
# -------------------------------------------------------------------------
select_cancer <- 'PAAD'
load(paste0('data/',select_cancer,'tcga_omics.rda'))
load('score_list.rda')

cohort <- names(score_list)[grepl('TCGA',names(score_list))]
if(length(cohort)==0){
  message('Missing TCGA dataset in the cohorts where you calculated the score')
}

exp <- score_list[[cohort]]%>%tibble::rownames_to_column('ID')
exp <- exp[order(exp$score),]
exp <- exp[substr(exp$ID,14,16)=='01A',]
exp$ID <- substr(exp$ID,1,12)
exp <- exp[exp$ID%in%colnames(tcga_multiomic),]

## 自定义名字
self_name <- 'TIS'

lzq_alter_landscape(rank = exp[,c('ID','score')],Input = self_name,alterdata = tcga_multiomic)
















