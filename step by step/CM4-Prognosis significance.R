rm(list = ls())
library(tidyverse)
library(ggsci)
library(cowplot)
library(ggpubr)
options(digits = 2)

select_cancer <- 'CRC'
load(paste0('data/',select_cancer,'/symbol.rda'))
rm(total_expr_list)
load('score_list.rda')
Sur_ids <- c('OS','OS.time','RFS','RFS.time','DFS','DFS.time',
             'PFS','PFS.time','DSS','DSS.time') ## 生存变量ID
total_clin_list <- total_clin_list[all_cohort_name]
names <- names(total_clin_list)

# -------------------------------------------------------------------------
### 汇总生存变量 ###
# -------------------------------------------------------------------------

clin_var_list <- lapply(total_clin_list,function(x){colnames(x)[-1]})
clin_vars <- data.frame()
for (i in names) {
  clin_vars <- rbind(clin_vars,data.frame(Cohort=i,Var=clin_var_list[[i]]))
}
clin_vars <- clin_vars[order(clin_vars$Var),]
Sur_vars <- clin_vars[clin_vars$Var%in%Sur_ids,]
Sur_ids <- unique(Sur_vars$Var)

# -------------------------------------------------------------------------
### Prepare score ###
# -------------------------------------------------------------------------

score_list <- lapply(score_list,function(x){
  x <- data.frame(ID=rownames(x),Score=x[,3])
  return(x)
})

# -------------------------------------------------------------------------
### Prognostic Significance ###
# -------------------------------------------------------------------------

Color_palette <- 'jco' ## 选择颜色面板
alpha <- 0.8 ## 选择颜色透明度

if(Color_palette=='npg'){
  cols <- pal_npg(alpha = alpha)(10)
}
if(Color_palette=='nejm'){
  cols <- pal_nejm(alpha = alpha)(8)
}
if(Color_palette=='jco'){
  cols <- pal_jco(alpha = alpha)(10)
}
if(Color_palette=='d3'){
  cols <- pal_d3(alpha = alpha)(10)
}
if(Color_palette=='lancet'){
  cols <- pal_lancet(alpha = alpha)(9)
}
if(Color_palette=='jama'){
  cols <- pal_jama(alpha = alpha)(7)
}

cutoff <- 'optimal' ## Please input cutoff method! The method include 'median','mean','quantile','optimal','custom', just pick one
percent <- 0.3 ## If users select the "custom" cutoff method, the percent option will appearand are set to a range of 0.1 to 0.9

source('lzq_survplot.R')

Sur_ids2 <- gsub('.time','',Sur_ids)%>%unique()
Sur_vars$Var <- gsub('.time','',Sur_vars$Var)
Sur_vars <- distinct(Sur_vars,Cohort,Var,.keep_all = T)

plotsur_list <- list()
coxr <- data.frame()
for (i in Sur_ids2) {
  for (j in Sur_vars$Cohort[Sur_vars$Var==i]) {
    tmp <- na.omit(merge(total_clin_list[[j]][,c('ID',i,paste0(i,'.time'))],score_list[[j]],by=1)[,-1])
    tmp <- as.data.frame(apply(tmp,2,as.numeric))
    fit <- summary(coxph(gsurv(select_survars = i),tmp))
    coxr <- rbind(coxr,data.frame(Sur_var=i,Cohort=j,HR=fit$conf.int[,1],
                                  HRL=fit$conf.int[,3],HRH=fit$conf.int[,4],P=fit$coefficients[,5]))
    plotsur_list[[paste0(i,'-',j)]] <- lzq_survplot2(Sur_ids = colnames(tmp),cohort = j,data = tmp,gene = 'Score',Input = self_name,
                                                     cols = cols,cutoff = cutoff,percent = percent)
  }
}

plot_grid(plotlist=plotsur_list,nrow = 2)

lzq_coxplot(coxr)
##保存森林图片注意高度要根据队列的y的数量来定 coxr的行数
