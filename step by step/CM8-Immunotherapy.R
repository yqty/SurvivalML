rm(list = ls())
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(ggsci)
library(cowplot)
library(pROC)
library(survminer)
library(survival)
options(digits = 2)

load('ICI_data-symbol.rda')
source('lzq_refit.R')
load('score_list.rda')

# -------------------------------------------------------------------------
### select or input your gene list (Symbol)###
# -------------------------------------------------------------------------

#Gene_Input <- c("CCL5", "CD27", "CD274", "CD276", "CD8A", "CMKLR1", "CXCL9", "HLA-DQA1", "HLA-DRB1", "HLA-E", "IDO1", "LAG3", "NKG7", "PDCD1LG2", "PSMB10", "STAT1", "TIGIT")
#self_name <- 'TIS' ## 自定义Score的名字 **通用**

# -------------------------------------------------------------------------
### 看看多少个队列有这些基因，少于一半输入的基因时，则不使用这个队列 ###
# -------------------------------------------------------------------------

index <- sapply(ICI_data,function(x){mean(Gene_Input%in%rownames(x$expr))==1})
select_cohorts <- ICI_data[index]
select_cohorts2 <- lapply(select_cohorts,function(x){as.data.frame(t(x$expr[Gene_Input,]))})

# -------------------------------------------------------------------------
### calculating risk score ###
# -------------------------------------------------------------------------

score_list <- lzq_refit(testdata = select_cohorts2)

# -------------------------------------------------------------------------
# plot boxplot
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

plotbox_list <- list()
for (i in names(select_cohorts)) {
  dd <- merge(ICI_data[[i]]$clin[,c('ID','Response')],score_list[[i]],by=1)%>%na.omit()
  test_type <- ifelse(shapiro.test(dd[,'score'])$p.val<0.05,'wilcox.test','t.test')
  plotbox_list[[i]] <- ggplot(dd,aes_string('Response','score'))+
    geom_jitter(shape = 21,size=2,width = 0.2,aes_string(fill='Response',color='Response'))+
    geom_boxplot(outlier.colour = NA,aes_string(fill='Response'),color='black',size=0.6,alpha=0.65)+
    geom_violin(alpha=0.5,aes_string(fill='Response'),color=NA,trim = T)+
    stat_compare_means(label.y = max(dd[,'score'])*1.05,method = test_type,color='black',size=5)+
    scale_fill_manual(values = cols)+
    scale_color_manual(values = cols)+
    expand_limits(y=max(dd[,'score'])*1.1)+
    theme_bw(base_rect_size = 2)+
    labs(y=self_name,title = paste0(gsub(' \\(.*\\)','',unlist(strsplit(ICI_data[[i]]$ann,'---'))[1]),'\n', gsub('\\(|\\)','',str_extract(unlist(strsplit(ICI_data[[i]]$ann,'---'))[1],'\\(.*\\)'))))+
    theme(axis.text.y = element_text(size = 11, colour = 'black'),
          axis.text.x = element_text(size = 14,colour = 'black'),
          axis.title.x = element_text(size = 15,colour = 'darkred',face='bold'),
          axis.title.y = element_text(size = 15,colour = 'darkred',face='bold'),
          panel.grid = element_blank(),
          panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
          panel.background = element_rect(fill='#f3f6f6'),
          plot.title = element_text(hjust = 0.5,size = 15,colour = 'darkred',face='bold'),
          legend.position = 'none',
          axis.ticks.x = element_blank())
}

plot_grid(plotlist=plotbox_list,nrow = 3)
plotbox_list$GSE100797

# -------------------------------------------------------------------------
# plot ROC

auc_cols <- c(pal_npg()(10),pal_d3()(10)[-c(3:6)],pal_jco()(3))[index]
names(auc_cols) <- names(score_list)

plotroc_list <- list()
for (i in names(select_cohorts)) {
  dd <- merge(ICI_data[[i]]$clin[,c('ID','Response')],score_list[[i]],by=1)%>%na.omit()
  fit <- roc(dd[,2],dd[,3],auc=T)
  rr <- data.frame(x=1-fit$specificities,y=fit$sensitivities)
  rr <- rr[order(rr$y,rr$x),]
  
  plotroc_list[[i]] <- ggplot(rr,aes(x,y))+
    geom_line(size=1,color=auc_cols[[i]])+
    labs(x='1-Specificity',y='Sensitivity',color=NULL)+
    theme_bw(base_rect_size = 2)+
    geom_abline(slope = 1,color='grey70')+
    scale_x_continuous(expand = c(0.01,0.01))+
    scale_y_continuous(expand = c(0.01,0.01))+
    labs(title = unlist(strsplit(ICI_data[[i]]$ann,'---'))[1],
         subtitle = unlist(strsplit(ICI_data[[i]]$ann,'---'))[2])+
    theme(axis.text.y = element_text(size = 11, colour = 'black'),
          axis.text.x = element_text(size = 11,colour = 'black'),
          axis.title.x = element_text(size = 15,colour = 'darkred',face='bold'),
          axis.title.y = element_text(size = 15,colour = 'darkred',face='bold'),
          panel.grid = element_blank(),
          panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
          panel.background = element_rect(fill='#f3f6f6'),
          plot.title = element_text(hjust = 0.5,size = 15,colour = 'darkred',face='bold'),
          plot.subtitle = element_text(hjust = 0.5,size = 11,colour = 'black'),
          legend.text = element_text(size=12),
          legend.key = element_blank(),
          legend.background = element_blank(),
          legend.position = c(0.995,0.03),
          legend.justification = c(1,0))+
    annotate('text',x = 0.75,y = 0.1,label=paste0('AUC = ',sprintf("%.3f",fit$auc)),size=4.8,color='black')
}

plot_grid(plotlist=plotroc_list,nrow = 3)
plotroc_list$GSE100797

# -------------------------------------------------------------------------
# plot survival
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

os_cohort <- os_cohort[os_cohort%in%names(select_cohorts)]
pfs_cohort <- pfs_cohort[pfs_cohort%in%names(select_cohorts)]

if(length(os_cohort)!=0){
  plotos_list <- list()
  for (i in os_cohort) {
    dd <- merge(ICI_data[[i]]$clin[,c('ID','OS','OS.time')],score_list[[i]],by=1)%>%na.omit()
    plotos_list[[i]] <- tryCatch(lzq_survplot2(Sur_ids = c('OS','OS.time'),cohort = unlist(strsplit(ICI_data[[i]]$ann,'---'))[1],
                                               data = dd,gene = 'score',Input = self_name,
                                               cols = cols,cutoff = cutoff,percent = percent),error=function(e)NA)
  }
}

if(length(pfs_cohort)!=0){
  plotpfs_list <- list()
  for (i in pfs_cohort) {
    dd <- merge(ICI_data[[i]]$clin[,c('ID','PFS','PFS.time')],score_list[[i]],by=1)%>%na.omit()
    plotpfs_list[[i]] <- tryCatch(lzq_survplot2(Sur_ids = c('PFS','PFS.time'),cohort = unlist(strsplit(ICI_data[[i]]$ann,'---'))[1],
                                                data = dd,gene = 'score',Input = self_name,
                                                cols = cols,cutoff = cutoff,percent = percent),error=function(e)NA)
    
  }
}

if(length(pfs_cohort)==0&length(os_cohort)==0){
  message('No available datasets with prognostic information')
}


