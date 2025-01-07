rm(list = ls())
library(tidyverse)
library(ggsci)
library(cowplot)
library(ggpubr)
options(digits = 2)

select_cancer <- 'PAAD'
load(paste0('data/',select_cancer,'/ensembl.rda'))
rm(total_expr_list)
load('score_list.rda')
Sur_ids <- c('OS','OS.time','RFS','RFS.time','DFS','DFS.time',
             'PFS','PFS.time','DSS','DSS.time') ## 生存变量ID
total_clin_list <- total_clin_list[all_cohort_name]
names <- names(total_clin_list)

# -------------------------------------------------------------------------
### 汇总临床变量 ###
# -------------------------------------------------------------------------

clin_var_list <- lapply(total_clin_list,function(x){colnames(x)[-1]})
clin_vars <- data.frame()
for (i in names) {
  clin_vars <- rbind(clin_vars,data.frame(Cohort=i,Var=clin_var_list[[i]]))
}
clin_vars <- clin_vars[order(clin_vars$Var),]
Sur_ids <- clin_vars[clin_vars$Var%in%Sur_ids,]
clin_vars <- clin_vars[!clin_vars$Var%in%Sur_ids$Var,]
clin_ids <- unique(clin_vars$Var)
clin_ids <- clin_ids[clin_ids!='Tissue']

# -------------------------------------------------------------------------
### Prepare score ###
# -------------------------------------------------------------------------

score_list <- lapply(score_list,function(x){
  x <- data.frame(ID=rownames(x),Score=x[,1])
  return(x)
})

# -------------------------------------------------------------------------
### Clinical association ###
# -------------------------------------------------------------------------

Color_palette <- 'npg' ## 选择颜色面板
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

library(showtext)
showtext_auto()
clin_plot_list <- list()
for (i in clin_ids) {
  plots <- list()
  for (j in clin_vars$Cohort[clin_vars$Var==i]) {
    tmp <- na.omit(merge(total_clin_list[[j]][,c('ID',i)],score_list[[j]],by=1)[,-1])
    tmp[,1] <- as.character(tmp[,1])
    width <- length(unique(tmp[,1]))
    
    if(width>1){
      if(shapiro.test(tmp[,2])$p.val<0.05){
        test_type <- ifelse(width>2,'kruskal.test','wilcox.test')
      }else{
        test_type <- ifelse(width>2,'anova','t.test')}
      
      plots[[j]] <- ggplot(tmp,aes_string(i,'Score'))+
        geom_jitter(shape = 21,size=2,width = 0.2,aes_string(fill=i,color=i))+
        geom_boxplot(outlier.colour = NA,aes_string(fill=i),color='black',size=0.6,alpha=0.65)+
        geom_violin(alpha=0.5,aes_string(fill=i),color=NA,trim = T)+
        stat_compare_means(label.y = max(tmp[,'Score'])*1.05,method = test_type,color='black',size=5)+
        scale_fill_manual(values = cols)+
        scale_color_manual(values = cols)+
        expand_limits(y=max(tmp[,'Score'])*1.1)+
        theme_bw(base_rect_size = 2)+
        labs(y=self_name,title = j)+
        theme(axis.text.y = element_text(size = 11, colour = 'black'),
              axis.text.x = element_text(size = 14,colour = 'black',angle = 50,hjust = 1),
              axis.title.x = element_text(size = 15,colour = 'darkred',face='bold',family = 'STKaiti'),
              axis.title.y = element_text(size = 15,colour = 'darkred',face='bold'),
              panel.grid = element_blank(),
              panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
              panel.background = element_rect(fill='#f3f6f6'),
              plot.title = element_text(hjust = 0.5,size = 15,colour = 'darkred',face='bold'),
              legend.position = 'none',
              axis.ticks.x = element_blank())
      plots[[j]]$width <- width
    }}
  clin_plot_list[[i]] <- plots
}

clin_plot_list <- clin_plot_list[sapply(clin_plot_list,length)!=0]
output_clin_plots <- list()
for (i in names(clin_plot_list)) {
  width_vector <- sapply(clin_plot_list[[i]],function(x){x$width})%>%as.numeric()
  width_vector <- ifelse(width_vector==2,1,ifelse(width_vector==3,1.4,ifelse(width_vector==4,1.8,2.2)))
  #nrow <- ifelse(length(width_vector)<=5,1,ifelse(length(width_vector)<=10,2,ifelse(length(width_vector)<=15,3,4)))
  output_clin_plots[[i]] <- plot_grid(plotlist=clin_plot_list[[i]],rel_widths = width_vector,nrow = 1)
}


output_clin_plots$ACT
output_clin_plots$Age
output_clin_plots$Microsatellite
