rm(list = ls())
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(ggsci)
library(ggcor)
library(RobustRankAggreg)
options(digits = 2)

select_cancer <- 'LUAD'
load(paste0('data/',select_cancer,'/Drugdata.rda'))

# -------------------------------------------------------------------------
### Prepare score ###
# -------------------------------------------------------------------------

load('score_list.rda')
score_list <- lapply(score_list,function(x){
  x <- data.frame(ID=rownames(x),Score=x[,3])
  return(x)
})

Drug_list <- lapply(Drug_list,function(x){x[names(score_list)]})

# -------------------------------------------------------------------------
##for example GDSC1
dd <- Drug_list$PRISM

cor_list <- map2(score_list,dd,function(x,y){
  ID <- intersect(x$ID,rownames(y))
  x <- x[match(ID,x$ID),]
  y <- y[ID,]
  l <- cor(y,x$Score)%>%as.data.frame()%>%tibble::rownames_to_column('ID')
  colnames(l)[2] <- 'Cor'
  return(l)
})

if(length(names(cor_list))>1){
  cor_data <- Reduce(function(x,y){merge(x,y,by=1,all = T)},cor_list)%>%
    tibble::column_to_rownames('ID')
}

if(length(names(cor_list))==1){
  cor_data <- cor_list[[1]]
  rownames(cor_data) <- NULL
  cor_data <- cor_data%>%tibble::column_to_rownames('ID')
}

colnames(cor_data) <- names(score_list)

# -------------------------------------------------------------------------

cordata <- cor_data
poslist <- lapply(as.list(cordata),function(x){
  x <- data.frame(ID=rownames(cordata),R=x)%>%na.omit()
  x <- x[x$R>0,]%>%arrange(desc(R))
  return(x$ID)
})
neglist <- lapply(as.list(cordata),function(x){
  x <- data.frame(ID=rownames(cordata),R=x)%>%na.omit()
  x <- x[x$R<0,]%>%arrange(R)
  return(x$ID)
})

pos_agg=aggregateRanks(poslist)
pos_agg <- pos_agg[pos_agg$Score<0.05,]
neg_agg=aggregateRanks(neglist)
neg_agg <- neg_agg[neg_agg$Score<0.05,]

max_num <- 30 ##可以自己设置
pos_num <- ifelse(nrow(pos_agg)<max_num,nrow(pos_agg),max_num)
neg_num <- ifelse(nrow(neg_agg)<max_num,nrow(neg_agg),max_num)

if(ncol(cordata)!=1){
  posdata <- cordata[pos_agg$Name[1:pos_num],]
  negdata <- cordata[neg_agg$Name[1:neg_num],]
}else{
  posdata <- cordata[pos_agg$Name[1:pos_num],]%>%as.data.frame()
  rownames(posdata) <- pos_agg$Name[1:pos_num]
  colnames(posdata) <- names(total_expr_list)[show_cohort_index]
  negdata <- cordata[neg_agg$Name[1:neg_num],]%>%as.data.frame()
  rownames(negdata) <- neg_agg$Name[1:neg_num]
  colnames(negdata) <- names(total_expr_list)[show_cohort_index]
}

# -------------------------------------------------------------------------

cols <- colorRamp2(c(-1,-0.5,0,0.5,1),c('#0077b6','#48cae4','white',alpha('#e3170a',0.5),'#e3170a'))

right_ann1 <- HeatmapAnnotation(type = rep('A',nrow(posdata)),which = 'row',
                                show_legend = F,show_annotation_name = F,
                                simple_anno_size =unit(2,'mm'),
                                col = list(type = c('A' = '#e9c46a','B' = '#2a9d8f')))
right_ann2 <- HeatmapAnnotation(type = rep('B',nrow(negdata)),which = 'row',
                                show_legend = F,show_annotation_name = F,
                                simple_anno_size =unit(2,'mm'),
                                col = list(type = c('A' = '#e9c46a','B' = '#2a9d8f')))
library(showtext)
showtext_auto(enable = TRUE)
p1 <- Heatmap(as.matrix(posdata),na_col = NULL,
              name = 'Correlation',border = F,
              col = cols,
              row_names_side = 'left',
              row_title_side = 'right',
              row_title = '← High expression indicates resistance →',
              row_title_gp = gpar(fontfamily='SimSun'),
              cluster_rows = T,
              right_annotation = right_ann1,
              width= ncol(posdata) * unit(5.5, "mm"),
              height= nrow(posdata) * unit(5, "mm"),
              show_row_dend = F,show_column_dend = F)

p2 <- Heatmap(as.matrix(negdata),na_col = NULL,
              name = 'Correlation',border = F,
              col = cols,
              row_names_side = 'left',
              row_title_side = 'right',
              row_title = '← High expression indicates sensitivity →',
              row_title_gp = gpar(fontfamily='SimSun'),
              cluster_rows = T,
              right_annotation = right_ann2,
              width= ncol(negdata) * unit(5.5, "mm"),
              height= nrow(negdata) * unit(5, "mm"),
              show_row_dend = F,show_column_dend = F)
draw(p1%v%p2)


## 鼠标箭头放上去有相关性系数大小，点击格子会出来相关性散点图

Database <- 'GDSC1'
Cohort <- 'TCGA_LUAD'
Drug_ID <- 'docetaxel'

cor_tmp <- merge(score_list[[Cohort]],cbind(rownames(Drug_list[[Database]][[Cohort]]),Drug_list[[Database]][[Cohort]][,Drug_ID]),by=1)
cor_tmp <- na.omit(cor_tmp)
cor_tmp$V2 <- as.numeric(cor_tmp$V2)

fit <- cor.test(cor_tmp[,2],cor_tmp[,3],method='spearman')

ggplot(cor_tmp,aes_string('Score','V2')) +
  geom_point(col='#6a4c93',alpha=0.8,size=2) +
  geom_smooth(method=lm, se=T,na.rm=T, fullrange=T,size=1.8,col="#1982c4") +
  geom_rug(col="#3d405b",size=1) +
  theme_bw(base_rect_size = 1.5) +
  labs(x=self_name,y=paste0('IC50 of ',Drug_ID,' (',Database,')'),title=Cohort,
       subtitle = paste0('Cor = ',sprintf("%.3f", fit$estimate),'; Pval = ',format(fit$p.value, scientific = T))) +
  theme(plot.title = element_text(hjust = 0.5,size = 14,colour = 'darkred',face='bold'),
        plot.subtitle = element_text(hjust = 0.5,size = 11,colour = 'black'),
        axis.text = element_text(size = 10,colour = 'black'),
        axis.title.x = element_text(size = 13,colour = 'darkred',face='bold'),
        axis.title.y = element_text(size = 13,colour = 'darkred',face='bold'),
        axis.ticks = element_line(size=1),
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.background = element_rect(fill='#f3f6f6'))

