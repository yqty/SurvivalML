rm(list = ls())
library(tidyverse)
library(ggcor)
library(circlize)
options(digits = 2)

select_cancer <- 'LUAD'
load(paste0('data/',select_cancer,'/ims_expression.rda'))

# -------------------------------------------------------------------------
### Prepare score ###
# -------------------------------------------------------------------------

load('score_list.rda')
score_list <- lapply(score_list,function(x){
  x <- data.frame(ID=rownames(x),Score=x[,1],row.names = rownames(x))
  return(x)
})

ims_expr_list <- ims_expr_list[names(score_list)]

cor_list <- map2(score_list,ims_expr_list,function(x,y){
  id <- intersect(rownames(x),colnames(y))
  x <- x[id,]
  y <- y[,id]
  l <- cor(t(y),x$Score)%>%as.data.frame()%>%tibble::rownames_to_column('ID')
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
### plot heatmap ###
# -------------------------------------------------------------------------

cols <- colorRamp2(c(-1,-0.5,0,0.5,1),c('#0077b6','#48cae4','white',alpha('#e3170a',0.5),'#e3170a'))

if(ncol(cor_data)==1){
  ims <- ims[ims$V3%in%rownames(cor_data),]
  tmp <- cor_data[match(ims$V3,rownames(cor_data)),]%>%as.data.frame()
  rownames(tmp) <- ims$V3
  colnames(tmp) <- colnames(cor_data)
}else{
  tmp <- cor_data[match(ims$V3,rownames(cor_data)),]
}

left_ann <- HeatmapAnnotation(Immunomodulators = ims$V2,which = 'row',
                              show_legend = T,show_annotation_name = F,
                              simple_anno_size =unit(2,'mm'),
                              col = list(Immunomodulators = c('Antigen presentation' = pal_nejm(alpha = 0.8)(8)[1],
                                                              'Immunoinhibitor' = pal_nejm(alpha = 0.8)(8)[2],
                                                              'Immunostimulator' = pal_nejm(alpha = 0.8)(8)[3],
                                                              'Chemokine' = pal_nejm(alpha = 0.8)(8)[4],
                                                              'Receptor' = pal_nejm(alpha = 0.8)(8)[5]
                              )))

tryCatch(print(Heatmap(as.matrix(tmp),na_col = NULL,
                       name = 'Correlation',border = F,
                       col = cols,
                       row_names_side = 'left',
                       row_split = rep(c(1,2,3,4,5),times=as.numeric(table(ims$V2))),
                       row_gap = unit(5,'mm'),
                       left_annotation = left_ann,
                       row_title = NULL,
                       width= ncol(tmp) * unit(5.5, "mm"),
                       height= nrow(tmp) * unit(5, "mm"),
                       cluster_rows = T,show_row_dend = F)),
         error=function(e){
           print(Heatmap(as.matrix(tmp),na_col = NULL,
                         name = 'Correlation',border = F,
                         col = cols,
                         row_names_side = 'left',
                         row_split = rep(c(1,2,3,4,5),times=as.numeric(table(ims$V2))),
                         row_gap = unit(5,'mm'),
                         left_annotation = left_ann,
                         row_title = NULL,
                         width= ncol(tmp) * unit(5.5, "mm"),
                         height= nrow(tmp) * unit(5, "mm"),
                         cluster_rows = F,show_row_dend = F))
         })

Cohort <- 'TCGA_LUAD'
IMS_ID <- 'CD274'
if(!IMS_ID%in%rownames(ims_expr_list[[Cohort]])){
  message(paste0(Cohort,' absents ',IMS_ID,', please change another gene with red or blue color!'))
}
IMS_ID2 <- IMS_ID

cor_tmp <- merge(score_list[[Cohort]],t(ims_expr_list[[Cohort]][IMS_ID,]),by.x=1,by.y=0)[,-1]
cor_tmp <- na.omit(cor_tmp)
colnames(cor_tmp) <- gsub('-','_',colnames(cor_tmp))
IMS_ID <- gsub('-','_',IMS_ID)

fit <- cor.test(cor_tmp[,1],cor_tmp[,2],method = 'spearman')

ggplot(cor_tmp,aes_string('Score',IMS_ID)) +
  geom_point(col='#ef476f',alpha=0.9,size=2) +
  geom_smooth(method=lm, se=T,na.rm=T, fullrange=T,size=1.8,col="#9a8c98") +
  geom_rug(col="#f9c74f",size=1) +
  theme_bw(base_rect_size = 1.5) +
  labs(x=self_name,y=paste0(IMS_ID2,' Expression (z-score)'),title=Cohort,
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









