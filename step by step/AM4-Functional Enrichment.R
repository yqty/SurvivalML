rm(list = ls())
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(RColorBrewer)
library(cowplot)
options(digits = 2)

select_cancer <- 'GBM'
load(paste0('data/',select_cancer,'/mRNA.rda'))

# -------------------------------------------------------------------------
### Prepare score ###
# -------------------------------------------------------------------------

load('score_list.rda')
score_list <- lapply(score_list,function(x){
  x <- data.frame(ID=rownames(x),Score=x[,1],row.names = rownames(x))
  return(x)
})

# -------------------------------------------------------------------------
### 计算该评分与其它基因在每个队列的相关性，取平均值 ###
# -------------------------------------------------------------------------

total_expr_list <- total_expr_list[names(score_list)]
cor_list <- map2(score_list,total_expr_list,function(x,y){
  y <- y[,x$ID]
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

cor_data$mean_cor <- rowMeans(cor_data,na.rm = T)
cor_data <- cor_data[order(cor_data$mean_cor,decreasing = T),]
cordata <- data.frame(ID=rownames(cor_data),mean_cor=cor_data$mean_cor)

ID <- bitr(cordata$ID,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)
cordata <- merge(ID,cordata,by=1)%>%dplyr::arrange(desc(mean_cor))
cordata <- cordata[cordata$mean_cor<2,]%>%na.omit()

# -------------------------------------------------------------------------
### Over Representation Analysis ###
# -------------------------------------------------------------------------
##可自定义选择数量 50-1000
select_pcor_ids <- cordata$ENTREZID[2:501]
select_ncor_ids <- cordata$ENTREZID[nrow(cordata):(nrow(cordata)-499)]

# GO analysis-------------------------------------------------------------------------

source('lzq_Enrichment.R')

## qvalueCutoff ont 可以自定义
GO_pos <- enrichGO(gene = select_pcor_ids,OrgDb = org.Hs.eg.db,keyType = "ENTREZID",ont = 'all',
                   pvalueCutoff = 0.05,qvalueCutoff = 0.05,minGSSize = 1,maxGSSize = 50000) 
GO_neg <- enrichGO(gene = select_ncor_ids,OrgDb = org.Hs.eg.db,keyType = "ENTREZID",ont = 'all',
                   pvalueCutoff = 0.05,qvalueCutoff = 0.05,minGSSize = 1,maxGSSize = 50000)

tmp <- lzq_GO_merge(GO_pos,GO_neg)
lzq_GO_plot(tmp,width = 0.35) ## width可以自己设置
# 显示一行注释 红色是正相关基因的富集通路，蓝色是负相关基因的富集通路

# KEGG analysis-------------------------------------------------------------------------

KEGG_pos <- enrichKEGG(gene = select_pcor_ids,organism = 'hsa',pvalueCutoff = 0.05,qvalueCutoff = 0.05,
                       minGSSize = 1,maxGSSize = 50000)
KEGG_neg <- enrichKEGG(gene = select_ncor_ids,organism = 'hsa',pvalueCutoff = 0.05,qvalueCutoff = 0.05,
                       minGSSize = 1,maxGSSize = 50000)

tmp <- lzq_KEGG_merge(KEGG_pos,KEGG_neg)
lzq_KEGG_plot(tmp,width = 0.3) ## width可以自己设置
# 显示一行注释 红色是正相关基因的富集通路，蓝色是负相关基因的富集通路

# -------------------------------------------------------------------------
### Gene set enrichment analysis ###
# -------------------------------------------------------------------------

load('GSEA_genesets.rda')
kegg_sets <- kegg_sets[kegg_sets$gene%in%cordata$ENTREZID,]
kegg_sets <- kegg_sets[kegg_sets$term%in%names(table(kegg_sets$term)>1),]
go_sets <- go_sets[go_sets$gene%in%cordata$ENTREZID,]
go_sets <- go_sets[go_sets$term%in%names(table(go_sets$term)>1),]
hall_sets <- hall_sets[hall_sets$gene%in%cordata$ENTREZID,]
hall_sets <- hall_sets[hall_sets$term%in%names(table(hall_sets$term)>1),]

gene_rank_list <- cordata$mean_cor
names(gene_rank_list) <- cordata$ENTREZID

## GSEA-GO
GSEA_GO <- GSEA(gene_rank_list,pvalueCutoff = 1,TERM2GENE = go_sets)
GSEA_GO@result <- lzq_GSEA_merge(GSEA_GO)[[1]]
cate <- lzq_GSEA_merge(GSEA_GO)[[2]]+lzq_GSEA_merge(GSEA_GO)[[3]]

ridgeplot(GSEA_GO,showCategory = cate,
          fill = "p.adjust",
          core_enrichment = TRUE,
          label_format = 80)+
  geom_vline(xintercept = 0,linetype=2,color='grey30',alpha=0.7)+
  theme_classic(base_rect_size = 2)+
  labs(x='Enrichment Score',y=NULL,title = 'GSEA-GO Analysis')+
  scale_x_continuous(expand = c(0.01,0.01))+
  scale_fill_gradient(low = alpha('#21b6af', 0.7),high = alpha('#eeba4d', 0.7))+
  theme(axis.text.x = element_text(size = 10,colour = 'black'),
        axis.text.y = element_text(size = 12,colour = 'black'),
        axis.title.x = element_text(size = 13,colour = 'darkred',face='bold'),
        axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 14,colour = 'darkred',face='bold'),
        legend.title = element_text(size=13,colour = 'black'),
        legend.text = element_text(size=12,colour = 'black'))


## GSEA-KEGG
GSEA_KEGG <- GSEA(gene_rank_list,pvalueCutoff = 1,TERM2GENE = kegg_sets)
GSEA_KEGG@result <- lzq_GSEA_merge(GSEA_KEGG)[[1]]
cate <- lzq_GSEA_merge(GSEA_KEGG)[[2]]+lzq_GSEA_merge(GSEA_KEGG)[[3]]

ridgeplot(GSEA_KEGG,showCategory = cate,
          fill = "p.adjust",
          core_enrichment = TRUE,
          label_format = 80)+
  geom_vline(xintercept = 0,linetype=2,color='grey30',alpha=0.7)+
  theme_classic(base_rect_size = 2)+
  labs(x='Enrichment Score',y=NULL,title = 'GSEA-KEGG Analysis')+
  scale_x_continuous(expand = c(0.01,0.01))+
  scale_fill_gradient(low = alpha('#21b6af', 0.7),high = alpha('#eeba4d', 0.7))+
  theme(axis.text.x = element_text(size = 10,colour = 'black'),
        axis.text.y = element_text(size = 12,colour = 'black'),
        axis.title.x = element_text(size = 13,colour = 'darkred',face='bold'),
        axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 14,colour = 'darkred',face='bold'),
        legend.title = element_text(size=13,colour = 'black'),
        legend.text = element_text(size=12,colour = 'black'))

## GSEA-Hallmark
GSEA_hall <- GSEA(gene_rank_list,pvalueCutoff = 1,TERM2GENE = hall_sets)
GSEA_hall@result <- lzq_GSEA_merge(GSEA_hall)[[1]]
cate <- lzq_GSEA_merge(GSEA_hall)[[2]]+lzq_GSEA_merge(GSEA_hall)[[3]]

ridgeplot(GSEA_hall,showCategory = cate,
          fill = "p.adjust",
          core_enrichment = TRUE,
          label_format = 80)+
  geom_vline(xintercept = 0,linetype=2,color='grey30',alpha=0.7)+
  theme_classic(base_rect_size = 2)+
  labs(x='Enrichment Score',y=NULL,title = 'GSEA-Hallmark Analysis')+
  scale_x_continuous(expand = c(0.01,0.01))+
  scale_fill_gradient(low = alpha('#21b6af', 0.7),high = alpha('#eeba4d', 0.7))+
  theme(axis.text.x = element_text(size = 10,colour = 'black'),
        axis.text.y = element_text(size = 12,colour = 'black'),
        axis.title.x = element_text(size = 13,colour = 'darkred',face='bold'),
        axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 14,colour = 'darkred',face='bold'),
        legend.title = element_text(size=13,colour = 'black'),
        legend.text = element_text(size=12,colour = 'black'))

# -------------------------------------------------------------------------
### Select pathway ID from enriched pathways ###
# -------------------------------------------------------------------------

Pathway_ID <- 'Dopamine transport'

source('lzq_gseaplot.R')
options(digits = 2)
lzq_gseaplot(GSEA_GO,Pathway_ID,'firebrick3')
ggsave(filename = 'gsea19.pdf',width = 4,height = 4)




