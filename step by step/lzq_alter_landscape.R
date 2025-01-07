lzq_alter_landscape <- function(rank,Input,alterdata){
  rank$Group <- ifelse(rank[,2]>median(rank[,2]),paste0('High ',Input),paste0('Low ',Input))
  
  tmp <- alterdata[,rank$ID]
  tmp[tmp==''] <- NA
  
  Top = HeatmapAnnotation(Group=anno_block(gp = gpar(fill = c('#B0997F','#ffc857')), height = unit(5.5,'mm'),
                                           labels = c(paste0('Low ',Input),paste0('High ',Input)), 
                                           labels_gp = gpar(cex = 0.85, col = "white",fontface='bold')),
                          border = T,
                          show_annotation_name = F)
  
  low <- tmp[,rank$ID[rank$Group==paste0('Low ',Input)]]
  low <- apply(low,1,function(x){sum(!is.na(x))})/ncol(tmp)
  high <- tmp[,rank$ID[rank$Group==paste0('High ',Input)]]
  high <- apply(high,1,function(x){sum(!is.na(x))})/ncol(tmp)
  pct1 <- data.frame(row.names = rownames(tmp),Low=low,High=high)
  
  ll  <- c()
  for (i in rownames(pct1)) {
    dd <- data.frame(x=c(pct1[i,1]*ncol(tmp),(1-pct1[i,1])*ncol(tmp)),
                     y=c(pct1[i,2]*ncol(tmp),(1-pct1[i,2])*ncol(tmp)))
    p <- fisher.test(dd)$p.value
    ll <- c(ll,ifelse(p<0.0001,'****',ifelse(p<0.001,'***',ifelse(p<0.01,'**',ifelse(p<0.05,'*','')))))
  }
  
  low <- tmp[,rank$ID[rank$Group==paste0('Low ',Input)]]
  low <- apply(low,1,function(x){sum(!is.na(x))})/apply(tmp,1,function(x){sum(!is.na(x))})
  high <- tmp[,rank$ID[rank$Group==paste0('High ',Input)]]
  high <- apply(high,1,function(x){sum(!is.na(x))})/apply(tmp,1,function(x){sum(!is.na(x))})
  pct2 <- data.frame(row.names = rownames(tmp),Low=low,High=high)
  
  right_anno <- anno_barplot(as.matrix(pct2),
                             which = "row",
                             border = F,
                             show_annotation_name = F,
                             gp = gpar(fill = c('#B0997F','#ffc857'),border=NA,lty="blank"), 
                             bar_width = 0.6,
                             width = unit(1.8, "cm"),
                             height = unit(1, "cm"))
  right <- rowAnnotation(Percent = right_anno,
                         annotation_name_side="top",
                         annotation_name_rot = 0,
                         annotation_name_gp = gpar(fontsize = 11),
                         ann=anno_text(ll))
  
  print(Heatmap(as.matrix(tmp),
                na_col = 'WhiteSmoke',
                border = T,
                name = ' ',
                col = c('Mut'='#67AB9F','Gain'=alpha('#ff477e',0.8),'Loss'=alpha('#0096c7',0.8)),
                top_annotation = Top,
                right_annotation = right,
                column_split = rev(as.numeric(factor(rank$Group))),
                row_split = factor(cna_ann,levels = c('Mutation','Gain','Loss')),
                cluster_rows = F,row_title = NULL,
                show_column_names = F,cluster_columns = F,column_title = NULL,
                show_row_names = T,row_names_side = 'left',
                row_names_gp = gpar(fontface='italic',fontsize=10),
                width= ncol(tmp) * unit(0.5, "mm"),
                height= nrow(tmp) * unit(4, "mm"),
                show_heatmap_legend = T))
}
