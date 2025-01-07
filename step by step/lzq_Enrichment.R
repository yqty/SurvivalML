lzq_GO_merge <- function(GO_pos,GO_neg){
  d1 <- GO_pos@result%>%filter(ONTOLOGY=='BP')
  d2 <- GO_neg@result%>%filter(ONTOLOGY=='BP')
  
  num1 <- num2 <- 8
  if(nrow(d1)<8){
    num1 <- nrow(d1)
    if(nrow(d2)< 16-num1){
      num2 <- nrow(d2)
    }else{num2 <- 16-num1}
  }else{
    if(nrow(d2)<8){
      num2 <- nrow(d2)
      if(nrow(d1) < 16-num2){
        num1 <- nrow(d1)
      }
    }
  }
  
  if(num1!=0){dd1 <- head(d1,num1);dd1$type <- 'Positive correlation'}else{dd1 <- d1;dd1$type <- NULL}
  if(num2!=0){dd2 <- head(d2,num2);dd2$type <- 'Negative correlation'}else{dd2 <- d2;dd2$type <- NULL}
  
  tmp1 <- rbind(dd1,dd2)%>%dplyr::select(ID,Description,GeneRatio,type,ONTOLOGY)%>%
    tidyr::separate(GeneRatio,into = c('Enriched','Detected'))%>%
    dplyr::rename(Ontology=ONTOLOGY)%>%
    dplyr::mutate(Enriched = as.numeric(Enriched),
                  Detected = as.numeric(Detected),
                  Ontology = ifelse(Ontology=='BP','Biological Process',
                                    ifelse(Ontology=='CC','Cellular Component','Molecular Function')))%>%
    dplyr::arrange(Ontology) %>%
    dplyr::mutate(ID = factor(ID, rev(unique(ID))))
  
  d1 <- GO_pos@result%>%filter(ONTOLOGY=='CC')
  d2 <- GO_neg@result%>%filter(ONTOLOGY=='CC')
  
  num1 <- num2 <- 8
  if(nrow(d1)<8){
    num1 <- nrow(d1)
    if(nrow(d2)< 16-num1){
      num2 <- nrow(d2)
    }else{num2 <- 16-num1}
  }else{
    if(nrow(d2)<8){
      num2 <- nrow(d2)
      if(nrow(d1) < 16-num2){
        num1 <- nrow(d1)
      }
    }
  }
  
  if(num1!=0){dd1 <- head(d1,num1);dd1$type <- 'Positive correlation'}else{dd1 <- d1;dd1$type <- NULL}
  if(num2!=0){dd2 <- head(d2,num2);dd2$type <- 'Negative correlation'}else{dd2 <- d2;dd2$type <- NULL}
  
  tmp2 <- rbind(dd1,dd2)%>%dplyr::select(ID,Description,GeneRatio,type,ONTOLOGY)%>%
    tidyr::separate(GeneRatio,into = c('Enriched','Detected'))%>%
    dplyr::rename(Ontology=ONTOLOGY)%>%
    dplyr::mutate(Enriched = as.numeric(Enriched),
                  Detected = as.numeric(Detected),
                  Ontology = ifelse(Ontology=='BP','Biological Process',
                                    ifelse(Ontology=='CC','Cellular Component','Molecular Function')))%>%
    dplyr::arrange(Ontology) %>%
    dplyr::mutate(ID = factor(ID, rev(unique(ID))))
  
  d1 <- GO_pos@result%>%filter(ONTOLOGY=='MF')
  d2 <- GO_neg@result%>%filter(ONTOLOGY=='MF')
  
  num1 <- num2 <- 8
  if(nrow(d1)<8){
    num1 <- nrow(d1)
    if(nrow(d2)< 16-num1){
      num2 <- nrow(d2)
    }else{num2 <- 16-num1}
  }else{
    if(nrow(d2)<8){
      num2 <- nrow(d2)
      if(nrow(d1) < 16-num2){
        num1 <- nrow(d1)
      }
    }
  }
  if(num1!=0){dd1 <- head(d1,num1);dd1$type <- 'Positive correlation'}else{dd1 <- d1;dd1$type <- NULL}
  if(num2!=0){dd2 <- head(d2,num2);dd2$type <- 'Negative correlation'}else{dd2 <- d2;dd2$type <- NULL}
  
  tmp3 <- rbind(dd1,dd2)%>%dplyr::select(ID,Description,GeneRatio,type,ONTOLOGY)%>%
    tidyr::separate(GeneRatio,into = c('Enriched','Detected'))%>%
    dplyr::rename(Ontology=ONTOLOGY)%>%
    dplyr::mutate(Enriched = as.numeric(Enriched),
                  Detected = as.numeric(Detected),
                  Ontology = ifelse(Ontology=='BP','Biological Process',
                                    ifelse(Ontology=='CC','Cellular Component','Molecular Function')))%>%
    dplyr::arrange(Ontology) %>%
    dplyr::distinct(Description,.keep_all = T)%>%
    dplyr::mutate(ID = factor(ID, rev(unique(ID))))
  tmp <- rbind(tmp1,tmp2,tmp3)%>%dplyr::distinct(ID,.keep_all = T)
  tmp$ID <- factor(tmp$ID,rev(tmp$ID))
  return(tmp)
}

lzq_GO_plot <- function(tmp,width){
  tmp_l1 <- data.frame(nrow(tmp) - (cumsum(table(tmp$Ontology)) - table(tmp$Ontology)/2))
  tmp_l1$start <- tmp_l1$Freq - table(tmp$Ontology)/2
  tmp_l1$end <- tmp_l1$Freq + table(tmp$Ontology)/2
  
  m1 <- ifelse(log(max(tmp$Detected),base = 20)*1000<1000,1000,log(max(tmp$Detected),base = 20)*1000)
  m2 <- ifelse(log(max(tmp$Enriched),base = 20)*1000<1000,1000,log(max(tmp$Enriched),base = 20)*1000)
  
  p1 <- ggplot(tmp) +
    geom_col(mapping = aes(ID, Detected),
             color = "black", fill = "#e6b8a2",
             width = 0.75, 
             show.legend = F) +
    geom_text(mapping = aes(ID, Detected, label = Detected),hjust=-0.3, size = 2.5) +
    scale_y_log10(limits = c(1, m1),expand = c(0,0),breaks=c(1,100),labels=c(1,100)) +
    coord_flip() + 
    theme_classic() +
    labs(x = NULL, y = NULL, title = "Detected Genes") +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          plot.margin = unit(c(1,2,1,0), "mm"),
          axis.text.x = element_text(size=7),
          plot.title = element_text(hjust = 0.2, size = 10,face = 'bold',color='black',vjust=-1))
  
  p2 <- ggplot(tmp) +
    geom_col(mapping = aes(ID, Enriched, fill = type),color = "black", width = 0.75, show.legend = F) +
    geom_text(mapping = aes(ID, Enriched, label = Enriched),hjust=-0.3, size = 2.5) +
    scale_y_log10(limits = c(1, m2),expand = expansion(),breaks=c(1,100,1000),labels=c(1,100,1000)) +
    scale_fill_manual(values = alpha(c('#00b4d8','#ff477e'),0.9)) +
    coord_flip() + 
    theme_classic() +
    labs(x = NULL, y = NULL, title = "Enriched Genes",fill=NULL) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          axis.text.x = element_text(size=7),
          plot.margin = unit(c(1,1,1,1), "mm"),
          plot.title = element_text(hjust = 0.2, size = 10,face = 'bold',vjust=-1))
  
  p0 <- ggplot(tmp) +
    geom_text(mapping = aes(ID, 0, label = ID, color = type),size = 3, show.legend = F, hjust = 0) +
    scale_color_manual(values = alpha(c('#00b4d8','#ff477e'),0.9)) +
    scale_y_continuous(expand = c(0,0),limits = c(0,0.02)) +
    coord_flip() + theme_void() +
    labs(x = NULL, y = NULL, title = "Identifiers") +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          plot.margin = unit(c(1,0,1,3), "mm"),
          plot.title = element_text(hjust = 0.01, size = 10,face = 'bold',vjust=1.5))
  
  p3 <- ggplot(tmp) +
    geom_text(mapping = aes(ID, 0, label = Description, color = type), 
              size = 3, show.legend = F, hjust = 0) +
    scale_color_manual(values = alpha(c('#00b4d8','#ff477e'),0.9)) +
    scale_y_continuous(expand = c(0,0),limits = c(0,0.1)) +
    coord_flip() + 
    theme_void() +
    labs(x = NULL, y = NULL, title = "Description") +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          plot.margin = unit(c(1,0.2,1,1), "mm"),
          plot.title = element_text(hjust = 0, size = 10,face = 'bold',vjust=1.5))
  
  p4 <- ggplot(tmp_l1) +
    geom_segment(mapping = aes(x = start+0.1, xend = end-0.1, y = -0.1, yend = -0.1), size = 2)+
    geom_text(mapping = aes(Freq, 0, label = Var1), size = 3, show.legend = F, hjust = 0) +
    scale_y_continuous(expand = expansion(), limits = c(-0.1,1)) +
    scale_x_continuous(expand = expansion()) +
    coord_flip() + 
    theme_void() +
    labs(x = NULL, y = NULL, title = "Ontoloty") +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          plot.margin = unit(c(1,1,1,0.2), "mm"),
          plot.title = element_text(hjust = 0.1,size = 10,face = 'bold',vjust=1.5))
  
  cowplot::plot_grid(p0,p1,p2,p3,p4,align = "h", nrow = 1, 
                     rel_widths = c(0.1,0.15,0.15,width,0.2))
}

lzq_KEGG_merge <- function(KEGG_pos,KEGG_neg){
  d1 <- KEGG_pos@result
  d2 <- KEGG_neg@result
  
  num1 <- num2 <- 25
  if(nrow(d1)<25){
    num1 <- nrow(d1)
    if(nrow(d2)< 50-num1){
      num2 <- nrow(d2)
    }else{num2 <- 50-num1}
  }else{
    if(nrow(d2)<25){
      num2 <- nrow(d2)
      if(nrow(d1) < 50-num2){
        num1 <- nrow(d1)
      }
    }
  }
  
  if(num1!=0){dd1 <- head(d1,num1);dd1$type <- 'Positive correlation'}else{dd1 <- d1;dd1$type <- NULL}
  if(num2!=0){dd2 <- head(d2,num2);dd2$type <- 'Negative correlation'}else{dd2 <- d2;dd2$type <- NULL}
  
  load('keggdatabase.rda')
  tmp <- rbind(dd1,dd2)%>%dplyr::select(ID,Description,GeneRatio,type)%>%
    dplyr::mutate(ID = gsub('hsa','ko',ID))%>%
    tidyr::separate(GeneRatio,into = c('Enriched','Detected'))%>%
    dplyr::mutate(Enriched = as.numeric(Enriched),
                  Detected = as.numeric(Detected),
                  level3 = kegg$level3[match(ID, kegg$id)],
                  level2 = kegg$level2[match(ID, kegg$id)],
                  level1 = kegg$level1[match(ID, kegg$id)])%>%
    dplyr::arrange(level1, level2, type, level3) %>%
    dplyr::distinct(Description,.keep_all = T)
  tmp <- tmp[order(tmp$level1,tmp$level2),]
  tmp$level2 <- factor(tmp$level2,levels = unique(tmp$level2))
  tmp$level1 <- factor(tmp$level1,levels = unique(tmp$level1))
  tmp$Description <- factor(tmp$Description,levels = rev(tmp$Description))
  tmp$ID <- factor(tmp$ID,levels = rev(tmp$ID))
  return(tmp)
}

lzq_KEGG_plot <- function(tmp,width){
  tmp_l1 <- data.frame(nrow(tmp) - (cumsum(table(tmp$level1)) - table(tmp$level1)/2))
  tmp_l1$start <- tmp_l1$Freq - table(tmp$level1)/2
  tmp_l1$end <- tmp_l1$Freq + table(tmp$level1)/2
  tmp_l2 <- data.frame(nrow(tmp) - (cumsum(table(tmp$level2)) - table(tmp$level2)/2))
  tmp_l2$start <- tmp_l2$Freq - table(tmp$level2)/2
  tmp_l2$end <- tmp_l2$Freq + table(tmp$level2)/2
  
  m1 <- ifelse(log(max(tmp$Detected),base = 20)*1000<1000,1000,log(max(tmp$Detected),base = 20)*1000)
  m2 <- ifelse(log(max(tmp$Enriched),base = 20)*1000<1000,1000,log(max(tmp$Enriched),base = 20)*1000)
  
  p1 <- ggplot(tmp) +
    geom_col(mapping = aes(ID, Detected),
             color = "black", fill = "#e6b8a2",
             width = 0.75, 
             show.legend = F) +
    geom_text(mapping = aes(ID, Detected, label = Detected),hjust=-0.3, size = 2.5) +
    scale_y_log10(limits = c(1, m1),expand = c(0,0),breaks=c(1,100),labels=c(1,100)) +
    coord_flip()+ 
    theme_classic() +
    labs(x = NULL, y = NULL, title = "Detected Genes") +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          plot.margin = unit(c(1,2,1,0), "mm"),
          axis.text.x = element_text(size=7),
          plot.title = element_text(hjust = 0.2, size = 10,face = 'bold',color='black',vjust=-1))
  
  p2 <- ggplot(tmp) +
    geom_col(mapping = aes(ID, Enriched, fill = type),color = "black", width = 0.75, show.legend = F) +
    geom_text(mapping = aes(ID, Enriched, label = Enriched),hjust=-0.3, size = 2.5) +
    scale_y_log10(limits = c(1, m2),expand = expansion()) +
    scale_fill_manual(values = alpha(c('#00b4d8','#ff477e'),0.9)) +
    coord_flip() + 
    theme_classic() +
    labs(x = NULL, y = NULL, title = "Enriched Genes",fill=NULL) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          axis.text.x = element_text(size=7),
          plot.margin = unit(c(1,1,1,1), "mm"),
          plot.title = element_text(hjust = 0.2, size = 10,face = 'bold',vjust=-1))
  
  p0 <- ggplot(tmp) +
    geom_text(mapping = aes(ID, 0, label = ID, color = type),size = 3, show.legend = F, hjust = 0) +
    scale_color_manual(values = alpha(c('#00b4d8','#ff477e'),0.9)) +
    scale_y_continuous(expand = c(0,0),limits = c(0,0.02)) +
    coord_flip() + theme_void() +     
    labs(x = NULL, y = NULL, title = "Terms") +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          plot.margin = unit(c(1,0,1,3), "mm"),
          plot.title = element_text(hjust = 0.01, size = 10,face = 'bold',vjust=1.5))
  
  p3 <- ggplot(tmp) +
    geom_text(mapping = aes(ID, 0, label = Description, color = type), 
              size = 3, show.legend = F, hjust = 0) +
    scale_color_manual(values = alpha(c('#00b4d8','#ff477e'),0.9)) +
    scale_y_continuous(expand = c(0,0),limits = c(0,0.1)) +
    coord_flip() + 
    theme_void() +
    labs(x = NULL, y = NULL, title = "Description (Level 3)") +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          plot.margin = unit(c(1,0.2,1,0.2), "mm"),
          plot.title = element_text(hjust = 0, size = 10,face = 'bold',vjust=1.5))
  
  p4 <- ggplot(tmp_l2) +
    geom_segment(mapping = aes(x = start+0.1, xend = end-0.1, y = -0.1, yend = -0.1), size = 2)+
    geom_text(mapping = aes(Freq, 0, label = Var1), 
              size = 3, show.legend = F, hjust = 0) +
    scale_y_continuous(expand = c(0,0), limits = c(-0.1,1)) +
    scale_x_continuous(expand = expansion()) +
    coord_flip() + 
    theme_void() +
    labs(x = NULL, y = NULL, title = "Category (Level 2)") +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          plot.margin = unit(c(1,0.2,1,0.2), "mm"),
          plot.title = element_text(hjust = 0.2, size = 10,face = 'bold',vjust=1.5))
  
  p5 <- ggplot(tmp_l1) +
    geom_segment(mapping = aes(x = start+0.1, xend = end-0.1, y = -0.1, yend = -0.1), size = 2)+
    geom_text(mapping = aes(Freq, 0, label = Var1), 
              size = 3, show.legend = F, hjust = 0) +
    scale_y_continuous(expand = c(0,0), limits = c(-0.1,1)) +
    scale_x_continuous(expand = expansion()) +
    coord_flip() + 
    theme_void() +
    labs(x = NULL, y = NULL, title = "Category (Level 1)") +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          plot.margin = unit(c(1,0.2,1,0), "mm"),
          plot.title = element_text(hjust = 0.2, size = 10,face = 'bold',vjust=1.5))
  cowplot::plot_grid(p0,p1, p2,p3,p4,p5, align = "h", nrow = 1, 
                     rel_widths = c(0.08,0.15,0.15,width, 0.25, 0.25))
}

lzq_GSEA_merge <- function(GSEA_res){
  GSEA_res@result <- GSEA_res@result[order(GSEA_res@result$NES,decreasing = T),]
  num1 <- num2 <- 15
  if(sum(GSEA_res@result$NES>0)<15){
    num1 <- sum(GSEA_res@result$NES>0)
    if(sum(GSEA_res@result$NES<0)< 30-num1){
      num2 <- sum(GSEA_res@result$NES<0)
    }else{num2 <- 30-num1}
  }else{
    if(sum(GSEA_res@result$NES<0)<15){
      num2 <- sum(GSEA_res@result$NES<0)
      if(sum(GSEA_res@result$NES>0) < 30-num2){
        num1 <- sum(GSEA_res@result$NES>0)
      }
    }
  }
  
  if(num1!=0){dd1 <- head(GSEA_res@result[GSEA_res@result$NES>0,],num1)}else{dd1 <- GSEA_res@result[GSEA_res@result$NES>0,]}
  if(num2!=0){dd2 <- tail(GSEA_res@result[GSEA_res@result$NES<0,],num2)}else{dd2 <- GSEA_res@result[GSEA_res@result$NES<0,]}
  
  tmp <- rbind(dd1,dd2)
  return(list(tmp,num1,num2))
}






