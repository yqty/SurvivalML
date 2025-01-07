lzq_gseaplot <- function(GSEA_result,Pathway_ID,color){
  gsdata <- enrichplot:::gsInfo(GSEA_result,Pathway_ID)
  label <- paste0('NES = ',sprintf("%.3f", GSEA_result@result$NES[GSEA_result@result$ID==Pathway_ID]),'; FDR = ',
                  format(GSEA_result@result$p.adjust[GSEA_result@result$ID==Pathway_ID], scientific = T))
  p1 <- ggplot(gsdata, aes(x)) + 
    geom_line(aes(y = runningScore), size=1, color=color) +
    scale_x_continuous(expand=c(0,0)) + 
    scale_y_continuous(expand=c(0,0.01)) +
    labs(x=NULL,y="Enrichment Score",title=Pathway_ID,subtitle = label) +
    theme_bw(base_rect_size = 2) +
    theme(plot.title = element_text(hjust = 0.5,size = 14,colour = 'darkred',face='bold'),
          plot.subtitle = element_text(face = 'italic',hjust = 0.5,size = 11,colour = 'black'),
          panel.grid = element_blank(),
          axis.text.x=element_blank(),
          axis.text.y = element_text(size = 10,colour = 'black'),
          axis.title.y = element_text(size = 13,colour = 'darkred',face='bold'),
          axis.ticks.x=element_blank(),
          axis.line.x=element_blank(),
          plot.margin=margin(t = 0.2, r = .2, b=-0.07, l=.2, unit="cm"))
  p2 <- ggplot(gsdata, aes(x)) + 
    geom_linerange(aes(ymin = ymin,ymax = ymax),color='grey30') + 
    xlab(NULL) + ylab(NULL) + 
    theme_bw(base_rect_size = 2) +
    theme(legend.position = "none", 
          panel.grid = element_blank(),
          plot.margin = margin(t = -0.1, b = 0, unit = "cm"), 
          axis.ticks = element_blank(), axis.text = element_blank(), 
          axis.line.x = element_blank()) + 
    scale_x_continuous(expand=c(0,0)) + 
    scale_y_continuous(expand=c(0,0))
  p2
  v <- seq(1, sum(gsdata$position), length.out = 9)
  inv <- findInterval(rev(cumsum(gsdata$position)), v)
  if (min(inv) == 0) {inv <- inv + 1}
  col <- c(rev(brewer.pal(5, "Blues")), brewer.pal(5,"Reds"))
  ymin <- min(p2$data$ymin)
  yy <- max(p2$data$ymax - p2$data$ymin) * 0.3
  xmin <- which(!duplicated(inv))
  xmax <- xmin + as.numeric(table(inv)[as.character(unique(inv))])
  d <- data.frame(ymin = ymin, ymax = yy, xmin = xmin, 
                  xmax = xmax, col = col[unique(inv)])
  p2 <- p2 + geom_rect(aes(xmin = xmin, xmax = xmax,ymin = ymin, ymax = 0, fill = I(col)),
                       data = d,alpha = 0.9, inherit.aes = FALSE)
  p2
  
  p3 <- ggplot(gsdata, aes(x)) + 
    labs(y='Correlation',x='Rank in Ordered Dataset')+
    scale_x_continuous(expand=c(0,0)) + 
    scale_y_continuous(expand=c(0,0)) +
    geom_segment(aes(x = x, xend = x,y = geneList, yend = 0, color = geneList))+
    scale_color_continuous(type = 'viridis')+
    theme_bw(base_rect_size = 2) +
    theme(plot.title = element_text(hjust = 0.5,size = 14,colour = 'darkred',face='bold'),
          plot.subtitle = element_text(face = 'italic',hjust = 0.5,size = 11,colour = 'black'),
          panel.grid = element_blank(),
          legend.position = 'none',
          axis.text.x=element_text(size = 10,colour = 'black'),
          axis.text.y = element_text(size = 10,colour = 'black'),
          axis.title = element_text(size = 13,colour = 'darkred',face='bold'),
          plot.margin=margin(t = -.17, r = .2, b=.2, l=.2, unit="cm"))
  p3
  
  require(cowplot)
  plotlist <- list(p1, p2, p3)
  n <- length(plotlist)
  plotlist[[n]] <- plotlist[[n]] +
    theme(axis.line.x = element_line(),
          axis.ticks.x = element_line(),
          axis.text.x = element_text())
  plot_grid(plotlist = plotlist, ncol = 1, align="v", rel_heights = c(1.5, .2, 1))
}
