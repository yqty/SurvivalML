gsurv <- function(select_survars){
  if(select_survars=='OS'){
    gsurv <-  as.formula(Surv(OS.time,OS)~.)
  }
  if(select_survars=='RFS'){
    gsurv <-  as.formula(Surv(RFS.time,RFS)~.)
  }
  if(select_survars=='DFS'){
    gsurv <-  as.formula(Surv(DFS.time,DFS)~.)
  }
  if(select_survars=='PFS'){
    gsurv <-  as.formula(Surv(PFS.time,PFS)~.)
  }
  if(select_survars=='DSS'){
    gsurv <-  as.formula(Surv(DSS.time,DSS)~.)
  }
  return(gsurv)
}

lzq_coxplot <- function(coxr){
  coxr$Sur_var <- factor(coxr$Sur_var,levels = c('OS','DFS','RFS','DSS','PFS'))
  coxr <- coxr[order(coxr$HR,decreasing = T),]
  coxr <- coxr[order(coxr$Sur_var,decreasing = F),]
  coxr$ll <- ifelse(coxr$P<0.0001,'****',ifelse(coxr$P<0.001,'***',ifelse(coxr$P<0.01,'**',ifelse(coxr$P<0.05,'*',''))))
  coxr$y <- factor(paste0(coxr$Sur_var,'_',coxr$Cohort),levels = rev(paste0(coxr$Sur_var,'_',coxr$Cohort)))
  cols2 <- rev(c('#B0997F','#FF6666','#ffc857','#93a8ac','#119da4'))
  print(ggplot(coxr,aes(HR,y))+
          geom_rect(xmin=log10(min(coxr$HRL*0.7)), xmax=log10(max(coxr$HRH*1.3)), ymin=cumsum(rev(table(coxr$Sur_var)))[4]+0.5, ymax=cumsum(rev(table(coxr$Sur_var)))[5]+0.5,fill=alpha(cols2[1],0.003))+
          geom_rect(xmin=log10(min(coxr$HRL*0.7)), xmax=log10(max(coxr$HRH*1.3)), ymin=cumsum(rev(table(coxr$Sur_var)))[3]+0.5, ymax=cumsum(rev(table(coxr$Sur_var)))[4]+0.5,fill=alpha(cols2[2],0.006))+
          geom_rect(xmin=log10(min(coxr$HRL*0.7)), xmax=log10(max(coxr$HRH*1.3)), ymin=cumsum(rev(table(coxr$Sur_var)))[2]+0.5, ymax=cumsum(rev(table(coxr$Sur_var)))[3]+0.5,fill=alpha(cols2[3],0.003))+
          geom_rect(xmin=log10(min(coxr$HRL*0.7)), xmax=log10(max(coxr$HRH*1.3)), ymin=cumsum(rev(table(coxr$Sur_var)))[1]+0.5, ymax=cumsum(rev(table(coxr$Sur_var)))[2]+0.5,fill=alpha(cols2[4],0.002))+
          geom_rect(xmin=log10(min(coxr$HRL*0.7)), xmax=log10(max(coxr$HRH*1.3)), ymin=0.5, ymax=cumsum(rev(table(coxr$Sur_var)))[1]+0.5,fill=alpha('#583101',0.006))+
          geom_vline(xintercept = 1,linetype=2,color='grey50')+
          geom_errorbar(aes(xmin=HRL,xmax=HRH),width=0.1,size=0.8)+
          geom_point(shape=15,size=3,aes(color=Sur_var))+
          scale_color_manual(values = cols2)+
          geom_text(aes(label=ll),vjust=1.6,size=4)+
          scale_x_log10()+
          labs(x='Hazard ratio',title = 'Cox regression anlaysis')+
          theme_bw(base_rect_size = 0)+
          geom_rect(xmin=log10(min(coxr$HRL*0.7)), xmax=log10(max(coxr$HRH*1.3)), ymin=cumsum(rev(table(coxr$Sur_var)))[4]+0.5, ymax=cumsum(rev(table(coxr$Sur_var)))[5]+0.5,fill=alpha('#B0997F',0.01))+
          theme(axis.text.y = element_text(size = 11,colour = 'black'),
                axis.text.x = element_text(size = 8,colour = 'black'),
                axis.title.x = element_text(size = 12,colour = 'darkred',face='bold'),
                axis.title.y = element_blank(),
                strip.background = element_blank(),
                strip.text = element_text(size = 13,colour = 'darkred',face='bold'),
                panel.grid = element_blank(),
                panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed",size = 0.2),
                panel.background = element_rect(fill='#f3f6f6'),
                plot.title = element_text(hjust = 0.5,size = 12,colour = 'darkred',face='bold'),
                legend.position = 'bottom',
                legend.title = element_blank(),
                legend.text = element_text(size = 12,colour = 'black'))+
          guides(color=guide_legend(nrow = 1)))
}
lzq_survplot <- function(Cohort = NULL,
                         sur_type = 'Overall survival',
                         data,
                         time = 'time',
                         event = 'status',
                         var,
                         color = c("#1F78B4", "#E31A1C", "#FDBF6F"),
                         cutoff = NULL,
                         percent = NULL,
                         main = 'CD274',
                         label = c('High', 'Low')) {
  tmp <- data[, c(time, event, var)]
  colnames(tmp) <- c('time', 'event', 'var')
  spf <- function(label = label) {
    customize_labels <- function (p,
                                  font.title = NULL,
                                  font.subtitle = NULL,
                                  font.caption = NULL,
                                  font.x = NULL,
                                  font.y = NULL,
                                  font.xtickslab = NULL,
                                  font.ytickslab = NULL)
    {
      original.p <- p
      if (is.ggplot(original.p))
        list.plots <- list(original.p)
      else if (is.list(original.p))
        list.plots <- original.p
      else
        stop("Can't handle an object of class ", class (original.p))
      .set_font <- function(font) {
        font <- ggpubr:::.parse_font(font)
        ggtext::element_markdown (
          size = font$size,
          face = font$face,
          colour = font$color
        )
      }
      for (i in 1:length(list.plots)) {
        p <- list.plots[[i]]
        if (is.ggplot(p)) {
          if (!is.null(font.title))
            p <- p + theme(plot.title = .set_font(font.title))
          if (!is.null(font.subtitle))
            p <- p + theme(plot.subtitle = .set_font(font.subtitle))
          if (!is.null(font.caption))
            p <- p + theme(plot.caption = .set_font(font.caption))
          if (!is.null(font.x))
            p <- p + theme(axis.title.x = .set_font(font.x))
          if (!is.null(font.y))
            p <- p + theme(axis.title.y = .set_font(font.y))
          if (!is.null(font.xtickslab))
            p <- p + theme(axis.text.x = .set_font(font.xtickslab))
          if (!is.null(font.ytickslab))
            p <- p + theme(axis.text.y = .set_font(font.ytickslab))
          list.plots[[i]] <- p
        }
      }
      if (is.ggplot(original.p))
        list.plots[[1]]
      else
        list.plots
    }
    pp <- ggsurvplot(
      fit,
      tmp,
      pval = TRUE,
      pval.method = T,
      ylab = NULL,
      xlab = 'Time in years',
      size = 1.3,
      conf.int = F,
      legend.title = main,
      legend.labs = label,
      legend = 'none',
      risk.table = TRUE,
      risk.table.pos = 'out',
      tables.col = "strata",
      risk.table.title = "Number at risk",
      risk.table.height = .3,
      risk.table.y.text.col = T,
      risk.table.y.text = T,
      risk.table.y.title = F,
      palette = color,
      font.main = 15,
      ggtheme = theme_bw(base_rect_size = 2)
    )
    pp$plot <- customize_labels(
      pp$plot,
      font.x        = c(14, "bold", "darkred"),
      font.y        = c(14, "bold", "darkred"),
      font.xtickslab = c(12, "plain", "black"),
      font.ytickslab = c(12, "plain", 'black')
    )
    pp$plot <- pp$plot + labs(y = sur_type, x = NULL, title = Cohort) +
      theme(
        plot.title = element_text(
          face = "bold",
          colour = "darkred",
          size = 18,
          hjust = 0.5
        ),
        panel.background = element_rect(fill = "#f3f6f6", color = NA),
        panel.grid.minor  = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "#cacfd2", linetype = "dashed"),
        axis.line = element_line(color = "#606F7B"),
        legend.position = 'none'
      )
    pp$table <- customize_labels(
      pp$table,
      font.title  = c(14, "bold", "darkgreen"),
      font.x        = c(15, "bold", "darkred"),
      font.y        = c(14, "bold", "darkred"),
      font.xtickslab = c(12, "plain", "black"),
      font.ytickslab = c(12, "bold")
    ) +
      theme(
        panel.background = element_rect(fill = "#f3f6f6", color = NA),
        panel.grid.minor  = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "#cacfd2", linetype = "dashed"),
        axis.line = element_line(color = "#606F7B"),
        legend.position = 'none'
      )
    
    # merge table and plot with patchwork
    library(patchwork)
    pp.out = pp$plot + pp$table + patchwork::plot_layout(nrow = 2,heights = c(3,1))
    return(pp.out)
  }
  
  if (is.null(cutoff)) {
    stop(
      "Please input cutoff method! The method include 'median','mean','quantile','optimal','custom', just pick one"
    )
  }
  if (is.character(cutoff)) {
    if (cutoff == 'median') {
      tmp$Group <- ifelse(tmp$var > median(tmp$var), 'High', 'Low')
      fit <- survfit(Surv(time, event) ~ Group, tmp)
      return(spf(label))
    }
    if (cutoff == 'optimal') {
      cut <- surv_cutpoint(tmp, 'time', 'event', 'var', minprop = 0.15)
      tmp$Group <-
        ifelse(tmp$var > cut$cutpoint$cutpoint, 'High', 'Low')
      fit <- survfit(Surv(time, event) ~ Group, tmp)
      return(spf(label))
    }
    if (cutoff == 'mean') {
      tmp$Group <- ifelse(tmp$var > mean(tmp$var), 'High', 'Low')
      fit <- survfit(Surv(time, event) ~ Group, tmp)
      return(spf(label))
    }
    if (cutoff == 'quantile') {
      v1 <- as.numeric(quantile(tmp$var)[2])
      v2 <- as.numeric(quantile(tmp$var)[4])
      tmp <- subset(tmp, var > v2 | var < v1)
      tmp$Group <- ifelse(tmp$var > v1, 'High', 'Low')
      fit <- survfit(Surv(time, event) ~ Group, tmp)
      return(spf(label))
    }
    if (cutoff == 'custom') {
      if (is.null(percent)) {
        stop('You must specify the percent value for the smaller group!')
      }
      if (is.numeric(percent)) {
        x <- as.numeric(quantile(tmp$var, percent))
        tmp$Group <- ifelse(tmp$var > x, 'Low', 'High')
        fit <- survfit(Surv(time, event) ~ Group, tmp)
        return(spf(label))
      }
    }
  }
  
}

lzq_survplot2 <-
  function(Sur_ids,
           cohort,
           data,
           gene,
           cols,
           cutoff,
           Input,
           percent) {
    if ('OS' %in% Sur_ids) {
      if (cutoff != 'custom') {
        p.out = lzq_survplot(
          Cohort = cohort,
          sur_type = 'Overall survival',
          data = data,
          time = 'OS.time',
          event = 'OS',
          var = gene,
          color = cols,
          cutoff = cutoff,
          main = Input
        )
        return(p.out)
      } else{
        p.out = lzq_survplot(
          Cohort = cohort,
          sur_type = 'Overall survival',
          data = data,
          time = 'OS.time',
          event = 'OS',
          var = gene,
          color = cols,
          cutoff = 'custom',
          percent = percent,
          main = Input
        )
        return(p.out)
      }
    }
    
    if ('RFS' %in% Sur_ids) {
      if (cutoff != 'custom') {
        p.out = lzq_survplot(
          Cohort = cohort,
          sur_type = 'Relapse-free survival',
          data = data,
          time = 'RFS.time',
          event = 'RFS',
          var = gene,
          color = cols,
          cutoff = cutoff,
          main = Input
        )
        return(p.out)
      } else{
        p.out = lzq_survplot(
          Cohort = cohort,
          sur_type = 'Relapse-free survival',
          data = data,
          time = 'RFS.time',
          event = 'RFS',
          var = gene,
          color = cols,
          cutoff = 'custom',
          percent = percent,
          main = Input
        )
        return(p.out)
      }
    }
    
    if ('DFS' %in% Sur_ids) {
      if (cutoff != 'custom') {
        p.out = lzq_survplot(
          Cohort = cohort,
          sur_type = 'Disease-free survival',
          data = data,
          time = 'DFS.time',
          event = 'DFS',
          var = gene,
          color = cols,
          cutoff = cutoff,
          main = Input
        )
        return(p.out)
      } else{
        p.out = lzq_survplot(
          Cohort = cohort,
          sur_type = 'Disease-free survival',
          data = data,
          time = 'DFS.time',
          event = 'DFS',
          var = gene,
          color = cols,
          cutoff = 'custom',
          percent = percent,
          main = Input
        )
        return(p.out)
      }
    }
    
    if ('PFS' %in% Sur_ids) {
      if (cutoff != 'custom') {
        p.out = lzq_survplot(
          Cohort = cohort,
          sur_type = 'Progress-free survival',
          data = data,
          time = 'PFS.time',
          event = 'PFS',
          var = gene,
          color = cols,
          cutoff = cutoff,
          main = Input
        )
        return(p.out)
      } else{
        p.out = lzq_survplot(
          Cohort = cohort,
          sur_type = 'Progress-free survival',
          data = data,
          time = 'PFS.time',
          event = 'PFS',
          var = gene,
          color = cols,
          cutoff = 'custom',
          percent = percent,
          main = Input
        )
        return(p.out)
      }
    }
    
    if ('DSS' %in% Sur_ids) {
      if (cutoff != 'custom') {
        p.out = lzq_survplot(
          Cohort = cohort,
          sur_type = 'Disease-specific survival',
          data = data,
          time = 'DSS.time',
          event = 'DSS',
          var = gene,
          color = cols,
          cutoff = cutoff,
          main = Input
        )
        return(p.out)
      } else{
        p.out = lzq_survplot(
          Cohort = cohort,
          sur_type = 'Disease-specific survival',
          data = data,
          time = 'DSS.time',
          event = 'DSS',
          var = gene,
          color = cols,
          cutoff = 'custom',
          percent = percent,
          main = Input
        )
        return(p.out)
      }
    }
  }