rm(list = ls())
library(tidyverse)
library(ggsci)
library(survival)
library(survminer)
library(patchwork)
library(timeROC)
library(rms)

source('lzq_survplot.R')
load('score_list.rda')

# 单因素与Cindex --------------------------------------------------------------

if(select_survars=='OS'){
  gsurv <-  as.formula(Surv(OS.time,OS)~score)
}
if(select_survars=='RFS'){
  gsurv <-  as.formula(Surv(RFS.time,RFS)~score)
}
if(select_survars=='DFS'){
  gsurv <-  as.formula(Surv(DFS.time,DFS)~score)
}
if(select_survars=='PFS'){
  gsurv <-  as.formula(Surv(PFS.time,PFS)~score)
}
if(select_survars=='DSS'){
  gsurv <-  as.formula(Surv(DSS.time,DSS)~score)
}

coxr <- lapply(score_list,function(x){
  x[,2] <- as.numeric(x[,2])
  x[,3] <- as.numeric(x[,3])
  fit <- summary(coxph(gsurv,x))
  Cindex <- fit$concordance[1]
  Cse <- fit$concordance[2]
  HR <- fit$coefficients[,2]
  HRL <- fit$conf.int[,3]
  HRR <- fit$conf.int[,4]
  P <- fit$coefficients[,5]
  return(data.frame(HR=HR,HRL=HRL,HRR=HRR,P=P,Cindex=Cindex,Cse=Cse))
})
coxr <- Reduce(rbind,coxr)
coxr$ID <- names(score_list)
coxr$ll <- ifelse(coxr$P<0.0001,'****',ifelse(coxr$P<0.001,'***',ifelse(coxr$P<0.01,'**',ifelse(coxr$P<0.05,'*',''))))

p1 <- ggplot(coxr,aes(HR,ID))+
  geom_vline(xintercept = 1,linetype=2,color='grey50')+
  geom_errorbar(aes(xmin=HRL,xmax=HRR),width=0.1,size=0.8)+
  geom_point(shape=15,size=3,color='#ffb627')+
  geom_text(aes(label=ll),vjust=1.8,size=4)+
  scale_x_log10()+
  labs(x='Hazard ratio',title = 'Cox regression')+
  theme_bw(base_rect_size = 0)+
  theme(axis.text.y = element_text(size = 12,colour = 'darkred',face='bold'),
        axis.text.x = element_text(size = 10,colour = 'black'),
        axis.title.x = element_text(size = 13,colour = 'darkred',face='bold'),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.background = element_rect(fill='#f3f6f6'),
        plot.title = element_text(hjust = 0.5,size = 14,colour = 'darkred',face='bold'),
        legend.position = 'none')

p2 <- ggplot(coxr,aes(Cindex,ID))+
  geom_errorbar(aes(xmin=Cindex-Cse,xmax=Cindex+Cse),width=0.1,size=0.8)+
  geom_bar(stat = 'identity',width = 0.6,color='black',fill='#2ec4b6',size=0.5)+
  geom_vline(xintercept = 0.5,linetype=2,color='grey50')+
  labs(x='Index of concordance',title = 'C-index')+
  theme_bw(base_rect_size = 0)+
  scale_x_continuous(expand = c(0,0))+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 10,colour = 'black'),
        axis.title.x = element_text(size = 13,colour = 'darkred',face='bold'),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.background = element_rect(fill='#f3f6f6'),
        plot.title = element_text(hjust = 0.5,size = 14,colour = 'darkred',face='bold'),
        legend.position = 'none')
p1+p2

#take "TCGA_LUAD" as an example
all_cohort_name

Cohort <- 'GSE37745' # GSE50081 GSE72094 GSE42127

## KM survival curve ##
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

cutoff <- 'median' ## Please input cutoff method! The method include 'median','mean','quantile','optimal','custom', just pick one
percent <- 0.3 ## If users select the "custom" cutoff method, the percent option will appearand are set to a range of 0.1 to 0.9

p1 <- lzq_survplot2(Sur_ids = select_survars,
              cohort = Cohort,
              data = score_list[[Cohort]],
              gene = 'score',
              Input = self_name,
              cols = cols,cutoff = cutoff,percent = percent)
p1
ggsave(plot = p1,filename = paste0(Cohort,'-KM.pdf'),width = 4.5,height = 5)


## ROC curve ##

tmp <- score_list[[Cohort]]

if(max(tmp[,1],na.rm = T)<1){
  time_cut <- c('3-Month AUC','6-Month AUC','9-Month AUC')
}
if(max(tmp[,1],na.rm = T)<2&max(tmp[,1],na.rm = T)>1){
  time_cut <- c('4-Month AUC','8-Month AUC','12-Month AUC')
}
if(max(tmp[,1],na.rm = T)<3&max(tmp[,1],na.rm = T)>2){
  time_cut <- c('6-Month AUC','12-Month AUC','24-Month AUC')
}
if(max(tmp[,1],na.rm = T)<5&max(tmp[,1],na.rm = T)>3){
  time_cut <- c('1-Year AUC','2-Year AUC','3-Year AUC')
}
if(max(tmp[,1],na.rm = T)>5){
  time_cut <- c('1-Year AUC','3-Year AUC','5-Year AUC')
}
get_times <- function(time_cut){
  if(grepl('Year',time_cut[1])){
    times <- sapply(time_cut,function(x){unlist(strsplit(x,'-'))[1]})%>%as.numeric()
  }
  if(grepl('Month',time_cut[1])){
    times <- (sapply(time_cut,function(x){unlist(strsplit(x,'-'))[1]})%>%as.numeric())/12
  }
  return(times)
}
times <- get_times(time_cut)

tt <- timeROC(tmp[,1],tmp[,2],tmp[,3],cause = 1,times = times,ROC = T,weighting = 'marginal')
tp <- tt$TP%>%as.data.frame()%>%pivot_longer(cols = 1:3,names_to = 'time',values_to = 'tp')
fp <- tt$FP%>%as.data.frame()%>%pivot_longer(cols = 1:3,names_to = 'time',values_to = 'fp')

dd <- tp
dd$fp <- fp$fp
dd$time <- ifelse(dd$time==unique(dd$time)[1],paste0(time_cut[1],' = ',sprintf("%.3f",tt$AUC[1])),
                  ifelse(dd$time==unique(dd$time)[2],paste0(time_cut[2],' = ',sprintf("%.3f",tt$AUC[2])),
                         paste0(time_cut[3],' = ',sprintf("%.3f",tt$AUC[3]))))

p2 <- ggplot(dd,aes(fp,tp,color=time))+
  geom_line(size=1)+
  labs(x='1-Specificity',y='Sensitivity',color=NULL)+
  theme_bw(base_rect_size = 1.5)+
  geom_abline(slope = 1,color='grey70')+
  ggtitle(Cohort)+
  theme(axis.text.y = element_text(size = 12,colour = 'black'),
        axis.text.x = element_text(size = 12,colour = 'black'),
        axis.title = element_text(size = 15,colour = 'darkred',face='bold'),
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.background = element_rect(fill='#f3f6f6'),
        plot.title = element_text(hjust = 0.5,size = 18,colour = 'darkred',face='bold'),
        legend.text = element_text(size=14),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.position = c(0.99,0.013),
        legend.justification = c(1,0))+
  scale_color_manual(values = c('#2a9d8f','#fcbf49','#eb5e28'))+
  scale_x_continuous(expand = c(0.01,0.01))+
  scale_y_continuous(expand = c(0.01,0.01))
p2

## Calibration curve ##

cox1 <- cph(gsurv,surv=T,x=T, y=T,data=tmp,time.inc = times[1]) 
cal1 <- tryCatch(calibrate(cox1, cmethod="KM", method="boot", u=times[1], m= ceiling(nrow(tmp)*0.95/5), B=1000),error=function(e){NA})
cox2 <- cph(gsurv,surv=T,x=T, y=T,data=tmp,time.inc = times[2]) 
cal2 <- tryCatch(calibrate(cox2, cmethod="KM", method="boot", u=times[2], m= ceiling(nrow(tmp)*0.95/5), B=1000),error=function(e){NA})
cox3 <- cph(gsurv,surv=T,x=T, y=T,data=tmp,time.inc = times[3]) 
cal3 <- tryCatch(calibrate(cox3, cmethod="KM", method="boot", u=times[3], m= ceiling(nrow(tmp)*0.95/5), B=1000),error=function(e){NA})

dd <- Reduce(rbind,list(tryCatch(cal1[,c('mean.predicted',"KM")],error=function(e)data.frame(mean.predicted=rep(NA,5),KM=rep(NA,5))),
                        tryCatch(cal2[,c('mean.predicted',"KM")],error=function(e)data.frame(mean.predicted=rep(NA,5),KM=rep(NA,5))),
                        tryCatch(cal3[,c('mean.predicted',"KM")],error=function(e)data.frame(mean.predicted=rep(NA,5),KM=rep(NA,5)))))%>%
  as.data.frame()
dd$time <- paste0(rep(gsub(' AUC','',time_cut),each=5),' Calibration')
dd <- na.omit(dd)
colnames(dd)[1:2] <- c('Predicted','Observed')

p3 <- ggplot(dd,aes(Predicted,Observed))+
  geom_abline(slope = 1,color='grey70')+
  geom_line(size=1,alpha=0.7,aes(color=time))+
  theme_bw(base_rect_size = 1.5)+
  geom_point(size=3,shape=21,aes(fill=time),stroke=1)+
  labs(x='Predicted survival (%)',y='Observed survival (%)',title = Cohort)+
  theme(axis.text.y = element_text(size = 12,colour = 'black'),
        axis.text.x = element_text(size = 12,colour = 'black'),
        axis.title = element_text(size = 15,colour = 'darkred',face='bold'),
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.background = element_rect(fill='#f3f6f6'),
        plot.title = element_text(hjust = 0.5,size = 18,colour = 'darkred',face='bold'),
        legend.text = element_text(size=14),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.99,0.013),
        legend.justification = c(1,0))+
  scale_color_manual(values = c('#2a9d8f','#fcbf49','#eb5e28'))+
  scale_fill_manual(values = c('#2a9d8f','#fcbf49','#eb5e28'))+
  scale_x_continuous(expand = c(0.03,0))+
  scale_y_continuous(expand = c(0,0.03))
p3

## DCA
tmp$score <- (tmp$score-min(tmp$score))/(max(tmp$score)-min(tmp$score))
library(dcurves)
##从times变量中选择 这就是 Timepoint for DCA的参数
select_time <- times[2]

dd <- tryCatch((dcurves::dca(gsurv, 
                             data = tmp, 
                             label = list(all='All',none='None',score=self_name),
                             time = select_time)),error=function(e){message('Calculating failed due to some error, please replace another timepoint!')}) 
my <- dd$dca%>%na.omit()
xmax <- max(my$threshold[my$variable=='score'])
dd <- tryCatch((dcurves::dca(gsurv, 
                             data = tmp, thresholds = seq(0,xmax,by=0.01),
                             label = list(all='All',none='None',score=self_name),
                             time = select_time)),error=function(e){message('Calculating failed due to some error, please replace another timepoint!')}) 

p4 <- dd%>%plot(smooth=T)+
  theme_bw(base_rect_size = 1.5)+
  ggtitle(Cohort)+
  theme(axis.text.y = element_text(size = 12,colour = 'black'),
        axis.text.x = element_text(size = 12,colour = 'black'),
        axis.title = element_text(size = 15,colour = 'darkred',face='bold'),
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.background = element_rect(fill='#f3f6f6'),
        plot.title = element_text(hjust = 0.5,size = 18,colour = 'darkred',face='bold'),
        legend.text = element_text(size=14),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.95,0.95),
        legend.justification = c(1,1))+
  scale_color_manual(values = c('#2a9d8f','#fcbf49','#eb5e28'))

p2+p3+p4
ggsave(filename = paste0(Cohort,'-ROC+Calibration+DCA.pdf'),width = 4.5*3,height = 4.6)


