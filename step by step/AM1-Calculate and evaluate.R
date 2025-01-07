rm(list = ls())
library(survival)
library(survminer)
library(timeROC)
library(rms)
library(tidyverse)
library(ggsci)
library(ComplexHeatmap)
library(circlize)
library(patchwork)
source('lzq_consensus_model.R')
source('lzq_survplot.R')
source('lzq_surmodel.R')

# -------------------------------------------------------------------------

select_cancer <- 'GBM'
load(paste0('data/',select_cancer,'/symbol.rda'))
names <- names(total_clin_list)

# -------------------------------------------------------------------------
## select event
Survar_Input <- 'Overall survival' ## 用户自己选择
select_survars <- ifelse(Survar_Input=='Overall survival','OS',
                         ifelse(Survar_Input=='Relapse-free survival','RFS',
                                ifelse(Survar_Input=='Disease-free survival','DFS',
                                       ifelse(Survar_Input=='Progression-free survival','PFS','DSS'))))
available_cohorts <- names[sapply(total_clin_list,function(x){select_survars%in%colnames(x)})]

# -------------------------------------------------------------------------
## select cohorts from available cohorts
##这里我选择了全部
select_cohorts <- c('GSE74187','GSE83300','CGGA_301','CGGA_325','CGGA_693')

# -------------------------------------------------------------------------

expdata <- total_expr_list[select_cohorts]

overlap_genes <- Reduce(intersect,lapply(expdata,rownames))
message(paste0('A total of ',length(overlap_genes),' genes overlap in your selected cohorts'))

expdata <- lapply(expdata,function(x){x <- x[overlap_genes,]%>%na.omit();return(x)})

alldata <- list()
for (i in names(expdata)) {
  alldata[[i]] <- merge(total_clin_list[[i]][,c('ID',select_survars,paste0(select_survars,'.time'))],t(expdata[[i]]),by.x=1,by.y=0)%>%
    column_to_rownames('ID')%>%na.omit()
}
alldata <- lapply(alldata,function(x){x[x[,2]>0,]})
sample_size <- sapply(alldata,nrow)
message(paste0('A total of ',sum(sample_size),' eligible samples present in your selected cohorts'))

# -------------------------------------------------------------------------
## select consensus pvalue
##这里我选择了0.05
consensus_pval <- 0.05

# -------------------------------------------------------------------------
## select training cohort, the other cohorts were naturally testing cohorts
##这里我选择了GSE103479
train_cohort_name <- "CGGA_693"

# -------------------------------------------------------------------------
## 通过设置的consensus_pval来找出在所有队列中稳定的预后相关基因

load(paste0('data/',select_cancer,'/cm_surlist.rda'))
unires2 <- surlist[[select_survars]][,select_cohorts]%>%na.omit()

PS <- apply(unires2,1,function(x){sum(abs(x)<consensus_pval&x<0)})
RS <- apply(unires2,1,function(x){sum(abs(x)<consensus_pval&x>0)})

unires2$PS <- PS
unires2$RS <- RS

# -------------------------------------------------------------------------
## select num_sig
## Whether to filter out genes with inconsistent prognosis
num_sig <- 5 ## mean a gene must be signifincant (p < consensus_pval) in more than eight cohorts
filter <- TRUE

if(filter){
  d1 <- unires2[unires2$PS!=0&unires2$RS==0,]
  d2 <- unires2[unires2$RS!=0&unires2$PS==0,]
  unires2 <- rbind(d1,d2)
}

select_ids <- c(rownames(unires2)[unires2$PS>=num_sig],rownames(unires2)[unires2$RS>=num_sig])

tmp <- unires2[select_ids,!colnames(unires2)%in%c('PS','RS')]
for (i in 1:ncol(tmp)) {
  tmp[,i] <- ifelse(tmp[,i]>0&abs(tmp[,i])<consensus_pval,'Risky',ifelse(tmp[,i]<0&abs(tmp[,i])<consensus_pval,'Protective','Not Significant'))
}

if(length(unique(c(t(tmp))))==3){
  cols <- c('#F2F2F2',alpha(c('#00b4d8','#ff477e'),0.9))
}
if(length(unique(c(t(tmp))))==2){
  if(sum(unique(c(t(tmp)))%in%c('Protective','Risky'))==2){
    cols <- alpha(c('#00b4d8','#ff477e'),0.9)
  }
  if(sum(unique(c(t(tmp)))%in%c('Not Significant','Risky'))==2){
    cols <- c('#F2F2F2',alpha('#ff477e',0.9))
  }
  if(sum(unique(c(t(tmp)))%in%c('Not Significant','Protective'))==2){
    cols <- c('#F2F2F2',alpha('#00b4d8',0.9))
  }
}
if(length(unique(c(t(tmp))))==1){
  if(unique(c(t(tmp)))=='Protective'){
    cols <- alpha('#00b4d8',0.9)
  }
  if(unique(c(t(tmp)))=='Risky'){
    cols <- alpha('#ff477e',0.9)
  }
}

##输出
draw(Heatmap(t(tmp),col = cols,
        name = 'Type',border = T,
        rect_gp = gpar(col='black',lwd=2),
        row_names_side = 'left',
        column_names_side = 'bottom',
        row_split = 1:ncol(tmp),
        row_title = NULL,
        heatmap_legend_param=list(labels_gp = gpar(fontsize = 12), border = T,nrow=1,title_position = 'leftcenter',
                                  title_gp = gpar(fontsize = 12, fontface = "bold"))), heatmap_legend_side = "top")

## 建立10个模型
alldata <- lapply(alldata,function(x){return(x[,c(paste0(select_survars,'.time'),select_survars,select_ids)])})
traindata <- alldata[[train_cohort_name]]

result <- lzq_consensus_model(select_survars = select_survars,
                              alldata = alldata,
                              RSF_nsplit = 5,RSF_nodesize = 10, # for gbm
                              traindata = traindata)

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

coxrl <- lapply(result,function(model){
  dd <- lapply(model,function(x){
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
  dd <- Reduce(rbind,dd)
  dd$ID <- names(model)
  dd$ll <- ifelse(dd$P<0.0001,'****',ifelse(dd$P<0.001,'***',ifelse(dd$P<0.01,'**',ifelse(dd$P<0.05,'*',''))))
  return(dd)
})

for (i in names(coxrl)) {coxrl[[i]]$model <- i}

tmp <- Reduce(rbind,coxrl)
tmp$ID <- factor(tmp$ID,levels = c(sort(select_cohorts[select_cohorts!=train_cohort_name],decreasing = T),train_cohort_name))

cols <- c(pal_npg(alpha = 1)(10),pal_jco(alpha = 1)(10),pal_d3(alpha = 1)(10),pal_nejm(alpha = 1)(8))
options(digits = 2)

p1 <- ggplot(tmp,aes(HR,ID))+
  geom_vline(xintercept = 1,linetype=2,color='grey50')+
  geom_errorbar(aes(xmin=HRL,xmax=HRR),width=0.1,size=0.8)+
  geom_point(shape=15,size=3,aes(color=model))+
  scale_color_manual(values = cols)+
  geom_text(aes(label=ll),vjust=1.8,size=4)+
  scale_x_log10()+
  labs(x='Hazard ratio')+
  theme_bw(base_rect_size = 0)+
  facet_wrap(~model,nrow = 1,scales = 'free_x')+
  theme(axis.text.y = element_text(size = 12,colour = 'darkred',face='bold'),
        axis.text.x = element_text(size = 8,colour = 'black'),
        axis.title.x = element_text(size = 13,colour = 'darkred',face='bold'),
        axis.title.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 13,colour = 'darkred',face='bold'),
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.background = element_rect(fill='#f3f6f6'),
        plot.title = element_text(hjust = 0.5,size = 14,colour = 'darkred',face='bold'),
        legend.position = 'none')
p1
tmp2 <- tmp[,c(5,6,7,9)]

# -------------------------------------------------------------------------
## 在选择最佳模型的时候是否纳入训练集作为评估
include_train <- F
if(include_train){
  mm <- group_by(tmp2,model)%>%summarise(Cindex=mean(Cindex))
  mm$Cse <- 0
  mm$ID <- 'Mean C-index'
  mm <- mm[,c('Cindex','Cse','ID','model')]
  tmp2 <- rbind(tmp2,mm)
}else{
  mm <- group_by(tmp2[tmp2$ID!=train_cohort_name,],model)%>%summarise(Cindex=mean(Cindex))
  mm$Cse <- 0
  mm$ID <- 'Mean C-index'
  mm <- mm[,c('Cindex','Cse','ID','model')]
  tmp2 <- rbind(tmp2,mm)
}
tmp2$ID <- factor(tmp2$ID,levels = c('Mean C-index',sort(select_cohorts[select_cohorts!=train_cohort_name],decreasing = T),train_cohort_name))
tmp2$ll2 <- sprintf('%.3f',tmp2$Cindex)

p2 <- ggplot(tmp2,aes(Cindex,ID))+
  geom_errorbar(aes(xmin=Cindex-Cse,xmax=Cindex+Cse),width=0.1,size=0.8)+
  geom_bar(stat = 'identity',width = 0.6,color='black',aes(fill=model),size=0.5)+
  geom_vline(xintercept = 0.05,linetype=2,color='grey50')+
  geom_text(aes(0.3,ID,label=ll2),tmp2[tmp2$ID=='Mean C-index',],color='white',fontface='bold',size=3.5)+
  scale_fill_manual(values = cols)+
  labs(x='Index of concordance')+
  theme_bw(base_rect_size = 0)+
  facet_wrap(~model,nrow = 1,scales = 'free_x')+
  theme(axis.text.y = element_text(size = 12,colour = 'darkred',face='bold'),
        axis.text.x = element_text(size = 8,colour = 'black'),
        axis.title.x = element_text(size = 13,colour = 'darkred',face='bold'),
        axis.title.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.background = element_rect(fill='#f3f6f6'),
        plot.title = element_text(hjust = 0.5,size = 14,colour = 'darkred',face='bold'),
        legend.position = 'none')
p1/p2

# -------------------------------------------------------------------------
#输出最佳模型
if(include_train){
  message('Based on the mean Cindex across all cohorts, the optimal model is ',mm$model[which.max(mm$Cindex)])
}else{
  message('Based on the mean Cindex across all validation cohorts, the optimal model is ',mm$model[which.max(mm$Cindex)])
}  
score_list <- result[[mm$model[which.max(mm$Cindex)]]]
save(score_list,file = 'score_list.rda')

method <- mm$model[which.max(mm$Cindex)]

results <- lzq_surmodel(select_survars=select_survars,
                        alldata=alldata[train_cohort_name],
                        method = method,##可选 StepCox, RSF, Lasso, Ridge, Enet, GBM, SVM, plsRcox, Coxboost, SuperPC
                        seed=1234,
                        StepCox_direction='backward', ##包括三个选项：forward backward both
                        RSF_nodesize=5, ## 3-30
                        RSF_nsplit=10, ## 2-20
                        RSF_splitrule="logrank", ##包括三个选项：logrank bs.gradient logrankscore
                        Lasso_lamda_rule='lambda.min', #包括两个选项 lambda.min lambda.1se
                        Ridge_lamda_rule='lambda.min', #包括两个选项 lambda.min lambda.1se
                        Enet_alpha=0.5, ##可选择0.1-0.9
                        Enet_lamda_rule='lambda.min', #包括两个选项 lambda.min lambda.1se
                        GBM_nodesize=5,
                        SVM_type= 'vanbelle1', ##包括四个选项：regression vanbelle1 vanbelle2 hybrid'
                        SVM_diffmeth= 'makediff3', ## 包括三个选项： makediff1 makediff2 and makediff3
                        SVM_optmeth= 'quadprog',   ##包括两个选项：quadprog or ipop
                        SVM_kernel='add_kernel', ##包括四个选项：lin_kernel add_kernel rbf_kernel poly_kernel
                        plsRcox_lambda_rule='lambda.min', #包括两个选项 lambda.min lambda.1se
                        Coxboost_type='verweij', #包括两个选项 verweij naive
                        SuperPC_ncomponents=1 ## Number of principal components to compute. Should be 1,2 or 3.
)
save(results,file = 'fit.rda')

# -------------------------------------------------------------------------
#可以评估每个模型
Cohort <- 'CGGA_693'
model <- 'Ridge'
##输出生存曲线，ROC calibration

## KM survival curve ##
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

lzq_survplot2(Sur_ids = select_survars,
              cohort = Cohort,
              data = result[[model]][[Cohort]],
              gene = 'score',
              Input = model,
              cols = cols,cutoff = cutoff,percent = percent)


## ROC curve ##

tmp <- result[[model]][[Cohort]]

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
time_cut <- c('1-Year AUC','2-Year AUC','3-Year AUC')
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

ggplot(dd,aes(fp,tp,color=time))+
  geom_line(size=1)+
  labs(x='1-Specificity',y='Sensitivity',color=NULL)+
  theme_bw(base_rect_size = 1.5)+
  geom_abline(slope = 1,color='grey70')+
  ggtitle(Cohort)+
  theme(axis.text.y = element_text(size = 10,colour = 'black'),
        axis.text.x = element_text(size = 10,colour = 'black'),
        axis.title = element_text(size = 13,colour = 'darkred',face='bold'),
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.background = element_rect(fill='#f3f6f6'),
        plot.title = element_text(hjust = 0.5,size = 14,colour = 'darkred',face='bold'),
        legend.text = element_text(size=12),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.position = c(0.995,0.012),
        legend.justification = c(1,0))+
  scale_color_manual(values = c('#2a9d8f','#fcbf49','#eb5e28'))+
  scale_x_continuous(expand = c(0.01,0.01))+
  scale_y_continuous(expand = c(0.01,0.01))

## Calibration curve ##
library(rms)
cox1 <- cph(gsurv,surv=T,x=T, y=T,data=tmp,time.inc = times[1]) 
cal1 <- tryCatch(calibrate(cox1, cmethod="KM", method="boot", u=times[1], m= ceiling(nrow(tmp)*0.95/5), B=1000),error=function(e){NA})
cox2 <- cph(gsurv,surv=T,x=T, y=T,data=tmp,time.inc = times[2]) 
cal2 <- tryCatch(calibrate(cox2, cmethod="KM", method="boot", u=times[2], m= ceiling(nrow(tmp)*0.95/5), B=1000),error=function(e){NA})
cox3 <- cph(gsurv,surv=T,x=T, y=T,data=tmp,time.inc = times[3]) 
cal3 <- tryCatch(calibrate(cox3, cmethod="KM", method="boot", u=times[3], m= ceiling(nrow(tmp)*0.95/5), B=1000),error=function(e){NA})

dd <- Reduce(rbind,list(tryCatch(cal1[,c('mean.predicted',"KM")],error=function(e)data.frame(mean.predicted=rep(NA,5),KM=rep(NA,5))),
                        tryCatch(cal2[,c('mean.predicted',"KM")],error=function(e)data.frame(mean.predicted=rep(NA,5),KM=rep(NA,5))),
                        tryCatch(cal3[,c('mean.predicted',"KM")],error=function(e)data.frame(mean.predicted=rep(NA,5),KM=rep(NA,5)))))%>%as.data.frame()
dd$time <- rep(gsub(' AUC','',time_cut),each=5)
dd <- na.omit(dd)
colnames(dd)[1:2] <- c('Predicted','Observed')

ggplot(dd,aes(Predicted,Observed))+
  geom_abline(slope = 1,color='grey70')+
  geom_line(size=1,alpha=0.7,aes(color=time))+
  theme_bw(base_rect_size = 1.5)+
  geom_point(size=3,shape=21,aes(fill=time),stroke=1)+
  labs(x='Predicted survival (%)',y='Observed survival (%)',title = Cohort)+
  theme(axis.text.y = element_text(size = 10,colour = 'black'),
        axis.text.x = element_text(size = 10,colour = 'black'),
        axis.title = element_text(size = 13,colour = 'darkred',face='bold'),
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.background = element_rect(fill='#f3f6f6'),
        plot.title = element_text(hjust = 0.5,size = 14,colour = 'darkred',face='bold'),
        legend.text = element_text(size=12),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.position = c(0.995,0.012),
        legend.justification = c(1,0))+
  scale_color_manual(values = c('#2a9d8f','#fcbf49','#eb5e28'))+
  scale_fill_manual(values = c('#2a9d8f','#fcbf49','#eb5e28'))+
  scale_x_continuous(expand = c(0.03,0))+
  scale_y_continuous(expand = c(0,0.03))


## DCA
tmp$score <- (tmp$score-min(tmp$score))/(max(tmp$score)-min(tmp$score))
library(dcurves)
##从times变量中选择 这就是 Timepoint for DCA的参数
select_time <- times[2]

dd <- tryCatch((dcurves::dca(gsurv, 
                             data = tmp, 
                             label = list(all='All',none='None',score=model),
                             time = select_time)),error=function(e){message('Calculating failed due to some error, please replace another timepoint!')}) 
my <- dd$dca%>%na.omit()
xmax <- max(my$threshold[my$variable=='score'])
dd <- tryCatch((dcurves::dca(gsurv, 
                             data = tmp, thresholds = seq(0,xmax,by=0.01),
                             label = list(all='All',none='None',score=model),
                             time = select_time)),error=function(e){message('Calculating failed due to some error, please replace another timepoint!')}) 

dd%>%plot(smooth=T)+
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




