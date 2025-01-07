rm(list = ls())
library(tidyverse)
library(ggsci)
source('lzq_surmodel.R')

select_cancer <- 'LUAD'
load(paste0('data/',select_cancer,'/symbol.rda'))
names <- names(total_clin_list)
Survar_Input <- 'Overall survival' ## 用户自己选择
select_survars <- ifelse(Survar_Input=='Overall survival','OS',
                         ifelse(Survar_Input=='Relapse-free survival','RFS',
                                ifelse(Survar_Input=='Disease-free survival','DFS',
                                       ifelse(Survar_Input=='Progression-free survival','PFS','DSS'))))
Select_cohorts <- names[sapply(total_clin_list,function(x){select_survars%in%colnames(x)})]

###Select_cohorts <- Select_cohorts[1:5] ###用户可以自定义 选择哪几个队列

# -------------------------------------------------------------------------
### select or input your gene list (Symbol)###
# -------------------------------------------------------------------------

Gene_Input <- c('ANLN','CCNB1','CCNB2','CDCA2','CDO1','CENPF','CENPW','CYP4B1',
                'DEPDC1','GJB3','GMPSP1','HMMR','IL22RA1','LIMCH1','MAMDC2','PARM1',
                'PLOD2','POC1A','PRKCE','RHOC','SEPTIN4','TEDC2','TPX2','TTK','UHRF1')
if(sum(Reduce(union,lapply(total_expr_list,rownames))%in%Gene_Input)==0){
  stop('No gene is in all data sets')
}

self_name <- 'Pindex' ## 自定义Score的名字 **通用**

# -------------------------------------------------------------------------
###看看基因在队列中的是否出现 ###
# -------------------------------------------------------------------------

tmp <- lapply(lapply(total_expr_list[Select_cohorts],rownames),function(x){Gene_Input%in%x})
bar <- data.frame(Freq=sapply(tmp,sum))%>%tibble::rownames_to_column('ID')

##网站输出
for (i in names(tmp)) {
  if(sum(tmp[[i]])!=length(Gene_Input)){
    cat(paste0(i,' absents ',paste0(Gene_Input[!tmp[[i]]],collapse = ', '),'\n'))
  }
}

##网站输出
cols <- c(pal_npg(alpha = 0.8)(10),pal_d3(alpha = 0.8)(10)[-c(3:6)],pal_jco(alpha = 0.8)(10),pal_nejm(alpha = 0.8)(8),pal_lancet(alpha = 0.8)(8),pal_jama(alpha = 0.8)(7))
ggplot(bar,aes(Freq,reorder(ID,Freq),fill=ID))+
  geom_bar(stat = 'identity',width=0.6,color='black')+
  geom_text(aes(label=Freq),hjust=1.5,color='white',fontface='bold')+
  labs(x = 'Frequence')+
  theme_bw(base_rect_size = 0)+
  scale_fill_manual(values = cols)+
  scale_x_continuous(expand = c(0,0))+
  theme(axis.text.y = element_text(size = 12,colour = 'black'),
        axis.text.x = element_text(size = 11,colour = 'black'),
        axis.title.x = element_text(size = 13,colour = 'black'),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.background = element_rect(fill='#f3f6f6'),
        legend.position = 'none',
        legend.text = element_text(colour = 'black'),
        legend.title = element_blank(),
        legend.background = element_rect(fill = NA),
        legend.key.width = unit(4,'mm'),
        legend.key.height = unit(4,'mm'))

# -------------------------------------------------------------------------
### 选择训练集和验证集 ###
# -------------------------------------------------------------------------

##可供选择的队列名（只考虑包含所有基因的队列）：

availble_cohorts <- bar$ID[bar$Freq==length(Gene_Input)]

train_cohort_name <- 'TCGA_LUAD'
# test_cohort_name <- availble_cohorts[availble_cohorts!=train_cohort_name]
test_cohort_name <- availble_cohorts[!availble_cohorts%in%c(train_cohort_name,'GSE19188')]
all_cohort_name <- c(train_cohort_name,test_cohort_name)
alldata <- lapply(all_cohort_name,function(x){
  merge(na.omit(total_clin_list[[x]][,c('ID',paste0(select_survars,'.time'),select_survars)]),na.omit(t(total_expr_list[[x]][Gene_Input,])),by.x=1,by.y=0)%>%tibble::column_to_rownames('ID')
})
alldata <- lapply(alldata,function(x){x[x[,1]>0,]})
names(alldata) <- all_cohort_name

# -------------------------------------------------------------------------
### 选择机器学习方法 ###
# -------------------------------------------------------------------------

method <- 'RSF'
##可选 StepCox, RSF, Lasso, Ridge, Enet, GBM, SVM, plsRcox, Coxboost, SuperPC
results <- lzq_surmodel(select_survars=select_survars,
                           alldata=alldata,
                           method = method,##可选 StepCox, RSF, Lasso, Ridge, Enet, GBM, SVM, plsRcox, Coxboost, SuperPC
                           seed=1234,
                           StepCox_direction='backward', ##包括三个选项：forward backward both
                           RSF_nodesize=10, ## 3-30
                           RSF_nsplit=10, ## 2-20
                           RSF_splitrule="logrank", ##包括三个选项：logrank bs.gradient logrankscore
                           Lasso_lamda_rule='lambda.min', #包括两个选项 lambda.min lambda.1se
                           Ridge_lamda_rule='lambda.min', #包括两个选项 lambda.min lambda.1se
                           Enet_alpha=0.5, ##可选择0.1-0.9
                           Enet_lamda_rule='lambda.min', #包括两个选项 lambda.min lambda.1se
                           GBM_nodesize=10,
                           SVM_type= 'vanbelle1', ##包括四个选项：regression vanbelle1 vanbelle2 hybrid'
                           SVM_diffmeth= 'makediff3', ## 包括三个选项： makediff1 makediff2 and makediff3
                           SVM_optmeth= 'quadprog',   ##包括两个选项：quadprog or ipop
                           SVM_kernel='add_kernel', ##包括四个选项：lin_kernel add_kernel rbf_kernel poly_kernel
                           plsRcox_lambda_rule='lambda.min', #包括两个选项 lambda.min lambda.1se
                           Coxboost_type='verweij', #包括两个选项 verweij naive
                           SuperPC_ncomponents=1 ## Number of principal components to compute. Should be 1,2 or 3.
)
score_list <- results$rs
save(Gene_Input,select_survars,all_cohort_name,self_name,results,score_list,file = 'score_list.rda')

load(paste0('data/',select_cancer,'/cm_surlist.rda'))
tmp <- surlist[[select_survars]][,all_cohort_name]%>%na.omit()
tmp <- tmp[Gene_Input,]
for (i in 1:ncol(tmp)) {
  tmp[,i] <- ifelse(tmp[,i]>0&abs(tmp[,i])<0.05,'Risky',ifelse(tmp[,i]<0&abs(tmp[,i])<0.05,'Protective','Not Significant'))
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
library(ComplexHeatmap)
draw(Heatmap(t(tmp),col = cols,
             name = 'Type',border = T,
             rect_gp = gpar(col='black',lwd=2),
             row_names_side = 'left',
             column_names_side = 'bottom',
             row_split = 1:ncol(tmp),
             row_title = NULL,
             width = ncol(t(tmp)) * unit(5,'mm'),
             height = nrow(t(tmp)) * unit(6,'mm'),
             heatmap_legend_param=list(labels_gp = gpar(fontsize = 12), border = T,nrow=1,title_position = 'leftcenter',
                                       title_gp = gpar(fontsize = 12, fontface = "bold"))), heatmap_legend_side = "top")
# dev.copy2pdf(file='PS-feature prognosis.pdf',width=10,height=10)

