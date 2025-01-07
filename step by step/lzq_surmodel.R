lzq_surmodel <- function(
  select_survars,
  alldata,
  method,
  seed=1234,
  StepCox_direction='backward', ##包括三个选项：forward backward both
  RSF_nodesize=5,
  RSF_nsplit=10,
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

){
  library(survival)
  library(randomForestSRC)
  library(glmnet)
  library(gbm)
  library(survivalsvm)
  library(plsRcox)
  library(superpc)
  library(CoxBoost)

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

  if(method=='StepCox'){
    fit <- tryCatch(step(coxph(gsurv,alldata[[1]]),direction = StepCox_direction),error=function(e)NA)
    rs <- tryCatch(lapply(alldata,function(x){cbind(x[,1:2],score=as.numeric(predict(fit,type = 'risk',newdata = x)))}),error=function(e)NA)
  }

  if(method=='RSF'){
    fit <- rfsrc(formula = gsurv,
                 data = alldata[[1]],
                 ntree = 10000,
                 nodesize = RSF_nodesize,
                 nsplit = RSF_nsplit,
                 splitrule = RSF_splitrule,
                 importance = T,
                 seed = seed)
    best <- which.min(fit$err.rate)
    fit <- rfsrc(formula = gsurv,
                 data = alldata[[1]],
                 ntree = best,
                 nodesize = RSF_nodesize,
                 nsplit = RSF_nsplit,
                 splitrule = RSF_splitrule,
                 importance = T,
                 seed = seed)
    rs <- lapply(alldata,function(x){cbind(x[,1:2],score=as.numeric(predict(fit,newdata = x)$predicted))})
  }

  if(method=='Lasso'){
    set.seed(seed)
    fit = cv.glmnet(as.matrix(alldata[[1]][,-c(1,2)]), as.matrix(Surv(alldata[[1]][,1],alldata[[1]][,2])),family = "cox",alpha=1,nfolds = 10)
    lambda <- ifelse(Lasso_lamda_rule=='lambda.min',fit$lambda.min,fit$lambda.1se)
    rs <- lapply(alldata,function(x){cbind(x[,1:2],score=as.numeric(predict(fit,type='link',newx=as.matrix(x[,-c(1,2)]),s=lambda)))})
  }

  if(method=='Ridge'){
    set.seed(seed)
    fit = cv.glmnet(as.matrix(alldata[[1]][,-c(1,2)]), as.matrix(Surv(alldata[[1]][,1],alldata[[1]][,2])),family = "cox",alpha=0,nfolds = 10)
    lambda <- ifelse(Ridge_lamda_rule=='lambda.min',fit$lambda.min,fit$lambda.1se)
    rs <- lapply(alldata,function(x){cbind(x[,1:2],score=as.numeric(predict(fit,type='link',newx=as.matrix(x[,-c(1,2)]),s=lambda)))})
  }

  if(method=='Enet'){
    set.seed(seed)
    fit = cv.glmnet(as.matrix(alldata[[1]][,-c(1,2)]), as.matrix(Surv(alldata[[1]][,1],alldata[[1]][,2])),family = "cox",alpha=Enet_alpha,nfolds = 10)
    lambda <- ifelse(Enet_lamda_rule=='lambda.min',fit$lambda.min,fit$lambda.1se)
    rs <- lapply(alldata,function(x){cbind(x[,1:2],score=as.numeric(predict(fit,type='link',newx=as.matrix(x[,-c(1,2)]),s=lambda)))})
  }

  if(method=='GBM'){
    set.seed(seed)
    fit <- tryCatch(gbm(formula = gsurv,
                        data = alldata[[1]],
                        distribution = 'coxph',
                        n.trees = 10000,
                        interaction.depth = 3,
                        n.minobsinnode = GBM_nodesize,
                        shrinkage = 0.001,
                        cv.folds = 10,n.cores = 6),error=function(e)NA)
    n_tree <- tryCatch(which.min(fit$cv.error),error=function(e)NA)
    rs <- tryCatch(lapply(alldata,function(x){cbind(x[,1:2],score=as.numeric(predict(fit,x,n.trees = n_tree,type = 'link')))}),error=function(e)NA)
  }

  if(method=='SVM'){
    fit = survivalsvm(formula = gsurv, data = alldata[[1]], gamma.mu = 0.1,
                      type = SVM_type, diff.meth = SVM_diffmeth,
                      opt.meth = SVM_optmeth, kernel = SVM_kernel,
                      sgf.sv = 5,sigf = 7, maxiter = 20, margin = 0.05,
                      bound = 10, eig.tol = 1e-06,
                      conv.tol = 1e-07, posd.tol = 1e-08)
    rs <- lapply(alldata,function(x){cbind(x[,1:2],score=as.numeric(predict(fit,x)$predicted))})
  }

  if(method=='plsRcox'){
    set.seed(seed)
    nfold <- ifelse(nrow(alldata[[1]])<=60,3,ifelse(nrow(alldata[[1]])<=100,5,10))
    cv.plsRcox.res=tryCatch(cv.plsRcox(list(x=alldata[[1]][,-c(1,2)],time=alldata[[1]][,1],status=alldata[[1]][,2]),
                                       nt=10,verbDFSe = FALSE,nfold = nfold),
                            error=function(e){
                              tryCatch(cv.plsRcox(list(x=alldata[[1]][,-c(1,2)],time=alldata[[1]][,1],status=alldata[[1]][,2]),
                                                  nt=1,verbDFSe = FALSE,nfold = 5),
                                       error=function(e){tryCatch(cv.plsRcox(list(x=alldata[[1]][,-c(1,2)],time=alldata[[1]][,1],status=alldata[[1]][,2]),
                                                                             nt=1,verbDFSe = FALSE,nfold = 3),
                                                                  error=function(e){cv.plsRcox(list(x=alldata[[1]][,-c(1,2)],time=alldata[[1]][,1],status=alldata[[1]][,2]),
                                                                                               nt=1,verbDFSe = FALSE,nfold = 2)})
                                       })})
    lambda <- ifelse(plsRcox_lambda_rule=='lambda.min',cv.plsRcox.res[[5]],cv.plsRcox.res[[6]])
    fit <- plsRcox(alldata[[1]][,-c(1,2)],time=alldata[[1]][,1],event=alldata[[1]][,2],nt=lambda)
    rs <- lapply(alldata,function(x){cbind(x[,1:2],score=as.numeric(predict(fit,type="lp",newdata=x[,-c(1,2)])))})
  }

  if(method=='Coxboost'){
    set.seed(seed)
    pen <- optimCoxBoostPenalty(alldata[[1]][,1],alldata[[1]][,2],as.matrix(alldata[[1]][,-c(1,2)]),
                                trace=TRUE,start.penalty=1000,parallel = T)
    cv.res <- cv.CoxBoost(alldata[[1]][,1],alldata[[1]][,2],as.matrix(alldata[[1]][,-c(1,2)]),
                          maxstepno=1000,K=10,type=Coxboost_type,penalty=pen$penalty)
    fit <- CoxBoost(alldata[[1]][,1],alldata[[1]][,2],as.matrix(alldata[[1]][,-c(1,2)]),
                    stepno=cv.res$optimal.step,penalty=pen$penalty)
    rs <- lapply(alldata,function(x){cbind(x[,1:2],score=as.numeric(predict(fit,newdata=x[,-c(1,2)], newtime=x[,1], newstatus=x[,2], type="lp")))})
  }

  if(method=='SuperPC'){
    data <- list(x=t(alldata[[1]][,-c(1,2)]),y=alldata[[1]][,1],censoring.status=alldata[[1]][,2],featurenames=colnames(alldata[[1]])[-c(1,2)])
    set.seed(seed)
    fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5)
    set.seed(seed)
    cv.fit <- superpc.cv(fit,data,n.threshold = 20,
                         n.fold = 10,
                         n.components=SuperPC_ncomponents,
                         min.features=0,
                         max.features=nrow(data$x),
                         compute.fullcv= TRUE,
                         compute.preval=TRUE)
    threshold <- cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])]
    rs <- lapply(alldata,function(w){
      test <- list(x=t(w[,-c(1,2)]),y=w[,1],censoring.status=w[,2],featurenames=colnames(w)[-c(1,2)])
      ff <- tryCatch(superpc.predict(fit,data,test,threshold = threshold,n.components = SuperPC_ncomponents),error=function(e){NA})
      rr <- tryCatch(as.numeric(ff$v.pred),error=function(e){rep(0,nrow(w))})
      rr2 <- cbind(w[,1:2],score=rr)
      return(rr2)})
  }
  rs <- lapply(rs,function(x){
    x$score[x$score%in%c(Inf,-Inf)] <- 0
    x$score <- scale(x$score)%>%as.numeric()
    return(x)})
  if(sum(is.na(rs))>=1){
    message('Fitting failed due to some error, please replace a method or genelist')
  }else{
    if(sum(sapply(rs,function(x){length(unique(x$score))==1}))>0){
      message('Fitting failed due to some error, please replace a method or genelist')
    }
  }
  return(list(rs=rs,
              select_survars=select_survars,
              method=method,
              fit=fit,
              lambda = tryCatch(lambda,error=function(e)NA),
              n_tree = tryCatch(n_tree,error=function(e)NA),
              threshold = tryCatch(threshold,error=function(e)NA),
              SuperPC_ncomponents=SuperPC_ncomponents,
              data = tryCatch(data,error=function(e)NA)
  ))
}
