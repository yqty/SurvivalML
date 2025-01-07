lzq_consensus_model <- function(
  select_survars,
  alldata,
  traindata,
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
  
  results <- list()
  
  fit <- tryCatch(step(coxph(gsurv,traindata),direction = StepCox_direction),error=function(e)NA)
  rs <- tryCatch(lapply(alldata,function(x){cbind(x[,1:2],score=as.numeric(scale(predict(fit,type = 'risk',newdata = x))))}),error=function(e)NA)
  results[['StepCox']] <- rs
  results <- results[!is.na(results)]
  
  fit <- rfsrc(formula = gsurv,
               data = traindata,
               ntree = 1000,
               nodesize = RSF_nodesize,
               nsplit = RSF_nsplit,
               splitrule = RSF_splitrule,
               importance = T,
               seed = seed)
  best <- which.min(fit$err.rate)
  fit <- rfsrc(formula = gsurv,
               data = traindata,
               ntree = best,
               nodesize = RSF_nodesize,
               nsplit = RSF_nsplit,
               splitrule = RSF_splitrule,
               importance = T,
               seed = seed)
  rs <- lapply(alldata,function(x){cbind(x[,1:2],score=as.numeric(scale(predict(fit,newdata = x)$predicted)))})
  results[['RSF']] <- rs
  
  
  set.seed(seed)
  fit = cv.glmnet(as.matrix(traindata[,-c(1,2)]), as.matrix(Surv(traindata[,1],traindata[,2])),family = "cox",alpha=1,nfolds = 10)
  lambda <- ifelse(Lasso_lamda_rule=='lambda.min',fit$lambda.min,fit$lambda.1se)
  rs <- lapply(alldata,function(x){cbind(x[,1:2],score=as.numeric(scale(predict(fit,type='link',newx=as.matrix(x[,-c(1,2)]),s=lambda))))})
  results[['Lasso']] <- rs
  
  set.seed(seed)
  fit = cv.glmnet(as.matrix(traindata[,-c(1,2)]), as.matrix(Surv(traindata[,1],traindata[,2])),family = "cox",alpha=0,nfolds = 10)
  lambda <- ifelse(Ridge_lamda_rule=='lambda.min',fit$lambda.min,fit$lambda.1se)
  rs <- lapply(alldata,function(x){cbind(x[,1:2],score=as.numeric(scale(predict(fit,type='link',newx=as.matrix(x[,-c(1,2)]),s=lambda))))})
  results[['Ridge']] <- rs
  
  set.seed(seed)
  fit = cv.glmnet(as.matrix(traindata[,-c(1,2)]), as.matrix(Surv(traindata[,1],traindata[,2])),family = "cox",alpha=Enet_alpha,nfolds = 10)
  lambda <- ifelse(Enet_lamda_rule=='lambda.min',fit$lambda.min,fit$lambda.1se)
  rs <- lapply(alldata,function(x){cbind(x[,1:2],score=as.numeric(scale(predict(fit,type='link',newx=as.matrix(x[,-c(1,2)]),s=lambda))))})
  results[['Enet']] <- rs
  
  set.seed(seed)
  fit <- gbm(formula = gsurv,
             data = traindata,
             distribution = 'coxph',
             n.trees = 1000,
             interaction.depth = 3,
             n.minobsinnode = GBM_nodesize,
             shrinkage = 0.001,
             cv.folds = 10,n.cores = 6)
  rs <- lapply(alldata,function(x){cbind(x[,1:2],score=as.numeric(scale(predict(fit,x,n.trees = which.min(fit$cv.error),type = 'link'))))})
  results[['GBM']] <- rs
  
  fit = survivalsvm(formula = gsurv, data = traindata, gamma.mu = 0.1,
                    type = SVM_type, diff.meth = SVM_diffmeth, 
                    opt.meth = SVM_optmeth, kernel = SVM_kernel,
                    sgf.sv = 5,sigf = 7, maxiter = 20, margin = 0.05, 
                    bound = 10, eig.tol = 1e-06,
                    conv.tol = 1e-07, posd.tol = 1e-08)
  rs <- lapply(alldata,function(x){cbind(x[,1:2],score=as.numeric(predict(fit,x)$predicted))})
  results[['SVM']] <- rs
  
  set.seed(seed)
  nfold <- ifelse(nrow(alldata[[1]])<=60,3,ifelse(nrow(alldata[[1]])<=100,5,10))
  cv.plsRcox.res=tryCatch(cv.plsRcox(list(x=traindata[,-c(1,2)],time=traindata[,1],status=traindata[,2]),
                                     nt=10,verbDFSe = FALSE,nfold = nfold),
                          error=function(e){
                            tryCatch(cv.plsRcox(list(x=traindata[,-c(1,2)],time=traindata[,1],status=traindata[,2]),
                                                nt=1,verbDFSe = FALSE,nfold = 5),
                                     error=function(e){tryCatch(cv.plsRcox(list(x=traindata[,-c(1,2)],time=traindata[,1],status=traindata[,2]),
                                                                           nt=1,verbDFSe = FALSE,nfold = 3),
                                                                error=function(e){cv.plsRcox(list(x=traindata[,-c(1,2)],time=traindata[,1],status=traindata[,2]),
                                                                                             nt=1,verbDFSe = FALSE,nfold = 2)})
                                     })})
  lambda <- ifelse(plsRcox_lambda_rule=='lambda.min',cv.plsRcox.res[[5]],cv.plsRcox.res[[6]])
  fit <- plsRcox(traindata[,-c(1,2)],time=traindata[,1],event=traindata[,2],nt=lambda)
  rs <- lapply(alldata,function(x){cbind(x[,1:2],score=as.numeric(scale(predict(fit,type="lp",newdata=x[,-c(1,2)]))))})
  results[['plsRcox']] <- rs
  
  set.seed(seed)
  pen <- optimCoxBoostPenalty(traindata[,1],traindata[,2],as.matrix(traindata[,-c(1,2)]),
                              trace=TRUE,start.penalty=1000,parallel = T)
  cv.res <- cv.CoxBoost(traindata[,1],traindata[,2],as.matrix(traindata[,-c(1,2)]),
                        maxstepno=1000,K=10,type=Coxboost_type,penalty=pen$penalty)
  fit <- CoxBoost(traindata[,1],traindata[,2],as.matrix(traindata[,-c(1,2)]),
                  stepno=cv.res$optimal.step,penalty=pen$penalty)
  rs <- lapply(alldata,function(x){cbind(x[,1:2],score=as.numeric(predict(fit,newdata=x[,-c(1,2)], newtime=x[,1], newstatus=x[,2], type="lp")))})
  results[['Coxboost']] <- rs
  
  data <- list(x=t(traindata[,-c(1,2)]),y=traindata[,1],censoring.status=traindata[,2],featurenames=colnames(traindata)[-c(1,2)])
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
  rs <- tryCatch(lapply(alldata,function(w){
    test <- list(x=t(w[,-c(1,2)]),y=w[,1],censoring.status=w[,2],featurenames=colnames(w)[-c(1,2)])
    ff <- superpc.predict(fit,data,test,threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = SuperPC_ncomponents)
    rr <- as.numeric(ff$v.pred)
    rr2 <- cbind(w[,1:2],score=rr)
    return(rr2)}),error=function(e){NA})
  results[['SuperPC']] <- rs
  results <- results[!is.na(results)]
  
  return(results)
}

