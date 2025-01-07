lzq_refit <- function(testdata){
  if(results$select_survars=='OS'){
    gsurv <-  as.formula(Surv(OS.time,OS)~.)
  }
  if(results$select_survars=='RFS'){
    gsurv <-  as.formula(Surv(RFS.time,RFS)~.)
  }
  if(results$select_survars=='DFS'){
    gsurv <-  as.formula(Surv(DFS.time,DFS)~.)
  }
  if(results$select_survars=='PFS'){
    gsurv <-  as.formula(Surv(PFS.time,PFS)~.)
  }
  if(results$select_survars=='DSS'){
    gsurv <-  as.formula(Surv(DSS.time,DSS)~.)
  }
  
  library(survival)
  library(randomForestSRC)
  library(glmnet)
  library(gbm)
  library(survivalsvm)
  library(plsRcox)
  library(superpc)
  library(CoxBoost)

  if(results$method=='StepCox'){
    rs <- lapply(testdata,function(x){data.frame(ID=rownames(x),score=as.numeric(predict(results$fit,type = 'risk',newdata = x)))})
  }
  
  if(results$method=='RSF'){
    rs <- lapply(testdata,function(x){data.frame(ID=rownames(x),score=as.numeric(predict(results$fit,newdata = x)$predicted))})
  }
  
  if(results$method=='Lasso'){
    rs <- lapply(testdata,function(x){data.frame(ID=rownames(x),score=as.numeric(predict(results$fit,type='link',newx=as.matrix(x),s=results$lambda)))})
  }
  
  if(results$method=='Ridge'){
    rs <- lapply(testdata,function(x){data.frame(ID=rownames(x),score=as.numeric(predict(results$fit,type='link',newx=as.matrix(x),s=results$lambda)))})
  }
  
  if(results$method=='Enet'){
    rs <- lapply(testdata,function(x){data.frame(ID=rownames(x),score=as.numeric(predict(results$fit,type='link',newx=as.matrix(x),s=results$lambda)))})
  }
  
  if(results$method=='GBM'){
    rs <- lapply(testdata,function(x){data.frame(ID=rownames(x),score=as.numeric(predict(results$fit,x,n.trees = results$n_tree,type = 'link')))})
  }
  
  if(results$method=='SVM'){
    rs <- lapply(testdata,function(x){data.frame(ID=rownames(x),score=as.numeric(predict(results$fit,x)$predicted))})
  }
  
  if(results$method=='plsRcox'){
    rs <- lapply(testdata,function(x){data.frame(ID=rownames(x),score=as.numeric(predict(results$fit,type="lp",newdata=x)))})
  }
  
  if(results$method=='Coxboost'){
    rs <- lapply(testdata,function(x){data.frame(ID=rownames(x),score=as.numeric(predict(results$fit,newdata=x, type="lp")))})
  }
  
  if(results$method=='SuperPC'){
    rs <- lapply(testdata,function(w){
      test <- list(x=t(w),featurenames=colnames(w))
      ff <- superpc.predict(results$fit,results$data,test,threshold = results$threshold,n.components = results$SuperPC_ncomponents)
      rr <- as.numeric(ff$v.pred)
      rr2 <- data.frame(ID=rownames(w),score=rr)
      return(rr2)})
  }
  return(rs)
}