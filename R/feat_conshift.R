#A is a binary vector
#X is a data frame of covariates
#Y is a data frame of outcomes/labels

#' @export
feat_conshift<-function(learners,loss.fun,A,X,Y,nfolds=5,
                   cond.loss.SL.control=list(family=gaussian(),SL.library=c("SL.glm","SL.gam","SL.randomForest")),
                   conf.level=.95){
    #check inputs
    if(any(!(A %in% c(0,1)))){
        stop("A must take values 0 or 1")
    }
    if(all(A==0) || all(A==1)){
        stop("A must contain both values 0 and 1")
    }
    n<-length(A)
    
    assert_that(is.data.frame(X))
    assert_that(nrow(X)==n)
    assert_that(is.data.frame(Y))
    assert_that(nrow(Y)==n)
    assert_that(is.number(conf.level) && conf.level>=0.5 && conf.level<1)
    
    assert_that(!anyNA(A))
    assert_that(!anyNA(X))
    assert_that(!anyNA(Y[A==0,]))
    
    if(any(c("Y","X","newX") %in% names(cond.loss.SL.control))){
        stop("Y, X and newX must not be specified in cond.loss.SL.control")
    }
    if(!("family" %in% names(cond.loss.SL.control))){
        cond.loss.SL.control<-c(list(family=gaussian()),cond.loss.SL.control)
    }
    
    if(is.function(learners)){
        learners<-list(learners)
    }else if(!is.list(learners)){
        stop("learners must be a function or a list of functions")
    }else if(!all(sapply(learners,is.function))){
        stop("Some elements of learners are not functions")
    }
    assert_that(is.function(loss.fun))
    
    folds<-create.folds(A,nfolds)
    fold.sizes<-sapply(folds,length)
    
    risk.est<-SE<-numeric(length(learners))
    for(learner.index in 1:length(learners)){
        learner<-learners[[learner.index]]
        loss<-as.numeric(loss.fun(X,Y,learner))
        cond.loss<-CV.target.est(A,X,as.matrix(loss),folds,cond.loss.SL.control)
        
        risk.v<-sigma2.v<-matrix(nrow=nfolds,ncol=1)
        for(v in 1:nfolds){
            A.test<-A[folds[[v]]]
            loss.test<-loss[folds[[v]]]
            cond.loss.test<-cond.loss[[v]][,1]
            gamma.test<-mean(A.test)
            
            risk.v[v]<-mean(ifelse(A.test==0,(loss.test-cond.loss.test)/(1-gamma.test),0)+cond.loss.test)
            IF<-ifelse(A.test==0,(loss.test-cond.loss.test)/(1-gamma.test),0)+cond.loss.test-risk.v[v]-(mean(cond.loss.test)-risk.v[v])/(1-gamma.test)*(A.test-gamma.test)
            sigma2.v[v]<-mean(IF^2)
        }
        
        risk.est[learner.index]<-t(fold.sizes)%*%risk.v/n
        sigma2<-t(fold.sizes)%*%sigma2.v/n
        SE[learner.index]<-sqrt(sigma2/n)
    }
    
    z<-qnorm((conf.level+1)/2)
    data.frame(risk.est=risk.est,SE=SE,CI.lower=risk.est-z*SE,CI.upper=risk.est+z*SE)
}








#A is a binary vector
#X is a data frame of covariates
#Y is a data frame of outcomes/labels with one column

#' @export
linMSE_feat_conshift<-function(learners,A,X,Y,nfolds=5,
                          cond.Y.SL.control=list(family=gaussian(),SL.library=c("SL.glm","SL.gam","SL.randomForest")),
                          conf.level=.95){
    #check inputs
    if(any(!(A %in% c(0,1)))){
        stop("A must take values 0 or 1")
    }
    if(all(A==0) || all(A==1)){
        stop("A must contain both values 0 and 1")
    }
    n<-length(A)
    
    assert_that(is.data.frame(X))
    assert_that(nrow(X)==n)
    assert_that(is.data.frame(X))
    assert_that(nrow(Y)==n)
    assert_that(is.number(conf.level) && conf.level>=0.5 && conf.level<1)
    
    assert_that(ncol(Y)==1)
    Y.vec<-Y[,1,drop=TRUE]
    
    assert_that(!anyNA(A))
    assert_that(!anyNA(X))
    assert_that(!anyNA(Y))
    
    if(any(c("Y","X","newX") %in% names(cond.Y.SL.control))){
        stop("Y, X and newX must not be specified in cond.Y.SL.control")
    }
    if(!("family" %in% names(cond.Y.SL.control))){
        cond.Y.SL.control<-c(list(family=gaussian()),cond.Y.SL.control)
    }
    
    if(is.function(learners)){
        learners<-list(learners)
    }else if(!is.list(learners)){
        stop("learners must be a function or a list of functions")
    }else if(!all(sapply(learners,is.function))){
        stop("Some elements of learners are not functions")
    }
    
    folds<-create.folds(A,nfolds)
    fold.sizes<-sapply(folds,length)
    
    cond.Y<-CV.target.est(A,X,as.matrix(Y.vec),folds,cond.Y.SL.control)
    fx<-do.call(cbind,lapply(learners,function(learner){
        learner(X)
    }))
    loss<-fx*(-2*matrix(Y.vec,nrow=n,ncol=length(learners),byrow=FALSE)+fx)
    
    risk.est<-SE<-numeric(length(learners))
    for(learner.index in 1:length(learners)){
        risk.v<-sigma2.v<-matrix(nrow=nfolds,ncol=1)
        for(v in 1:nfolds){
            A.test<-A[folds[[v]]]
            cond.Y.test<-cond.Y[[v]][,1]
            gamma.test<-mean(A.test)
            fx.test<-fx[folds[[v]],learner.index]
            loss.test<-loss[folds[[v]],learner.index]
            cond.loss.test<-fx.test*(-2*cond.Y.test+fx.test)
            
            risk.v[v]<-mean(ifelse(A.test==0,(loss.test-cond.loss.test)/(1-gamma.test),0)+cond.loss.test)
            IF<-ifelse(A.test==0,(loss.test-cond.loss.test)/(1-gamma.test),0)+cond.loss.test-risk.v[v]-(mean(cond.loss.test)-risk.v[v])/(1-gamma.test)*(A.test-gamma.test)
            sigma2.v[v]<-mean(IF^2)
        }
        
        risk.est[learner.index]<-t(fold.sizes)%*%risk.v/n
        sigma2<-t(fold.sizes)%*%sigma2.v/n
        SE[learner.index]<-sqrt(sigma2/n)
    }
    
    z<-qnorm((conf.level+1)/2)
    data.frame(risk.est=risk.est,SE=SE,CI.lower=risk.est-z*SE,CI.upper=risk.est+z*SE)
}








#A is a binary vector
#X is a data frame of covariates
#Y is a data frame of outcomes/labels with one column

#' @export
cross_entropy_feat_conshift<-function(learners,A,X,Y,nfolds=5,
                                      cond.Y.SL.control=list(family=binomial(),SL.library=c("SL.glm","SL.gam","SL.randomForest")),
                                      conf.level=.95){
    #check inputs
    if(any(!(A %in% c(0,1)))){
        stop("A must take values 0 or 1")
    }
    if(all(A==0) || all(A==1)){
        stop("A must contain both values 0 and 1")
    }
    n<-length(A)
    
    if(any(!(Y %in% c(0,1)))){
        stop("Y must take values 0 or 1")
    }
    
    assert_that(is.data.frame(X))
    assert_that(nrow(X)==n)
    assert_that(is.data.frame(Y))
    assert_that(nrow(Y)==n)
    assert_that(is.number(conf.level) && conf.level>=0.5 && conf.level<1)
    
    assert_that(ncol(Y)==1)
    Y.vec<-Y[,1,drop=TRUE]
    assert_that(all(Y.vec %in% c(0,1)))
    
    assert_that(!anyNA(A))
    assert_that(!anyNA(X))
    assert_that(!anyNA(Y))
    
    if(any(c("Y","X","newX") %in% names(cond.Y.SL.control))){
        stop("Y, X and newX must not be specified in cond.Y.SL.control")
    }
    if(!("family" %in% names(cond.Y.SL.control))){
        cond.Y.SL.control<-c(list(family=gaussian()),cond.Y.SL.control)
    }
    if(any(c("Y","X","newX") %in% names(prop.score.SL.control))){
        stop("Y, X and newX must not be specified in prop.score.SL.control")
    }
    
    if(is.function(learners)){
        learners<-list(learners)
    }else if(!is.list(learners)){
        stop("learners must be a function or a list of functions")
    }else if(!all(sapply(learners,is.function))){
        stop("Some elements of learners are not functions")
    }
    
    folds<-create.folds(A,nfolds)
    fold.sizes<-sapply(folds,length)
    
    cond.Y<-CV.target.est(A,X,as.matrix(Y.vec),folds,cond.Y.SL.control)
    fx<-do.call(cbind,lapply(learners,function(learner){
        learner(X)
    }))
    loss<--matrix(Y.vec,nrow=n,ncol=length(learners),byrow=FALSE)*logit(fx)-log(1-fx)
    
    risk.est<-SE<-numeric(length(learners))
    for(learner.index in 1:length(learners)){
        risk.v<-sigma2.v<-matrix(nrow=nfolds,ncol=1)
        for(v in 1:nfolds){
            A.test<-A[folds[[v]]]
            cond.Y.test<-cond.Y[[v]][,1]
            gamma.test<-mean(A.test)
            fx.test<-fx[folds[[v]],learner.index]
            loss.test<-loss[folds[[v]],learner.index]
            cond.loss.test<-fx.test*(-cond.Y.test*logit(fx.test)-log(1-fx.test))
            
            risk.v[v]<-mean(ifelse(A.test==0,(loss.test-cond.loss.test)/(1-gamma.test),0)+cond.loss.test)
            IF<-ifelse(A.test==0,(loss.test-cond.loss.test)/(1-gamma.test),0)+cond.loss.test-risk.v[v]-(mean(cond.loss.test)-risk.v[v])/(1-gamma.test)*(A.test-gamma.test)
            sigma2.v[v]<-mean(IF^2)
        }
        
        risk.est[learner.index]<-t(fold.sizes)%*%risk.v/n
        sigma2<-t(fold.sizes)%*%sigma2.v/n
        SE[learner.index]<-sqrt(sigma2/n)
    }
    
    z<-qnorm((conf.level+1)/2)
    data.frame(risk.est=risk.est,SE=SE,CI.lower=risk.est-z*SE,CI.upper=risk.est+z*SE)
}
