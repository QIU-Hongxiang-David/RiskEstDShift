#' @title Risk estimation and inference under concept shift in feature/covariate
#' @name feat_conshift
#' @export
#' @description Using data from both source population (label/outcome may be missing) and target population (fully observed), estimate the risk of a list of learners in the target population under concept shift (in feature/covariate) assumption (that is, feature/covariate distributions are identical in both populations but the distributions of label/outcome given feature/covariate might differ). The estimator is semiparametrically efficient when the nuisance function is consistently estimated and is asymptotically normal even if the nuisance function is inconsistently estimated.
#' @param learners a list of learners (functions) that can take covariate `X` as input; can also be a function if only one learner is of interest
#' @param loss.fun a user specified loss function depending on covariate `X`, outcome `Y` and the learner that defines the risk (expected loss in the target population for the learner). Must take in three arguments `X`, `Y` and `learner` (an element of `learners`) and output a numeric vector with length equal to the total sample size (may be missing for data points from source population)
#' @param A a vector of binary indicator of data point coming from source population (encoded by `1`) or target population (encoded by `0`)
#' @param X data frame of features/covariates
#' @param Y data frame of labels/outcomes
#' @param nfolds number of folds used for the cross-fit estimator; default to 5
#' @param cond.loss.SL.control a named list containing options passed to \code{\link[SuperLearner:SuperLearner]{SuperLearner::SuperLearner}} to estimate the loss conditional on feature/covariate `X` (a nuisance function). Must not specify `Y`, `X` or `newX`. Default to `list(family=gaussian(),SL.library=c("SL.glm","SL.gam","SL.randomForest"))`
#' @param conf.level confidence level in the output (asymptotic) confidence interval
#' @param losses alternative input in place of `learners`, `loss.fun` and `Y`; a matrix of loss function evaluated at each data point and for each learner (row index is for data point; column index is for learner); default to `NULL`
#' @return a list with the following elements:
#' * `estimates`: a data frame with each row being the output for one learner in `learners`
#' * `empirical.IF`: a matrix of influence functions evaluated at empirical data poins. Each column corresponds to a learner in `learners` and each row corresponds to a data point
#' @section Using `empirical.IF` to infer contrasts of risks:
#' See \link{covshift}
feat_conshift<-function(learners,loss.fun,A,X,Y,nfolds=5,
                   cond.loss.SL.control=list(family=gaussian(),SL.library=c("SL.glm","SL.gam","SL.randomForest")),
                   conf.level=.95,
                   losses=NULL){
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
    assert_that(is.number(conf.level) && conf.level>=0.5 && conf.level<1)
    
    assert_that(!anyNA(A))
    assert_that(!anyNA(X))
    
    if(any(c("Y","X","newX") %in% names(cond.loss.SL.control))){
        stop("Y, X and newX must not be specified in cond.loss.SL.control")
    }
    if(!("family" %in% names(cond.loss.SL.control))){
        cond.loss.SL.control<-c(list(family=gaussian()),cond.loss.SL.control)
    }
    
    if(is.null(losses)){
        assert_that(is.data.frame(Y))
        assert_that(nrow(Y)==n)
        assert_that(!anyNA(Y[A==0,]))
        
        if(is.function(learners)){
            learners<-list(learners)
        }else if(!is.list(learners)){
            stop("learners must be a function or a list of functions")
        }else if(!all(sapply(learners,is.function))){
            stop("Some elements of learners are not functions")
        }
        assert_that(is.function(loss.fun))
        nlearners<-length(learners)
    }else{
        assert_that(is.matrix(losses),nrow(losses)==n,!anyNA(losses[A==0,]))
        nlearners<-ncol(losses)
    }
    
    folds<-create.folds(A,nfolds)
    fold.sizes<-sapply(folds,length)
    
    risk.est<-SE<-numeric(nlearners)
    empirical.IF<-matrix(nrow=n,ncol=nlearners)
    for(learner.index in 1:nlearners){
        if(is.null(losses)){
            learner<-learners[[learner.index]]
            loss<-as.numeric(loss.fun(X,Y,learner))
        }else{
            loss<-losses[,learner.index]
        }
        cond.loss<-CV.target.est(A,X,as.matrix(loss),folds,cond.loss.SL.control)
        
        risk.v<-sigma2.v<-matrix(nrow=nfolds,ncol=1)
        IF.v<-list()
        for(v in 1:nfolds){
            # A.train<-A[-folds[[v]]]
            A.test<-A[folds[[v]]]
            loss.test<-loss[folds[[v]]]
            cond.loss.test<-cond.loss[[v]][,1]
            # pi.train<-1-mean(A.train)
            pi.test<-1-mean(A.test)
            
            risk.v[v]<-mean(ifelse(A.test==0,(loss.test-cond.loss.test)/pi.test,0)+cond.loss.test)
            IF<-ifelse(A.test==0,(loss.test-cond.loss.test)/pi.test,0)+cond.loss.test-risk.v[v]+(mean(cond.loss.test)-risk.v[v])/pi.test*(1-A.test-pi.test)
            sigma2.v[v]<-mean(IF^2)
            IF.v[[v]]<-data.frame(id=folds[[v]],IF=IF)
        }
        
        risk.est[learner.index]<-t(fold.sizes)%*%risk.v/n
        sigma2<-t(fold.sizes)%*%sigma2.v/n
        SE[learner.index]<-sqrt(sigma2/n)
        
        all.IF<-do.call(rbind,IF.v)
        all.IF<-all.IF[order(all.IF$id),]
        empirical.IF[,learner.index]<-all.IF$IF
    }
    
    z<-qnorm((conf.level+1)/2)
    estimates<-data.frame(risk.est=risk.est,SE=SE,CI.lower=risk.est-z*SE,CI.upper=risk.est+z*SE)
    list(estimates=estimates,empirical.IF=empirical.IF)
}








#' @title Linear MSE estimation and inference under concept shift in feature/covariate
#' @name linMSE_feat_conshift
#' @export
#' @description Using data from both source population (outcome may be missing) and target population (fully observed), estimate the linear MSE of a list of learners in the target population under concept shift (in feature/covariate) assumption (that is, feature/covariate distributions are identical in both populations but the distributions of label/outcome given feature/covariate might differ). The estimator is semiparametrically efficient when the nuisance function is consistently estimated and is asymptotically normal even if the nuisance function is inconsistently estimated. This function is more computationally efficient than using the more general function \code{\link{feat_conshift}} when multiple learners are of interest because the structure of the loss function is exploited.
#' @param learners a list of learners (functions) that can take covariate `X` as input; can also be a function if only one learner is of interest
#' @param A a vector of binary indicator of data point coming from source population (encoded by `1`) or target population (encoded by `0`)
#' @param X data frame of features/covariates
#' @param Y data frame of labels/outcomes with one column
#' @param nfolds number of folds used for the cross-fit estimator; default to 5
#' @param cond.Y.SL.control a named list containing options passed to \code{\link[SuperLearner:SuperLearner]{SuperLearner::SuperLearner}} to estimate the conditional mean of label/outcome `Y` given feature/covariate `X` (a nuisance function). Must not specify `Y`, `X` or `newX`. Default to `list(family=gaussian(),SL.library=c("SL.glm","SL.gam","SL.randomForest"))`
#' @param conf.level confidence level in the output (asymptotic) confidence interval
#' @param Yhat alternative input in place of `learners`; a matrix of predicted P(Y=1) at each data point and for each learner (row index is for data point; column index is for learner); default to `NULL`
#' @return a list with the following elements:
#' * `estimates`: a data frame with each row being the output for one learner in `learners`
#' * `empirical.IF`: a matrix of influence functions evaluated at empirical data poins. Each column corresponds to a learner in `learners` and each row corresponds to a data point
#' @section Using `empirical.IF` to infer contrasts of risks:
#' See \link{covshift}
linMSE_feat_conshift<-function(learners,A,X,Y,nfolds=5,
                          cond.Y.SL.control=list(family=gaussian(),SL.library=c("SL.glm","SL.gam","SL.randomForest")),
                          conf.level=.95,
                          Yhat=NULL){
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
    assert_that(!anyNA(Y.vec[A==0]))
    
    if(any(c("Y","X","newX") %in% names(cond.Y.SL.control))){
        stop("Y, X and newX must not be specified in cond.Y.SL.control")
    }
    if(!("family" %in% names(cond.Y.SL.control))){
        cond.Y.SL.control<-c(list(family=gaussian()),cond.Y.SL.control)
    }
    
    if(is.null(Yhat)){
        if(is.function(learners)){
            learners<-list(learners)
        }else if(!is.list(learners)){
            stop("learners must be a function or a list of functions")
        }else if(!all(sapply(learners,is.function))){
            stop("Some elements of learners are not functions")
        }
        nlearners<-length(learners)
    }else{
        assert_that(is.matrix(Yhat),nrow(Yhat)==n,!anyNA(Yhat[A==0,]))
        nlearners<-ncol(Yhat)
    }
    
    folds<-create.folds(A,nfolds)
    fold.sizes<-sapply(folds,length)
    
    cond.Y<-CV.target.est(A,X,as.matrix(Y.vec),folds,cond.Y.SL.control)
    if(is.null(Yhat)){
        fx<-do.call(cbind,lapply(learners,function(learner){
            learner(X)
        }))
    }else{
        fx<-Yhat
    }
    loss<-fx*(-2*matrix(Y.vec,nrow=n,ncol=nlearners,byrow=FALSE)+fx)
    
    risk.est<-SE<-numeric(nlearners)
    empirical.IF<-matrix(nrow=n,ncol=nlearners)
    for(learner.index in 1:nlearners){
        risk.v<-sigma2.v<-matrix(nrow=nfolds,ncol=1)
        IF.v<-list()
        for(v in 1:nfolds){
            # A.train<-A[-folds[[v]]]
            A.test<-A[folds[[v]]]
            cond.Y.test<-cond.Y[[v]][,1]
            # pi.train<-1-mean(A.train)
            pi.test<-1-mean(A.test)
            fx.test<-fx[folds[[v]],learner.index]
            loss.test<-loss[folds[[v]],learner.index]
            cond.loss.test<-fx.test*(-2*cond.Y.test+fx.test)
            
            risk.v[v]<-mean(ifelse(A.test==0,(loss.test-cond.loss.test)/pi.test,0)+cond.loss.test)
            IF<-ifelse(A.test==0,(loss.test-cond.loss.test)/pi.test,0)+cond.loss.test-risk.v[v]+(mean(cond.loss.test)-risk.v[v])/pi.test*(1-A.test-pi.test)
            sigma2.v[v]<-mean(IF^2)
            IF.v[[v]]<-data.frame(id=folds[[v]],IF=IF)
        }
        
        risk.est[learner.index]<-t(fold.sizes)%*%risk.v/n
        sigma2<-t(fold.sizes)%*%sigma2.v/n
        SE[learner.index]<-sqrt(sigma2/n)
        
        all.IF<-do.call(rbind,IF.v)
        all.IF<-all.IF[order(all.IF$id),]
        empirical.IF[,learner.index]<-all.IF$IF
    }
    
    z<-qnorm((conf.level+1)/2)
    estimates<-data.frame(risk.est=risk.est,SE=SE,CI.lower=risk.est-z*SE,CI.upper=risk.est+z*SE)
    list(estimates=estimates,empirical.IF=empirical.IF)
}








#' @title Cross-entropy risk estimation and inference under concept shift in feature/covariate
#' @name cross_entropy_feat_conshift
#' @export
#' @description Using data from both source population (outcome may be missing) and target population (fully observed), estimate the cross-entropy risk of a list of learners in the target population under concept shift (in feature/covariate) assumption (that is, feature/covariate distributions are identical in both populations but the distributions of label/outcome given feature/covariate might differ). The estimator is semiparametrically efficient when the nuisance function is consistently estimated and is asymptotically normal even if the nuisance function is inconsistently estimated. This function is more computationally efficient than using the more general function \code{\link{feat_conshift}} when multiple learners are of interest because the structure of the loss function is exploited.
#' @param learners a list of learners (functions) that can take covariate `X` as input; can also be a function if only one learner is of interest
#' @param A a vector of binary indicator of data point coming from source population (encoded by `1`) or target population (encoded by `0`)
#' @param X data frame of features/covariates
#' @param Y data frame of binary labels/outcomes (0 or 1) with one column
#' @param nfolds number of folds used for the cross-fit estimator; default to 5
#' @param cond.Y.SL.control a named list containing options passed to \code{\link[SuperLearner:SuperLearner]{SuperLearner::SuperLearner}} to estimate the conditional mean of label/outcome `Y` given feature/covariate `X` (a nuisance function). Must not specify `Y`, `X` or `newX`. Default to `list(family=binomial(),SL.library=c("SL.glm","SL.gam","SL.randomForest"))`
#' @param conf.level confidence level in the output (asymptotic) confidence interval
#' @param Yhat alternative input in place of `learners`; a matrix of predicted P(Y=1) at each data point and for each learner (row index is for data point; column index is for learner); default to `NULL`
#' @return a list with the following elements:
#' * `estimates`: a data frame with each row being the output for one learner in `learners`
#' * `empirical.IF`: a matrix of influence functions evaluated at empirical data poins. Each column corresponds to a learner in `learners` and each row corresponds to a data point
#' @section Using `empirical.IF` to infer contrasts of risks:
#' See \link{covshift}
cross_entropy_feat_conshift<-function(learners,A,X,Y,nfolds=5,
                                      cond.Y.SL.control=list(family=binomial(),SL.library=c("SL.glm","SL.gam","SL.randomForest")),
                                      conf.level=.95,
                                      Yhat=NULL){
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
    
    assert_that(ncol(Y)==1)
    Y.vec<-Y[,1,drop=TRUE]
    if(any(!(Y.vec[A==0] %in% c(0,1)))){
        stop("Y must take values 0 or 1")
    }
    
    assert_that(!anyNA(A))
    assert_that(!anyNA(X))
    assert_that(!anyNA(Y.vec[A==0]))
    
    if(any(c("Y","X","newX") %in% names(cond.Y.SL.control))){
        stop("Y, X and newX must not be specified in cond.Y.SL.control")
    }
    if(!("family" %in% names(cond.Y.SL.control))){
        cond.Y.SL.control<-c(list(family=gaussian()),cond.Y.SL.control)
    }
    
    if(is.null(Yhat)){
        if(is.function(learners)){
            learners<-list(learners)
        }else if(!is.list(learners)){
            stop("learners must be a function or a list of functions")
        }else if(!all(sapply(learners,is.function))){
            stop("Some elements of learners are not functions")
        }
        nlearners<-length(learners)
    }else{
        assert_that(is.matrix(Yhat),nrow(Yhat)==n,!anyNA(Yhat[A==0,]))
        nlearners<-ncol(Yhat)
    }
    
    folds<-create.folds(A,nfolds)
    fold.sizes<-sapply(folds,length)
    
    cond.Y<-CV.target.est(A,X,as.matrix(Y.vec),folds,cond.Y.SL.control)
    if(is.null(Yhat)){
        fx<-do.call(cbind,lapply(learners,function(learner){
            learner(X)
        }))
    }else{
        fx<-Yhat
    }
    loss<--matrix(Y.vec,nrow=n,ncol=nlearners,byrow=FALSE)*logit(fx)-log(1-fx)
    
    risk.est<-SE<-numeric(nlearners)
    empirical.IF<-matrix(nrow=n,ncol=nlearners)
    for(learner.index in 1:nlearners){
        risk.v<-sigma2.v<-matrix(nrow=nfolds,ncol=1)
        IF.v<-list()
        for(v in 1:nfolds){
            # A.train<-A[-folds[[v]]]
            A.test<-A[folds[[v]]]
            cond.Y.test<-cond.Y[[v]][,1]
            # pi.train<-1-mean(A.train)
            pi.test<-1-mean(A.test)
            fx.test<-fx[folds[[v]],learner.index]
            loss.test<-loss[folds[[v]],learner.index]
            cond.loss.test<-fx.test*(-cond.Y.test*logit(fx.test)-log(1-fx.test))
            
            risk.v[v]<-mean(ifelse(A.test==0,(loss.test-cond.loss.test)/pi.test,0)+cond.loss.test)
            IF<-ifelse(A.test==0,(loss.test-cond.loss.test)/pi.test,0)+cond.loss.test-risk.v[v]+(mean(cond.loss.test)-risk.v[v])/pi.test*(1-A.test-pi.test)
            sigma2.v[v]<-mean(IF^2)
            IF.v[[v]]<-data.frame(id=folds[[v]],IF=IF)
        }
        
        risk.est[learner.index]<-t(fold.sizes)%*%risk.v/n
        sigma2<-t(fold.sizes)%*%sigma2.v/n
        SE[learner.index]<-sqrt(sigma2/n)
        
        all.IF<-do.call(rbind,IF.v)
        all.IF<-all.IF[order(all.IF$id),]
        empirical.IF[,learner.index]<-all.IF$IF
    }
    
    z<-qnorm((conf.level+1)/2)
    estimates<-data.frame(risk.est=risk.est,SE=SE,CI.lower=risk.est-z*SE,CI.upper=risk.est+z*SE)
    list(estimates=estimates,empirical.IF=empirical.IF)
}







feat_conshift_naive<-function(learners,loss.fun,A,X,Y,nfolds=5,
                        cond.loss.SL.control=list(family=gaussian(),SL.library=c("SL.glm","SL.gam","SL.randomForest")),
                        conf.level=.95,
                        losses=NULL){
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
    assert_that(is.number(conf.level) && conf.level>=0.5 && conf.level<1)
    
    assert_that(!anyNA(A))
    assert_that(!anyNA(X))
    
    if(any(c("Y","X","newX") %in% names(cond.loss.SL.control))){
        stop("Y, X and newX must not be specified in cond.loss.SL.control")
    }
    if(!("family" %in% names(cond.loss.SL.control))){
        cond.loss.SL.control<-c(list(family=gaussian()),cond.loss.SL.control)
    }
    
    if(is.null(losses)){
        assert_that(is.data.frame(Y))
        assert_that(nrow(Y)==n)
        assert_that(!anyNA(Y[A==0,]))
        
        if(is.function(learners)){
            learners<-list(learners)
        }else if(!is.list(learners)){
            stop("learners must be a function or a list of functions")
        }else if(!all(sapply(learners,is.function))){
            stop("Some elements of learners are not functions")
        }
        assert_that(is.function(loss.fun))
        nlearners<-length(learners)
    }else{
        assert_that(is.matrix(losses),nrow(losses)==n,!anyNA(losses[A==0,]))
        nlearners<-ncol(losses)
    }
    
    folds<-create.folds(A,nfolds)
    fold.sizes<-sapply(folds,length)
    
    risk.est<-SE<-numeric(nlearners)
    empirical.IF<-matrix(nrow=n,ncol=nlearners)
    for(learner.index in 1:nlearners){
        if(is.null(losses)){
            learner<-learners[[learner.index]]
            loss<-as.numeric(loss.fun(X,Y,learner))
        }else{
            loss<-losses[,learner.index]
        }
        cond.loss<-CV.target.est(A,X,as.matrix(loss),folds,cond.loss.SL.control)
        
        risk.v<-sigma2.v<-matrix(nrow=nfolds,ncol=1)
        IF.v<-list()
        for(v in 1:nfolds){
            # A.train<-A[-folds[[v]]]
            A.test<-A[folds[[v]]]
            loss.test<-loss[folds[[v]]]
            cond.loss.test<-cond.loss[[v]][,1]
            # pi.train<-1-mean(A.train)
            pi.test<-1-mean(A.test)
            
            risk.v[v]<-mean(cond.loss.test)
            IF<-ifelse(A.test==0,(loss.test-cond.loss.test)/pi.test,0)+cond.loss.test-risk.v[v]+(mean(cond.loss.test)-risk.v[v])/pi.test*(1-A.test-pi.test)
            sigma2.v[v]<-mean(IF^2)
            IF.v[[v]]<-data.frame(id=folds[[v]],IF=IF)
        }
        
        risk.est[learner.index]<-t(fold.sizes)%*%risk.v/n
        sigma2<-t(fold.sizes)%*%sigma2.v/n
        SE[learner.index]<-sqrt(sigma2/n)
        
        all.IF<-do.call(rbind,IF.v)
        all.IF<-all.IF[order(all.IF$id),]
        empirical.IF[,learner.index]<-all.IF$IF
    }
    
    z<-qnorm((conf.level+1)/2)
    estimates<-data.frame(risk.est=risk.est,SE=SE,CI.lower=risk.est-z*SE,CI.upper=risk.est+z*SE)
    list(estimates=estimates,empirical.IF=empirical.IF)
}
