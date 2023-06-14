#A is a binary vector
#X is a data frame of covariates
#Y is a data frame of outcomes/labels
np<-function(learners,loss.fun,A,X,Y,conf.level=.95,losses=NULL){
    #check inputs
    if(any(!(A %in% c(0,1)))){
        stop("A must take values 0 or 1")
    }
    if(all(A==0) || all(A==1)){
        stop("A must contain both values 0 and 1")
    }
    n<-length(A)
    
    
    assert_that(is.number(conf.level) && conf.level>=0.5 && conf.level<1)
    
    assert_that(!anyNA(A))
    
    if(is.null(losses)){
        assert_that(is.data.frame(X))
        assert_that(nrow(X)==n)
        assert_that(is.data.frame(Y))
        assert_that(nrow(Y)==n)
        assert_that(!anyNA(X[A==0,]))
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
    
    risk.est<-SE<-numeric(nlearners)
    empirical.IF<-matrix(nrow=n,ncol=nlearners)
    for(learner.index in 1:nlearners){
        if(is.null(losses)){
            learner<-learners[[learner.index]]
            loss<-as.numeric(loss.fun(X,Y,learner))
        }else{
            loss<-losses[,learner.index]
        }
        
        risk.est[learner.index]<-mean(loss[A==0])
        empirical.IF[,learner.index]<-ifelse(A==0,loss-risk.est[learner.index],0)/mean(A==0)
        SE[learner.index]<-sqrt(mean(empirical.IF[,learner.index]^2)/n)
    }
    
    z<-qnorm((conf.level+1)/2)
    estimates<-data.frame(risk.est=risk.est,SE=SE,CI.lower=risk.est-z*SE,CI.upper=risk.est+z*SE)
    list(estimates=estimates,empirical.IF=empirical.IF)
}









#A is a binary vector
#X is a data frame of covariates
#Y is a data frame of outcomes/labels
linMSE_np<-function(learners,A,X,Y,conf.level=.95,Yhat=NULL){
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
    assert_that(!anyNA(X[A==0,]))
    assert_that(!anyNA(Y[A==0,]))
    
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
    
    if(is.null(Yhat)){
        fx<-do.call(cbind,lapply(learners,function(learner){
            learner(X)
        }))
    }
    else{
        fx<-Yhat
    }
    loss<-fx*(-2*matrix(Y.vec,nrow=n,ncol=nlearners,byrow=FALSE)+fx)
    
    risk.est<-SE<-numeric(nlearners)
    empirical.IF<-matrix(nrow=n,ncol=nlearners)
    for(learner.index in 1:nlearners){
        risk.est[learner.index]<-mean(loss[A==0,learner.index])
        empirical.IF[,learner.index]<-ifelse(A==0,loss-risk.est[learner.index],0)/mean(A==0)
        SE[learner.index]<-sqrt(mean(empirical.IF[,learner.index]^2)/n)
    }
    
    z<-qnorm((conf.level+1)/2)
    estimates<-data.frame(risk.est=risk.est,SE=SE,CI.lower=risk.est-z*SE,CI.upper=risk.est+z*SE)
    list(estimates=estimates,empirical.IF=empirical.IF)
}




#A is a binary vector
#X is a data frame of covariates
#Y is a data frame of outcomes/labels with one column
cross_entropy_np<-function(learners,A,X,Y,conf.level=.95,Yhat=NULL){
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
    assert_that(!anyNA(X[A==0,]))
    assert_that(!anyNA(Y[A==0,]))
    
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
        risk.est[learner.index]<-mean(loss[A==0,learner.index])
        empirical.IF[,learner.index]<-ifelse(A==0,loss-risk.est[learner.index],0)/mean(A==0)
        SE[learner.index]<-sqrt(mean(empirical.IF[,learner.index]^2)/n)
    }
    
    z<-qnorm((conf.level+1)/2)
    estimates<-data.frame(risk.est=risk.est,SE=SE,CI.lower=risk.est-z*SE,CI.upper=risk.est+z*SE)
    list(estimates=estimates,empirical.IF=empirical.IF)
}
