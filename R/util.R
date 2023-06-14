logit<-binomial()$linkfun
expit<-binomial()$linkinv

#create k folds of a vector id with A distributed roughly evenly in each fold
create.folds<-function(A,k){
    if(sum(A)<k || sum(1-A)<k){
        warning("Too few observations from source or target population for sample splitting! Devision by zero might occur!")
    }
    id<-1:length(A)
    id1<-id[A==1]
    id0<-id[A==0]
    
    order1<-sample.int(length(id1))
    d1<-suppressWarnings(data.frame(cbind(id1[order1],1:k)))
    names(d1)<-c("id","fold.id")
    
    order0<-sample.int(length(id0))
    d0<-suppressWarnings(data.frame(cbind(id0[order0],k:1)))
    names(d0)<-c("id","fold.id")
    
    d<-rbind(d0,d1)
    lapply(tapply(d$id,d$fold.id,identity,simplify=FALSE),sort)
}

#initial SL estimate of nuisance functions
#X is a data frame
#Y is a matrix with each column being an outcome
CV.est<-function(X,Y,folds,SL.control){
    n.Y<-ncol(Y)
    lapply(1:length(folds),function(v){
        X.train<-X[-folds[[v]],,drop=FALSE]
        X.test<-X[folds[[v]],,drop=FALSE]
        
        out<-matrix(nrow=length(folds[[v]]),ncol=n.Y)
        
        for(Y.index in 1:n.Y){
            Y.train<-Y[-folds[[v]],Y.index]
            if(all(Y.train==0)){
                out[,Y.index]<-rep(0,length(folds[[v]]))
            }else if(all(Y.train==1)){
                out[,Y.index]<-rep(1,length(folds[[v]]))
            }else{
                args<-c(list(Y=Y.train,X=X.train,newX=X.test),SL.control)
                SL.model<-do.call(SuperLearner,args)
                out[,Y.index]<-predict(SL.model)$pred
            }
        }
        out
    })
}

#initial SL estimate of nuisance functions using taarget data (A=0) only
#X is a data frame
#Y is a matrix with each column being an outcome
CV.target.est<-function(A,X,Y,folds,SL.control){
    n.Y<-ncol(Y)
    lapply(1:length(folds),function(v){
        A.train<-A[-folds[[v]]]
        X.train<-X[-folds[[v]],,drop=FALSE]
        X.test<-X[folds[[v]],,drop=FALSE]
        
        out<-matrix(nrow=length(folds[[v]]),ncol=n.Y)
        
        for(Y.index in 1:n.Y){
            Y.train<-Y[-folds[[v]],Y.index]
            if(all(Y.train[A.train==0]==0)){
                out[,Y.index]<-rep(0,length(folds[[v]]))
            }else if(all(Y.train[A.train==0]==1)){
                out[,Y.index]<-rep(1,length(folds[[v]]))
            }else{
                args<-c(list(Y=Y.train[A.train==0],X=X.train[A.train==0,,drop=FALSE],newX=X.test),SL.control)
                SL.model<-do.call(SuperLearner,args)
                out[,Y.index]<-predict(SL.model)$pred
            }
        }
        out
    })
}

#find the index i of x such that x[j] (j<=i) are all TRUE
find.max.true<-function(x){
    for(i in 1:length(x)){
        if(!x[i]){
            return(i-1)
        }
    }
    length(x)
}
