##' @name cvtest
##' @description General cross-validatad error tester
##' @param data is a data frame containing predictors and response. We presume all predictors are used and data has dimensions n x m.
##' @param response is a response vector of length n
##' @param proc is a functional on data, Y which returns a function of a dataframe analogous to data 
##' @param fold is the number of folds, default 10
##' @param seed random seed
##' @param ... passed to proc
##' @returns list of predicted responses, folds, and root mean square error
cvtest=function(data,response,proc,nfold=10,seed=Sys.time(),...) {
  
	n=dim(data)[1]; m=dim(data)[2]
	set.seed(seed) # replicably randomise fold choice
	fold=rep(1:nfold,1+floor(n/nfold))[1:n][order(runif(n))]
	predicted=rep(0,n)
	for (i in 1:nfold) {
		foldi=which(fold==i); foldx=which(fold != i)
		model=proc(data[foldx,],response[foldx],...)
		predicted_response=model(data[foldi,])
		predicted[foldi]=predicted_response
	}
	return(list(predicted=predicted,fold=fold,rms=sqrt(mean((predicted-response)^2,na.rm=T))))
}



##' @name findbestmodel
##' @description Finds the best model to use on a given dataset using cross-validation
##' @param X matrix of predictors
##' @param Y vector of responses
##' @param fold fold assignments
##' @param models vector of names of functions; each function fits a predictive model to training data (Xtest,Ytest) using cross-validation to tune hyperparameters if required, and returns a function applicable to new data X;
##' @param comp function comparing Y to a prediction for Y; assume lower is better
##' @return function corresponding to best model
findbestmodel=function(X,Y,fold,models,comp) {
nfold=length(unique(fold))
cvtop=rep(0,nfold)
predictions=matrix(0,length(Y),length(models))
for (i in 1:nfold) {
	inf=which(fold==i); outf=which(fold!=i)
	Xtest=X[inf,]; Ytest=Y[inf]
	Xtrain=X[outf,]; Ytrain=Y[outf]
  y_error=rep(0,length(models))	
	for (j in 1:length(models)) {
		Yj=get(models[j])(Xtrain,Ytrain)(Xtest)
		predictions[inf,j]=pred
		y_error[j]=get(comp)(Yj,Ytest)
	}
  cvtop[i]=which.min(y_error)
}
# mode?!
uniqv=unique(cvtop)
im=uniqv[which.max(tabulate(match(cvtop, uniqv)))]
return(list(best=models[im],predictions=predictions))
}



##' @name findbestmodel_surv
##' @description Finds the best model to use on a given dataset using cross-validation, for survival data
##' @param X matrix of predictors
##' @param Y vector of responses
##' @param Yt survival times
##' @param fold fold assignments
##' @param models vector of names of functions; each function fits a predictive model to training data (Xtest,Ytest) using cross-validation to tune hyperparameters if required, and returns a function applicable to new data X;
##' @param comp function comparing Y to a prediction for Y; assume lower is better
##' @return function corresponding to best model
findbestmodel_surv=function(X,Y,Yt,fold,models,comp) {
nfold=length(unique(fold))
cvtop=rep(0,nfold)
predictions=matrix(0,length(Y),length(models))
for (i in 1:nfold) {
	inf=which(fold==i); outf=which(fold!=i)
	Xtest=X[inf,]; Ytest=Y[inf]; Ytestt=Yt[inf]
	Xtrain=X[outf,]; Ytrain=Y[outf]; Ytraint=Yt[outf]
  y_error=rep(0,length(models))	
	for (j in 1:length(models)) {
		pred=get(models[j])(Xtrain,Ytrain,Ytraint)(Xtest)
		predictions[inf,j]=pred
		y_error[j]=get(comp)(pred,Ytest,Ytestt)
	}
  cvtop[i]=which.min(y_error)
}
# mode?!
uniqv=unique(cvtop)
im=uniqv[which.max(tabulate(match(cvtop, uniqv)))]
return(list(best=models[im],predictions=predictions))
}





##' @name evaluateperformance
##' @description evaluate performance of a model using cross-validation
##' @param X matrix of predictors
##' @param Y vector of responses
##' @param model name of model to use; assumed that the model
##' @param model name of model to use; we assume that this function fits a predictive model to training data (Xtest,Ytest) using given parameters and returns a function applicable to new data X;
##' @param pars parameters passed to model
##' @return cross-validated outcome vector of same length as Y
evaluateperformance=function(X,Y,fold=NULL,model,pars) {

if (is.null(fold)) fold=rep(1:10,length(Y))[1:length(Y)][order(runif(length(Y)))]
nfold=length(unique(fold))
Yout=rep(0,nfold)
for (i in 1:nfold) {
	inf=which(fold==i); outf=which(fold!=i)
	Xtest=X[inf,]; Ytest=Y[inf]
	Xtrain=X[outf,]; Ytrain=Y[outf]
  Ypred=get(model)(Xtrain,Ytrain,optimise=F,pars=pars)(Xtest)
	Yout[inf]=Ypred
}

return(Yout)
}


##' @name evaluateperformance_surv
##' @description evaluate performance of a model using cross-validation, on survival data
##' @param X matrix of predictors
##' @param Y vector of responses
##' @param Yt vector of survival times
##' @param model name of model to use; assumed that the model
##' @param model name of model to use; we assume that this function fits a predictive model to training data (Xtest,Ytest) using given parameters and returns a function applicable to new data X;
##' @param pars parameters passed to model
##' @return cross-validated outcome vector of same length as Y
evaluateperformance_surv=function(X,Y,Yt,fold=NULL,model,pars) {
if (is.null(fold)) fold=rep(1:10,length(Y))[1:length(Y)][order(runif(length(Y)))]

nfold=length(unique(fold))
Yout=rep(0,nfold)
for (i in 1:nfold) {
	inf=which(fold==i); outf=which(fold!=i)
	Xtest=X[inf,]; Ytest=Y[inf]; Ytestt=Yt[inf]
	Xtrain=X[outf,]; Ytrain=Y[outf]; Ytraint=Yt[outf]
  Ypred=get(model)(Xtrain,Ytrain,Ytraint,optimise=F,pars=pars)(Xtest)
	Yout[inf]=Ypred
}

return(Yout)
}




##' @name lrproc
##' @description non-penalised logistic regression function for cvtest
##' @param data is a data frame containing predictors and response. We presume all predictors are used and data has dimensions n x m.
##' @param Y outcomes for training data
##' @param optimise set to TRUE to optimise hyperparameters by cross-validation. No hyperparameters in this case; ignored
##' @param pars list of hyperparameters to consider. Ignored in this case.
##' @param logit return probabilities of Y=1 rather than logit-transformed
##' @param return_hyperparameters set to TRUE to return optimal hyperparameters instead of function
##' @param ... passed to glm
##' @returns either function mapping data to predicted values, or best hyperparameters
lrproc=function(data,Y,optimise=TRUE, pars=NULL, return_hyperparameters=FALSE,logit=FALSE,...) {
  if (length(unique(Y))==2)
  	gx=glm(Y~.,data=data,family = binomial(link = "logit"),...)
  else 
  	gx=glm(Y~.,data=data,family = "gaussian",...)
	cx=gx$coefficients; cxi=cx[1]; cxbeta=cx[2:length(cx)]; w=which(is.finite(cxbeta))
	if (return_hyperparameters) return(1) else 
		if (logit) return(function(X) logit(cxi+ (as.matrix(X[,w]) %*% cxbeta[w]))) else
			return(function(X) cxi+ (as.matrix(X[,w]) %*% cxbeta[w]))
}


##' @name lrproc_surv
##' @description non-penalised logistic regression function for cvtest where outcome is a survival time
##' @param data is a data frame containing predictors and response. We presume all predictors are used and data has dimensions n x m.
##' @param Ys event or censorship for training data
##' @param Yt time-to-event or time-to-censorship for training data
##' @param optimise set to TRUE to optimise hyperparameters by cross-validation. No hyperparameters in this case; ignored
##' @param pars list of hyperparameters to consider. Ignored in this case.
##' @param exp if TRUE, return hazard ratios, otherwise return log hazard ratio.
##' @param return_hyperparameters set to TRUE to return optimal hyperparameters instead of function
##' @param ... passed to glm
##' @returns either function mapping data to predicted values, or best hyperparameters
lrproc_surv=function(data,Ys,Yt,optimise=TRUE, pars=NULL, return_hyperparameters=FALSE,exp=FALSE,logit=F,...) {
	gx=coxph(Surv(time=Yt,event=Ys)~.,data=data,...)
	if (return_hyperparameters) return(1) else 
		if (logit) return(function(X) predict(gx,newdata=X,type="risk")) else 
			return(function(X) predict(gx,newdata=X,type="lp"))
}





##' @name lassoproc
##' @description L1-penalised logistic regression function for cvtest
##' @param data is a data frame containing predictors and response. We presume all predictors are used and data has dimensions n x m.
##' @param Y outcomes for training data
##' @param optimise set to TRUE to optimise hyperparameters (lambda) by cross-validation. 
##' @param pars list of hyperparameters to consider.
##' @param logit return probabilities of Y=1 rather than logit-transformed
##' @param return_hyperparameters set to TRUE to return optimal hyperparameters instead of function
##' @returns either function mapping data to predicted values, or best hyperparameters
lassoproc=function(data,Y,optimise=T,pars=NULL,logit=FALSE,return_hyperparameters=FALSE) {
	if (length(unique(Y))==2) famx="binomial" else famx="gaussian"
	if (optimise) {
		cxx=cv.glmnet(as.matrix(data),Y,family = famx,lambda=pars)
  	cx=glmnet(as.matrix(data),Y,family=famx,lambda=cxx$lambda.min); 
	} else cx=glmnet(as.matrix(data),Y,family=famx,lambda=pars); 
	cxi=cx$a0; cxbeta=as.numeric(cx$beta); w=which(is.finite(cxbeta))
	if (return_hyperparameters) return(cxx$lambda.min) else 
		if (logit) return(function(X) logit(cxi+ (as.matrix(X[,w]) %*% cxbeta[w]))) else 
			return(function(X) cxi+ (as.matrix(X[,w]) %*% cxbeta[w]))
}



##' @name lassoproc_surv
##' @description L1-penalised logistic regression function with survival time
##' @param data is a data frame containing predictors and response. We presume all predictors are used and data has dimensions n x m.
##' @param Y outcomes for training data
##' @param Yt survival times
##' @param optimise set to TRUE to optimise hyperparameters (lambda) by cross-validation. 
##' @param pars list of hyperparameters to consider.
##' @param logit return probabilities of Y=1 rather than logit-transformed
##' @param return_hyperparameters set to TRUE to return optimal hyperparameters instead of function
##' @returns either function mapping data to predicted values, or best hyperparameters
lassoproc_surv=function(data,Y,Yt,optimise=T,pars=NULL,logit=FALSE,return_hyperparameters=FALSE) {
	if (optimise) {
		cxx=cv.glmnet(as.matrix(data),Surv(time=Yt,event=Y),family="cox")
  	cx=glmnet(as.matrix(data),Surv(time=Yt,event=Y),family="cox",lambda=cxx$lambda.min)
	} else cx=glmnet(as.matrix(data),Surv(time=Yt,event=Y),family="cox",lambda=pars)
	cxi=0; cxbeta=as.numeric(cx$beta); w=which(is.finite(cxbeta))
	if (return_hyperparameters) return(cxx$lambda.min) else 
		if (logit) return(function(X) logit(cxi+ (as.matrix(X[,w]) %*% cxbeta[w]))) else 
			return(function(X) cxi+ (as.matrix(X[,w]) %*% cxbeta[w]))
}




##' @name rfproc
##' @description random forest classifier for cvtest
##' @param data is a data frame containing predictors and response. We presume all predictors are used and data has dimensions n x m.
##' @param Y outcome
##' @param return_hyperparameters set to TRUE to return optimal hyperparameters instead of function
##' @param optimise set to TRUE to optimise hyperparameters by cross-validation. No hyperparameters in this case; ignored
##' @param pars list of hyperparameters to consider. Ignored in this case.
##' @param ... passed to glm
##' @returns either function mapping data to predicted values, or best hyperparameters
rfproc=function(data,Y,optimise=T,pars=NULL,return_hyperparameters=FALSE,...) {
	cp=randomForest(Y~.,data=data)
	if (return_hyperparameters) return(1) else return(function(X) predict(cp,X))
}


##' @name rfproc_surv
##' @description random forest classifier for cvtest, working on survival data
##' @param data is a data frame containing predictors and response. We presume all predictors are used and data has dimensions n x m.
##' @param Y outcomes for training data
##' @param Yt survival time
##' @param return_hyperparameters set to TRUE to return optimal hyperparameters instead of function
##' @param optimise set to TRUE to optimise hyperparameters by cross-validation. No hyperparameters in this case; ignored
##' @param pars list of hyperparameters to consider. Ignored in this case.
##' @param ... passed to glm
##' @returns either function mapping data to predicted values, or best hyperparameters
rfproc_surv=function(data,Y,Yt,optimise=T,pars=NULL,return_hyperparameters=FALSE,...) {
	nobs=dim(data)[1]
	X=as.matrix(data)
	cxx=rfsrc(Surv(time = Yt, event=Y)~.,data=cbind(data,Y,Yt), block.size = 1)
	if (return_hyperparameters) return(1) else return(function(X) predict(cxx,newdata=X)$predicted)
}


##' @name auroc
##' @description finds best fit between a predictor for Y and Y based on area under the ROC. 
##' @param predictor predicted Y
##' @param Y vector of responses
##' @return area under ROC curve
auroc=function(predictor,Y) {
	n=length(predictor); sp=sort(predictor)
	sens=rep(0,n); spec=rep(0,n); 
	n0=length(which(Y==0)); n1=length(which(Y==1))
	for (i in 1:n) {
		sens[i]=length(which(predictor > sp[i] & Y==1))/n1
		spec[i]=length(which(predictor <= sp[i] & Y==0))/n0
	}
	return(-trapz(1-spec,sens))
}

##' @name draw_roc
##' @description draws an ROC curve
##' @param predictor predictor of response
##' @param Y response
##' @param ... passed to plot
##' @returns plot
draw_roc=function(predictor,Y,add=F,...) {
	n=length(predictor); sp=sort(predictor)
	sens=rep(0,n); spec=rep(0,n); 
	n0=length(which(Y==0)); n1=length(which(Y==1))
	for (i in 1:n) {
		sens[i]=length(which(predictor > sp[i] & Y==1))/n1
		spec[i]=length(which(predictor <= sp[i] & Y==0))/n0
	}
	if (!add) {
		plot(1-spec,sens,xlim=c(0,1),ylim=c(0,1),type="l",xlab="1-Specificity",ylab="Sensitivity",xaxs="i",yaxs="i",...); 
		abline(0,1,col="red")
	} else {
		lines(1-spec,sens,...); 
	}
	return(-trapz(1-spec,sens))
}

# Trapezoid rule
trapz=function(x,y) {
	idx = 2:length(x) 
	return (as.double( (x[idx] - x[idx-1]) %*% (y[idx] + y[idx-1])) / 2)
}

# Logit and inverse logit
logit=function(x) 1/(1+exp(-x))
ilogit=function(y) -log((1/y)-1)
	
# Format output of ci.cvAUC
formatci=function(cc,digits=3) 
 paste0(signif(cc$cvAUC,digits=digits)," (",
 signif(cc$ci[1],digits=digits),",",
 signif(cc$ci[2],digits=digits),")")




##' @name draw_roc_surv
##' @description draws an ROC curve at a given time point for survival data, averaging over cross-validation folds
##' @param yp predictor of response
##' @param Y response event
##' @param Yt response time
##' @param fold vector of cross-validation fold assignments
##' @param time timepoint to draw curve at
##' @param add set to T to overplot
##' @param ... passed to plot
##' @returns plot
draw_roc_surv=function(yp,Y,Yt,fold,time,add=T,...) {

mx=sort(unique(yp))
nfold=length(unique(fold))

fp=matrix(0,length(mx),nfold)
tp=fp

if (!add) plot(0,type="n",xlim=c(0,1),ylim=c(0,1),xaxs="i",yaxs="i");
for (i in 1:nfold) {
w=which(fold==i)
r1=risksetROC(Stime=Yt[w],status=Y[w],marker=yp[w],predict.time=time,plot=F)
fp[,i]=approx(r1$marker,r1$FP[2:(length(r1$FP)-1)],mx,rule=2)$y
tp[,i]=approx(r1$marker,r1$TP[2:(length(r1$FP)-1)],mx,rule=2)$y
}
lines(rowMeans(fp),rowMeans(tp),...)
if (!add) abline(0,1,col="red")
}



##' @name draw_roc_cv
##' @description draws an ROC curve averaging over cross-validation folds
##' @param yp predictor of response
##' @param Y response event
##' @param fold vector of cross-validation fold assignments
##' @param add set to T to overplot
##' @param ... passed to plot
##' @returns plot
draw_roc_cv=function(yp,Y,fold,add=F,...) {
n=length(yp); sp=sort(yp)
sens=rep(0,n); spec=rep(0,n);
for (i in 1:n) {
 sens[i]=mean(sapply(1:max(fold),function(x) 
    length(which(yp[which(fold==x)]>sp[i] & Y[which(fold==x)]==1))/length(which(Y[which(fold==x)]==1))))
 spec[i]=mean(sapply(1:max(fold),function(x) 
    length(which(yp[which(fold==x)]<=sp[i] & Y[which(fold==x)]==0))/length(which(Y[which(fold==x)]==0))))
} 
if (!add) {
  plot(1-spec,sens,xlim=c(0,1),ylim=c(0,1),type="l",xlab="1-Specificity",ylab="Sensitivity",xaxs="i",yaxs="i",...); 
  abline(0,1,col="red")
} else {
 lines(1-spec,sens,...)
}
}



##' @name draw_calibration
##' @description draws a calibration curve
##' @param Yp probabilistic predictor of response
##' @param Y response event
##' @param nc split predictor into nc-iles
##' @param ... passed to plot
##' @returns plot
draw_calibration=function(Yp,Y,nc=10,add=T,...) {
cal=rep(0,nc)
for (i in 1:nc) {
w=which(Yp>= (i-1)/nc  & Yp<i/nc)
cal[i]=mean(Y[w])
}
if (add) points( (1:nc)/nc - (1/(2*nc)),cal,...) else plot((1:nc)/(1+nc),cal,...)
}



