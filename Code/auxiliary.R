##' Logit and logistic
logistic=function(x) 1/(1+exp(-x))
logit=function(x) log(x/(1-x))


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
##' @param logistic return probabilities of Y=1 rather than logistic-transformed
##' @param return_hyperparameters set to TRUE to return optimal hyperparameters instead of function
##' @param ... passed to glm
##' @returns either function mapping data to predicted values, or best hyperparameters
lrproc=function(data,Y,optimise=TRUE, pars=NULL, return_hyperparameters=FALSE,logistic=FALSE,...) {
  if (length(unique(Y))==2)
  	gx=glm(Y~.,data=data,family = binomial(link = "logit"),...)
  else 
  	gx=glm(Y~.,data=data,family = "gaussian",...)
	cx=gx$coefficients; cxi=cx[1]; cxbeta=cx[2:length(cx)]; w=which(is.finite(cxbeta))
	if (return_hyperparameters) return(1) else 
		if (logistic) return(function(X) logistic(cxi+ (as.matrix(X[,w]) %*% cxbeta[w]))) else
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
lrproc_surv=function(data,Ys,Yt,optimise=TRUE, pars=NULL, return_hyperparameters=FALSE,exp=FALSE,logistic=F,...) {
	gx=coxph(Surv(time=Yt,event=Ys)~.,data=data,...)
	if (return_hyperparameters) return(1) else 
		if (logistic) return(function(X) predict(gx,newdata=X,type="risk")) else 
			return(function(X) predict(gx,newdata=X,type="lp"))
}





##' @name lassoproc
##' @description L1-penalised logistic regression function for cvtest
##' @param data is a data frame containing predictors and response. We presume all predictors are used and data has dimensions n x m.
##' @param Y outcomes for training data
##' @param optimise set to TRUE to optimise hyperparameters (lambda) by cross-validation. 
##' @param pars list of hyperparameters to consider.
##' @param logistic return probabilities of Y=1 rather than logistic-transformed
##' @param return_hyperparameters set to TRUE to return optimal hyperparameters instead of function
##' @returns either function mapping data to predicted values, or best hyperparameters
lassoproc=function(data,Y,optimise=T,pars=NULL,logistic=FALSE,return_hyperparameters=FALSE) {
	if (length(unique(Y))==2) famx="binomial" else famx="gaussian"
	if (optimise) {
		cxx=cv.glmnet(as.matrix(data),Y,family = famx,lambda=pars)
  	cx=glmnet(as.matrix(data),Y,family=famx,lambda=cxx$lambda.min); 
	} else cx=glmnet(as.matrix(data),Y,family=famx,lambda=pars); 
	cxi=cx$a0; cxbeta=as.numeric(cx$beta); w=which(is.finite(cxbeta))
	if (return_hyperparameters) return(cxx$lambda.min) else 
		if (logistic) return(function(X) logistic(cxi+ (as.matrix(X[,w]) %*% cxbeta[w]))) else 
			return(function(X) cxi+ (as.matrix(X[,w]) %*% cxbeta[w]))
}



##' @name lassoproc_surv
##' @description L1-penalised logistic regression function with survival time
##' @param data is a data frame containing predictors and response. We presume all predictors are used and data has dimensions n x m.
##' @param Y outcomes for training data
##' @param Yt survival times
##' @param optimise set to TRUE to optimise hyperparameters (lambda) by cross-validation. 
##' @param pars list of hyperparameters to consider.
##' @param logistic return probabilities of Y=1 rather than logistic-transformed
##' @param return_hyperparameters set to TRUE to return optimal hyperparameters instead of function
##' @returns either function mapping data to predicted values, or best hyperparameters
lassoproc_surv=function(data,Y,Yt,optimise=T,pars=NULL,logistic=FALSE,return_hyperparameters=FALSE) {
	if (optimise) {
		cxx=cv.glmnet(as.matrix(data),Surv(time=Yt,event=Y),family="cox")
  	cx=glmnet(as.matrix(data),Surv(time=Yt,event=Y),family="cox",lambda=cxx$lambda.min)
	} else cx=glmnet(as.matrix(data),Surv(time=Yt,event=Y),family="cox",lambda=pars)
	cxi=0; cxbeta=as.numeric(cx$beta); w=which(is.finite(cxbeta))
	if (return_hyperparameters) return(cxx$lambda.min) else 
		if (logistic) return(function(X) logistic(cxi+ (as.matrix(X[,w]) %*% cxbeta[w]))) else 
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


##' @name aucse
##' @description standard error of AUROC
##' @param n1 number of 1's/0's
##' @param n2 number of 0's/1's
##' @param auc AUROC
##' @returns estimate of standard error
aucse=function(n1,n2,auc) {
  q1=auc/(2-auc); q2=2*(auc^2)/(1+auc)
  num=auc*(1-auc) + (n1-1)*(q1- (auc^2)) + (n2-1)*(q2-(auc^2))
  return(sqrt(num/(n1*n2)))
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
	sens=c(1,sens,0); spec=c(0,spec,1)
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
draw_roc_surv=function(yp,Y,Yt,fold=NULL,time,add=T,...) {
  
  
  if (!add) plot(0,type="n",xlim=c(0,1),ylim=c(0,1),xaxs="i",yaxs="i");
  if (!is.null(fold)) {
    nfold=length(unique(fold))
    mx=sort(unique(yp))
    
    fp=matrix(0,length(mx),nfold)
    tp=fp
    for (i in 1:nfold) {
      w=which(fold==i)
      r1=risksetROC(Stime=Yt[w],status=Y[w],marker=yp[w],predict.time=time,plot=F)
      fp[,i]=approx(r1$marker,r1$FP[2:(length(r1$FP)-1)],mx,rule=2)$y
      tp[,i]=approx(r1$marker,r1$TP[2:(length(r1$FP)-1)],mx,rule=2)$y
    }
    lines(rowMeans(fp),rowMeans(tp),...)
  } else {
    r1=risksetROC(Stime=Yt,status=Y,marker=yp,predict.time=time,plot=F)
    lines(r1$FP,r1$TP,...)
  }
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




##' Find optimal threshold according to Youden index, bootstrapped, and return a contingency table and accuracy
##' @name opt_threshold
##' @param yp vector of predicted yes/no
##' @param y vector of yes/no
##' @param cv vector of cross-validation folds
##' @param nboot either NULL to estimate Youden's index directly, or number of reps to use for bootstrap resampling
##' @param cutoffs vector of cutoffs to consider
##' @param printh set to a value to print results with the given header
##' @return invisibly: list with optimal threshold, contingency table, accuracy, and 95% CI.
opt_threshold=function(yp,y,cv,boot=NULL,cutoffs=quantile(yp,(1:500)/501),printh=NULL) {
  if (!is.null(boot)) {
    ycx=rep(NA,boot)
    for (i in 1:boot) {
      ir=sample(1:length(y),rep=TRUE); 
      yr=y[ir]; ypr=yp[ir]; cvr=cv[ir]
      tx=table(cvr,yr); if (min(tx[,2])>2) {
        aucx=quickroc(yr,ypr,cvr,res=cutoffs)
        obj=colMeans(aucx$sens) + colMeans(aucx$spec)
        ycx[i]=cutoffs[which.max(obj)]
      } else ycx[i]=NA
    }
    
    # Optimal thresholld
    opt=mean(ycx,na.rm=T)
  } else {
    aucx=quickroc(y,yp,cv,res=cutoffs)
    obj=colMeans(aucx$sens) + colMeans(aucx$spec)
    opt=cutoffs[which.max(obj)]
  }

  # Contingency table at optimal threshold
  cont=table(y,yp>opt)
  cont=rbind(cont,colSums(cont))
  cont=cbind(cont,rowSums(cont))
  colnames(cont)=c("Pred. no","Pred. yes","Total")
  rownames(cont)=c("No","Yes","Total")
  
  
  # Sensitivity and specificity
  sens=cont[2,2]/(cont[2,2] + cont[2,1])
  ci_sens=prop.test(cont[2,2],cont[2,2] + cont[2,1])$conf.int[1:2]
  
  spec=cont[1,1]/(cont[1,2] + cont[1,1])
  ci_spec=prop.test(cont[1,1],cont[1,2] + cont[1,1])$conf.int[1:2]
  
  
  # Accuracy and 95% CI
  acc=(cont[1,1]+cont[2,2])/cont[3,3]
  ci=prop.test(cont[1,1]+cont[2,2],cont[3,3])$conf.int[1:2]
  
  if (!is.null(printh)) {
    cat(paste0("Optimal threshold, contingency table and best accuracy for predictor ",printh,"\n\n"))
    cat(paste0("Value of optimal threshold: ",signif(opt,digits=3),"\n\n"))
    cat(paste0("Contingency table: \n"))
    print(cont)
    cat("\n\n")
    cat(paste0("Sensitivity: ",signif(sens,digits=3)," (95% CI [",signif(ci_sens[1],digits=3),",",signif(ci_sens[2],digits=3),"])"))
    cat("\n\n")
    cat(paste0("Specificity: ",signif(spec,digits=3)," (95% CI [",signif(ci_spec[1],digits=3),",",signif(ci_spec[2],digits=3),"])"))
    cat("\n\n")
    cat(paste0("Accuracy: ",signif(acc,digits=3)," (95% CI [",signif(ci[1],digits=3),",",signif(ci[2],digits=3),"])"))
    cat("\n\n\n")
  }
  
  return(invisible(list(threshold=opt,cont=cont,acc=acc,ci=ci)))
  
}


##' quickroc() 
##' Comprehensive plotting function for receiver-operator characteristic curve. Also calculates AUROC and standard error. 
##' 
##' Rather than returning points corresponding to every cutoff, only returns output at prespecified cutoffs
##'
##' SE of AUROC with no CV structure is from Hanley and McNeil 1982. SE of AUROC with CV folds is from LeDell et al 2012
##'
##' @param y class labels, 0/1 or logical
##' @param ypred predictions Pr(Y=1), numeric vector
##' @param cv cross-validation fold assignments, if relevant. Changes estimate of standard error.
##' @param res resolution. Returns this many equally-spaced points along the curve. Set res to null to return all points. Set res to a vector to use specific cutoffs.
##' @export 
##' @return list containing: spec, specificity for res points in every cv fold; sens, sensitivity for res points in every cv fold; auc, areas under the curve for each fold and average (note length is 1 greater than number of CV folds); se, standard error for AUC in each fold and standard error for average auc (note length is 1 greater than number of CV folds)
##' @examples 
quickroc=function(y,ypred,cv=NULL,res=NULL) {
  if (is.null(cv)) cv=rep(1,length(y))
  if (!(length(y)==length(ypred))) stop("Parameters y and ypred should have the same length")
  
  sens=c(); spec=c(); auc=c(); se=c(); cutoffs=c();
  for (i in 1:max(cv)) {
    w=which(cv==i)
    y0=y[w]; 
    ypred0=ypred[w]
    
    yt=sum(y0); yl=length(y0)
    opred=order(ypred0)
    #ipred=order(opred) # can use ipred to reorder in the order of original ypred
    
    sy=y0[opred]; sp=ypred0[opred]
    
    # Cutoffs and number of samples
    cutoffs0=sp
    
    csy=cumsum(sy); 
    csy1=cumsum(1-sy); 
    
    sens0=1- (csy/yt)
    spec0= csy1/(yl-yt)
    
    auc0=integral(sens0,spec0)
    se0=aucse(as.numeric(yt),as.numeric(yl-yt),auc0)
    
    if (!is.null(res)) {
      if (length(res)==1) {
        ds=cumsum(sqrt((spec0[1:(yl-1)]-spec0[2:yl])^2 + (sens0[1:(yl-1)]-sens0[2:yl])^2))
        ds=ds/ds[yl-1]
        lsp=(1:(yl-1))/yl
        sub=round(yl*approx(ds,lsp,n=res)$y)
        sens0=sens0[sub]
        spec0=spec0[sub]
        cutoffs0=cutoffs0[sub]
      } else {
        sens0=suppressWarnings(approx(cutoffs0,sens0,xout=res,rule=2)$y)
        spec0=suppressWarnings(approx(cutoffs0,spec0,xout=res,rule=2)$y)
        cutoffs0=res
      }
    }
    
    auc=c(auc,auc0)
    se=c(se,se0)
    spec=rbind(spec,spec0)
    sens=rbind(sens,sens0)
    cutoffs=rbind(cutoffs,cutoffs0)
  }
  
  if (length(auc)>1) {
    auc=c(auc,mean(auc))
    se=c(se,ci.cvAUC(ypred,y,folds=cv)$se)
  }
  
  out=list(sens=sens,spec=spec,cutoffs=cutoffs,auc=auc,se=se)
  class(out)="xROC"
  return(out)
}


##' integral() 
##' Quick form for trapezoidal integration over range of x
##'
##' @param x x co-ordinates, or nx2 matrix of points 
##' @param y y co-ordinates
##' @return trapezoidal estimate of integral of y[x] over range of x.
integral=function(x,y=NULL) {
  if (is.null(y)) {
    y=x[,2]; x=x[,1]
  }
  ox=order(x); xs=x[ox]; ys=y[ox]
  sum((xs[-1]-xs[-length(xs)])*(ys[-1]+ys[-length(ys)]))/2
}


# Internal function to compute SE of AUC
aucse=function(n1,n2,auc) {
  q1=auc/(2-auc); q2=2*(auc^2)/(1+auc)
  num=auc*(1-auc) + (n1-1)*(q1- (auc^2)) + (n2-1)*(q2-(auc^2))
  return(sqrt(num/(n1*n2)))
}




##' Find optimal threshold for survival ROC according to Youden index, bootstrapped, and return a contingency table and accuracy
##' @name opt_threshold_surv
##' @param yp vector of predicted yes/no
##' @param y vector of outcomes; event or censored
##' @param yt vector of survival times
##' @param tx time at which to consider ROC
##' @param cv vector of cross-validation folds
##' @param nboot either NULL to estimate Youden's index directly, or number of reps to use for bootstrap resampling
##' @param cutoffs vector of cutoffs to consider
##' @param printh set to a value to print results with the given header
##' @return invisibly: list with optimal threshold, contingency table, accuracy, and 95% CI.
opt_threshold_surv=function(yp,y,yt,tx,cv,boot=NULL,cutoffs=quantile(yp,(1:500)/501),printh=NULL) {
  sub=which(yt>=tx | y==1) # non-censored, or survived past time
  y=yt[sub]<tx; yp=yp[sub]; cv=cv[sub]
  
  if (!is.null(boot)) {
    ycx=rep(NA,boot)
    for (i in 1:boot) {
      ir=sample(1:length(y),rep=TRUE); 
      yr=y[ir]; ypr=yp[ir]; cvr=cv[ir]
      tx=table(cvr,yr); if (min(tx[,2])>2) {
        aucx=quickroc(yr,ypr,cvr,res=cutoffs)
        obj=colMeans(aucx$sens) + colMeans(aucx$spec)
        ycx[i]=cutoffs[which.max(obj)]
      } else ycx[i]=NA
    }
    
    # Optimal thresholld
    opt=mean(ycx,na.rm=T)
  } else {
    aucx=quickroc(y,yp,cv,res=cutoffs)
    obj=colMeans(aucx$sens) + colMeans(aucx$spec)
    opt=cutoffs[which.max(obj)]
  }
  
  # Contingency table at optimal threshold
  cont=table(y,yp>opt)
  cont=rbind(cont,colSums(cont))
  cont=cbind(cont,rowSums(cont))
  colnames(cont)=c("Pred. no","Pred. yes","Total")
  rownames(cont)=c("No","Yes","Total")
  
  # Sensitivity and specificity
  sens=cont[2,2]/(cont[2,2] + cont[2,1])
  ci_sens=prop.test(cont[2,2],cont[2,2] + cont[2,1])$conf.int[1:2]
  
  spec=cont[1,1]/(cont[1,2] + cont[1,1])
  ci_spec=prop.test(cont[1,1],cont[1,2] + cont[1,1])$conf.int[1:2]
  
  # Accuracy and 95% CI
  acc=(cont[1,1]+cont[2,2])/cont[3,3]
  ci=prop.test(cont[1,1]+cont[2,2],cont[3,3])$conf.int[1:2]
  
  if (!is.null(printh)) {
    cat(paste0("Optimal threshold, contingency table and best accuracy for predictor ",printh," at time ",tx," \n\n"))
    cat(paste0("Value of optimal threshold: ",signif(opt,digits=3),"\n\n"))
    cat(paste0("Contingency table: \n"))
    print(cont)
    cat("\n\n")
    cat(paste0("Sensitivity: ",signif(sens,digits=3)," (95% CI [",signif(ci_sens[1],digits=3),",",signif(ci_sens[2],digits=3),"])"))
    cat("\n\n")
    cat(paste0("Specificity: ",signif(spec,digits=3)," (95% CI [",signif(ci_spec[1],digits=3),",",signif(ci_spec[2],digits=3),"])"))
    cat("\n\n")
    cat(paste0("Accuracy: ",signif(acc,digits=3)," (95% CI [",signif(ci[1],digits=3),",",signif(ci[2],digits=3),"])"))
    cat("\n\n\n")
  }
  
  return(invisible(list(threshold=opt,cont=cont,acc=acc,ci=ci)))
  
}

##' Bootstrap resampler for ROC curves for survival. 
##' 
##' @name se_auc_surv
##' @param Y event indicator
##' @param Yt time indicator
##' @param px predictor
##' @param ntrial number of trials
##' @param time time at which to compare
##' @return ntrial copies of AUC
se_auc_surv=function(Y,Yt,px,ntrial,time=tx) {
  xauc=rep(0,ntrial);
  for (i in 1:ntrial) {
    s=sample(1:length(Y),rep=T)
    xauc[i]=risksetROC(Stime=Yt[s],status=Y[s],marker=px[s],predict.time=tx,plot=F)$AUC
  }
  return(xauc)
}
  
  
##' Bootstrap comparison of ROC curves for survival
##' 
##' @name roc_test_surv
##' @param v1 AUC for predictor 1
##' @param v2 AUC for predictor 2
##' @param s1 bootstrap resamples for predictor 1
##' @param s2 bootstrap resamples for predictor 1
##' @param time time at which to compare
##' @return p-value, invisibly
roc_test_surv=function(v1,v2,s1,s2,time) {
  ntrial=length(s1)
  pval=t.test(s1,s2,var.equal=FALSE)$p.value
  mx1=length(which(s1>v2)); lx1=length(which(s1<=v2))
  p1=2*min(mx1/ntrial,lx1/ntrial)
  mx2=length(which(s1>v2)); lx2=length(which(s1<=v2))
  p2=2*min(mx2/ntrial,lx2/ntrial)
  pval=max(p1,p2)

  cat("Bootstrap test for comparison of two survival ROC curves (two-sided)\n\n")
  cat(paste0("Time: ",time,"\n\n"))
  cat(paste0("Number of trials: ",ntrial,"\n\n"))
  cat(paste0("P-value: ",signif(pval,digits=3),"\n\n"))
  return(invisible(pval))
}




##' Basic concordance
##' @name conc
##' @param y target
##' @param yp prediction
##' @return list of two values: concordance and SE
conc=function(y,yp) {
  s1=which(y==1); s0=which(y==0)
  n=length(y)
  tot=0
  for (i in 1:length(s1)) tot=tot + length(which(yp[s1[i]]>yp[s0]))
  denom=length(s1)*length(s0)
  xconc=tot/denom;
  concse=aucse(length(s1),length(s0),xconc)
  return(list(concordance=xconc,se=concse))
}


##' Auxiliary function to generate table 1
##' @param tab data frame which should contain all fields.
make_tab1=function(tab) {
  
  # Auxiliary functions
  miqr=function(x,dg=0) {
    if (dg==0) out=paste0(round(median(x,na.rm=T))," ± ",
                          round((quantile(x,0.75,na.rm=T)-quantile(x,0.25,na.rm=T))/2))
    if (dg>0) out=paste0(round(median(x,na.rm=T),digits=dg)," ± ",
                         round((quantile(x,0.75,na.rm=T)-quantile(x,0.25,na.rm=T))/2,digits=dg))
    return(out)
  }
  nperc=function(x) {
    paste0(sum(x,na.rm=T)," (",
           round(mean(100*x,na.rm=T)),")")
  }
  pform=function(x,y) {
    x1=x[which(is.finite(x))]
    y1=y[which(is.finite(y))]
    p=t.test(x1,y1)$p.value
    if (p<0.001) out="<0.001" else out=as.character(signif(p,digits=2))
  }
  
  
  tab1=c()
  
  # Total n
  rx=c("Total n",  
       dim(tab)[1],  
       "",   
       length(which(tab$Dead_before_followup==0)),
       "",
       "")
  tab1=rbind(tab1,rx)
  
  # Age at PEA
  rx=c("Age at PEA, yrs",
       length(which(is.finite(tab$age_pea))),
       miqr(tab$age_pea),
       "",
       "",
       "")
  tab1=rbind(tab1,rx)
  
  # Male sex
  rx=c("Male sex, n (%)",
       length(which(!is.na(tab$sex))),
       nperc(tab$sex),
       "",
       "",
       "")
  tab1=rbind(tab1,rx)
  
  # BMI
  rx=c("BMI, kg/m²",
       length(which(!is.na(tab$BMI))),
       miqr(tab$BMI),
       "",
       "",
       "")
  tab1=rbind(tab1,rx)
  
  # FEV1/FVC
  rx=c("FEV1/FVC, %",
       length(which(!is.na(tab$pft_bl_fev1.fvc_pc))),
       miqr(tab$pft_bl_fev1.fvc_pc),
       "",
       "",
       "")
  tab1=rbind(tab1,rx)
  
  # Smoker
  rx=c("Smoker*, n (%)",
       length(which(!is.na(tab$smoking_status))),
       nperc(tab$smoking_status>1),
       "",
       "",
       "")
  tab1=rbind(tab1,rx)
  
  # Comorbidities
  tab1=rbind(tab1,c("Comorbidities, n (%)","","","","",""))
  rx=c("   Atrial arrhythmia",
       length(which(!is.na(tab$comorbid_af))),
       nperc(tab$comorbid_af),
       "",
       "",
       "")
  tab1=rbind(tab1,rx)
  rx=c("   Systemic hypertension",
       length(which(!is.na(tab$comorbid_htn))),
       nperc(tab$comorbid_htn),
       "",
       "",
       "")
  tab1=rbind(tab1,rx)
  rx=c("   Diabetes Mellitis",
       length(which(!is.na(tab$comorbid_dm))),
       nperc(tab$comorbid_dm),
       "",
       "",
       "")
  tab1=rbind(tab1,rx)
  rx=c("   Chronic renal disease", 
       length(which(!is.na(tab$comorbid_renal))),
       nperc(tab$comorbid_renal),
       "",
       "",
       "")
  tab1=rbind(tab1,rx)
  rx=c("   Ischaemic heart disease†",
       length(which(!is.na(tab$comorbid_ihd))),
       nperc(tab$comorbid_ihd),
       "",
       "",
       "")
  tab1=rbind(tab1,rx)
  rx=c("   History of malignancy",
       length(which(!is.na(tab$comorbid_malig_solid_haem))),
       nperc(tab$comorbid_malig_solid_haem),
       "",
       "",
       "")
  tab1=rbind(tab1,rx)
  rx=c("   Thrombophilia",  
       length(which(!is.na(tab$comorbid_thrombophilia))),
       nperc(tab$comorbid_thrombophilia),
       "",
       "",
       "")
  tab1=rbind(tab1,rx)
  rx=c("   Thyroid dysfunction",
       length(which(!is.na(tab$comorbid_thyroid))),
       nperc(tab$comorbid_thyroid),
       "",
       "",
       "")
  tab1=rbind(tab1,rx)
  
  # Vasodilators
  rx=c("Pulmonary vasodilator, n (%)",
       length(which(!is.na(tab$med_vasodilator))),
       nperc(tab$med_vasodilator),
       "",
       "",
       "")
  tab1=rbind(tab1,rx)
  
  # Haemodynamics
  tab1=rbind(tab1,c("Haemodynamics","","","","",""))
  rx=c("   Mean PAP, mmHg",
       length(which(!is.na(tab$mpap_bl))),
       miqr(tab$mpap_bl),
       length(which(!is.na(tab$mpap_fu1))),
       miqr(tab$mpap_fu1),
       pform(tab$mpap_bl,tab$mpap_fu1))
  tab1=rbind(tab1,rx)
  rx=c("   PVR, dynes cm s-5",
       length(which(!is.na(tab$pvr_bl))),
       miqr(tab$pvr_bl),
       length(which(!is.na(tab$pvr_fu1))),
       miqr(tab$pvr_fu1),
       pform(tab$pvr_bl,tab$pvr_fu1))
  tab1=rbind(tab1,rx)
  rx=c("   PAWP, mmHg",
       length(which(!is.na(tab$pcwp_bl))),
       miqr(tab$pcwp_bl),
       length(which(!is.na(tab$pcwp_fu1))),
       miqr(tab$pcwp_fu1),
       pform(tab$pcwp_bl,tab$pcwp_fu1))
  tab1=rbind(tab1,rx)
  rx=c("   CI, l/min/m²",
       length(which(!is.na(tab$ci_bl))),
       miqr(tab$ci_bl,dg=1),
       length(which(!is.na(tab$ci_fu1))),
       miqr(tab$ci_fu1,dg=1),
       pform(tab$ci_bl,tab$ci_fu1))
  tab1=rbind(tab1,rx)
  
  
  
  # Functional status
  tab1=rbind(tab1,c("Functional status","","","","",""))
  rx=c("   NYHA, 1/2/3/4 %",
       length(which(!is.na(tab$nyha_bl))),
       miqr(tab$nyha_bl),
       length(which(!is.na(tab$nyha_fu1))),
       miqr(tab$nyha_fu1),
       pform(tab$nyha_bl,tab$nyha_fu1))
  tab1=rbind(tab1,rx)
  rx=c("   6MWD ‡, metres",
       length(which(!is.na(tab$sixmwt_bl))),
       miqr(tab$sixmwt_bl),
       length(which(!is.na(tab$sixmwt_fu1))),
       miqr(tab$sixmwt_fu1),
       pform(tab$sixmwt_bl,tab$sixmwt_fu1))
  tab1=rbind(tab1,rx)
  
  # CAMPHOR
  tab1=rbind(tab1,c("CAMPHOR","","","","",""))
  rx=c("   Symptoms",
       length(which(!is.na(tab$BL.Symptom))),
       miqr(tab$BL.Symptom),
       length(which(!is.na(tab$FU.symptom))),
       miqr(tab$FU.symptom),
       pform(tab$BL.Symptom,tab$FU.symptom))
  tab1=rbind(tab1,rx)
  rx=c("   Activity",
       length(which(!is.na(tab$BL.Activity))),
       miqr(tab$BL.Activity),
       length(which(!is.na(tab$FU.activity))),
       miqr(tab$FU.activity),
       pform(tab$BL.Activity,tab$FU.activity))
  tab1=rbind(tab1,rx)
  rx=c("   Quality of Life",
       length(which(!is.na(tab$BL.QoL))),
       miqr(tab$BL.QoL),
       length(which(!is.na(tab$FU.qol))),
       miqr(tab$FU.qol),
       pform(tab$BL.QoL,tab$FU.qol))
  tab1=rbind(tab1,rx)
  
  # Intra-op
  tab1=rbind(tab1,c("Intra-operative","","","","",""))
  rx=c("   CPB time, mins",
       length(which(!is.na(tab$bypass_mins))),
       miqr(tab$bypass_mins),
       "",
       "",
       "")
  tab1=rbind(tab1,rx)
  rx=c("   DHCA time, mins",
       length(which(!is.na(tab$dhca_min_total))),
       miqr(tab$dhca_min_total),
       "",
       "",
       "")
  tab1=rbind(tab1,rx)
  
  # Other surgery
  tab1=rbind(tab1,c("Other surgery, n (%)","","","","",""))
  rx=c("   CABG",
       length(which(!is.na(tab$cabg))),
       nperc(tab$cabg),
       "",
       "",
       "")
  tab1=rbind(tab1,rx)
  rx=c("   AVR",
       length(which(!is.na(tab$avr))),
       nperc(tab$avr),
       "",
       "",
       "")
  tab1=rbind(tab1,rx)
  rx=c("   MVR",
       length(which(!is.na(tab$mvr))),
       nperc(tab$mvr),
       "",
       "",
       "")
  tab1=rbind(tab1,rx)
  rx=c("   ASD/PFO closure",
       length(which(!is.na(tab$asd))),
       nperc(tab$asd),
       "",
       "",
       "")
  tab1=rbind(tab1,rx)
  
  # Complications
  tab1=rbind(tab1,c("Complications, n (%)","","","","",""))
  rx=c("   CPAP",
       length(which(!is.na(tab$cpap))),
       nperc(tab$cpap),
       "",
       "",
       "")
  tab1=rbind(tab1,rx)
  rx=c("   Haemofiltration",
       length(which(!is.na(tab$renal_cvvd))),
       nperc(tab$renal_cvvd),
       "",
       "",
       "")
  tab1=rbind(tab1,rx)
  rx=c("   ECMO",
       length(which(!is.na(tab$ecmo))),
       nperc(tab$ecmo),
       "",
       "",
       "")
  tab1=rbind(tab1,rx)
  rx=c("   Pneumonia",
       length(which(!is.na(tab$chest_infection))),
       nperc(tab$chest_infection),
       "",
       "",
       "")
  tab1=rbind(tab1,rx)
  rx=c("   Return to theatre",
       length(which(!is.na(tab$return_theatre))),
       nperc(tab$return_theatre),
       "",
       "",
       "")
  tab1=rbind(tab1,rx)
  rx=c("   Reperfusion injury",
       length(which(!is.na(tab$comp_reperfusion))),
       nperc(tab$comp_reperfusion),
       "",
       "",
       "")
  tab1=rbind(tab1,rx)
  
  # Intubation
  rx=c("Intubation, days",
       length(which(!is.na(tab$time_to_extubation))),
       nperc(tab$time_to_extubation),
       "",
       "",
       "")
  tab1=rbind(tab1,rx)
  
  # ICU stay
  rx=c("ICU stay, days",
       length(which(!is.na(tab$icu_days))),
       miqr(tab$icu_days),
       "",
       "",
       "")
  tab1=rbind(tab1,rx)
  
  # Total inpatient stay
  rx=c("Total inpatient stay, days",
       length(which(!is.na(tab$hosp_days))),
       miqr(tab$hosp_days),
       "",
       "",
       "")
  tab1=rbind(tab1,rx)
  
  # Inpatient death
  rx=c("Inpatient death, n (%)",
       length(which(!is.na(tab$death_in_hosp))),
       nperc(tab$death_in_hosp),
       "",
       "",
       "")
  tab1=rbind(tab1,rx)
  
  return(tab1)
}