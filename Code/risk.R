##**********************************************************************
##  Risk prediction in PEA                                          ####
##  Supporting code for main analyses                               
##  James Liley                                                     
##  9/3/19                                                          
##**********************************************************************


##**********************************************************************
## Switches                                                         ####
##**********************************************************************

do_assoc_test=TRUE # set to TRUE to do association tests
normalise=F

# force redo
force_redo_dm=FALSE
force_redo_5m=FALSE
force_redo_dq=FALSE
force_redo_all=FALSE

# missingness
low_missingness=FALSE

# random seed
set.seed(1)


##**********************************************************************
## Set data directories                          ####
##**********************************************************************

## D?
if (file.exists("~/hpc")) loc="med_sshfs" 
if (file.exists("/rds/project/who1000-1/rds-who1000-cbrc/user/sg550/")) loc="hpc"
if (file.exists("~/Old/CTEPH/Risk/")) loc="home"


if (loc=="med_sshfs") setwd("~/hpc/CTEPH/Risk/")
if (loc=="hpc") setwd("~/CTEPH/Risk/")
if (loc=="home") setwd("~/Old/CTEPH/Risk/")

datadir="data/" # Data location

# Prefix all saved filenames with this according to missingness
if (low_missingness) prefix="nmiss_" else prefix=""



##**********************************************************************
## Packages and scripts                                             ####
##**********************************************************************

library(glmnet)
library(randomForest)
library(snpStats)
library(survival)
library(randomForestSRC)
library(pec)
library(penalized)
library(cvAUC)
library(risksetROC)
library(lmtest)

source("code_new/process_raw.R") # Process raw data into usable format
source("code_new/auxiliary.R") # Auxiliary functions


##**********************************************************************
## Read data                                                        ####
##**********************************************************************

load("data/maintable.RData") # contains main table and predictor details


# Restrict to only <50% missing
if (low_missingness) {
 missingness=rep(0,dim(tab)[2])
 for (i in 1:dim(tab)[2]) missingness[i]=length(which(is.na(tab[,i])))+length(which(as.character(tab[,i])=="NA"))
 missingness=missingness/dim(tab)[1]
 tab=tab[which(missingness<0.5)]
}



##**********************************************************************
## Exclusions and summary                                           ####
##**********************************************************************


## Exclusions
e1=which(tab$excluded==1)
e2=which(tab$mpap_bl<25)
e3=which(!is.finite(tab$HospNo))

# Excluded because of mpap<25 mmHg at baseline
print(paste("Excluded due to mPap<25 mmHg at baseline:", length(e2)))

# Excluded for other reason
print(paste("Excluded due to other reason:", length(setdiff(e1,e2))))

print(paste("Included total:", dim(tab)[1]-length(unique(c(e1,e2,e3)))))



exclude=unique(c(e1,e2,e3))
tab=tab[-exclude,]
tab_raw=tab_raw[-exclude,]

## Variable exclusions:

var_exclude=c(); thresh=0.1*dim(tab)[1]
for (i in 1:dim(tab)[2]) if (length(which(!is.na(tab[,i])))<thresh) var_exclude=c(var_exclude,colnames(tab)[i])
tab=tab[,setdiff(colnames(tab),var_exclude)]
predictors=intersect(predictors,colnames(tab))
noninv_predictors=intersect(noninv_predictors,colnames(tab))
preop_predictors=intersect(preop_predictors,colnames(tab))
discharge_predictors=intersect(discharge_predictors,colnames(tab))



# Outcome summaries

# Mortality
print(paste("Deaths before followup:", sum(tab$Dead_before_followup,na.rm=T)))
print(paste("Deaths within 5y:", sum(tab$surv_days_max5yrs_census,na.rm=T)))
print(paste("Survived 5y:", length(which(tab$surv_days_max5yrs>= 1826 & tab$surv_days_max5yrs_census==0))))
print(paste("Censored before 5y:", length(which(tab$surv_days_max5yrs < 1826 & tab$surv_days_max5yrs_census==0))))

# Morbidity
print(paste("CAMPHOR available at baseline:", length(which(is.finite(tab$BL.Symptom+tab$BL.Activity+tab$BL.QoL)))))
print(paste("CAMPHOR available at followup:", length(which(is.finite(tab$FU.symptom+tab$FU.activity+tab$FU.qol)))))
print(paste("Mean CAMPHOR at baseline:", mean(tab$BL.Symptom+tab$BL.Activity+tab$BL.QoL,na.rm=T),
  "(",mean(tab$BL.Symptom,na.rm=T),mean(tab$BL.Activity,na.rm=T),
  mean(tab$BL.QoL,na.rm=T),")"))
print(paste("Mean CAMPHOR at followup:", mean(tab$FU.symptom+tab$FU.activity+tab$FU.qol,na.rm=T),
  "(",mean(tab$FU.symptom,na.rm=T),mean(tab$FU.activity,na.rm=T),
  mean(tab$FU.qol,na.rm=T),")"))
print(paste("Patients with CAMPHOR at baseline but not FU, survived to FU:",
length(which(!is.na(tab$BL.Activity + tab$BL.Symptom + tab$BL.QoL) & is.na(tab$FU.activity  + tab$FU.symptom + tab$FU.qol) & !(tab$Dead_before_followup>0)))))
print(paste("Patients with CAMPHOR at baseline but not FU, did not survive to FU:",
length(which(!is.na(tab$BL.Activity + tab$BL.Symptom + tab$BL.QoL) & is.na(tab$FU.activity  + tab$FU.symptom + tab$FU.qol) & (tab$Dead_before_followup>0)))))




##**********************************************************************
## Predictor matrix                                                 ####
##**********************************************************************

Xall=tab[,predictors]

# Remove correlated variables
#Xallcor=cor(Xall); m=dim(Xall)[2]
#ntrim=colnames(Xall)
#for (i in 1:(m-1)) {
#    w=colnames(Xall)[(1+i):m][which(abs(Xallcor[(1+i):m,i])>0.9)]
#    ntrim=setdiff(ntrim,w)
#}
#Xall=Xall[,ntrim]




##**********************************************************************
## Trait to be predicted (target)                                   ####
##**********************************************************************

# Randomise to check correctness/internal consistency of algorithm
randomise=F

# Traits to be predicted
YDM=tab$Dead_before_followup # More powerful than in-hospital death/30 day mortality
Y5M=tab$surv_days_max5yrs_census; Y5M_time=tab$surv_days_max5yrs
YDC=(tab$BL.Activity + tab$BL.Symptom + tab$BL.QoL) - (tab$FU.activity  +tab$FU.symptom + tab$FU.qol); wm=which(tab$Dead_before_followup>0); YDC[wm]=(tab$BL.Activity + tab$BL.Symptom + tab$BL.QoL)[wm] - max(tab$FU.activity  + tab$FU.symptom + tab$FU.qol)

if (randomise) {
   YDM=YDM[order(runif(length(YDM)))]
   Y5M=Y5M[order(runif(length(Y5M)))]
   YDC=YDC[order(runif(length(YDC)))]
}


##**********************************************************************
## Association tests                                                ####
##**********************************************************************

if (do_assoc_test) {
   
outstat=c()

vars=intersect(colnames(Xall),preop_predictors)

for (i in 1:length(vars)) {
    xx=as.numeric(tab[,vars[i]])
    
  missingness=length(which(is.na(xx)))/dim(tab)[1]
  nm=which(!is.na(xx))
  
  if (length(nm)> 0.1*dim(tab)[1]) {
  testpm=t.test(xx[intersect(nm,which(YDM==1))],xx[intersect(nm,which(YDM==0))])
  ppm=testpm$p.value
  xpm=testpm$statistic # t-value
  dpm=testpm$parameter
  
  test5m=summary(coxph(Surv(time = Y5M_time, event=Y5M)~ xx,subset=nm))
  p5m=test5m$logtest['pvalue']
  x5m=test5m$logtest['test']
  d5m=test5m$logtest['df']
  
  if (vars[i] %in% c("FU.symptom","FU.activity","FU.qol")) {
     pdc=NA; xdc=NA; ddc=NA
  } else {
    nmdq=intersect(nm,which(!is.na(YDC)))
    testdc=cor.test(xx[nmdq],YDC[nmdq])
    pdc=testdc$p.value
    xdc=testdc$statistic
    ddc=testdc$parameter
  }
  
  outstat=rbind(outstat,c(missingness,ppm,xpm,dpm,p5m,x5m,d5m,pdc,xdc,ddc))
  } else outstat=rbind(outstat,c(missingness,rep(NA,9)))
}
rownames(outstat)=vars
colnames(outstat)=c("missing","p.pm","t.pm","df.pm","p.5m","lr.5m","df.5m","p.dc","t.dc","df.dc")


# Determine significant values
fdr=0.1 # FDR control level

pval=as.vector(outstat[,2:7]); pval=pval[which(!is.na(pval))]
pcut=sort(pval)[max(rank(pval)[which(pval*length(pval)/rank(pval)<fdr)])]

outx=outstat[which(rowMins(outstat[,c(2,5,8)])<pcut),]

# Print
fullnames=read.table("data/raw_input_datasets/predictor_details.txt",stringsAs=F,sep=",") # details of predictors

processnum=function(x,digits=2,sthresh=1e-3,pthresh=pcut) {
  if (abs(x)<pthresh) suff="$^*$" else suff=""
                         if (abs(x)>sthresh | x==0) line=paste0(signif(x,digits=digits),suff)
                          else line=paste0(formatC(x,format="e",digits=digits),suff)
                line
}
printx=function(x,names,digits=2,sthresh=1e-3,pthresh=pcut) {
for (i in 1:dim(x)[1]) {
      xx=x[i,]
    line0=paste0(fullnames[match(rownames(x)[i],names),2]," & ",signif(xx[1],digits=2))
    line1=paste0(processnum(xx[2])," & ",processnum(xx[3])," (",round(xx[4]),")")
    line2=paste0(processnum(xx[5])," & ",processnum(xx[6])) #," (",round(xx[7]),")")
    line3=paste0(processnum(xx[8])," & ",processnum(xx[9])," (",round(xx[10]),")")
    line=paste0(line0," & ",line1," & ",line2," & ",line3," \\")
    line=gsub("e-([0-9][0-9])","$ \times 10^\\{-\\1\\} $",line,fixed=F)
    line=gsub("\\{-0([0-9])}\\}","\\{-\\1\\} $",line,fixed=F)
#   line=gsub("%","\%",fixed=T,line)
    print(line,quote=F)
  }
} 
printx(outx[order(outx[,2]),],names=fullnames[,1])

mx=match(rownames(outstat),fullnames[,1])
outwrite=cbind(fullnames[mx,3],fullnames[mx,2],outstat)
colnames(outwrite)[1:2]=c("Full_variable_name","Short_variable_name")
write.csv(outwrite,file=paste0("data/",prefix,"variable_table.csv"),row.names=F)
}





##**********************************************************************
## Imputation (singular, mean value)                                ####
##**********************************************************************

for (i in 1:dim(Xall)[2]) Xall[which(!is.finite(Xall[,i])),i]=mean(Xall[which(is.finite(Xall[,i])),i])
if (normalise) for (i in 1:dim(Xall)[2]) Xall[,i]=(Xall[,i]-mean(Xall[,i]))/sd(Xall[,i])
rmXall=c(); for (i in 1:dim(Xall)[2]) if (!(sd(Xall[,i])>0)) rmXall=c(rmXall,colnames(Xall)[i])
Xall=Xall[,setdiff(colnames(Xall),rmXall)]




##**********************************************************************
## Overall assessment of predictability                             ####
##**********************************************************************



Y=YDM

predset=c("noninv","preop","discharge","all")
fit_pm=rep(1,4); names(fit_pm) = paste0("pm_",predset)

for (i in 1:4) {
 if (i==1) X=Xall[,intersect(colnames(Xall),noninv_predictors)]
 if (i==2) X=Xall[,intersect(colnames(Xall),preop_predictors)]
 if (i==3) X=Xall[,intersect(colnames(Xall),discharge_predictors)]
 if (i==4) X=Xall[,intersect(colnames(Xall),predictors)]
 X=X[which(!is.na(Y)),]
 Y=Y[which(!is.na(Y))]

 l1=lm(Y~.,data=X); l0=lm(Y~1,data=X)
 fit_pm[i]=lrtest(l1,l0)[5][2,1]
}




Y=Y5M; Yt=Y5M_time
w1=which(!is.na(Y+Yt))
Y=Y[w1]; Yt=Yt[w1]

predset=c("noninv","preop","discharge","all")
fit_5m=rep(1,4); names(fit_5m) = paste0("x5m_",predset)

for (i in 1:4) {
 if (i==1) X=Xall[,intersect(colnames(Xall),noninv_predictors)]
 if (i==2) X=Xall[,intersect(colnames(Xall),preop_predictors)]
 if (i==3) X=Xall[,intersect(colnames(Xall),discharge_predictors)]
 if (i==4) X=Xall[,intersect(colnames(Xall),predictors)]

 X=X[w1,]


 l1=coxph(Surv(event=Y,time=Yt)~.,data=X); 
 l0=coxph(Surv(event=Y,time=Yt)~1,data=X)
 fit_5m[i]=summary(l1)[[9]][3]
}





Y=YDC
w1=which(!is.na(Y))
Y=Y[w1]; 

predset=c("noninv","preop","discharge","all")
fit_dq=rep(1,4); names(fit_dq) = paste0("dq_",predset)

morbidity_pars_fu=c("FU.activity","FU.symptom","FU.qol")

for (i in 1:4) {
 if (i==1) X=Xall[,intersect(colnames(Xall),noninv_predictors)]
 if (i==2) X=Xall[,intersect(colnames(Xall),preop_predictors)]
 if (i==3) X=Xall[,intersect(colnames(Xall),discharge_predictors)]
 if (i==4) X=Xall[,intersect(colnames(Xall),predictors)]
 
 # Remove morbidity predictors and NA observations
 X=X[w1,setdiff(colnames(X),morbidity_pars_fu)]

 l1=lm(Y~.,data=X); 
 l0=lm(Y~1,data=X)
 fit_dq[i]=lrtest(l1,l0)[5][2,1]
}


print(paste0("Overall assessment of fit: "))
print(c(fit_pm, fit_5m, fit_dq))


##**********************************************************************
## Predictive analysis for YDM                                      ####
##**********************************************************************

Y=YDM
w1=which(!is.na(Y))
Y=Y[w1]

# Split observations into k=10 equally sized folds
nfold=10
fold1=rep(1:nfold,1+floor(length(Y)/nfold))[1:length(Y)]; fold1=fold1[order(runif(length(fold1)))]
fold2=rep(1:nfold,1+floor(length(Y)/nfold))[1:length(Y)]; fold2=fold2[order(runif(length(fold2)))]



ydm_name=paste0("./data/",prefix,"dm_outputs")

predset=c("noninv","preop","discharge","all")

if (!file.exists(ydm_name) | (force_redo_dm|force_redo_all)) {
for (i in 1:4) {

pred=predset[i]
if (i==1) X=Xall[,intersect(colnames(Xall),noninv_predictors)]
if (i==2) X=Xall[,intersect(colnames(Xall),preop_predictors)]
if (i==3) X=Xall[,intersect(colnames(Xall),discharge_predictors)]
if (i==4) X=Xall[,intersect(colnames(Xall),predictors)]
X=X[w1,]

# Function to assess accuracy of predictions
ydif=function(pred,Y) 1-auroc(pred,Y) # 'Difference' between predicted and expected Y

models=c("lrproc","lassoproc","rfproc")
best=findbestmodel(X,Y,fold1,models,"ydif")$best # best model
assign(paste0("best_",predset[i],"_dm"),best)
for (j in 1:length(models)) assign(paste0("bestpars",j),get(models[j])(X,Y,return_hyperparameters=T)) # best parameters for each model 
for (j in 1:length(models)) assign(paste0("Ypred_",predset[i],"_dm_",models[j]),evaluateperformance(X,Y,fold=fold2,models[j],get(paste0("bestpars",j))))
print(i)
}

save(Y,Ypred_noninv_dm_lrproc,Ypred_noninv_dm_lassoproc,Ypred_noninv_dm_rfproc,Ypred_preop_dm_lrproc,Ypred_preop_dm_lassoproc,Ypred_preop_dm_rfproc,Ypred_discharge_dm_lrproc,Ypred_discharge_dm_lassoproc,Ypred_discharge_dm_rfproc,Ypred_all_dm_lrproc,Ypred_all_dm_lassoproc,Ypred_all_dm_rfproc,best_noninv_dm,best_preop_dm,best_discharge_dm,best_all_dm,fold1,fold2,ydif,file=ydm_name)
} else load(ydm_name)



# plot each outcome as ROC curve
# LM
pdf(paste0("./plots/",prefix,"dm_lr.pdf"),width=4,height=4)
plot(0,type="n",xlim=c(0,1),ylim=c(0,1),xlab="1-Specificity",ylab="Sensitivity",xaxs="i",yaxs="i")
draw_roc_cv(Ypred_noninv_dm_lrproc,Y,fold=fold2,add=T,col=2)
draw_roc_cv(Ypred_preop_dm_lrproc,Y,fold=fold2,add=T,col=3)
draw_roc_cv(Ypred_discharge_dm_lrproc,Y,fold=fold2,add=T,col=4)
#draw_roc_cv(Ypred_all_dm_lrproc,Y,fold=fold2,add=T,col=5)
abline(0,1)
#legend("bottomright",c("Noninv.","Preop.","Discharge","All"),lty=1,col=2:5)
legend("bottomright",c("Noninv.","Preop.","Discharge"),lty=1,col=2:4)
dev.off()



# Lasso
pdf(paste0("./plots/",prefix,"dm_lasso.pdf"),width=4,height=4)
plot(0,type="n",xlim=c(0,1),ylim=c(0,1),xlab="1-Specificity",ylab="Sensitivity",xaxs="i",yaxs="i")
draw_roc_cv(Ypred_noninv_dm_lassoproc,Y,fold=fold2,add=T,col=2)
draw_roc_cv(Ypred_preop_dm_lassoproc,Y,fold=fold2,add=T,col=3)
draw_roc_cv(Ypred_discharge_dm_lassoproc,Y,fold=fold2,add=T,col=4)
#draw_roc_cv(Ypred_all_dm_lassoproc,Y,fold=fold2,add=T,col=5)
abline(0,1)
#legend("bottomright",c("Noninv.","Preop.","Discharge","All"),lty=1,col=2:5)
legend("bottomright",c("Noninv.","Preop.","Discharge"),lty=1,col=2:4)
dev.off()

# RF
pdf(paste0("./plots/",prefix,"dm_rf.pdf"),width=4,height=4)
plot(0,type="n",xlim=c(0,1),ylim=c(0,1),xlab="1-Specificity",ylab="Sensitivity",xaxs="i",yaxs="i")
draw_roc_cv(Ypred_noninv_dm_rfproc,Y,fold=fold2,add=T,col=2)
draw_roc_cv(Ypred_preop_dm_rfproc,Y,fold=fold2,add=T,col=3)
draw_roc_cv(Ypred_discharge_dm_rfproc,Y,fold=fold2,add=T,col=4)
#draw_roc_cv(Ypred_all_dm_rfproc,Y,fold=fold2,add=T,col=5)
abline(0,1)
#legend("bottomright",c("Noninv.","Preop.","Discharge","All"),lty=1,col=2:5)
legend("bottomright",c("Noninv.","Preop.","Discharge"),lty=1,col=2:4)
dev.off()





# plot calibrations
# LM
pdf(paste0("./plots/",prefix,"dm_lr_calibration.pdf"),width=4,height=4)
plot(0,type="n",xlim=c(0,1.1),ylim=c(0,1.1),xlab="Predicted risk",ylab="Observed risk",xaxs="i",yaxs="i")
draw_calibration(logit(Ypred_noninv_dm_lrproc),Y,type="l",add=T,col=2)
draw_calibration(logit(Ypred_preop_dm_lrproc),Y,type="l",add=T,col=3)
draw_calibration(logit(Ypred_discharge_dm_lrproc),Y,type="l",add=T,col=4)
#draw_calibration(logit(Ypred_all_dm_lrproc),Y,type="l",add=T,col=5)
abline(0,1)
#legend("bottomright",c("Noninv.","Preop.","Discharge","All"),lty=1,col=2:5)
legend("bottomright",c("Noninv.","Preop.","Discharge"),lty=1,col=2:4)
dev.off()

# Lasso
pdf(paste0("./plots/",prefix,"dm_lasso_calibration.pdf"),width=4,height=4)
plot(0,type="n",xlim=c(0,1.1),ylim=c(0,1.1),xlab="Predicted risk",ylab="Observed risk",xaxs="i",yaxs="i")
draw_calibration(logit(Ypred_noninv_dm_lassoproc),Y,type="l",add=T,col=2)
draw_calibration(logit(Ypred_preop_dm_lassoproc),Y,type="l",add=T,col=3)
draw_calibration(logit(Ypred_discharge_dm_lassoproc),Y,type="l",add=T,col=4)
#draw_calibration(logit(Ypred_all_dm_lassoproc),Y,type="l",add=T,col=5)
abline(0,1)
#legend("bottomright",c("Noninv.","Preop.","Discharge","All"),lty=1,col=2:5)
legend("bottomright",c("Noninv.","Preop.","Discharge"),lty=1,col=2:4)
dev.off()

# RF
pdf(paste0("./plots/",prefix,"dm_rf_calibration.pdf"),width=4,height=4)
plot(0,type="n",xlim=c(0,1.1),ylim=c(0,1.1),xlab="Predicted risk",ylab="Observed risk",xaxs="i",yaxs="i")
draw_calibration(Ypred_noninv_dm_rfproc,Y,type="l",add=T,col=2)
draw_calibration(Ypred_preop_dm_rfproc,Y,type="l",add=T,col=3)
draw_calibration(Ypred_discharge_dm_rfproc,Y,type="l",add=T,col=4)
#draw_calibration(Ypred_all_dm_rfproc,Y,type="l",add=T,col=5)
abline(0,1)
#legend("bottomright",c("Noninv.","Preop.","Discharge","All"),lty=1,col=2:5)
legend("bottomright",c("Noninv.","Preop.","Discharge"),lty=1,col=2:4)
dev.off()







# Save AUCs and confidence intervals as table

tab_pm=c()

ccn=formatci(ci.cvAUC(Ypred_noninv_dm_lrproc,Y,fold=fold2))
ccp=formatci(ci.cvAUC(Ypred_preop_dm_lrproc,Y,fold=fold2))
ccd=formatci(ci.cvAUC(Ypred_discharge_dm_lrproc,Y,fold=fold2))
cca=formatci(ci.cvAUC(Ypred_all_dm_lrproc,Y,fold=fold2))

tab_pm=cbind(tab_pm,c(ccn,ccp,ccd,cca))

ccn=formatci(ci.cvAUC(Ypred_noninv_dm_lassoproc,Y,fold=fold2))
ccp=formatci(ci.cvAUC(Ypred_preop_dm_lassoproc,Y,fold=fold2))
ccd=formatci(ci.cvAUC(Ypred_discharge_dm_lassoproc,Y,fold=fold2))
cca=formatci(ci.cvAUC(Ypred_all_dm_lassoproc,Y,fold=fold2))

tab_pm=cbind(tab_pm,c(ccn,ccp,ccd,cca))

ccn=formatci(ci.cvAUC(Ypred_noninv_dm_rfproc,Y,fold=fold2))
ccp=formatci(ci.cvAUC(Ypred_preop_dm_rfproc,Y,fold=fold2))
ccd=formatci(ci.cvAUC(Ypred_discharge_dm_rfproc,Y,fold=fold2))
cca=formatci(ci.cvAUC(Ypred_all_dm_rfproc,Y,fold=fold2))

tab_pm=cbind(tab_pm,c(ccn,ccp,ccd,cca))

rownames(tab_pm)=c("NONINV","PREOP","DISCHARGE","ALL")
colnames(tab_pm)=c("Linear","Lasso","RF")

save(tab_pm,file=paste0("data/",prefix,"outputs_pm.RData"))







##**********************************************************************
## Predictive analysis for Y5M                                      ####
##**********************************************************************

Y=Y5M
Yt=Y5M_time
w1=which(!is.na(Y+Yt))
Y=Y[w1]; Yt=Yt[w1]

# Split observations into k=10 equally sized folds
nfold=10
fold1=rep(1:nfold,1+floor(length(Y)/nfold))[1:length(Y)]; fold1=fold1[order(runif(length(fold1)))]
fold2=rep(1:nfold,1+floor(length(Y)/nfold))[1:length(Y)]; fold2=fold2[order(runif(length(fold2)))]


nse=30 # run this many iterations to estimate bootstrapped standard error

predset=c("noninv","preop","discharge","all")

y5m_name=paste0("./data/",prefix,"y5m_outputs")

if (!file.exists(y5m_name) | (force_redo_5m|force_redo_all)) {

for (i in 1:4) {

pred=predset[i]
if (i==1) X=Xall[,intersect(colnames(Xall),noninv_predictors)]
if (i==2) X=Xall[,intersect(colnames(Xall),preop_predictors)]
if (i==3) X=Xall[,intersect(colnames(Xall),discharge_predictors)]
if (i==4) X=Xall[,intersect(colnames(Xall),predictors)]
X=X[w1,]

# Partial likeilhood-based measure of discretion
ydif=function(pred,Y,Yt) { # Harrell's concordance
      1-survConcordance(Surv(Yt,Y)~pred)$concordance
}

models=c("lrproc_surv","lassoproc_surv","rfproc_surv")
best=findbestmodel_surv(X,Y,Yt,fold1,models,"ydif") # best model
assign(paste0("best_",predset[i],"_5m"),best$best)
for (j in 1:length(models)) assign(paste0("bestpars",j),get(models[j])(X,Y,Yt,return_hyperparameters=T)) # best parameters for each model 

for (j in 1:length(models)) {
 yp=evaluateperformance_surv(X,Y,Yt,fold=fold2,models[j],get(paste0("bestpars",j)))
 assign(paste0("Ypred_",predset[i],"_5m_",models[j]),yp)
 
 outc=c()
 
 # Run nse bootstrap samples in order to estimate standard error
 for (k in 1:nse) {
  boot=sample(length(Y),length(Y),replace=T)
  yboot= evaluateperformance_surv(X[boot,],Y[boot],Yt[boot],fold=fold2,model=models[j],get(paste0("bestpars",j)))
  outc=c(outc,mean(sapply(1:max(fold2),
    function(x) ydif(yboot[which(fold2==x)],
                Y[boot][which(fold2==x)],
                Yt[boot][which(fold2==x)]))))
  }
  assign(paste0("sepred_",predset[i],"_5m_",models[j]),outc)
  print(j)
 }
print(i)

}

save(Ypred_noninv_5m_lrproc_surv, Ypred_noninv_5m_lassoproc_surv, Ypred_noninv_5m_rfproc_surv, Ypred_preop_5m_lrproc_surv, Ypred_preop_5m_lassoproc_surv, Ypred_preop_5m_rfproc_surv, Ypred_discharge_5m_lrproc_surv, Ypred_discharge_5m_lassoproc_surv, Ypred_discharge_5m_rfproc_surv, Ypred_all_5m_lrproc_surv, Ypred_all_5m_lassoproc_surv, Ypred_all_5m_rfproc_surv,sepred_noninv_5m_lrproc_surv, sepred_noninv_5m_lassoproc_surv, sepred_noninv_5m_rfproc_surv, sepred_preop_5m_lrproc_surv, sepred_preop_5m_lassoproc_surv, sepred_preop_5m_rfproc_surv, sepred_discharge_5m_lrproc_surv, sepred_discharge_5m_lassoproc_surv, sepred_discharge_5m_rfproc_surv, sepred_all_5m_lrproc_surv, sepred_all_5m_lassoproc_surv, sepred_all_5m_rfproc_surv,best_noninv_5m,best_preop_5m,best_discharge_5m,best_all_5m,fold=fold2,ydif,file=y5m_name)

} else load(y5m_name)


# plot each outcome as time-dependent ROC curve

tx=1826 # one-year survival 

pdf(paste0("./plots/",prefix,"x5m_lr.pdf"),width=4,height=4)
plot(0,type="n",xlim=c(0,1),ylim=c(0,1),xlab="1-Specificity",ylab="Sensitivity",xaxs="i",yaxs="i")
draw_roc_surv(Ypred_noninv_5m_lrproc_surv,Y,Yt,fold2,tx,add=T,col=2)
draw_roc_surv(Ypred_preop_5m_lrproc_surv,Y,Yt,fold2,tx,add=T,col=3)
draw_roc_surv(Ypred_discharge_5m_lrproc_surv,Y,Yt,fold2,tx,add=T,col=4)
draw_roc_surv(Ypred_all_5m_lrproc_surv,Y,Yt,fold2,tx,add=T,col=5)
abline(0,1)
legend("bottomright",c("Noninv.","Preop.","Discharge","All"),lty=1,col=2:5)
dev.off()

pdf(paste0("./plots/",prefix,"x5m_lasso.pdf"),width=4,height=4)
plot(0,type="n",xlim=c(0,1),ylim=c(0,1),xlab="1-Specificity",ylab="Sensitivity",xaxs="i",yaxs="i")
draw_roc_surv(Ypred_noninv_5m_lassoproc_surv,Y,Yt,fold2,tx,add=T,col=2)
draw_roc_surv(Ypred_preop_5m_lassoproc_surv,Y,Yt,fold2,tx,add=T,col=3)
draw_roc_surv(Ypred_discharge_5m_lassoproc_surv,Y,Yt,fold2,tx,add=T,col=4)
draw_roc_surv(Ypred_all_5m_lassoproc_surv,Y,Yt,fold2,tx,add=T,col=5)
abline(0,1)
legend("bottomright",c("Noninv.","Preop.","Discharge","All"),lty=1,col=2:5)
dev.off()

pdf(paste0("./plots/",prefix,"x5m_rf.pdf"),width=4,height=4)
plot(0,type="n",xlim=c(0,1),ylim=c(0,1),xlab="1-Specificity",ylab="Sensitivity",xaxs="i",yaxs="i")
draw_roc_surv(Ypred_noninv_5m_rfproc_surv,Y,Yt,fold2,tx,add=T,col=2)
draw_roc_surv(Ypred_preop_5m_rfproc_surv,Y,Yt,fold2,tx,add=T,col=3)
draw_roc_surv(Ypred_discharge_5m_rfproc_surv,Y,Yt,fold2,tx,add=T,col=4)
draw_roc_surv(Ypred_all_5m_rfproc_surv,Y,Yt,fold2,tx,add=T,col=5)
abline(0,1)
legend("bottomright",c("Noninv.","Preop.","Discharge","All"),lty=1,col=2:5)
dev.off()




# Save AUCs and confidence intervals as table

sets=c("noninv","preop","discharge","all")
mods=c("lrproc_surv","lassoproc_surv","rfproc_surv")
for (i in 1:3) for (j in 1:4) {
 pred=get(paste0("Ypred_",sets[j],"_5m_",mods[i]))
 ym=mean(sapply(1:max(fold2),function(x) ydif(pred[which(fold2==x)],Y[which(fold2==x)],Yt[which(fold2==x)])))
 yse=sd(get(paste0("sepred_",sets[j],"_5m_",mods[i])))
 ci=list(cvAUC=1-ym,ci=1-(ym - c(-1,1)*qnorm(0.05/2)*yse))
 assign(paste0("cipred_",sets[j],"_5m_",mods[i]),ci)
}
 
 
tab_5m=c()

ccn=formatci(cipred_noninv_5m_lrproc_surv)
ccp=formatci(cipred_preop_5m_lrproc_surv)
ccd=formatci(cipred_discharge_5m_lrproc_surv)
cca=formatci(cipred_all_5m_lrproc_surv)

tab_5m=cbind(tab_5m,c(ccn,ccp,ccd,cca))

ccn=formatci(cipred_noninv_5m_lassoproc_surv)
ccp=formatci(cipred_preop_5m_lassoproc_surv)
ccd=formatci(cipred_discharge_5m_lassoproc_surv)
cca=formatci(cipred_all_5m_lassoproc_surv)

tab_5m=cbind(tab_5m,c(ccn,ccp,ccd,cca))

ccn=formatci(cipred_noninv_5m_rfproc_surv)
ccp=formatci(cipred_preop_5m_rfproc_surv)
ccd=formatci(cipred_discharge_5m_rfproc_surv)
cca=formatci(cipred_all_5m_rfproc_surv)

tab_5m=cbind(tab_5m,c(ccn,ccp,ccd,cca))

rownames(tab_5m)=c("NONINV","PREOP","DISCHARGE","ALL")
colnames(tab_5m)=c("Linear","Lasso","RF")

save(tab_5m,file=paste0("data/",prefix,"outputs_5m.RData"))






##**********************************************************************
## Predictive analysis for YDQ                                      ####
##**********************************************************************

Y=YDC
w=which(is.finite(YDC))
Y=Y[w]

# Split observations into k=10 equally sized folds
nfold=10
fold1=rep(1:nfold,1+floor(length(Y)/nfold))[1:length(Y)]; fold1=fold1[order(runif(length(fold1)))]
fold2=rep(1:nfold,1+floor(length(Y)/nfold))[1:length(Y)]; fold2=fold2[order(runif(length(fold2)))]


morbidity_pars_fu=c("FU.activity","FU.symptom","FU.qol")
nse=30

predset=c("noninv","preop","discharge","all")

ydc_name=paste0("./data/",prefix,"dc_outputs")

if (!file.exists(ydc_name) | (force_redo_dq|force_redo_all)) {

for (i in 1:4) {

pred=predset[i]
if (i==1) X=Xall[w,intersect(colnames(Xall),noninv_predictors)]
if (i==2) X=Xall[w,intersect(colnames(Xall),preop_predictors)]
if (i==3) X=Xall[w,intersect(colnames(Xall),discharge_predictors)]
if (i==4) X=Xall[w,intersect(colnames(Xall),predictors)]

# Remove morbidity predictors
X=X[,setdiff(colnames(X),morbidity_pars_fu)]

# Function to determine accuracy of prediction
ydif=function(pred,Y) 1-cor(pred,Y,use="complete",method="spearman") # 'Difference' between predicted and expected Y

models=c("lrproc","lassoproc","rfproc")
best=findbestmodel(X,Y,fold1,models,"ydif")$best # best model
assign(paste0("best_",predset[i],"_dc"),best)
for (j in 1:length(models)) assign(paste0("bestpars",j),get(models[j])(X,Y,return_hyperparameters=T)) # best parameters for each model 

for (j in 1:length(models)) {
 yp=evaluateperformance(X,Y,fold=fold2,models[j],get(paste0("bestpars",j)))
 assign(paste0("Ypred_",predset[i],"_dc_",models[j]),yp)
 
 outc=c()
 
 # Run nse bootstrap samples in order to estimate standard error
 for (k in 1:nse) {
  boot=sample(length(Y),length(Y),replace=T)
  yboot= evaluateperformance(X[boot,],Y[boot],fold=fold2,model=models[j],get(paste0("bestpars",j)))
  outc=c(outc,mean(sapply(1:max(fold2),
    function(x) ydif(yboot[which(fold2==x)],
                Y[boot][which(fold2==x)]))))
  }
  assign(paste0("sepred_",predset[i],"_dc_",models[j]),outc)
  print(j)
 }
print(i)

}

save(Ypred_noninv_dc_lrproc,Ypred_noninv_dc_lassoproc,Ypred_noninv_dc_rfproc,Ypred_preop_dc_lrproc,Ypred_preop_dc_lassoproc,Ypred_preop_dc_rfproc,Ypred_discharge_dc_lrproc,Ypred_discharge_dc_lassoproc,Ypred_discharge_dc_rfproc,Ypred_all_dc_lrproc,Ypred_all_dc_lassoproc,Ypred_all_dc_rfproc,sepred_noninv_dc_lrproc,sepred_noninv_dc_lassoproc,sepred_noninv_dc_rfproc,sepred_preop_dc_lrproc,sepred_preop_dc_lassoproc,sepred_preop_dc_rfproc,sepred_discharge_dc_lrproc,sepred_discharge_dc_lassoproc,sepred_discharge_dc_rfproc,sepred_all_dc_lrproc,sepred_all_dc_lassoproc,sepred_all_dc_rfproc,best_noninv_dc,best_preop_dc,best_discharge_dc,best_all_dc,fold1,fold2,ydif, file=ydc_name)

} else load(ydc_name)


# plot each outcome as ROC curve
# LM
pdf(paste0("./plots/",prefix,"dq_lr.pdf"),width=4,height=4)
plot(0,type="n",xlim=c(0,1),ylim=c(0,1),xlab="1-Specificity",ylab="Sensitivity",xaxs="i",yaxs="i")
draw_roc_cv(Ypred_noninv_dc_lrproc,Y>0,fold=fold2,add=T,col=2)
draw_roc_cv(Ypred_preop_dc_lrproc,Y>0,fold=fold2,add=T,col=3)
draw_roc_cv(Ypred_discharge_dc_lrproc,Y>0,fold=fold2,add=T,col=4)
#draw_roc_cv(Ypred_all_dc_lrproc,Y>0,fold=fold2,add=T,col=5)
abline(0,1)
#legend("bottomright",c("Noninv.","Preop.","Discharge","All"),lty=1,col=2:5)
legend("bottomright",c("Noninv.","Preop.","Discharge"),lty=1,col=2:4)
dev.off()


# Lasso
pdf(paste0("./plots/",prefix,"dq_lasso.pdf"),width=4,height=4)
plot(0,type="n",xlim=c(0,1),ylim=c(0,1),xlab="1-Specificity",ylab="Sensitivity",xaxs="i",yaxs="i")
draw_roc_cv(Ypred_noninv_dc_lassoproc,Y>0,fold=fold2,add=T,col=2)
draw_roc_cv(Ypred_preop_dc_lassoproc,Y>0,fold=fold2,add=T,col=3)
draw_roc_cv(Ypred_discharge_dc_lassoproc,Y>0,fold=fold2,add=T,col=4)
#draw_roc_cv(Ypred_all_dc_lassoproc,Y>0,fold=fold2,add=T,col=5)
abline(0,1)
#legend("bottomright",c("Noninv.","Preop.","Discharge","All"),lty=1,col=2:5)
legend("bottomright",c("Noninv.","Preop.","Discharge"),lty=1,col=2:4)
dev.off()


# RF
pdf(paste0("./plots/",prefix,"dq_rf.pdf"),width=4,height=4)
plot(0,type="n",xlim=c(0,1),ylim=c(0,1),xlab="1-Specificity",ylab="Sensitivity",xaxs="i",yaxs="i")
draw_roc_cv(Ypred_noninv_dc_rfproc,Y>0,fold=fold2,add=T,col=2)
draw_roc_cv(Ypred_preop_dc_rfproc,Y>0,fold=fold2,add=T,col=3)
draw_roc_cv(Ypred_discharge_dc_rfproc,Y>0,fold=fold2,add=T,col=4)
#draw_roc_cv(Ypred_all_dc_rfproc,Y>0,fold=fold2,add=T,col=5)
abline(0,1)
#legend("bottomright",c("Noninv.","Preop.","Discharge","All"),lty=1,col=2:5)
legend("bottomright",c("Noninv.","Preop.","Discharge"),lty=1,col=2:4)
dev.off()





# Save AUCs and confidence intervals as table

sets=c("noninv","preop","discharge","all")
mods=c("lrproc","lassoproc","rfproc")
for (i in 1:3) for (j in 1:4) {
 pred=get(paste0("Ypred_",sets[j],"_dc_",mods[i]))
 ym=mean(sapply(1:max(fold2),function(x) ydif(pred[which(fold2==x)],Y[which(fold2==x)])))
 yse=sd(get(paste0("sepred_",sets[j],"_dc_",mods[i])))
 ci=list(cvAUC=1-ym,ci=1-(ym - c(-1,1)*qnorm(0.05/2)*yse))
 assign(paste0("cipred_",sets[j],"_dc_",mods[i]),ci)
}
 
 
tab_dc=c()

ccn=formatci(cipred_noninv_dc_lrproc)
ccp=formatci(cipred_preop_dc_lrproc)
ccd=formatci(cipred_discharge_dc_lrproc)
cca=formatci(cipred_all_dc_lrproc)

tab_dc=cbind(tab_dc,c(ccn,ccp,ccd,cca))

ccn=formatci(cipred_noninv_dc_lassoproc)
ccp=formatci(cipred_preop_dc_lassoproc)
ccd=formatci(cipred_discharge_dc_lassoproc)
cca=formatci(cipred_all_dc_lassoproc)

tab_dc=cbind(tab_dc,c(ccn,ccp,ccd,cca))

ccn=formatci(cipred_noninv_dc_rfproc)
ccp=formatci(cipred_preop_dc_rfproc)
ccd=formatci(cipred_discharge_dc_rfproc)
cca=formatci(cipred_all_dc_rfproc)

tab_dc=cbind(tab_dc,c(ccn,ccp,ccd,cca))

rownames(tab_dc)=c("NONINV","PREOP","DISCHARGE","ALL")
colnames(tab_dc)=c("Linear","Lasso","RF")

save(tab_dc,file=paste0("data/",prefix,"outputs_dc.RData"))





##**********************************************************************
## Final models and most important variables in RF                  ####
##**********************************************************************

set.seed(1)

predset=c("noninv","preop","discharge","all")
morbidity_pars_fu=c("FU.activity","FU.symptom","FU.qol")

for (i in 1:4) {
 if (i==1) X=Xall[,intersect(colnames(Xall),noninv_predictors)]
 if (i==2) X=Xall[,intersect(colnames(Xall),preop_predictors)]
 if (i==3) X=Xall[,intersect(colnames(Xall),discharge_predictors)]
 if (i==4) X=Xall[,intersect(colnames(Xall),predictors)]

 assign(paste0("mod_dm_",predset[i]),randomForest(factor(YDM)~.,data=X,subset=which(is.finite(YDM))))
 assign(paste0("mod_5m_",predset[i]),rfsrc(Surv(time = Y5M_time, event=Y5M)~.,data=cbind(X,Y5M,Y5M_time),subset=which(is.finite(Y5M+Y5M_time)), block.size = 1,importance=T))
 assign(paste0("mod_dq_",predset[i]),randomForest(YDC~.,data=X[,setdiff(colnames(X),morbidity_pars_fu)],subset=which(is.finite(YDC))))

 print(i)
}

save(mod_dm_noninv,mod_dm_preop,mod_dm_discharge,mod_dm_all,
	   mod_5m_noninv,mod_5m_preop,mod_5m_discharge,mod_5m_all,
     mod_dq_noninv,mod_dq_preop,mod_dq_discharge,mod_dq_all,
        file=paste0("data/",prefix,"final_models.RData"))

nn=20
m1=mod_dm_noninv$importance; m1=m1[order(-m1),]
m2=mod_dm_preop$importance; m2=m2[order(-m2),]
m3=mod_dm_discharge$importance; m3=m3[order(-m3),]
m4=mod_dm_all$importance; m4=m4[order(-m4),]
top10_dm=cbind(names(m1)[1:nn],m1[1:nn],names(m2)[1:nn],m2[1:nn],names(m3)[1:nn],m3[1:nn],names(m4)[1:nn],m4[1:nn])


m1=mod_5m_noninv$importance; m1=m1[order(-m1)]
m2=mod_5m_preop$importance; m2=m2[order(-m2)]
m3=mod_5m_discharge$importance; m3=m3[order(-m3)]
m4=mod_5m_all$importance; m4=m4[order(-m4)]
top10_5m=cbind(names(m1)[1:nn],m1[1:nn],names(m2)[1:nn],m2[1:nn],names(m3)[1:nn],m3[1:nn],names(m4)[1:nn],m4[1:nn])


m1=mod_dq_noninv$importance; m1=m1[order(-m1),]
m2=mod_dq_preop$importance; m2=m2[order(-m2),]
m3=mod_dq_discharge$importance; m3=m3[order(-m3),]
m4=mod_dq_all$importance; m4=m4[order(-m4),]
top10_dc=cbind(names(m1)[1:nn],m1[1:nn],names(m2)[1:nn],m2[1:nn],names(m3)[1:nn],m3[1:nn],names(m4)[1:nn],m4[1:nn])

rownames(top10_dm)=NULL
rownames(top10_5m)=NULL
rownames(top10_dc)=NULL

for (i in c(1,3,5,7)) {
top10_dm[,i]=fullnames[match(top10_dm[,i],fullnames[,1]),2]
top10_5m[,i]=fullnames[match(top10_5m[,i],fullnames[,1]),2]
top10_dc[,i]=fullnames[match(top10_dc[,i],fullnames[,1]),2]

top10_dm[,i+1]=signif(as.numeric(top10_dm[,i+1]),digits=3)
top10_5m[,i+1]=signif(as.numeric(top10_5m[,i+1]),digits=3)
top10_dc[,i+1]=signif(as.numeric(top10_dc[,i+1]),digits=3)
}

top10_dm

top10_5m

top10_dc


##**********************************************************************
## Comparison with EUROSCORE2                                       ####
##**********************************************************************

tab_raw=tab_raw[-exclude]

N=dim(tab_raw)[1]
cm=tab_raw$comorbid_other


# Coefficients
coef=read.csv("data/euroscore/coefficients.txt",sep=",",stringsAs=F,header=F)



#####################################
#### Infer EUROSCORE2 parameters ####
#####################################

## Simple parameters
e_age=tab_raw$age
nyha=tab_raw$nyha_bl; nyha[which(is.na(nyha))]=3
e_nyha2=(nyha==2)
e_nyha3=(nyha==3)
e_nyha4=(nyha==4)
e_ccs4=rep(0,N)
e_gender=(tab_raw$sex=="F")



## Creatinine clearance
bsa=tab_raw$bsa; bmi=tab_raw$BMI

# bsa=c1*(w^x1)*((h*100)^y1); BMI=w/(h^2)
c1=0.007184; x1=0.425; y1=0.725

ht=(bsa/(c1* (bmi^x1) * (100^y1)))^(1/(2*x1 + y1))
wt=bmi*(ht^2)
crcl=(140-e_age)*wt* (c(1,0.85)[1+e_gender])/(0.72*tab_raw$preop_creat)
# too many missing values at present

e_kd_50_85=rep(0,N) # Approximate those with renal impairment as 'poor' rather than 'moderate'
e_kd_dial=rep(0,N)
e_kd_50=rep(0,N)
e_kd_50[c(grep("kd",cm),grep("renal imp",cm))]=1



## LV function
lvi=intersect(union(grep("lv",cm),grep("LV",cm)),
              union(grep("impairment",cm),grep("dysfunction",cm)))
e_lv_mod=rep(0,N); e_lv_mod[lvi]=1
e_lv_poor=rep(0,N);
e_lv_vpoor=rep(0,N);
# LV function-- we do have data on this specifically from echo but only in a subset of patients. The mean was 62+/- 5. Almost everyone has a normal LVEF because it would be a contraindication to operating so imputation would be fine in this case as they should all be above 50


## Extracardiac arteriopathy
e_exa=rep(0,N)
exa=unique(c(grep("pvd",cm),grep("PVD",cm),grep("al vasc",cm)))
e_exa[exa]=1


## Recent MI
e_mi=rep(0,N)


## Poor mobility
e_mob=rep(0,N)
## 'severe impairment of mobility secondary to musculoskeletal or neurological dysfunction'


## Pulmonary hypertension
e_pah_31_55=rep(0,N)
e_pah_55=rep(0,N)
spap=tab_raw$spap_bl
e_pah_31_55[which(spap>31 & spap<55)]=1
e_pah_55[which(spap>=55)]=1
e_pah_55[which(is.na(spap))]=1 # impute


## Previous cardiac surgery
e_prevcs=rep(0,N)
prevcs=unique(c(setdiff(grep("valv",cm),grep("ndocard",cm)),grep("cabg",cm),grep("CABG",cm)))
e_prevcs[prevcs]=1

## Chronic lung disease
e_cld=rep(0,N)
cld=unique(c(grep("sthm",cm),grep("copd",cm),grep("COPD",cm)))
e_cld[cld]=1
# Chronic lung disease- we should have this in co-morbidities and PH does not count as "chronic lung disease"

## Emergency operation
e_urg=rep(0,N)
e_urg[which(tab_raw$emergency=="y")]=1
e_emerg=rep(0,N) # not done
e_salv=rep(0,N)


## Active endocarditis
e_endo=rep(0,N)
# Active endocarditis- another absolute contraindication so all will be no

## Weight of intervention
w3=(tab_raw$additional_cardiac>0.5)
e_weight_no_cabg=rep(0,N)
e_weight2=1-w3
e_weight3=w3

## Critical preoperative state
e_crit_preop=rep(0,N)

## Surgery on thoracic aorta
e_thor=rep(0,N)
e_thor[which(tab_raw$avr=="Y")]=1

## Diabetes
e_diab=rep(0,N)
e_diab[which(tab_raw$comorbid_dm=="y")]=1

## Constant term
e_const=rep(1,N)



#####################################
#### Compute score ##################
#####################################

euroscore=rep(0,N)
for (i in 1:dim(coef)[1]) {
    vals=as.numeric(get(coef[i,1]))
    vals[which(is.na(vals))]=mean(vals[which(!is.na(vals))])
    euroscore=euroscore + vals*as.numeric(coef[i,2])
}

#####################################
#### Draw ROC #######################
#####################################

# DM
pdf(paste0("./plots/",prefix,"euroscore2_roc.pdf"),width=4,height=4)
draw_roc(euroscore,YDM)
dev.off()

# 5M
pdf(paste0("./plots/",prefix,"euroscore2_roc_5m.pdf"),width=4,height=4)
Y=Y5M; Yt=Y5M_time; w1=which(!is.na(Y+Yt)); Y=Y[w1]; Yt=Yt[w1]; ex=euroscore[w1]; tx=1826 # one-year survival 
plot(0,type="n",xlim=c(0,1),ylim=c(0,1),xlab="1-Specificity",ylab="Sensitivity",xaxs="i",yaxs="i")
draw_roc_surv(euroscore,Y,Yt,rep(1,length(Y)),tx,add=T,col=1); abline(0,1,col="red")
dev.off()

# Lasso
pdf(paste0("./plots/",prefix,"euroscore2_roc_dq.pdf"),width=4,height=4)
draw_roc(euroscore,YDC>0)
dev.off()



#####################################
#### Save outputs for Shiny##########
#####################################

# Preop data (preliminary)
Xp=Xall[,intersect(colnames(Xall),preop_predictors)]


# Mean, best 20%, worst 20% survival 
# Preliminaries
survfit=predict(mod_5m_preop,Xp)
survint=rowSums(survfit$survival)
i_median=which.min(abs(survint-median(survint)))
i_20=which.min(abs(survint-quantile(survint,0.2)))
i_80=which.min(abs(survint-quantile(survint,0.8)))

survtime=survfit$time.interest
survmean=survfit$survival[i_median,]
surv20=survfit$survival[i_20,]
surv80=survfit$survival[i_80,]

# Mean value of all variables
mtab=colMeans(Xp)

# Logistic regression coefficients
cx=coxph(Surv(time=Y5M_time,event=Y5M)~.,data=Xp)
lr_coefficients=cx$coefficients

# Camphor scores
cph_fu=Xall$FU.symptom + Xall$FU.qol + Xall$FU.activity
cph_bl=Xall$BL.Symptom + Xall$FU.QoL + Xall$FU.Activity


save(mtab,survtime,survmean,surv20,surv80,mod_dm_preop,mod_5m_preop,mod_dq_preop,
  lr_coefficients,cph_fu,cph_bl,file=paste0("./data/data.RData"))





#####################################
#### Validation #####################
#####################################


set.seed(327462)
vtab=vtab0
Yc=vtab0$surgeon_predicted_mortality_mean

w=which((vtab$excluded==1))
vtab=vtab[-w,]
Yc=Yc[-w]

ds=(as.Date(vtab$Death_date,format="%d/%m/%Y") - as.Date(vtab$pea_date,format="%d/%m/%Y"))

nm=apply(vtab,2,function(x) length(which(is.na(x))))
subm=colnames(vtab)#[which(nm<250)]

xcol=intersect(intersect(subm,colnames(Xall)),preop_predictors)


ds[which(!is.finite(ds))]=500
Yv=(ds<35)
for (i in 1:ncol(vtab)) {
  if (colnames(vtab)[i] %in% names(mtab)) {
    vtab[which(is.na(vtab[,i])),i]=mtab[colnames(vtab)[i]]
  }
}

vtab=vtab[,xcol]
Xsub=Xall[,xcol]

train_dat=cbind(Xsub,Y=as.factor(Y))
nt=round(dim(train_dat)[1]/2)
train=1:nt; test=(nt+1):dim(train_dat)[1]
m1=randomForest(Y~.,data=train_dat[train,])
m2=randomForest(Y~.,data=train_dat)
m3=glm(Y~.,data=train_dat,family=binomial(link="logit"))

par(mfrow=c(1,2))

Yt=Y[test]
Ypt=predict(m1,train_dat[test,],type="prob")[,1]
rt=roc(Yt,Ypt)
plot(rt)

Ypv1=predict(m1,vtab,type="prob")[,1]
Ypv2=predict(m2,vtab,type="prob")[,1]
Ypv3=predict(m3,vtab,type="response")
rv1=roc(Yv,Ypv1)
rv2=roc(Yv,Ypv2)
rv3=roc(Yv,Ypv3)
rvc=roc(Yv,Yc)
plot(rv2)
lines(rv3,col="blue")
lines(rvc,col="red")
rv$auc

print(wilcox.test(Ypv2[which(Yv==1)],Ypv2[which(Yv==0)]))
print(wilcox.test(Ypv3[which(Yv==1)],Ypv3[which(Yv==0)]))
print(wilcox.test(Yc[which(Yv==1)],Yc[which(Yv==0)]))



#####################################
#### Save image #####################
#####################################

save.image(paste0("./data/",prefix,"all_outputs_us.RData"))

