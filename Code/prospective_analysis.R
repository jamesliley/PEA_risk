##**********************************************************************
##  Risk prediction in PEA                                          ####
##  Supporting code to analyse prospective data table                      
##  James Liley                                                         
##  05/10/2022                                                              
##**********************************************************************




##**********************************************************************
## Set data and output directories                                  ####
##**********************************************************************

## Data location
datadir="../../Data/" 

# Output location
outdir="Outputs/"

##**********************************************************************
## Packages and scripts                                             ####
##**********************************************************************

library(glmnet)
library(randomForest)
library(pROC)

source("Code/auxiliary.R") # Auxiliary functions


outfile=paste0(outdir,"summary_prospective.txt")
sink(outfile)


##**********************************************************************
## 3. PROSPECTIVE VALIDATION                                        ####
## Read datasets                                                    ####
##**********************************************************************

# Prospective data
load(paste0(datadir,"Temp_data/Design/prospective_table.RData"))
vtab0=vtab

# Discovery
load(paste0(datadir,"Temp_data/Workspaces/discovery.RData"))

##**********************************************************************
## Sample size calculation setup                                    ####
##**********************************************************************

# Training data
Y=YDM
train_dat=cbind(Xall,Y=as.factor(Y))


nrs=c(50,500) # Maximum and minimum sizes to test
nsim=100 # Simulate this many times

# Store simulation data in a list
sim_data=list()


##**********************************************************************
## Sample size calculation                                          ####
##**********************************************************************

# Random seed
set.seed(3625362)

# Function to run test for one potential sample size
test_n=function(nr,verbose=TRUE) {
  
  npoll=rep(NA,nsim)
  for (ss in 1:nsim) {
    
    # Take bootstrap sample 
    sboot=sample(1:dim(train_dat)[1],nr,rep=TRUE);
    cboot=setdiff(1:dim(train_dat)[1],sboot)
    
    # Simulated validation data and discovery data
    sim_val=train_dat[sboot,]
    sim_disc=train_dat[cboot,]
    
    # Variables with missingness in simulated validation sample <0.5
    nm_val=apply(sim_val,2,function(x) length(which(is.na(x))))/dim(sim_val)[1]
    subm_val=colnames(sim_val)[which(nm_val < 0.5)]
    
    # Intersection with pre-operative predictors
    xcol=intersect(intersect(subm_val,colnames(Xall)),preop_predictors)
    
    # Mean-value impute (according to original dataset)
    for (i in 1:ncol(sim_val)) {
      if (colnames(sim_val)[i] %in% names(mtab)) {
        sim_val[which(is.na(sim_val[,i])),i]=mtab[colnames(sim_val)[i]]
        sim_disc[which(is.na(sim_disc[,i])),i]=mtab[colnames(sim_disc)[i]]
      }
    }
    
    # Restrict to available predictors
    sim_val=sim_val[,c(xcol,"Y")]
    sim_disc=sim_disc[,c(xcol,"Y")]
    
    # Fit RF
    nt=round(dim(sim_disc)[1]/2)
    m_sim=randomForest(Y~.,data=sim_disc)
    
    # Predictions on simulated validation sample
    Ypred_sim=predict(m_sim,sim_val,type="prob")[,1]
    
    # Wilcoxon test - does the model perform better than randomly?
    if (length(which(sim_val$Y==1))>5) {
      wt_sim=wilcox.test(Ypred_sim[which(sim_val$Y==1)],Ypred_sim[which(sim_val$Y==0)])
      npoll[ss]=wt_sim$p.value
    } else npoll[ss]=1
    
    if (verbose) print(paste0("Completed ",ss," of ",nsim," trials."))
    
  }
  if (verbose) print(paste0("Power estimated at sample size: ",nr," is ",mean(npoll<0.05)))
  return(npoll)
}



# Run bisection method to find sample size for 90% power
sim_file="Reference/simdata_pea_risk.RData"

if (!file.exists(sim_file)) {
  power_threshold=0.9
  lower=nrs[1]
  upper=nrs[2]
  p_lower=test_n(lower)
  p_upper=test_n(upper)
  nx=c(lower,upper)
  sim_data[[1]]=p_lower
  sim_data[[2]]=p_upper
  s_index=3
  while(upper-lower > 1) {
    xmid=round((upper + lower)/2)
    pmid=test_n(xmid)
    xp=mean(pmid<0.05)
    if (xp < power_threshold) lower=xmid else upper=xmid
    sim_data[[s_index]]=pmid
    nx[s_index]=xmid
    s_index=s_index+1
  }
  
  # Save
  names(sim_data)=nx
  save(sim_data,file=sim_file)
  
} else load(sim_file)

# Print
cat("\n\n\n")
print (paste0("A sample size of ",names(sim_data)[length(sim_data)]," is sufficient to have ",
              round(power_threshold*100), "% power to reject the null hypothesis that prediction ",
              "of YDM from the discovery dataset is no better than random, assuming that individuals ",
              "in a prospective validation sample had identical distributions of covariates"))
cat("\n\n\n")


##**********************************************************************
## Set up dataset                                                   ####
##**********************************************************************

set.seed(327462)
Yc=vtab$surgeon_predicted_mortality_mean

# Remove exclusions
w=which((vtab$excluded==1|is.na(Yc)))
vtab=vtab[-w,]
Yc=Yc[-w]

# Establish time since PEA, if died.
ds=(as.Date(vtab$Death_date,format="%d/%m/%Y") - as.Date(vtab$pea_date,format="%d/%m/%Y"))
ds[which(!is.finite(ds))]=500
Yv=(ds<35) # Died <35 days after PEA

# Find number of missing values per column, and remove variables with  missingness > 0.5
nm=apply(vtab,2,function(x) length(which(is.na(x))))/dim(vtab)[1]
subm=colnames(vtab)[which(nm < 0.5)]

# Pre-operative predictors
xcol=intersect(intersect(subm,colnames(Xall)),preop_predictors)

# Mean-value impute (according to original dataset)
for (i in 1:ncol(vtab)) {
  if (colnames(vtab)[i] %in% names(mtab)) {
    vtab[which(is.na(vtab[,i])),i]=mtab[colnames(vtab)[i]]
  }
}

# Restrict to available predictors
vtab=vtab[,xcol]
Xsub=Xall[,xcol]




##**********************************************************************
## Train model to original data using available variables           ####
##**********************************************************************

Y=YDM
train_dat=cbind(Xsub,Y=as.factor(Y))
nt=round(dim(train_dat)[1]/2)
train=1:nt; test=(nt+1):dim(train_dat)[1]
m1=randomForest(Y~.,data=train_dat[train,])
m2=randomForest(Y~.,data=train_dat)

Yt=Y[test]
Ypt=predict(m1,train_dat[test,],type="prob")[,1]

# Internal ROC for model prediction
rt=roc(Yt,Ypt)



##**********************************************************************
## Compute and compare AUCs of model vs surgeon prediction          ####
##**********************************************************************

# Predictions on prospective data
Ypv=predict(m2,vtab,type="prob")[,2]


# ROC for model on prospective data
rv=roc(Yv,Ypv)


# ROC for surgical prediction on prospective data
rvc=roc(Yv,Yc)


# Comparison of whether prediction is better than random in each case
cat("\n\n","Wilcoxon rank-sum test comparing model-derived predictions for individuals who died post-PEA with predictions for individuals who did not\n")
print(wilcox.test(Ypv[which(Yv==1)],Ypv[which(Yv==0)]))

cat("\n\n","Wilcoxon rank-sum test comparing surgeon-estimated predictions for individuals who died post-PEA with predictions for individuals who did not\n")
print(wilcox.test(Yc[which(Yv==1)],Yc[which(Yv==0)]))

cat("\n\n","Permutation test comparing area-under-ROC curve using model-derived predictions and using surgeon-derived predictions\n")
print(roc.test(rv,rvc))

## AUCs, standard errors and confidence intervals
auc_mod=rv$auc
se_mod=aucse(sum(Yv),sum(1-Yv),rv$auc)
ci_mod=auc_mod - qnorm(0.05/2)*se_mod*c(-1,1)

auc_surg=rvc$auc
se_surg=aucse(sum(Yv),sum(1-Yv),rvc$auc)
ci_surg=auc_surg - qnorm(0.05/2)*se_surg*c(-1,1)


##**********************************************************************
## Plot ROC and calibration curves                                  ####
##**********************************************************************


# ROC
pdf(paste0(outdir,"Plots/",prefix,"prospective_roc_surg_vs_mod.pdf"),width=4,height=4)

plot(0,type="n",xlim=c(0,1),ylim=c(0,1),
     xlab="Specificity",ylab="Sensitivity")
abline(0,1,lty=2)
draw_roc(Ypv,Yv,add=T,col="black",lwd=2)
draw_roc(Yc,Yv,add=T,col="red",lwd=2)
legend("bottomright",
       c(paste0("Model-derived: ",signif(auc_mod,digits=2),
                " (",signif(ci_mod[1],digits=2),", ", signif(ci_mod[2],digits=2),")"),
         paste0("Surgeon-predicted: ",signif(auc_surg,digits=2),
                  " (",signif(ci_surg[1],digits=2),", ", signif(ci_surg[2],digits=2),")")),
       col=c("black","red"),lwd=2,bty="n",cex=0.5)

dev.off()

# Calibration

pdf(paste0(outdir,"Plots/prospective_calibration_surg_vs_mod.pdf"),width=4,height=4)
plot(0,type="n",xlim=c(0,0.5),ylim=c(0,1.1),xlab="Predicted risk",ylab="Observed risk",xaxs="i",yaxs="i")
abline(0,1)
draw_calibration(Ypv,Yv,nc=20,type="l",add=T,col="black",lwd=2)
draw_calibration(Yc/100,Yv,nc=20,type="l",add=T,col="red",lwd=2)
legend("topright",c("Model_derived","Surgeon-predicted"),lwd=2,col=c("black","red"),cex=0.5)
dev.off()

sink()


##**********************************************************************
## Save image                                                       ####
##**********************************************************************

save.image(file=paste0(datadir,"Temp_data/Workspaces/prospective.RData"))
