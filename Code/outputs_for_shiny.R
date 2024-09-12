##**********************************************************************
##  Risk prediction in PEA                                          ####
##  Supporting code for main analyses                               
##  James Liley                                                     
##  9/3/19                                                          
##**********************************************************************


##**********************************************************************
## Set data and output directories                                  ####
##**********************************************************************

## Data location
datadir="../../Data/" 

# Output location
outdir="Outputs/"


##**********************************************************************
## Load data from main analysis                                     ####
##**********************************************************************

load(paste0(datadir,"Temp_data/Workspaces/discovery.RData"))


##**********************************************************************
## Process outputs for Shiny                                        ####
##**********************************************************************

# Preop data (preliminary)
Xp=Xall[,intersect(colnames(Xall),preop_predictors)]


# Mean, best 20%, worst 20% survival 
survfit=predict(mod_5m_preop,Xp)
survtime=survfit$time.interest
survmean=colMedians(survfit$survival) 
surv20=apply(survfit$survival,2,function(x) quantile(x,0.2)) 
surv80=apply(survfit$survival,2,function(x) quantile(x,0.8)) 

# 'Typical' patient; this is not necessarily median (who is healthier than typical) or mean (who is less healthy). 
dm=rowSums((survfit$survival-survmean)^2)
ix=order(dm)[1:10]
set.seed(543643) # This is not from key-pounding - it was searched for to generate a hypothetical patient whose survival profile matched the median.
ttab=apply(Xp[ix,],2,function(x) sample(x,1))

# Mean values
mtab=colMeans(Xp,na.rm=T)

# Logistic regression coefficients
cx=coxph(Surv(time=Y5M_time,event=Y5M)~.,data=Xp)
lr_coefficients=cx$coefficients

# Camphor scores
cph_fu=Xall$FU.symptom + Xall$FU.qol + Xall$FU.activity
cph_bl=Xall$BL.Symptom + Xall$BL.QoL + Xall$BL.Activity


##**********************************************************************
## Save outputs                                                     ####
##**********************************************************************

shiny_data=paste0("../../Shiny/data.RData")
save(mtab,ttab,survtime,survmean,surv20,surv80,mod_dm_preop,mod_5m_preop,mod_dq_preop,
     lr_coefficients,cph_fu,cph_bl,file=shiny_data)



