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


##**********************************************************************
## Save outputs                                                     ####
##**********************************************************************

shiny_data=paste0("../../Shiny/data.RData")
save(mtab,survtime,survmean,surv20,surv80,mod_dm_preop,mod_5m_preop,mod_dq_preop,
     lr_coefficients,cph_fu,cph_bl,file=shiny_data)



