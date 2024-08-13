##**********************************************************************
##  Risk prediction in PEA                                          ####
##  Supporting code to generate data table                              
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
## Packages and scripts                                             ####
##**********************************************************************

library(snpStats)
library(matrixStats)
source("Code/predictor_lists.R") # Lists of predictors
source("Code/auxiliary.R") # Auxiliary functions


##**********************************************************************
## Read raw data and pre-format (discovery)                         ####
##**********************************************************************

# Main dataset
tab=read.csv(paste0(datadir,"Original/Discovery/LHDmaster17.04.19.csv"),sep=",",stringsAsFactors=FALSE) # raw data

# US data
t2=read.table(paste0(datadir,"Original/Discovery/PEA_echos_james.csv"),header=T,stringsAs=F,sep=",")
tab=cbind(tab,t2[,4:(dim(t2)[2]-1)])

# Details of predictors
predictor_details=read.table(paste0("Reference/predictor_details.txt"),stringsAs=F,sep=",") # details of predictors


##**********************************************************************
## Basic formatting                                                 ####
##**********************************************************************

# Fix NAs
tab[which(tab=="NA")]=NA

# Format dates
datex=colnames(tab)[which(grepl("date",colnames(tab)))]
for (i in 1:length(datex)) tab[,datex[i]]=as.Date(tab[,datex[i]],"%d/%m/%Y")
tab$echo_bl=as.Date(tab$echo_bl,"%d/%m/%Y")
tab$echo_fu1=as.Date(tab$echo_fu1,"%d/%m/%Y")

# Remove useless variables
tab=tab[,setdiff(colnames(tab),c("preop_billirubin","comorbid_malig_comments","X","X.1"))]

# Raw version with text as is
tab_raw=tab

# Text codings
tab$sex=as.numeric(tab$sex=="M")
tab[,"comorbid_other"]=as.numeric(!is.na(tab[,"comorbid_other"]))
cmty=c("cabg","pfo","avr","mvr","asd",colnames(tab)[grep("comorbid",colnames(tab))])
tab[,cmty]=as.numeric(tab[,cmty]=="y"|tab[,cmty]=="Y")
tab$smoking_status=4-as.numeric(as.factor(tab$smoking_status))
tab$comp_reperfusion=(tab$comp_reperfusion=="yes")
tab$death_in_hosp=(tab$death_in_hosp=="Y")
tab$la_bl_dilated=(tab$la_bl_dilated!="normal")

abl=c("ra_area_bl","ra_area_fu1","la_4c_area_bl","la_4c_area_fu1"); 
vbl=c("ra_vol_bl","ra_vol_fu1","la_4c_vol_bl","la_4c_vol_fu1")

for (i in 1:4) {
  
  ax=abl[i]; vx=vbl[i]
  
  nm=which(tab[,ax]=="normal")
  if (i %in% 1:2) tab[,ax][nm]=c(16.2,15.2)[1+(tab$sex)][nm]
  if (i %in% 3:4) tab[,ax][nm]=c(22,20)[1+(tab$sex)][nm]
  suppressWarnings(tab[,ax]<-as.numeric(tab[,ax]))
  
  # Substitute right atrial areas
  yra=tab[,ax]; xra=tab[,vx]; 
  w1=which(!is.na(xra+yra)); gxx=lm(yra~xra,subset=w1); 
  Fra=function(x) gxx$coefficients[1] + gxx$coefficients[2]*x
  wxra=which(is.na(yra) & !is.na(xra))
  tab[wxra,ax]=Fra(tab[wxra,vx])
  
}

# Ultrasound
tab$afrhythm=as.numeric((tab$rhythmn=="af") + (1-(tab$rhythmn %in% c("sinus","sinu", "sinus "))))
tab$pacedrhythm=as.numeric((tab$rhythmn=="paced")+ (1-(tab$rhythmn %in% c("sinus","sinu", "sinus "))))

tab$rv_dilatation=as.numeric(c(2,2,3,1,1,4,4)[as.numeric(factor(t2$rv_dilatation,levels=c("mild","mildly","moderate","normal","normal ","secere","severe")))])
tab$collapse...50.=as.numeric(tab$collapse...50.=="yes")
tab$effusion=as.numeric(tab$effusion=="yes")
tab$rv_systolic_function=as.numeric(c(2,3,1,1,4)[as.numeric(factor(t2$rv_systolic_function,levels=c("mild","moderate","noraml","normal","severe")))])
tab$IVC.diam=suppressWarnings(as.numeric(tab$IVC.diam))
tab$MR=c(2,3,1,2)[as.numeric(factor(t2$MR,levels=c("mild","moderate","none","trace")))]
tab$AS=c(2,3,3,1,1)[as.numeric(factor(t2$AS,levels=c("mild","moderare","moderate","none","noen")))]
tab$AR=c(2,3,1,1)[as.numeric(factor(t2$AR,levels=c("mild","moderate","none","none ")))]



##**********************************************************************
## Affix medication data                                            ####
##**********************************************************************


med_data_file=paste0(datadir,"Temp_data/Working/meds_by_patient.RData")
if (!file.exists(med_data_file)) {
  
  rx=read.csv(paste0(datadir,"Original/Discovery/pea_meds.csv"),header=T,stringsAs=F)
  
  dt=read.csv(paste0("Reference/drug_type.txt"),stringsAs=F)
  dt[,2]=trimws(dt[,2])
  dtype=unique(dt[,2]); 
  for (i in 1:length(dtype)) if (nchar(trimws(dtype[i]))>0) assign(trimws(dtype[i]),dt[which(dt[,2]==dtype[i]),1])
  dtype=dtype[which(nchar(dtype)>0)]
  
  outx=data.frame(HospNo=rx$HospNo,med_diuretic=0,med_vasodilator=0,
                  med_beta_block=0,med_acei=0,med_ca_chan_block=0,med_digoxin=0,med_amiodarone=0)
  
  for (i in 1:dim(outx)[1]) {
    dlist=rx[i,2]
    if (is.na(dlist)) outx[i,2:dim(outx)[2]]=NA else
      outx[i,2:dim(outx)[2]] = c(
        grepl(paste(diuretic,collapse="|"),dlist,fixed=F),
        grepl(paste(vasodilator,collapse="|"),dlist,fixed=F),
        grepl(paste(beta,collapse="|"),dlist,fixed=F),
        grepl(paste(acei,collapse="|"),dlist,fixed=F),
        grepl(paste(cacb,collapse="|"),dlist,fixed=F),
        grepl(paste(digoxin,collapse="|"),dlist,fixed=F),	 
        grepl(paste(amiodarone,collapse="|"),dlist,fixed=F))
  }
  
  # Additional vasodilator details
  vaso=read.csv(paste0(datadir,"Original/Discovery/pea_vaso_meds.csv"),header=T,stringsAs=F)
  vaso$med_vasodilator=rep(NA,dim(vaso)[1])
  vaso$med_vasodilator[which(!is.na(vaso[,2]) & (vaso[,2]!="n"))]=1
  vaso$med_vasodilator[which((vaso[,2]=="n"))]=0
  
  outx[match(vaso$HospNo,outx$HospNo),]$med_vasodilator=vaso$med_vasodilator
  
  save(outx,file=med_data_file)
  
} else load(med_data_file)

# affix totab
tab=cbind(tab,med_diuretic=NA,med_vasodilator=NA,
          med_beta_block=NA,med_acei=NA,med_ca_chan_block=NA,med_digoxin=NA,med_amiodarone=NA)
ms=intersect(outx$HospNo,tab$HospNo); mx=match(ms,tab$HospNo); my=match(ms,outx$HospNo)
for (tn in colnames(outx)[2:ncol(outx)]) tab[mx,tn]=outx[my,tn]



##**********************************************************************
## Combination variables                                            ####
##**********************************************************************

any_comorbidity=as.numeric(rowMaxs(as.matrix(tab[,cmty]),na.rm=T)>0)

any_medication=as.numeric(rowMaxs(as.matrix(tab[,colnames(tab)[grep("med",colnames(tab))]]),na.rm=T)>0)

tab$any_comorbidity=any_comorbidity
tab$any_medication=any_medication


##**********************************************************************
## Genetics. Requires plink software on system                      ####
##**********************************************************************

# PRSs from http://www.broadcvdi.org/informational/data

prs_file=paste0(datadir,"Temp_data/PRS/prs_scores.RData")
if (!file.exists(prs_file)) {
  
  # Genotype location
  input_loc=paste0(datadir,"Original/Discovery/Genetics/2017-01-batch123_imputed_update.info_0.5.maf_0.01")
  
  # Read fam and bim files
  bimx=read.table(paste0(input_loc,".bim"),stringsAs=F,header=F)
  famx=read.table(paste0(input_loc,".fam"),stringsAs=F,header=F)
  
  
  # Get list of samples with genotype ID
  lk=read.csv(paste0(datadir,"Original/Discovery/pheno_all_6.csv"),stringsAs=F, header=T)
  iid=rep("",dim(tab)[1])
  ww=which(tab$HospNo %in% lk$sample_id_4)
  iid[ww]=lk[match(tab$HospNo[ww],lk$sample_id_4),]$IID
  pts=iid[ww]
  
  # Read GRS coefficients
  afib=read.table(paste0(datadir,"PRS/afib.txt"),skip=10,stringsAs=F,header=T)
  cad=read.table(paste0(datadir,"PRS/cad.txt"),skip=10,stringsAs=F,header=T)
  t2d=read.table(paste0(datadir,"PRS/t2d.txt"),skip=10,stringsAs=F,header=T)
  ibd=read.table(paste0(datadir,"PRS/ibd.txt"),skip=10,stringsAs=F,header=T)
  
  # Indices for GRS coefficients
  afib_id=paste0(afib$chr,":",afib$position); rownames(afib)=afib_id
  cad_id=paste0(cad$chr,":",cad$position); rownames(cad)=cad_id
  t2d_id=paste0(t2d$chr,":",t2d$position); rownames(t2d)=t2d_id
  ibd_id=paste0(ibd$chr,":",ibd$position); rownames(ibd)=ibd_id
  
  # Top 
  ntop=10000
  afibx=order(-afib$effect_weight)[1:ntop]; afibs=paste0(afib$chr[afibx],":",afib$position[afibx])
  cadx=order(-cad$effect_weight)[1:ntop]; cads=paste0(cad$chr[cadx],":",cad$position[cadx])
  t2dx=order(-t2d$effect_weight)[1:ntop]; t2ds=paste0(t2d$chr[t2dx],":",t2d$position[t2dx])
  ibdx=order(-ibd$effect_weight)[1:ntop]; ibds=paste0(ibd$chr[ibdx],":",ibd$position[ibdx])
  grsnp=unique(c(afibs,cads,t2ds,ibds))
  
  
  # Extract relevant part of genotype file
  subsnp=intersect(grsnp,bimx[,2])
  subpt_iid=intersect(pts,famx[,2])
  subpt_fid=famx[match(subpt_iid,famx[,2]),1]
  prs_snps_file=paste0(datadir,"Temp_data/PRS/prs_snps.txt"); 
  write(subsnp,file=prs_snps_file)
  
  prs_pts_file=paste0(datadir,"Temp_data/PRS/prs_pts.txt")
  write.table(cbind(subpt_fid,subpt_iid),file=prs_pts_file,quote=F,row.names=F,col.names=F)
  system(paste0(
    "plink --noweb ",
    " --bfile ",input_loc,
    " --extract ",prs_snps_file," --keep ",prs_pts_file,
    " --recodeA"))
  temp_genotypes_file=paste0(datadir,"Temp_data/PRS/prs_genotypes.txt")
  system(paste0("mv plink.raw ",temp_genotypes_file))

  # Read back in
  snpx=read.table(paste0(datadir,"prs_genotypes.txt"),stringsAs=F,header=T)
  ptx=snpx$IID
  snpx=snpx[,7:dim(snpx)[2]]
  
  rsnp=gsub("X","",colnames(snpx)[1:dim(snpx)[2]])
  rsnp=gsub(".",":",rsnp,fixed=T)
  rsnp=gsub("_[ATCG]","",rsnp)
  
  ra1=unlist(lapply(strsplit(colnames(snpx),"_"),function(x) x[2]))
  
  
  # Normalise genotype matrix
  csnp=colMeans(snpx,na.rm=T)
  csd=sqrt(csnp*(1-csnp/2))
  snpx=t(t(snpx)/csd)
  
  
  
  # Work out PRSs
  ix_afib=intersect(afibs,rsnp); mix_afib=match(ix_afib,rsnp)
  coef_afib=afib[ix_afib,]$effect_weight
  afa1=afib[ix_afib,]$A1; afa2=afib[ix_afib,]$A2
  sign_afib=rep(0,length(ix_afib));
  sign_afib[which(afa1==ra1[mix_afib])]=1; sign_afib[which(afa2==ra1[mix_afib])]=-1
  afib_score= as.matrix(snpx[,mix_afib]) %*% (coef_afib*sign_afib); names(afib_score)=ptx
  
  
  ix_cad=intersect(cads,rsnp); mix_cad=match(ix_cad,rsnp)
  coef_cad=cad[ix_cad,]$effect_weight
  afa1=cad[ix_cad,]$A1; afa2=cad[ix_cad,]$A2
  sign_cad=rep(0,length(ix_cad));
  sign_cad[which(afa1==ra1[mix_cad])]=1; sign_cad[which(afa2==ra1[mix_cad])]=-1
  cad_score= as.matrix(snpx[,mix_cad]) %*% (coef_cad*sign_cad); names(cad_score)=ptx
  
  
  
  ix_t2d=intersect(t2ds,rsnp); mix_t2d=match(ix_t2d,rsnp)
  coef_t2d=t2d[ix_t2d,]$effect_weight
  afa1=t2d[ix_t2d,]$A1; afa2=t2d[ix_t2d,]$A2
  sign_t2d=rep(0,length(ix_t2d));
  sign_t2d[which(afa1==ra1[mix_t2d])]=1; sign_t2d[which(afa2==ra1[mix_t2d])]=-1
  t2d_score= as.matrix(snpx[,mix_t2d]) %*% (coef_t2d*sign_t2d); names(t2d_score)=ptx
  
  
  
  ix_ibd=intersect(ibds,rsnp); mix_ibd=match(ix_ibd,rsnp)
  coef_ibd=ibd[ix_ibd,]$effect_weight
  afa1=ibd[ix_ibd,]$A1; afa2=ibd[ix_ibd,]$A2
  sign_ibd=rep(0,length(ix_ibd));
  sign_ibd[which(afa1==ra1[mix_ibd])]=1; sign_ibd[which(afa2==ra1[mix_ibd])]=-1
  ibd_score= as.matrix(snpx[,mix_ibd]) %*% (coef_ibd*sign_ibd); names(ibd_score)=ptx
  
  
  prs_details_file=paste0(datadir,"Temp_data/PRS/prs_details.RData")
  save(snpx,ptx,rsnp,ra1,
       coef_afib,sign_afib,
       coef_cad,sign_cad,
       coef_t2d,sign_t2d,
       coef_ibd,sign_ibd,
       ix_afib,ix_cad,ix_t2d,ix_ibd,
       afib_score,cad_score,t2d_score,ibd_score,
       file=prs_details_file)
  
  lk=read.csv(paste0(datadir,"Original/Discovery/pheno_all_6.csv"),stringsAs=F, header=T)
  iid=rep("",dim(tab)[1])
  ww=which(tab$HospNo %in% lk$sample_id_4)
  iid[ww]=lk[match(tab$HospNo[ww],lk$sample_id_4),]$IID
  
  prs_afib=rep(NA,dim(tab)[1])
  prs_cad=rep(NA,dim(tab)[1])
  prs_t2d=rep(NA,dim(tab)[1])
  prs_ibd=rep(NA,dim(tab)[1])
  
  wxx=which(ptx %in% iid)
  mm=match(ptx[wxx],iid)
  prs_afib[mm]=afib_score[wxx]
  prs_cad[mm]=cad_score[wxx]
  prs_t2d[mm]=t2d_score[wxx]
  prs_ibd[mm]=ibd_score[wxx]
  
  save(iid,prs_afib,prs_cad,prs_t2d,prs_ibd,file=prs_file)
  
} else load(prs_file)

tab$iid=iid
tab$prs_afib=prs_afib
tab$prs_cad=prs_cad
tab$prs_t2d=prs_t2d
tab$prs_ibd=prs_ibd




##**********************************************************************
## Final formatting                                                 ####
##**********************************************************************


# Find predictors
predictors=intersect(predictors,colnames(tab))

# Remove variables with zero variance
zerovar=predictors[which(colVars(as.matrix(tab[,predictors]),na.rm=T)==0)]

# Remove highly correlated predictors, accounting for missingness
nzp=setdiff(predictors,zerovar)
ntrim=c(); corx=c()
for (i in 1:(length(nzp)-1)) {
  for (j in (i+1):length(nzp)) {
    x=tab[,nzp[i]]; y=tab[,nzp[j]]
    wx=length(which(is.finite(x))); wy=length(which(is.finite(y)))
    wxy=length(which(is.finite(x+y))); 
    cxy=suppressWarnings(cor(x,y,use="pairwise")); 
    if (!is.finite(cxy)) cxy=0
    if (cxy>0.95 & (max(wxy/wx,wxy/wy)>0.9)) {
      ntrim=c(ntrim,nzp[c(i,j)][which.min(c(wx,wy))])
      corx=rbind(corx,c(nzp[i],nzp[j],wx,wy,wxy,cxy))
    }
  }
}

var_remove=unique(c(zerovar,ntrim))
tab=tab[,setdiff(colnames(tab),var_remove)]


predictors=intersect(predictors,colnames(tab))
preop_predictors=intersect(preop_predictors,colnames(tab))
noninv_predictors=intersect(noninv_predictors,colnames(tab))
discharge_predictors=intersect(discharge_predictors,colnames(tab))


##**********************************************************************
## Predictor details                                                ####
##**********************************************************************

rownames(predictor_details)=predictor_details[,1]
colnames(predictor_details)=c("name","short_name","long_name","log_transform")


##**********************************************************************
## Save                                                             ####
##**********************************************************************

main_table_file=paste0(datadir,"Temp_data/Design/maintable.RData")
save(tab,tab_raw,predictor_details,predictors,preop_predictors,
     noninv_predictors,discharge_predictors,
     file=main_table_file)


