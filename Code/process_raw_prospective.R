##**********************************************************************
##  Risk prediction in PEA                                          ####
##  Supporting code to generate prospective data table                      
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
## Read raw data and pre-format (validation)                        ####
##**********************************************************************

prospective_file=paste0(datadir,"Original/Prospective/pea_prospective_final_anon.csv")
vtab=read.csv(prospective_file,sep=",",stringsAsFactors=FALSE) # raw data


##**********************************************************************
## Basic formatting: must match 'predictors'                        ####
##**********************************************************************

# Fix NAs
vtab[which(vtab=="NA")]=NA

variable_differences=rbind(
  c("body_surface_area","bsa"),
  c("body_mass_index","BMI"),
  c("splenectomy_yn","comorbid_splenectomy"),
  c("diabetes_yn","comorbid_dm"),
  c("hypertension_yn","comorbid_htn"),
  c("af_yn","comorbid_af"),
  c("ihd_yn","comorbid_ihd"),
  c("dyslipidaemia_yn","comorbid_chol"),
  c("thyroid_dysfunction_yn","comorbid_thyroid"),
  c("preop_hct","preop_hematocrit"),
  c("preop_platelets","preop_platelet"),
  c("Ca_channelblocker_yn","med_ca_chan_block"),
  c("Pre_op_vasodilator","med_vasodilator"),
  c("amiodarone_yn","med_amiodarone"),
  c("digoxin_yn","med_digoxin"),
  c("acei_yn","med_acei"),
  c("b_blocker_yn","med_beta_block"),
  c("preop_bnp","bnp_bl"),
  c("preop_nyha","nyha_bl"),
  c("prep_6mwd","sixmwt_bl"),
  c("preop_spap","spap"),
  c("perop_ci","ci_bl"),
  c("preop_co","co_bl"),
  c("preop_dpap","dpap_bl"),
  c("preop_mpap","mpap_bl"),
  c("preop_spap","spap_bl"),
  c("preop_pcwp","pcwp_bl"),
  c("pre_op_pvr","pvr_bl"),
  c("postop_spap","spap_fu1"),
  c("postop_dpap","dpap_fu1"),
  c("postop_mpap","mpap_fu1"),
  c("postop_pcwp","pcwp_fu1"),
  c("postop_co","co_fu1"),
  c("postop_ci","ci_fu1"),
  c("postop_pvr","pvr_fu1"),
  c("postop_6mwd","sixmwt_fu1"),
  c("preop_symptoms","BL.Symptom"),
  c("preop_activity","BL.Activity"),
  c("preop_qol","BL.QoL"),
  c("postop_symptoms","FU.symptom"),
  c("postop_activity","FU.activity"),
  c("postop_qol","FU.qol"))

variable_allnas=c(
  "comorbid_malig_solid_haem","comorbid_malig_solid","haem_malig_yn",
  "any_comorbidity","preop_mpv","preop_mcv","preop_wbc",
  "preop_abs_lym","preop_abs_mono","preop_sodium","preop_potassium","preop_creat",
  "preop_alk_phos", "preop_alb", "preop_bilirubin", "preop_alt", "preop_uric_acid", 
  "preop_tsh","med_diuretic",  "any_medication","smoking_status", "any_medication", "smoking_status", "ever_smoked", 
  "bypass_mins", "dhca_min_total", "disease_type_L", "disease_type_R", 
  "additional_cardiac_procedure", "cabg","pfo", "avr", "mvr", "asd", "icu_days", 
  "hosp_days", "comp_reperfusion","pft_bl_fev1", "pft_bl_fvc", "pft_bl_fev1pc", 
  "pft_bl_fvcpc","pft_bl_fev1.fvc_pc", "pa_sats", "max_spo2_bl", "min_spo2_bl", 
  "tlco_bl_pred.", "afrhythm", "pacedrhythm", "la_bl_dilated", 
  "rv_dilatation", "rv_systolic_function", "rv_basal_diam", "pa_diam", 
  "dEI", "sEI", "TAPSE", "TR_max_vel", "FAC", "IVC.diam", "collapse...50.", 
  "effusion", "LVEF", "E", "A", "E.A", "E.e.av", "E.e.lat", "Lat.e", 
  "DT", "PAT", "MR", "AS", "AR", "spap_bl", "la_4c_area_bl", "ra_area_bl", 
  "nyha_fu1", "la_4c_area_fu1", "ra_area_fu1", "bnp_fu1", "sixmwt_delta", 
  "prs_cad", "prs_t2d", "prs_afib", "prs_ibd"
)

variable_rm=c(
  "preop_blood_date","dob","preop_monocytes","preop_6mwd_date","preop_rhc_date",
  "other_significant_comorbidity_details","Pre.op.pulm.vasodilators","postop_date",
  "preop_camphor_date","postop_camphor_date"
)

variable_temprm=c("X","exclusion.reason","surgeon_predicted_mortality_mean",
                  "Death_date","cause.of.death","pea_date","excluded","centre")

cx=colnames(vtab)
for (i in 1:dim(variable_differences)[1]) {
  cx[which(cx==variable_differences[i,1])]=variable_differences[i,2]
}
colnames(vtab)=cx

for (i in 1:length(variable_allnas)) {
  vtab=cbind(vtab,rep(NA,dim(vtab)[1]))
  colnames(vtab)[length(colnames(vtab))]=variable_allnas[i]
}

vtab=vtab[,setdiff(colnames(vtab),variable_rm)]
#vtab=vtab[,setdiff(colnames(vtab),variable_temprm)]


##**********************************************************************
## Sets of variables to recode as yes/no or numeric                 ####
##**********************************************************************


recode_yn=c("med_vasodilator", "comorbid_splenectomy", 
            "comorbid_dm", "comorbid_htn", "comorbid_chol", "comorbid_af", 
            "comorbid_thyroid", "comorbid_ihd", "med_amiodarone", "med_digoxin", 
            "med_ca_chan_block", "med_acei", "med_beta_block",
            "comorbid_malig_solid_haem", "comorbid_malig_solid", 
            "haem_malig_yn", "any_comorbidity","med_diuretic",
            "any_medication")

to_numeric=c("preop_platelet", 
             "bnp_bl", "nyha_bl", "sixmwt_bl", "spap", "dpap_bl", "mpap_bl", 
             "pcwp_bl", "spap_fu1", "dpap_fu1", "mpap_fu1", "pcwp_fu1", "sixmwt_fu1", 
             "BL.Symptom", "BL.Activity", "BL.QoL", "FU.symptom", "FU.activity", 
             "FU.qol",  "preop_mpv", "preop_mcv", 
             "preop_wbc", "preop_abs_lym", "preop_abs_mono", "preop_sodium", 
             "preop_potassium", "preop_alb", "preop_bilirubin", "preop_alt", 
             "any_medication", "smoking_status", "ever_smoked", 
             "bypass_mins", "dhca_min_total", "disease_type_L", "disease_type_R", 
             "additional_cardiac_procedure", "cabg", "pfo", "avr", "mvr", 
             "asd", "icu_days", "hosp_days", "comp_reperfusion", "pft_bl_fev1", 
             "pft_bl_fvc", "pft_bl_fev1pc", "pft_bl_fvcpc", "pft_bl_fev1.fvc_pc", 
             "pa_sats", "max_spo2_bl", "min_spo2_bl", "tlco_bl_pred.", "afrhythm", 
             "pacedrhythm", "la_bl_dilated", "rv_dilatation", "rv_systolic_function", 
             "rv_basal_diam", "dEI", "sEI", "TAPSE", "IVC.diam", "collapse...50.", 
             "effusion", "LVEF", "E", "A", "E.A", "E.e.av", "E.e.lat", "Lat.e", 
             "DT", "PAT", "MR", "AS", "AR", "spap_bl", "la_4c_area_bl", "ra_area_bl", 
             "nyha_fu1", "la_4c_area_fu1", "ra_area_fu1", "bnp_fu1", "sixmwt_delta", 
             "prs_cad", "prs_t2d", "prs_afib", "prs_ibd")
vtab$sex=as.numeric(vtab$sex=="M")
for (v in 1:length(recode_yn)) 
  vtab[[recode_yn[v]]]=suppressWarnings(as.numeric(vtab[[recode_yn[v]]] %in% c("y","Y")))
for (v in 1:length(to_numeric)) 
  vtab[[to_numeric[v]]]=suppressWarnings(as.numeric(vtab[[to_numeric[v]]]))



##**********************************************************************
## Save table                                                       ####
##**********************************************************************

prospective_table_file=paste0(datadir,"Temp_data/Design/prospective_table.RData")
save(vtab,file=prospective_table_file)




