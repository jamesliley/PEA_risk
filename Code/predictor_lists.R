##**********************************************************************
## Lists of potential predictors by category                        #### 
##                                                                    
##**********************************************************************



### Key

# 1) PFO, AVR, MVR, ASD
#     These relate to other procedures performed at the time of endarterectomy and stand for:
#     PFO - patent foramen ovale closure
#     AVR - aortic valve replacement
#     MVR - mitral valve replacement
#     ASD - atrial septal defect closure
#     CABG - coronary bypass grafting
#     The column before  (additional_cardiac_procedure) codes in a binary form if any of the above were performed (0 being no and 1 being yes)
# 
# 2) co_fu1
#    This denotes cardiac output on right heart catheter at first follow-up within one year of PEA
# 
# 3) ci_fu1
#     This denotes cardiac index (cardiac output/body surface area) on right heart catheter at first follow-
#     up within one year of PEA. There may be some discrepancy between the ci_fu1 value and that
#    derived from co_fu1/bsa (body surface area) as the bsa is pulled from different source data.
#
# 4) preop_rdw
#    This was pulled automatically from the blood dataset I have and is the red cell distribution width
#
# 5) comp_reperfusion
#    This is coded binary if a patient experienced reperfusion lung injury in the perioperative period.
#    Comp stands for complication.
#
# 6) TPG and DPG
#     TPG - Transpulmonary pressure gradient (mean PAP - pcwp)
#     DPG - Diastolic pressure gradient (diastolic PAP - pcwp)
#     Both can be used as markers of pre-capillary involvement in left heart disease. The latter is thought to be less sensitive to changes in loading conditions/cardiac output.
# 
# 7) PVR: pulmonary vascular resistance
#    BSA: body surface area
#
# Other variables
# 



# Predictors which are likely to have a bearing on outcomes, but are not of predictive interest
nuisance_predictors=c("centre","pea_date")

# All potentially useful predictors
# Independent predictors only
predictors=c(
	
	# Anthropometry
	"age_pea", 
	"bsa", "BMI", 
	"sex", 
  
  # Preoperative condition
  "comorbid_splenectomy", "comorbid_ihd", "comorbid_malig_solid_haem", "comorbid_malig_solid", "comorbid_malig_type1", "comorbid_malig_type2", "comorbid_haem_group", "haem_malig_yn", "comorbid_haem_diag", "comorbid_af", "comorbid_thyroid", "comorbid_htn", "comorbid_dm", "comorbid_chol", "comorbid_other", "any_comorbidity",
  "preop_hb", "preop_hematocrit", "preop_mcv", "preop_rdw", "preop_platelet", "preop_mpv", "preop_wbc", "preop_abs_neut", "preop_abs_lym", "preop_abs_mono", "preop_sodium", "preop_potassium", "preop_urea", "preop_creat", "preop_alk_phos", "preop_alb", "preop_bilirubin", "preop_alt", "preop_bnp", "preop_uric_acid", "preop_crp", "preop_tsh", 
  "med_diuretic", "med_vasodilator", "med_beta_block", "med_acei", "med_ca_chan_block", "med_digoxin", "med_amiodarone", "any_medication",
  "smoking_status", "ever_smoked", 
  
  # Peri- and post- operative factors
  "bypass_mins", "dhca_min_total", 
  "disease_type_L", "disease_type_R", 
  "additional_cardiac_procedure", "cabg", "pfo", "avr", "mvr", "asd", 
  "icu_days", "hosp_days", 
  "comp_reperfusion", 
  
  # Baseline physiology and function not repeated
    "pft_bl_fev1", "pft_bl_fvc",  "pft_bl_fev1pc", "pft_bl_fvcpc",  "pft_bl_fev1.fvc_pc", #"pft_bl_fev1.fvc",
    "pa_sats","max_spo2_bl",  "min_spo2_bl" , "tlco_bl_pred.",
  
  "afrhythm", "pacedrhythm",  
  
  # Baseline ultrasound
    "la_bl_dilated",  "afrhythm", "pacedrhythm",  "rv_dilatation", "rv_systolic_function", "rv_basal_diam", "pa_diam",
    "spap", "dEI", "sEI", "TAPSE", "TR_max_vel", "FAC", "IVC.diam", "collapse...50.", "effusion", "LVEF",
    "E", "A", "E.A", "E.e.av", "E.e.lat", "Lat.e", "DT", "PAT", "MR", "AS", "AR", 
    
  # Baseline physiology and function repeated at follow up
    "pcwp_bl", "spap_bl", "dpap_bl", "mpap_bl",	 #"tpg_bl", "dpg_bl", 
    "co_bl", "ci_bl", "pvr_bl", 
    "sixmwt_bl", "BL.Symptom", "BL.Activity", "BL.QoL", "nyha_bl", 
    "la_4c_area_bl","ra_area_bl",
  "bnp_bl", 
  
  # Follow-up physiology and function
  "pcwp_fu1", "spap_fu1", "dpap_fu1", "mpap_fu1", # "tpg_fu1", "dpg_fu1",
  "co_fu1", "ci_fu1", "pvr_fu1", 
  "sixmwt_fu1", "FU.symptom", "FU.activity", "FU.qol", "nyha_fu1", 
  "la_4c_area_fu1","ra_area_fu1",
  "bnp_fu1",
  
        # Changes	
	"sixmwt_delta", 
	# "pvr_change",

        # Genetic risk scores
        "prs_cad","prs_t2d","prs_afib","prs_ibd"
)


# Potentially useful predictors available at discharge from hospital post-PEA
discharge_predictors=c(
	# Anthropometry
	"age_pea", 
	"bsa", "BMI", 
	"sex", 
  
  # Preoperative condition
  "comorbid_splenectomy", "comorbid_ihd", "comorbid_malig_solid_haem", "comorbid_malig_solid", "comorbid_malig_type1", "comorbid_malig_type2", "comorbid_haem_group", "haem_malig_yn", "comorbid_haem_diag", "comorbid_af", "comorbid_thyroid", "comorbid_htn", "comorbid_dm", "comorbid_chol", "comorbid_other", "any_comorbidity",
  "preop_hb", "preop_hematocrit", "preop_mcv", "preop_rdw", "preop_platelet", "preop_mpv", "preop_wbc", "preop_abs_neut", "preop_abs_lym", "preop_abs_mono", "preop_sodium", "preop_potassium", "preop_urea", "preop_creat", "preop_alk_phos", "preop_alb", "preop_bilirubin", "preop_alt", "preop_bnp", "preop_uric_acid", "preop_crp", "preop_tsh", 
  "med_diuretic", "med_vasodilator", "med_beta_block", "med_acei", "med_ca_chan_block", "med_digoxin", "med_amiodarone", "any_medication",
  "smoking_status", "ever_smoked", 
  
  # Peri- and post- operative factors
  "bypass_mins", "dhca_min_total", 
  "disease_type_L", "disease_type_R", 
  "additional_cardiac_procedure", "cabg", "pfo", "avr", "mvr", "asd", 
  "icu_days", "hosp_days", 
  "comp_reperfusion", 
  
  # Baseline physiology and function not repeated
    "pft_bl_fev1", "pft_bl_fvc",  "pft_bl_fev1pc", "pft_bl_fvcpc",  "pft_bl_fev1.fvc_pc", #"pft_bl_fev1.fvc",
    "pa_sats","max_spo2_bl",  "min_spo2_bl" , "tlco_bl_pred.",

    # Baseline ultrasound
      "la_bl_dilated",  "afrhythm", "pacedrhythm",  "rv_dilatation", "rv_systolic_function", "rv_basal_diam", "pa_diam",
      "spap", "dEI", "sEI", "TAPSE", "TR_max_vel", "FAC", "IVC.diam", "collapse...50.", "effusion", "LVEF",
      "E", "A", "E.A", "E.e.av", "E.e.lat", "Lat.e", "DT", "PAT", "MR", "AS", "AR", 
      
    # Baseline physiology and function repeated at follow up
    "pcwp_bl", "spap_bl", "dpap_bl", "mpap_bl",	#"tpg_bl", "dpg_bl", 
    "co_bl", "ci_bl", "pvr_bl", 
    "sixmwt_bl", "BL.Symptom", "BL.Activity", "BL.QoL", "nyha_bl", 
    "la_4c_area_bl","ra_area_bl",
  "bnp_bl", 
  
  # Genetic risk scores
 "prs_cad","prs_t2d","prs_afib","prs_ibd"
)







# Potentially useful predictors available immediately pre-operatively
preop_predictors=c(
	# Anthropometry
	"age_pea", 
	"bsa", "BMI", 
	"sex", 
  
  # Preoperative condition
  "comorbid_splenectomy", "comorbid_ihd", "comorbid_malig_solid_haem", "comorbid_malig_solid", "comorbid_malig_type1", "comorbid_malig_type2", "comorbid_haem_group", "haem_malig_yn", "comorbid_haem_diag", "comorbid_af", "comorbid_thyroid", "comorbid_htn", "comorbid_dm", "comorbid_chol", "comorbid_other", "any_comorbidity",
  "preop_hb", "preop_hematocrit", "preop_mcv", "preop_rdw", "preop_platelet", "preop_mpv", "preop_wbc", "preop_abs_neut", "preop_abs_lym", "preop_abs_mono", "preop_sodium", "preop_potassium", "preop_urea", "preop_creat", "preop_alk_phos", "preop_alb", "preop_bilirubin", "preop_alt", "preop_bnp", "preop_uric_acid", "preop_crp", "preop_tsh", 
  "med_diuretic", "med_vasodilator", "med_beta_block", "med_acei", "med_ca_chan_block", "med_digoxin", "med_amiodarone", "any_medication",
  "smoking_status", "ever_smoked", 
  

  # Baseline physiology and function not repeated
    "pft_bl_fev1", "pft_bl_fvc",  "pft_bl_fev1pc", "pft_bl_fvcpc",  "pft_bl_fev1.fvc_pc", #"pft_bl_fev1.fvc",
    "pa_sats","max_spo2_bl",  "min_spo2_bl" , "tlco_bl_pred.",

    # Baseline ultrasound
      "la_bl_dilated", "afrhythm", "pacedrhythm",  "rv_dilatation", "rv_systolic_function", "rv_basal_diam", "pa_diam",
      "spap", "dEI", "sEI", "TAPSE", "TR_max_vel", "FAC", "IVC.diam", "collapse...50.", "effusion", "LVEF",
      "E", "A", "E.A", "E.e.av", "E.e.lat", "Lat.e", "DT", "PAT", "MR", "AS", "AR", 
          
    # Baseline physiology and function repeated at follow up
    "pcwp_bl", "spap_bl", "dpap_bl", "mpap_bl",	#"tpg_bl", "dpg_bl", 
    "co_bl", "ci_bl", "pvr_bl", 
    "sixmwt_bl", "BL.Symptom", "BL.Activity", "BL.QoL", "nyha_bl", 
    "la_4c_area_bl","ra_area_bl",
  "bnp_bl", 
  

# Genetic risk scores
"prs_cad","prs_t2d","prs_afib","prs_ibd"
)







# Potentially useful non-invasive pre-operative predictors (include U/S, PFTs and bloods)
noninv_predictors=c(
	# Anthropometry
	"age_pea", 
	"bsa", "BMI", 
	"sex", 
  
  # Preoperative condition
  "comorbid_splenectomy", "comorbid_ihd", "comorbid_malig_solid_haem", "comorbid_malig_solid", "comorbid_malig_type1", "comorbid_malig_type2", "comorbid_haem_group", "haem_malig_yn", "comorbid_haem_diag", "comorbid_af", "comorbid_thyroid", "comorbid_htn", "comorbid_dm", "comorbid_chol", "comorbid_other", "any_comorbidity",
  "preop_hb", "preop_hematocrit", "preop_mcv", "preop_rdw", "preop_platelet", "preop_mpv", "preop_wbc", "preop_abs_neut", "preop_abs_lym", "preop_abs_mono", "preop_sodium", "preop_potassium", "preop_urea", "preop_creat", "preop_alk_phos", "preop_alb", "preop_bilirubin", "preop_alt", "preop_bnp", "preop_uric_acid", "preop_crp", "preop_tsh", 
  "med_diuretic", "med_vasodilator", "med_beta_block", "med_acei", "med_ca_chan_block", "med_digoxin", "med_amiodarone", "any_medication",
  "smoking_status", "ever_smoked", 
  
  # Baseline physiology and function not repeated
    "pft_bl_fev1", "pft_bl_fvc",  "pft_bl_fev1pc", "pft_bl_fvcpc",  "pft_bl_fev1.fvc_pc", #"pft_bl_fev1.fvc",
    "pa_sats","max_spo2_bl",  "min_spo2_bl" , "tlco_bl_pred.",
      
      # Baseline ultrasound
        "la_bl_dilated", "afrhythm", "pacedrhythm",  "rv_dilatation", "rv_systolic_function", "rv_basal_diam", "pa_diam",
	"spap", "dEI", "sEI", "TAPSE", "TR_max_vel", "FAC", "IVC.diam", "collapse...50.", "effusion", "LVEF",
	"E", "A", "E.A", "E.e.av", "E.e.lat", "Lat.e", "DT", "PAT", "MR", "AS", "AR", 
    
    # Baseline physiology and function repeated at follow up
    "sixmwt_bl", "BL.Symptom", "BL.Activity", "BL.QoL", "nyha_bl", 
    "la_4c_area_bl","ra_area_bl",
    "bnp_bl", 

     # Genetic risk scores
   "prs_cad","prs_t2d","prs_afib","prs_ibd"
)





discharge_predictors=intersect(discharge_predictors,predictors)
noninv_predictors=intersect(noninv_predictors,predictors)
preop_predictors=intersect(preop_predictors,predictors)
