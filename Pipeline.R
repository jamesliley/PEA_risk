##**********************************************************************
##  Risk prediction in PEA                                          ####
##  Overall pipeline                               
##  James Liley                                                     
##  14/07/23                                                          
##**********************************************************************

##**********************************************************************
## 1. DISCOVERY ANALYSIS                                            ####
##**********************************************************************

## Implement necessary functions                                    ####
source("Code/auxiliary.R") 

## Process raw data                                                 ####
source("Code/process_raw.R")

## Discovery analysis                                               ####
source("Code/discovery_analysis.R") 

## Discovery analysis (low missingness)                             ####
low_missingness <<- TRUE
source("Code/discovery_analysis.R") 
low_missingness <<- FALSE

# Output to Shiny
source("Code/outputs_for_shiny.R") 


##**********************************************************************
## 2. PROSPECTIVE VALIDATION                                        ####
##**********************************************************************

## Process raw data                                                 ####
source("Code/process_raw_prospective.R")

## Discovery analysis                                               ####
source("Code/prospective_analysis.R") 

