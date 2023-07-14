#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#  http://shiny.rstudio.com/
#
#
#
## General strategy
## Trim tab to only useful variables and import
## Inputs on left, default to average values
## Output: 
## Survival curve (cox-PH model); expected value and error in 5-year survival
##  For a typical patient: (mean)
##  For a patient with these characteristics: 

## Packages and scripts
library(shiny)
library(shinyBS)
library(randomForest)
library(randomForestSRC)
library(survival)
logit=function(x) 1/(1+exp(-x))
inverse_logit=function(x) -log((1/x)-1)


## Data
load("data.RData")

## Window specifications
left=30; top=80; w1=520; w2=900; div=20

ui <- fluidPage(
 align="center", 

 # Application title
 absolutePanel(left=left,top=0,width=w1+w2+ 2*left +div,h1("Risk in PEA procedure. Not for clinical use.")),

 absolutePanel(left=left,top=top,width=w1, 
   h3("Single patient"),
  
   actionButton("dembutton", "Demographics"),
    bsCollapse(id = "demographics", open = "Panel 3",
   bsCollapsePanel(title=NULL,value="dem_panel", 
   splitLayout(cellWidths = rep("50%",2),
    numericInput("age_pea","Age:",value = signif(mean(mtab['age_pea'],na.rm=T),digits=3)),
   radioButtons("sex","Sex: ",inline=T,choiceNames=c("Unk.","M","F"),choiceValues=c(NA,1,0))),
    splitLayout(cellWidths = rep("50%",2),
     numericInput("bsa","Body surface area: ",value= signif(mean(mtab['bsa'],na.rm=T),digits=3)),
     numericInput("BMI","Body mass index: ",value= signif(mean(mtab['BMI'],na.rm=T),digits=3)))
   )),

  actionButton("combutton", "Comorbidities"),
    bsCollapse(id = "comorbidities", open = "Panel 3",
   bsCollapsePanel(title=NULL,value="com_panel", 
   splitLayout(cellWidths = rep("50%",2),
   radioButtons("comorbid_splenectomy","Splenectomy",inline=T,choiceNames=c("Unk.","Y","N"),choiceValues=c(NA,1,0)), 
   radioButtons("comorbid_af","Atr. fibr.",inline=T,choiceNames=c("Unk.","Y","N"),choiceValues=c(NA,1,0))), 
   splitLayout(cellWidths = rep("50%",2),
   radioButtons("comorbid_dm","Diabetes",inline=T,choiceNames=c("Unk.","Y","N"),choiceValues=c(NA,1,0)), 
   radioButtons("comorbid_thyroid","Hypo/Hyper-thyroidism",inline=T,choiceNames=c("Unk.","Y","N"),choiceValues=c(NA,1,0))), 
   splitLayout(cellWidths = rep("50%",2),
   radioButtons("comorbid_htn","Hypertension",inline=T,choiceNames=c("Unk.","Y","N"),choiceValues=c(NA,1,0)), 
   radioButtons("comorbid_ihd","Isch. heart disease",inline=T,choiceNames=c("Unk.","Y","N"),choiceValues=c(NA,1,0))),   
   splitLayout(cellWidths = rep("50%",2),
   radioButtons("comorbid_chol","Dyslipidaemia",inline=T,choiceNames=c("Unk.","Y","N"),choiceValues=c(NA,1,0)), 
   radioButtons("comorbid_other","Other significant",inline=T,choiceNames=c("Unk.","Y","N"),choiceValues=c(NA,1,0))),   
   	style = "default")
   ),

  
  
   actionButton("medbutton", "Medications"),
    bsCollapse(id = "medications", open = "Panel 2",
   bsCollapsePanel(title=NULL,value="med_panel", 
   splitLayout(cellWidths = rep("50%",2),
    radioButtons("med_amiodarone","Amiodarone: ",inline=T,choiceNames=c("Unk.","Y","N"),choiceValues=c(NA,1,0)),
   radioButtons("med_digoxin","Digoxin: ",inline=T,choiceNames=c("Unk.","Y","N"),choiceValues=c(NA,1,0))
   ),
   splitLayout(cellWidths = rep("50%",2),
    radioButtons("med_ca_chan_block","Ca channel blocker: ",inline=T,choiceNames=c("Unk.","Y","N"),choiceValues=c(NA,1,0)),
   radioButtons("med_vasodilator","Vasodilator: ",inline=T,choiceNames=c("Unk.","Y","N"),choiceValues=c(NA,1,0))
   ),
   splitLayout(cellWidths = rep("50%",2),
    radioButtons("med_acei","ACE inhibitor: ",inline=T,choiceNames=c("Unk.","Y","N"),choiceValues=c(NA,1,0)),
    radioButtons("med_beta_block","Beta blocker: ",inline=T,choiceNames=c("Unk.","Y","N"),choiceValues=c(NA,1,0))
   ),
   style = "default")
   ),
  
   actionButton("bloodbutton", "Blood tests"),
    bsCollapse(id = "bloods", open = "Panel 3",   
    bsCollapsePanel(title=NULL,value="blood_panel",
    splitLayout(cellWidths = rep("50%",2),
   numericInput("preop_hematocrit","Hct:",value = signif(mean(mtab['preop_hematocrit'],na.rm=T),digits=3)),
   numericInput("preop_rdw","RDW:",value = signif(mean(mtab['preop_rdw'],na.rm=T),digits=3))),
   splitLayout(cellWidths = rep("50%",2),
   numericInput("preop_platelet","Plts:",value = signif(mean(mtab['preop_platelet'],na.rm=T),digits=3)),
   numericInput("preop_abs_mono","Abs. mono.:",value = signif(mean(mtab['preop_abs_mono'],na.rm=T),digits=3))),
   splitLayout(cellWidths = rep("50%",2),
   numericInput("bnp_bl","BNP:",value = signif(mean(mtab['bnp_bl'],na.rm=T),digits=3)),
   numericInput("preop_urea","Urea:",value = signif(mean(mtab['preop_urea'],na.rm=T),digits=3))),
   style = "default")
   ),

  
   actionButton("functionbutton", "Baseline function"),
    bsCollapse(id = "function", open = "Panel 3",   
    bsCollapsePanel(title=NULL,value="function_panel",
   splitLayout(cellWidths = rep("50%",2),
   numericInput("nyha_bl","NYHA class:",value = signif(mean(mtab['nyha_bl'],na.rm=T),digits=3)),
   numericInput("sixmwt_bl","6mt walk:",value = signif(mean(mtab['sixmwt_bl'],na.rm=T),digits=3))),
     splitLayout(cellWidths = rep("33%",3),
   numericInput("BL.QoL","QoL:",value = signif(mean(mtab['BL.QoL'],na.rm=T),digits=3)),
     numericInput("BL.Activity","Act.:",value = signif(mean(mtab['BL.Activity'],na.rm=T),digits=3)),
     numericInput("BL.Symptom","Sympt.:",value = signif(mean(mtab['BL.Symptom'],na.rm=T),digits=3))),
   style = "default")
   ),
  
   hr(),
  
   h3("Multiple patients"),
  
   actionButton("multiplebutton", "Table"),
    bsCollapse(id = "multiple", open = "Panel 3",   
    bsCollapsePanel(title=NULL,value="multiple_panel",
   textAreaInput("multiple_input",NULL,value = "",rows=5),
   style = "default")
   )
  
 ),
  
   # Show a plot of the generated distribution
 absolutePanel(left=left+w1+div,top=top,width=w2,
 tabsetPanel(
  tabPanel("Morbidity",h4("Expected change in CAMPHOR score"),plotOutput("camphplot",width=850,height=400),textOutput("camphortext")),
  tabPanel("Survival",h4("Survival curve"),plotOutput("survival",width=850,height=400)) #,textOutput("survivaltext"))
 ),
 br(), br(),
 hr(),
 splitLayout(cellWidths = rep("50%",2),
 actionButton("update","Update and recompute"),
 downloadButton("downloadData", "Download"))
 )
)


server <- function(input, output,session) {
  
  predvec=mtab
  # Evaulate baesline logistic regression-based score for 5-year survival probability
  nx=intersect(names(lr_coefficients),names(predvec))
  lx=as.numeric(lr_coefficients[nx])

  observeEvent(input$dembutton, (updateCollapse(session, "demographics", open = "dem_panel",close="dem_panel")))
  observeEvent(input$medbutton, (updateCollapse(session, "medications", open = "med_panel",close="med_panel")))
  observeEvent(input$combutton, (updateCollapse(session, "comorbidities", open = "com_panel",close="com_panel")))
  observeEvent(input$bloodbutton, (updateCollapse(session, "bloods", open = "blood_panel",close="blood_panel")))
  observeEvent(input$functionbutton, (updateCollapse(session, "function", open = "function_panel",close="function_panel")))  
  observeEvent(input$multiplebutton, (updateCollapse(session, "multiple", open = "multiple_panel",close="multiple_panel")))  
   
  # on button click
  go=eventReactive(input$update,{
  
  if (nchar(trimws(input$multiple_input))<1) {
  # Read in inputs
   nmi=names(input); nmp=intersect(nmi,names(predvec))
   for (i in 1:length(nmp)) eval(parse(text=paste0("predvec['",nmp[i],"']=as.numeric(input$",nmp[i],")")))
   w=which(is.na(predvec))
   predvec[w]=mtab[w]
   np=names(predvec)
   Xp=matrix(predvec,1,length(predvec)); colnames(Xp)=np
  } else {
   predx=read.table(text=input$multiple_input,header=T); dx=dim(predx)
   nmi=colnames(predx)[2:dx[2]]; nmp=intersect(nmi, names(predvec))
   Xp=t(matrix(mtab,length(mtab),dx[1])); colnames(Xp)=names(mtab)
   for (i in 1:length(nmp)) Xp[,nmp[i]]=predx[,nmp[i]]
  }
 
  mortprob=predict(mod_dm_preop,Xp,type="prob")[,2]

  survfit=predict(mod_5m_preop,data.frame(Xp))
  camphor=predict(mod_dq_preop,Xp)
  
  download=cbind(Xp,surv_score=survfit$predicted,camph_predicted=camphor)
  
  xscore=mortprob  
  return(list(xscore=xscore,mortprob=mortprob,predvec=Xp,survfit=survfit,cfit=camphor,download=download))
  })
  
  
  output$morttext=renderText({
  getv=go()
  with(getv,paste0("5-year mortality probability: ",signif(mortprob,digits=2)))
  })
 

  output$survival <- renderPlot({
   getv=go()
   with (getv,{
   	xx=survtime/365.25
    plot(xx,survmean, xlab = "Years", ylab="Survival",type="l",ylim=c(0,1))
    lines(xx,surv20,lty=2)
    lines(xx,surv80,lty=2)
    for (i in 1:dim(survfit$survival[,,drop=FALSE])[1]) lines(xx,t(survfit$survival[i, , drop = FALSE]),col="red")
   })
   legend("bottomleft",c("Average","Best and worst 20%","Specific patient(s)"),
    lty=c(1,2,1),
    col=c("black","black","red"))
   })
  

  output$camphplot <- renderPlot({
   getv=go()
   with (getv,{
    pred=cfit
   qbase=predvec[,"BL.QoL"] + predvec[,"BL.Activity"] + predvec[,"BL.Symptom"]
   qfu=qbase-pred

   dd1=density(cph_bl[which(is.finite(cph_bl))]); 
   dd2=density(cph_fu[which(is.finite(cph_fu))]); 
   sc=5; padx=0.1; gl=0.2

   rx=c(dd1$x,dd2$x);
   plot(0,type="n",ylim=range(rx),xlim=c(-padx,1+padx),ylab=expression("Camphor score"),xlab="Timepoint",xaxt="n")
   axis(1,at=0:1,labels=c("Baseline","Follow-up"))
   segments(rep(0,length(cph_bl)),cph_bl,rep(1,length(cph_fu)),cph_fu,col=gray(1-gl + gl*runif(length(cph_fu))))
   segments(0,mean(cph_bl,na.rm=T),1,mean(cph_fu,na.rm=T),lwd=4,col="darkgray")
   points(rep(0,length(cph_bl)),cph_bl,pch=16,col="gray"); points(rep(1,length(cph_fu)),cph_fu,pch=16,col="gray")
   lines(dd1$y*sc,dd1$x); lines(-dd1$y*sc,dd1$x); 
   lines(1+dd2$y*sc,dd2$x); lines(1-dd2$y*sc,dd2$x)
   segments(0,qbase,1,qfu,col="red",lwd=2)
    })
  })
 
 output$camphortext=renderText({
  getv=go()
   with(getv, {
    pred=cfit
   qbase=predvec["BL.QoL"] + predvec["BL.Activity"] + predvec["BL.Symptom"]
   qfu=qbase-pred
   paste0("Expected change in Camphor score: ",round(mean(pred)),". Gray lines show Camphor score changes for population; black lines show distribution of baseline and follow-up Camphor scores. Red line shows projected Camphor score change for this individual.")
   })
  })

  output$downloadData <- downloadHandler(
    filename = "pea_risk_predictions.csv",
    content = function(file) {
      getv=go()
      #with(getv, return(go$download))
    	write.csv(getv$download, file, row.names = FALSE)
    }
    )
}

# Run the application 
shinyApp(ui = ui, server = server)

