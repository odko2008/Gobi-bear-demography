## Gobi bear population abundance and survival estimate ###
## Capture mark-recapture ###
## Before running the CMR models, program Mark must be installed. 

library(RMark)

## Import RData ## 
load("./GobiBear_CMR.RData")
str(build.4_5_6.inp6)

############## Start analysis in RMark ###########

############## time intervals #####################
###### 3 primary sessions (2009 [4 secondary], 2013 [5 secondary], 2017 [6 secondary])
time.intervals.4_5_6 = c(0,0,0,4,0,0,0,0,4,0,0,0,0,0)
nrow(build.4_5_6.inp6) # 51 individuals 


####################### Subset of models that Have identifiable and justifiable parameters ###############
################################################## START ##################################################
##############################################B. MODELS SET 1. #############################################
##################################### Huggin's Full p and c with Heterogeneity; model=RDHFHhet #############
############################################### Sex as group ###############################################
#~~~~~~~~~~~~  Start WITH Sex ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######

## Process data for Huggin's p and c with Het RD, w/ sex

Gobi.RDpc1.sex.456.proc = process.data(data = build.4_5_6.inp6, 
                                       model = "RDHFHet", 
                                       time.intervals = time.intervals.4_5_6, 
                                       groups = c("Sex"))# "Year", "Captured"

## Make DDL ##
Gobi.RDpc1.sex.456.ddl = make.design.data(Gobi.RDpc1.sex.456.proc)

######## Fx to run models ###########
RDpc.sub.sex.456.mods = function()
{
  
  ############ Name and create formulas for each variable ########
  ######################## Survival (S) ##########################
  # Null
  S.dot = list(formula =~1)
  # by sex
  S.sex = list(formula =~Sex)
  # diffrent between open intervals
  S.time =list(formula=~time)
  # by finite mixture
  S.mix = list(formula=~mixture)
  ## by time and sex
  S.time.sex =list(formula=~time+Sex)
  ## Sex by year interaction  
  S.sex.session =list(formula=~Sex*session)
  
  #################################################################
  
  ########################## Detection ############################
  
  ######################## Detection (p=c) ##########################
  
  ## Null
  pc.dot = list(formula=~1, share = T)
  ## by sex
  pc.sex = list(formula=~Sex, share = T)
  ## by session
  pc.session =list(formula=~session, share = T)
  ## By finite mixture 
  pc.mix = list(formula =~ mixture, share = T)
  ## by session with mixture
  pc.session.mix =list(formula=~session+mixture, share = T)
  ## by sex with mixture
  pc.sex.mix = list(formula=~Sex+mixture, share = T)
  ## by session and sex
  pc.session.sex =list(formula=~session+Sex, share = T)
  
  ####################################################################
  ## Gamma single prime ##
  # Fixed to 1 (once gone, always gone)
  gammaP.fixed = list(formula =~1, fixed=1)
  
  ## Gamma double prime ##
  # Fixed to 0 (once in, always in[until death])
  gammaPP.fixed = list(formula=~1, fixed=0)
  
  ## Pi (mixture) ##
  # No mixture
  pi.fixed = list(formula=~1,fixed=1)
  # Null (but not fixed)
  pi.dot = list(formula=~1)
  # By sex
  pi.sex = list(formula=~Sex)
  
   ######################################## MODELS ####################################################

  ############################################## Set 1 ################################################
  ########################### Testing to find best predictors of detection ##############################
  ################# Survival dot, p=c (variable), gamma'=1, gamma''=0, pi.dot ##########################
  pc.dot.subset1 = mark(data = build.4_5_6.inp6, 
                        model = "RDHFHet", 
                        time.intervals = time.intervals.4_5_6, 
                        groups = c( "Sex"), model.parameters=list(S=S.dot, GammaPrime=gammaP.fixed,
                                                                  GammaDoublePrime=gammaPP.fixed,
                                                                  p=pc.dot, pi=pi.dot), threads = 2, output = F)
  
  pc.sex.subset1 = mark(data = build.4_5_6.inp6, 
                        model = "RDHFHet", 
                        time.intervals = time.intervals.4_5_6, 
                        groups = c( "Sex"), model.parameters=list(S=S.dot, GammaPrime=gammaP.fixed,
                                                                  GammaDoublePrime=gammaPP.fixed,
                                                                  p=pc.sex, pi=pi.dot), threads = 2, output = F)
  
  pc.session.set1 = mark(data = build.4_5_6.inp6, 
                         model = "RDHFHet", 
                         time.intervals = time.intervals.4_5_6, 
                         groups = c( "Sex"), model.parameters=list(S=S.dot, GammaPrime=gammaP.fixed,
                                                                   GammaDoublePrime=gammaPP.fixed,
                                                                   p=pc.session, pi=pi.dot), threads = 2, output = F)
  
  pc.mix.subset1 = mark(data = build.4_5_6.inp6, 
                        model = "RDHFHet", 
                        time.intervals = time.intervals.4_5_6, 
                        groups = c( "Sex"), model.parameters=list(S=S.dot, GammaPrime=gammaP.fixed,
                                                                  GammaDoublePrime=gammaPP.fixed,
                                                                  p=pc.mix, pi=pi.dot), threads = 2, output = F)
  
  
  pc.session.mix.subset1 = mark(data = build.4_5_6.inp6, 
                                model = "RDHFHet", 
                                time.intervals = time.intervals.4_5_6, 
                                groups = c( "Sex"), model.parameters=list(S=S.dot, GammaPrime=gammaP.fixed,
                                                                          GammaDoublePrime=gammaPP.fixed,
                                                                          p=pc.session.mix, pi=pi.dot), threads = 2, output = F)
  pc.sex.mix.subset1 = mark(data = build.4_5_6.inp6, 
                            model = "RDHFHet", 
                            time.intervals = time.intervals.4_5_6, 
                            groups = c( "Sex"), model.parameters=list(S=S.dot, GammaPrime=gammaP.fixed,
                                                                      GammaDoublePrime=gammaPP.fixed,
                                                                      p=pc.sex.mix, pi=pi.dot), threads = 2, output = F)
  
  pc.session.sex.subset1 = mark(data = build.4_5_6.inp6, 
                                model = "RDHFHet", 
                                time.intervals = time.intervals.4_5_6, 
                                groups = c( "Sex"), model.parameters=list(S=S.dot, GammaPrime=gammaP.fixed,
                                                                          GammaDoublePrime=gammaPP.fixed,
                                                                          p=pc.session.sex, pi=pi.dot), threads = 2, output = F)
  
  
  ######################################### STOP Subset 1 ######################################################
  
  #################################### START subset 2 ##########################################################
  ############################ Best detection model and covariate for survival ####################################
  
  S.sex_pc.mix.subset2 = mark(data = build.4_5_6.inp6, 
                              model = "RDHFHet", 
                              time.intervals = time.intervals.4_5_6, 
                              groups = c( "Sex"), model.parameters=list(S=S.sex, GammaPrime=gammaP.fixed,
                                                                        GammaDoublePrime=gammaPP.fixed,
                                                                        p=pc.mix, pi=pi.dot), threads = 2, output = F)
  
  S.time_pc.mix.subset2 = mark(data = build.4_5_6.inp6, 
                               model = "RDHFHet", 
                               time.intervals = time.intervals.4_5_6, 
                               groups = c( "Sex"), model.parameters=list(S=S.time, GammaPrime=gammaP.fixed,
                                                                         GammaDoublePrime=gammaPP.fixed,
                                                                         p=pc.mix, pi=pi.dot), threads = 2, output = F)
  
  S.time.sex_pc.mix.subset2 = mark(data = build.4_5_6.inp6, 
                                   model = "RDHFHet", 
                                   time.intervals = time.intervals.4_5_6, 
                                   groups = c( "Sex"), model.parameters=list(S=S.time.sex, GammaPrime=gammaP.fixed,
                                                                             GammaDoublePrime=gammaPP.fixed,
                                                                             p=pc.mix, pi=pi.dot), threads = 2, output = F)
  
   S.time.sex_pc.session.mix.subset2 = mark(data = build.4_5_6.inp6, 
                                           model = "RDHFHet", 
                                           time.intervals = time.intervals.4_5_6, 
                                           groups = c( "Sex"), model.parameters=list(S=S.sex, GammaPrime=gammaP.fixed,
                                                                                     GammaDoublePrime=gammaPP.fixed,
                                                                                     p=pc.session.mix, pi=pi.dot), threads = 2, output = F)
  
  S.time.sex_pc.session.mix.subset2 = mark(data = build.4_5_6.inp6, 
                                           model = "RDHFHet", 
                                           time.intervals = time.intervals.4_5_6, 
                                           groups = c( "Sex"), model.parameters=list(S=S.time, GammaPrime=gammaP.fixed,
                                                                                     GammaDoublePrime=gammaPP.fixed,
                                                                                     p=pc.session.mix, pi=pi.dot), threads = 2, output = F)
  
  S.time.sex_pc.session.mix.subset2 = mark(data = build.4_5_6.inp6, 
                                           model = "RDHFHet", 
                                           time.intervals = time.intervals.4_5_6, 
                                           groups = c( "Sex"), model.parameters=list(S=S.time.sex, GammaPrime=gammaP.fixed,
                                                                                     GammaDoublePrime=gammaPP.fixed,
                                                                                     p=pc.session.mix, pi=pi.dot), threads = 2, output = F)
  
  
  
  #################################### Stop Subset 2 ###########################################################
  
  ########################################### START subset 3 ####################################################
  ############### Best survival (S = S.dot) and detection (pc=mixture) model wit  testing pi covarates ############
  
  S.dot_pc.mix_pi.sex.subset3 = mark(data = build.4_5_6.inp6, 
                                     model = "RDHFHet", 
                                     time.intervals = time.intervals.4_5_6, 
                                     groups = c( "Sex"), model.parameters=list(S=S.dot, GammaPrime=gammaP.fixed,
                                                                               GammaDoublePrime=gammaPP.fixed,
                                                                               p=pc.mix, pi=pi.sex), threads = 2, output = F)
  
  S.dot_pc.dot_pi.sex.subset3 = mark(data = build.4_5_6.inp6, 
                                     model = "RDHFHet", 
                                     time.intervals = time.intervals.4_5_6, 
                                     groups = c("Sex"), model.parameters=list(S=S.dot, GammaPrime=gammaP.fixed,
                                                                              GammaDoublePrime=gammaPP.fixed,
                                                                              p=pc.dot, pi=pi.sex), threads = 2, output = F)
  
  ################################## Stop subset 3 ####################################################
  
  #################################################################  
  ############################################################
  ###############################################################  
  return(collect.models())
}
## Run models

RDpc.sub.sex.456.results = RDpc.sub.sex.456.mods()

## View model list
RDpc.sub.sex.456.results


####### pruned model list ################ 

## Keep models wanted for AIC table
Prune.sex.mods.V1.0 = remove.mark(RDpc.sub.sex.456.results,c(9,6,10,13))

## View AIC table
Prune.sex.mods.V1.0 ## Table 2 in the manuscript 

## View specific models in text doc
RDpc.sub.sex.456.results$S.dot_pc.dot_pi.sex.subset3 ## top model

RDpc.sub.sex.456.results$S.time.sex_pc.mix.subset2 ## second top model

