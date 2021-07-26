# Script to calculate estimated cases due to new exposure, calculate observed sensitivity 
# and calculate adjusted sensitivity accounting for new exposure

# uses bootstrap approach to calculate CI on observed sensitivity
# for each bootstrap calculates number of cases due to new exposure (drawn from bernoulli trials)
# calculates adjusted sensitivity if cases due to new exposure are removed
# generates plots

# load libraries
library(ggplot2)
library(reshape2)
library(gridExtra)

# set path to data
setwd("C:/Users/eidetsum/Filr/My Files/sync/COR/Post_trial_modelling/Public_dataset")

# set the number of samples to draw - set to 20,000 for results in paper
n_samp <- 20          

# These are parameters to explore different values for
# ARI scaling - use to vary the ARI in the model
# Current estimate is based on simple relationship with LTBI prevalence and average age of chort
ARI_scale_v <- c(0.5,1,2)
# Relative protection of prior infection (0.41 from Vynnycky, 0.81 from Andrews, 0 assumes no protection)
Rel_prot_v <- c(0.41,0.81)
# set the sample type to use - 1 or 2 
samp_type_v <- c("ts","osts")

# Load the dataset
dat <- read.csv("ADSL_version2.csv")

#############################################################################
# Data on prevalent cases and other withdrawl - used to check cohort sizes  #
# Total enrolled = 2923                                                     #
#                                                                           # 
#                   RISK11+/3HP+  RISK11+/3HP-  RISK11-                     #
# Enrolled          375           764           1784                        #
# Prevalent(ts)     25            22            14                          #
# Prevalent(osts)   29            27            18                          #
# Withdrawn         12            5             7                           # 
#                                                                           #
# mITT(ts) = 2923-25-22-14-12-5-7 = 2838                                    #
# 3hp-(ts) = 764+1784-22-14-5-7 = 2500                                      #
# mITT(osts) = 2923-29-27-18-12-5-7 = 2825                                  #
# 3hp-(osts) = 764+1784-27-18-5-7 = 2491                                    #
#############################################################################

########################################################################################################################


#########################################################################################################################
# Arrays to store results
par_runs <- length(ARI_scale_v)*length(Rel_prot_v)*length(samp_type_v)

# Outputs by month (1 to 15), sample (1 to n_samp), CORTIS threshold (60 or 26), sensitvity analysis combinations (ARI scaling, relative protection, endpoint definition)
ccp <- array(dim=c(15,n_samp,2,par_runs))    # Cumulative cases in RISK11+
ccn <- ccp                                   # Cumulative cases in RISK11-
ccpw <- ccp                                  # Weighted cumulative cases in RISK11+
ccnw <- ccp                                  # weighted cumulative cases in RISK11-
sens <- ccp                                  # observed sensitivity 
sens_reinf_n <- ccp                          # adjusted sensitivity - removing RISK11- cases due to new exposure
sens_reinf_all <- ccp                        # adjusted sensitivity - removing all cases due to new exposure
sens_ratio <- ccp                            # ratio of adjusted sensitivity to observed sensitivity    
# these are also by IGRA status (all, IGRA+ ,IGRA-)
rcp <- array(dim=c(15,n_samp,2,3,par_runs))  # reinfection cumulative cases in RISK11+
rcn <- rcp                                   # reinfection cumulative cases in RISK11-
rcpw <- rcp                                  # reinfection weighted cumulative cases in RISK11+
rcnw <- rcp                                  # reinfection weighted cumulative cases in RISK11-

# outputs by sample (1 to n_samp), CORTIS threshold and results (60+, 60-, 26+, 26-), sensitvity analysis combinations (ARI scaling, relative protection, endpoint definition)
prev_IGRA <- array(dim=c(n_samp,4,par_runs)) # Overall IGRA prevalence
prev_L <- prev_IGRA                          # prevalence of LTBI (corrected for IGRA performance)
ARI <- prev_IGRA                             # Annual risks of infection used in the model
n_inf <- prev_IGRA                           # number infected

# Summaries of outputs by month (1 to 15), (t, 50%, 2.5%, 97.5%), RISK11 threshold (60, 26),  sensitvity analysis combinations (ARI scaling, relative protection, endpoint definition)
range_sens <- array(dim=c(15,4,2,par_runs))
range_sens_reinf_n <- array(dim=c(15,4,2,par_runs))
range_sens_reinf_all <- array(dim=c(15,4,2,par_runs))
range_sens_ratio <- array(dim=c(15,4,2,par_runs))
range_cases_p <- array(dim=c(15,4,2,par_runs))
range_cases_n <- array(dim=c(15,4,2,par_runs))
# these are also by IGRA status (all, IGRA+ ,IGRA-)
range_cases_reinf_p <- array(dim=c(15,4,2,3,par_runs))
range_cases_reinf_n <- array(dim=c(15,4,2,3,par_runs))


# sensitivity analysis combination parameters
pars <- mat.or.vec(par_runs,3)

#########################################################################################################################
# Run the analysis
zz <- 1

for (aa in 1:length(ARI_scale_v)){       # start loop on ARI scaling
  
  for (rr in 1:length(Rel_prot_v)){      # start loop on relative protection
    
    for (ss in 1:length(samp_type_v)){   # start loop on endpoint definition
    
      # set the ARI scaling
      ARI_scale <- ARI_scale_v[aa]
      # set the relative protevtion 
      Rel_prot <- Rel_prot_v[rr]
      # set the sample type
      samp_type <- samp_type_v[ss]
      # store these in a useful way
      pars[zz,] <- c(ARI_scale,Rel_prot,ss)
      
      # create variables used to subset dataset based on sample type
      mitt_type <- paste("mitt_",samp_type,sep="")
      tevent_type <- paste("tevent_",samp_type,sep="")
      endpoint_type <- paste("endpoint_",samp_type,sep="")
    
      # Extract the data
      # First get the enrolled population
      enrolled <- dat[!is.na(dat$group), ]  
      enrolled[is.na(enrolled$IGRAstatus),"IGRAstatus"]="IGRA-" # set indeterminate IGRA to IGRA-
      # Then get the mITT group (i.e. drop prevalent cases)
      mitt <- enrolled[enrolled[,mitt_type]==1,]  
      # Then get the treatment negative
      mitt_Tn <- mitt[mitt$group!="RISK11+/3HP+",] 
      N_mitt_Tn <- dim(mitt_Tn)[1]
    
      for (i in 1:n_samp){     # start loop on random sample 
  
        # sample the parameters
        # sensitivity of IGRA (gold plus)
        SE_L_IGRA <- rbeta(1,306.46652,42.10844)  # beta distribution fitted to 88 (95% CI 84-91) https://www.atsjournals.org/doi/abs/10.1164/ajrccm-conference.2020.201.1_MeetingAbstracts.A5438 
        # specificty of IGRA (gold plus)
        SP_L_IGRA <- rbeta(1,356.98447,23.11142)  # beta distribution fitted to 88 (95% CI 84-91) https://www.atsjournals.org/doi/abs/10.1164/ajrccm-conference.2020.201.1_MeetingAbstracts.A5438 
        # annual risk of disease in first year after infection (from Vynnycky and Fine 1997)
        R_TB_p <- rbeta(1,947.9745,10000)
        # convert to monthly risks
        R_TB_p_mnth <- R_TB_p/12      # monthly risk 
        R_TB_r <- (1-Rel_prot)*R_TB_p # calculate risk of disease in first year after re-infection 
        R_TB_r_mnth <- R_TB_r/12      # monthly risk
        # average age of cohort
        a_age <- 28.4                 # average age of population 
        
        # sample the data 
        tt <- sample(N_mitt_Tn,replace=TRUE)   
        mitt_Tn_resamp <- mitt_Tn[tt,]
  
        # RISK11+(60) ##########################################################################
        mitt_Tn_p_resamp_60 <- mitt_Tn_resamp[mitt_Tn_resamp$risk11_score>=60,]
        N_p_60 <- dim(mitt_Tn_p_resamp_60)[1]

        # use the prevalence of IGRA and average age to calculate the ARI
        prev_IGRA[i,1,zz] <- dim(mitt_Tn_p_resamp_60[mitt_Tn_p_resamp_60$IGRAstatus=="IGRA+",])[1]/
        dim(mitt_Tn_p_resamp_60)[1] 
        prev_L[i,1,zz] <- (prev_IGRA[i,1,zz] - 1 + SP_L_IGRA)/(SE_L_IGRA - 1 + SP_L_IGRA)
        ARI[i,1,zz] <- ARI_scale*(-log(1-prev_L[i,1,zz])/a_age)
        RI_mnth <- ARI[i,1,zz]/12
        
        # vectors to store month of infection and disease for each participant
        m_inf <- rep(0,N_p_60)
        m_TB <- m_inf
  
        # for each participant
        for (cc in 1:N_p_60){
    
          # set the risk of disease, if IGRA+ use risk in previously infected 
          d_risk <- R_TB_p_mnth
          if (mitt_Tn_p_resamp_60[cc,"IGRAstatus"]=="IGRA+") d_risk <- R_TB_r_mnth
    
          # bernoulli trial for infection for each month that they are in the trial
          temp <- rbinom(round(mitt_Tn_p_resamp_60[cc,tevent_type]),1,RI_mnth)
  
          if (sum(temp)>0){                     # if they are infected at any point
            m_inf[cc] <- min(which(temp==1))    # find the first month infection happened
            # bernoulli trial for disease for each month after infection 
            temp1 <- rbinom(round(mitt_Tn_p_resamp_60[cc,tevent_type])-m_inf[cc],1,d_risk)  
            if (sum(temp1)>0){                            # if they do develop TB
              m_TB[cc] <- min(which(temp1==1))+m_inf[cc]  # find the first month they do and add time of infection
            }
          }
        }
        mitt_Tn_p_resamp_60 <- cbind(mitt_Tn_p_resamp_60,m_inf,m_TB)
        n_inf[i,1,zz] <- sum(m_inf>0) 
  
        # RISK11+(26) ##########################################################################
        mitt_Tn_p_resamp_26 <- mitt_Tn_resamp[mitt_Tn_resamp$risk11_score>=26,]
        N_p_26 <- dim(mitt_Tn_p_resamp_26)[1]
  
        prev_IGRA[i,2,zz] <- dim(mitt_Tn_p_resamp_26[mitt_Tn_p_resamp_26$IGRAstatus=="IGRA+",])[1]/
          dim(mitt_Tn_p_resamp_26)[1]   
        prev_L[i,2,zz] <- (prev_IGRA[i,2,zz] - 1 + SP_L_IGRA)/(SE_L_IGRA - 1 + SP_L_IGRA)
        ARI[i,2,zz] <- ARI_scale*(-log(1-prev_L[i,2,zz])/a_age)
        RI_mnth <- ARI[i,2,zz]/12
        
        m_inf <- rep(0,N_p_26)
        m_TB <- m_inf
    
        for (cc in 1:N_p_26){
    
          d_risk <- R_TB_p_mnth
          if (mitt_Tn_p_resamp_26[cc,"IGRAstatus"]=="IGRA+") d_risk <- R_TB_r_mnth
    
          temp <- rbinom(round(mitt_Tn_p_resamp_26[cc,tevent_type]),1,RI_mnth)
    
          if (sum(temp)>0){
            m_inf[cc] <- min(which(temp==1))
            temp1 <- rbinom(round(mitt_Tn_p_resamp_26[cc,tevent_type])-m_inf[cc],1,d_risk)
            if (sum(temp1)>0){
              m_TB[cc] <- min(which(temp1==1))+m_inf[cc]
            }
          }
        }
        mitt_Tn_p_resamp_26 <- cbind(mitt_Tn_p_resamp_26,m_inf,m_TB)
        n_inf[i,2,zz] <- sum(m_inf>0) 
  
        # RISK11-(60) ##########################################################################
        mitt_Tn_n_resamp_60 <- mitt_Tn_resamp[mitt_Tn_resamp$risk11_score<60,]
        N_n_60 <- dim(mitt_Tn_n_resamp_60)[1]
  
        prev_IGRA[i,3,zz] <- dim(mitt_Tn_n_resamp_60[mitt_Tn_n_resamp_60$IGRAstatus=="IGRA+",])[1]/
          dim(mitt_Tn_n_resamp_60)[1] 
        prev_L[i,3,zz] <- (prev_IGRA[i,3,zz] - 1 + SP_L_IGRA)/(SE_L_IGRA - 1 + SP_L_IGRA)
        ARI[i,3,zz] <- ARI_scale*(-log(1-prev_L[i,3,zz])/a_age)
        RI_mnth <- ARI[i,3,zz]/12
        
        m_inf <- rep(0,N_n_60)
        m_TB <- m_inf
  
        for (cc in 1:N_n_60){
    
          d_risk <- R_TB_p_mnth
          if (mitt_Tn_n_resamp_60[cc,"IGRAstatus"]=="IGRA+") d_risk <- R_TB_r_mnth
    
          temp <- rbinom(round(mitt_Tn_n_resamp_60[cc,tevent_type]),1,RI_mnth)
    
          if (sum(temp)>0){
            m_inf[cc] <- min(which(temp==1))
            temp1 <- rbinom(round(mitt_Tn_n_resamp_60[cc,tevent_type])-m_inf[cc],1,d_risk)
            if (sum(temp1)>0){
              m_TB[cc] <- min(which(temp1==1))+m_inf[cc]
            }
          }
    
        }
        mitt_Tn_n_resamp_60 <- cbind(mitt_Tn_n_resamp_60,m_inf,m_TB)
        n_inf[i,3,zz] <- sum(m_inf>0) 
  
        # RISK11-(26) ##########################################################################
        mitt_Tn_n_resamp_26 <- mitt_Tn_resamp[mitt_Tn_resamp$risk11_score<26,]
        N_n_26 <- dim(mitt_Tn_n_resamp_26)[1]
  
        prev_IGRA[i,4,zz] <- dim(mitt_Tn_n_resamp_26[mitt_Tn_n_resamp_26$IGRAstatus=="IGRA+",])[1]/
          dim(mitt_Tn_n_resamp_26)[1]   
        prev_L[i,4,zz] <- (prev_IGRA[i,4,zz] - 1 + SP_L_IGRA)/(SE_L_IGRA - 1 + SP_L_IGRA)
        ARI[i,4,zz] <- ARI_scale*(-log(1-prev_L[i,4,zz])/a_age)
        RI_mnth <- ARI[i,4,zz]/12
        
        m_inf <- rep(0,N_n_26)
        m_TB <- m_inf
  
        for (cc in 1:N_n_26){
    
          d_risk <- R_TB_p_mnth
          if (mitt_Tn_n_resamp_26[cc,"IGRAstatus"]==1) d_risk <- R_TB_r_mnth
    
          temp <- rbinom(round(mitt_Tn_n_resamp_26[cc,tevent_type]),1,RI_mnth)
    
          if (sum(temp)>0){
            m_inf[cc] <- min(which(temp==1))
            temp1 <- rbinom(round(mitt_Tn_n_resamp_26[cc,tevent_type])-m_inf[cc],1,d_risk)
            if (sum(temp1)>0){
              m_TB[cc] <- min(which(temp1==1))+m_inf[cc]
            }
          }
        } 
        mitt_Tn_n_resamp_26 <- cbind(mitt_Tn_n_resamp_26,m_inf,m_TB)
        n_inf[i,4,zz] <- sum(m_inf>0) 
  
        # now calculate outputs
  
        for (t in 1:15){                 # start loop on month
          
          # Cumulative cases in RISK11+
          ccp[t,i,1,zz] <- dim(mitt_Tn_p_resamp_60[mitt_Tn_p_resamp_60[,tevent_type]<=t&mitt_Tn_p_resamp_60[,endpoint_type]==1,])[1]
          ccp[t,i,2,zz] <- dim(mitt_Tn_p_resamp_26[mitt_Tn_p_resamp_26[,tevent_type]<=t&mitt_Tn_p_resamp_26[,endpoint_type]==1,])[1]
          # Cumulative cases in RISK11-
          ccn[t,i,1,zz] <- dim(mitt_Tn_n_resamp_60[mitt_Tn_n_resamp_60[,tevent_type]<=t&mitt_Tn_n_resamp_60[,endpoint_type]==1,])[1]
          ccn[t,i,2,zz] <- dim(mitt_Tn_n_resamp_26[mitt_Tn_n_resamp_26[,tevent_type]<=t&mitt_Tn_n_resamp_26[,endpoint_type]==1,])[1]
          # Weighted cumulative cases in RISK11+
          ccpw[t,i,1,zz] <- sum(mitt_Tn_p_resamp_60[mitt_Tn_p_resamp_60[,tevent_type]<=t&mitt_Tn_p_resamp_60[,endpoint_type]==1,"group_weights"])
          ccpw[t,i,2,zz] <- sum(mitt_Tn_p_resamp_26[mitt_Tn_p_resamp_26[,tevent_type]<=t&mitt_Tn_p_resamp_26[,endpoint_type]==1,"group_weights"])
          # Weighted cumulative cases in RISK11-
          ccnw[t,i,1,zz] <- sum(mitt_Tn_n_resamp_60[mitt_Tn_n_resamp_60[,tevent_type]<=t&mitt_Tn_n_resamp_60[,endpoint_type]==1,"group_weights"])
          ccnw[t,i,2,zz] <- sum(mitt_Tn_n_resamp_26[mitt_Tn_n_resamp_26[,tevent_type]<=t&mitt_Tn_n_resamp_26[,endpoint_type]==1,"group_weights"])
    
          # reinfection cases
          rcp[t,i,1,1,zz] <- dim(mitt_Tn_p_resamp_60[mitt_Tn_p_resamp_60$m_TB<=t&mitt_Tn_p_resamp_60$m_TB!=0,])[1]
          rcp[t,i,2,1,zz] <- dim(mitt_Tn_p_resamp_26[mitt_Tn_p_resamp_26$m_TB<=t&mitt_Tn_p_resamp_26$m_TB!=0,])[1]
          rcn[t,i,1,1,zz] <- dim(mitt_Tn_n_resamp_60[mitt_Tn_n_resamp_60$m_TB<=t&mitt_Tn_n_resamp_60$m_TB!=0,])[1]
          rcn[t,i,2,1,zz] <- dim(mitt_Tn_n_resamp_26[mitt_Tn_n_resamp_26$m_TB<=t&mitt_Tn_n_resamp_26$m_TB!=0,])[1]
          rcpw[t,i,1,1,zz] <- sum(mitt_Tn_p_resamp_60[mitt_Tn_p_resamp_60$m_TB<=t&mitt_Tn_p_resamp_60$m_TB!=0,"group_weights"])
          rcpw[t,i,2,1,zz] <- sum(mitt_Tn_p_resamp_26[mitt_Tn_p_resamp_26$m_TB<=t&mitt_Tn_p_resamp_26$m_TB!=0,"group_weights"])
          rcnw[t,i,1,1,zz] <- sum(mitt_Tn_n_resamp_60[mitt_Tn_n_resamp_60$m_TB<=t&mitt_Tn_n_resamp_60$m_TB!=0,"group_weights"])
          rcnw[t,i,2,1,zz] <- sum(mitt_Tn_n_resamp_26[mitt_Tn_n_resamp_26$m_TB<=t&mitt_Tn_n_resamp_26$m_TB!=0,"group_weights"])
    
          # Also count reinfection cases by initial IGRA status
          rcp[t,i,1,2,zz] <- dim(mitt_Tn_p_resamp_60[mitt_Tn_p_resamp_60$m_TB<=t&mitt_Tn_p_resamp_60$m_TB!=0&mitt_Tn_p_resamp_60$IGRAstatus=="IGRA+",])[1]
          rcp[t,i,2,2,zz] <- dim(mitt_Tn_p_resamp_26[mitt_Tn_p_resamp_26$m_TB<=t&mitt_Tn_p_resamp_26$m_TB!=0&mitt_Tn_p_resamp_26$IGRAstatus=="IGRA+",])[1]
          rcn[t,i,1,2,zz] <- dim(mitt_Tn_n_resamp_60[mitt_Tn_n_resamp_60$m_TB<=t&mitt_Tn_n_resamp_60$m_TB!=0&mitt_Tn_n_resamp_60$IGRAstatus=="IGRA+",])[1]
          rcn[t,i,2,2,zz] <- dim(mitt_Tn_n_resamp_26[mitt_Tn_n_resamp_26$m_TB<=t&mitt_Tn_n_resamp_26$m_TB!=0&mitt_Tn_n_resamp_26$IGRAstatus=="IGRA+",])[1]
          rcpw[t,i,1,2,zz] <- sum(mitt_Tn_p_resamp_60[mitt_Tn_p_resamp_60$m_TB<=t&mitt_Tn_p_resamp_60$m_TB!=0&mitt_Tn_p_resamp_60$IGRAstatus=="IGRA+","group_weights"])
          rcpw[t,i,2,2,zz] <- sum(mitt_Tn_p_resamp_26[mitt_Tn_p_resamp_26$m_TB<=t&mitt_Tn_p_resamp_26$m_TB!=0&mitt_Tn_p_resamp_26$IGRAstatus=="IGRA+","group_weights"])
          rcnw[t,i,1,2,zz] <- sum(mitt_Tn_n_resamp_60[mitt_Tn_n_resamp_60$m_TB<=t&mitt_Tn_n_resamp_60$m_TB!=0&mitt_Tn_n_resamp_60$IGRAstatus=="IGRA+","group_weights"])
          rcnw[t,i,2,2,zz] <- sum(mitt_Tn_n_resamp_26[mitt_Tn_n_resamp_26$m_TB<=t&mitt_Tn_n_resamp_26$m_TB!=0&mitt_Tn_n_resamp_26$IGRAstatus=="IGRA+","group_weights"])
    
          rcp[t,i,1,3,zz] <- dim(mitt_Tn_p_resamp_60[mitt_Tn_p_resamp_60$m_TB<=t&mitt_Tn_p_resamp_60$m_TB!=0&mitt_Tn_p_resamp_60$IGRAstatus=="IGRA-",])[1]
          rcp[t,i,2,3,zz] <- dim(mitt_Tn_p_resamp_26[mitt_Tn_p_resamp_26$m_TB<=t&mitt_Tn_p_resamp_26$m_TB!=0&mitt_Tn_p_resamp_26$IGRAstatus=="IGRA-",])[1]
          rcn[t,i,1,3,zz] <- dim(mitt_Tn_n_resamp_60[mitt_Tn_n_resamp_60$m_TB<=t&mitt_Tn_n_resamp_60$m_TB!=0&mitt_Tn_n_resamp_60$IGRAstatus=="IGRA-",])[1]
          rcn[t,i,2,3,zz] <- dim(mitt_Tn_n_resamp_26[mitt_Tn_n_resamp_26$m_TB<=t&mitt_Tn_n_resamp_26$m_TB!=0&mitt_Tn_n_resamp_26$IGRAstatus=="IGRA-",])[1]
          rcpw[t,i,1,3,zz] <- sum(mitt_Tn_p_resamp_60[mitt_Tn_p_resamp_60$m_TB<=t&mitt_Tn_p_resamp_60$m_TB!=0&mitt_Tn_p_resamp_60$IGRAstatus=="IGRA-","group_weights"])
          rcpw[t,i,2,3,zz] <- sum(mitt_Tn_p_resamp_26[mitt_Tn_p_resamp_26$m_TB<=t&mitt_Tn_p_resamp_26$m_TB!=0&mitt_Tn_p_resamp_26$IGRAstatus=="IGRA-","group_weights"])
          rcnw[t,i,1,3,zz] <- sum(mitt_Tn_n_resamp_60[mitt_Tn_n_resamp_60$m_TB<=t&mitt_Tn_n_resamp_60$m_TB!=0&mitt_Tn_n_resamp_60$IGRAstatus=="IGRA-","group_weights"])
          rcnw[t,i,2,3,zz] <- sum(mitt_Tn_n_resamp_26[mitt_Tn_n_resamp_26$m_TB<=t&mitt_Tn_n_resamp_26$m_TB!=0&mitt_Tn_n_resamp_26$IGRAstatus=="IGRA-","group_weights"])
    
          # Sensitivity 
          sens[t,i,1,zz] <- ccpw[t,i,1,zz]/(ccpw[t,i,1,zz]+ccnw[t,i,1,zz])
          sens[t,i,2,zz] <- ccpw[t,i,2,zz]/(ccpw[t,i,2,zz]+ccnw[t,i,2,zz])
    
          n_case_adj <- max(0,ccnw[t,i,1,zz]-rcnw[t,i,1,1,zz])
          sens_reinf_n[t,i,1,zz] <- (ccpw[t,i,1,zz])/(ccpw[t,i,1,zz]+n_case_adj)
          n_case_adj <- max(0,ccnw[t,i,2,zz]-rcnw[t,i,2,1,zz])
          sens_reinf_n[t,i,2,zz] <- (ccpw[t,i,2,zz])/(ccpw[t,i,2,zz]+n_case_adj)
    
          p_case_adj <- max(0,ccpw[t,i,1,zz]-rcpw[t,i,1,1,zz])
          n_case_adj <- max(0,ccnw[t,i,1,zz]-rcnw[t,i,1,1,zz])
          sens_reinf_all[t,i,1,zz] <- p_case_adj/(p_case_adj+n_case_adj)
          p_case_adj <- max(0,ccpw[t,i,2,zz]-rcpw[t,i,2,1,zz])
          n_case_adj <- max(0,ccnw[t,i,2,zz]-rcnw[t,i,2,1,zz])
          sens_reinf_all[t,i,2,zz] <- p_case_adj/(p_case_adj+n_case_adj)
        
          # ratio of adjusted sensitivty to observed sensitivity
          sens_ratio[t,i,1,zz] <- sens_reinf_all[t,i,1,zz]/sens[t,i,1,zz]
          sens_ratio[t,i,2,zz] <- sens_reinf_all[t,i,2,zz]/sens[t,i,2,zz]
          
        } # end loop on t

      } # end loop on sample (i)
      
      # get quantiles across samples
      
      for (t in 1:15){    # start loop on t
        
        range_sens[t,,1,zz] <- c(t,quantile(sens[t,,1,zz],probs=c(0.5,0.025,0.975),na.rm=TRUE))
        range_sens[t,,2,zz] <- c(t,quantile(sens[t,,2,zz],probs=c(0.5,0.025,0.975),na.rm=TRUE))
  
        range_sens_reinf_n[t,,1,zz] <- c(t,quantile(sens_reinf_n[t,,1,zz],probs=c(0.5,0.025,0.975),na.rm=TRUE))
        range_sens_reinf_n[t,,2,zz] <- c(t,quantile(sens_reinf_n[t,,2,zz],probs=c(0.5,0.025,0.975),na.rm=TRUE))
  
        range_sens_reinf_all[t,,1,zz] <- c(t,quantile(sens_reinf_all[t,,1,zz],probs=c(0.5,0.025,0.975),na.rm=TRUE))
        range_sens_reinf_all[t,,2,zz] <- c(t,quantile(sens_reinf_all[t,,2,zz],probs=c(0.5,0.025,0.975),na.rm=TRUE))
  
        range_sens_ratio[t,,1,zz] <- c(t,quantile(sens_ratio[t,,1,zz],probs=c(0.5,0.025,0.975),na.rm=TRUE))
        range_sens_ratio[t,,2,zz] <- c(t,quantile(sens_ratio[t,,2,zz],probs=c(0.5,0.025,0.975),na.rm=TRUE))
        
        range_cases_p[t,,1,zz] <- c(t,quantile(ccp[t,,1,zz],probs=c(0.5,0.025,0.975),na.rm=TRUE))
        range_cases_p[t,,2,zz] <- c(t,quantile(ccp[t,,2,zz],probs=c(0.5,0.025,0.975),na.rm=TRUE))
  
        range_cases_n[t,,1,zz] <- c(t,quantile(ccn[t,,1,zz],probs=c(0.5,0.025,0.975),na.rm=TRUE))
        range_cases_n[t,,2,zz] <- c(t,quantile(ccn[t,,2,zz],probs=c(0.5,0.025,0.975),na.rm=TRUE))
  
        range_cases_reinf_p[t,,1,1,zz] <- c(t,quantile(rcp[t,,1,1,zz],probs=c(0.5,0.025,0.975),na.rm=TRUE))
        range_cases_reinf_p[t,,2,1,zz] <- c(t,quantile(rcp[t,,2,1,zz],probs=c(0.5,0.025,0.975),na.rm=TRUE))
  
        range_cases_reinf_n[t,,1,1,zz] <- c(t,quantile(rcn[t,,1,1,zz],probs=c(0.5,0.025,0.975),na.rm=TRUE))
        range_cases_reinf_n[t,,2,1,zz] <- c(t,quantile(rcn[t,,2,1,zz],probs=c(0.5,0.025,0.975),na.rm=TRUE))
  
        range_cases_reinf_p[t,,1,2,zz] <- c(t,quantile(rcp[t,,1,2,zz],probs=c(0.5,0.025,0.975),na.rm=TRUE))
        range_cases_reinf_p[t,,2,2,zz] <- c(t,quantile(rcp[t,,2,2,zz],probs=c(0.5,0.025,0.975),na.rm=TRUE))
  
        range_cases_reinf_n[t,,1,2,zz] <- c(t,quantile(rcn[t,,1,2,zz],probs=c(0.5,0.025,0.975),na.rm=TRUE))
        range_cases_reinf_n[t,,2,2,zz] <- c(t,quantile(rcn[t,,2,2,zz],probs=c(0.5,0.025,0.975),na.rm=TRUE))
  
        range_cases_reinf_p[t,,1,3,zz] <- c(t,quantile(rcp[t,,1,3,zz],probs=c(0.5,0.025,0.975),na.rm=TRUE))
        range_cases_reinf_p[t,,2,3,zz] <- c(t,quantile(rcp[t,,2,3,zz],probs=c(0.5,0.025,0.975),na.rm=TRUE))
  
        range_cases_reinf_n[t,,1,3,zz] <- c(t,quantile(rcn[t,,1,3,zz],probs=c(0.5,0.025,0.975),na.rm=TRUE))
        range_cases_reinf_n[t,,2,3,zz] <- c(t,quantile(rcn[t,,2,3,zz],probs=c(0.5,0.025,0.975),na.rm=TRUE))
  
      } # end loop on t
      zz <- zz+1
    
    } # end loop on endpoint definition
  }   # end loop on relative protection
}     # end loop on ARI scaling

##########################################################################################################################
# CODE TO PULL OUT RESULTS IN PAPER ##########################################

# Report the ARI ##############
ARI_RISK_p_60 <- quantile(ARI[,1,5],probs=c(0.5,0.025,0.975))
ARI_RISK_p_26 <- quantile(ARI[,2,5],probs=c(0.5,0.025,0.975))
ARI_RISK_n_60 <- quantile(ARI[,3,5],probs=c(0.5,0.025,0.975))
ARI_RISK_n_26 <- quantile(ARI[,4,5],probs=c(0.5,0.025,0.975))

# Report the prevalence of IGRA+ #####################
IGRA_p_60 <- quantile(prev_IGRA[,1,5],probs=c(0.5,0.025,0.975))
IGRA_p_26 <- quantile(prev_IGRA[,2,5],probs=c(0.5,0.025,0.975))
IGRA_n_60 <- quantile(prev_IGRA[,3,5],probs=c(0.5,0.025,0.975))
IGRA_n_26 <- quantile(prev_IGRA[,4,5],probs=c(0.5,0.025,0.975))

# Cases ##################################################
cases_to_plot <- as.data.frame(rbind(cbind(range_cases_p[,,1,5],60,"RISK11+","Observed"),cbind(range_cases_p[,,2,5],26,"RISK11+","Observed"),
                                     cbind(range_cases_n[,,1,5],60,"RISK11-","Observed"),cbind(range_cases_n[,,1,5],26,"RISK11-","Observed"),
                                     cbind(range_cases_reinf_p[,,1,1,5],60,"RISK11+","New exposure"),cbind(range_cases_reinf_p[,,2,1,5],26,"RISK11+","New exposure"),
                                     cbind(range_cases_reinf_n[,,1,1,5],60,"RISK11-","New exposure"),cbind(range_cases_reinf_n[,,2,1,5],26,"RISK11-","New exposure")))
colnames(cases_to_plot) <- c("Month","Median","Low","High","Threshold","RISK11_status","Type")


cases_to_plot$Threshold_f = factor(cases_to_plot$Threshold, levels=c('60','26'))
levels(cases_to_plot$Threshold_f) <- c("RISK11 threshold = 60", "RISK11 threshold = 26")
cases_to_plot$Type_f = factor(cases_to_plot$Type, levels=c('Observed','New exposure'))

cases_15m <- cases_to_plot[cases_to_plot$Month==15,] ### THIS IS THE NUMBER OF CASES IN USEFUL FORMAT

# range of sensitivity ###########################################

# observed in trial 
observed_sens_60 <- range_sens[,,1,5]
observed_sens_26 <- range_sens[,,2,5]

range_sens_out <- rbind(cbind(range_sens[,,1,5],60,1,pars[5,1],pars[5,2],pars[5,3]),
                        cbind(range_sens[,,2,5],26,1,pars[5,1],pars[5,2],pars[5,3]),
                        cbind(range_sens_reinf_all[,,1,1],60,2,pars[1,1],pars[1,2],pars[1,3]),
                        cbind(range_sens_reinf_all[,,2,1],26,2,pars[1,1],pars[1,2],pars[1,3]))
for (i in 2:par_runs){
  range_sens_out <- rbind(range_sens_out,cbind(range_sens_reinf_all[,,1,i],60,2,pars[i,1],pars[i,2],pars[i,3]),
                      cbind(range_sens_reinf_all[,,2,i],26,2,pars[i,1],pars[i,2],pars[i,3]))
  
}
range_sens_out <- as.data.frame(range_sens_out)
colnames(range_sens_out) <- c("Month","Med","L","H","Threshold","Type","ARIscale","Protection","Sample")

# Main sensititivy results over time 
primary_adjusted_sens <- range_sens_out[range_sens_out$Type==2&
                                        range_sens_out$ARIscale==1&
                                        range_sens_out$Protection==0.41&
                                        range_sens_out$Sample==1,]
# Sensitivity at t15 for different assumptions
t_15_adjusted_sens <- range_sens_out[range_sens_out$Type==2&range_sens_out$Month==15,]


## CODE TO REPRODUCE FIGURE 1 FROM PAPER ##################################################################################
# Arrange data for plottin
sens_t <- rbind(cbind(seq(1,15),sens[,,1,1],60,1,ARI_plot[1],pars[1,2],pars[1,3]),
                cbind(seq(1,15),sens[,,2,1],26,1,ARI_plot[1],pars[1,2],pars[1,3]),
                cbind(seq(1,15),sens_reinf_all[,,1,1],60,2,ARI_plot[1],pars[1,2],pars[1,3]),
                cbind(seq(1,15),sens_reinf_all[,,2,1],26,2,ARI_plot[1],pars[1,2],pars[1,3]),
                cbind(seq(1,15),sens_reinf_n[,,1,1],60,3,ARI_plot[1],pars[1,2],pars[1,3]),
                cbind(seq(1,15),sens_reinf_n[,,2,1],26,3,ARI_plot[1],pars[1,2],pars[1,3]))
for (i in 2:par_runs){
  
  sens_t <- rbind(sens_t,cbind(seq(1,15),sens[,,1,i],60,1,ARI_plot[i],pars[i,2],pars[i,3]),
                  cbind(seq(1,15),sens[,,2,i],26,1,ARI_plot[i],pars[i,2],pars[i,3]),
                  cbind(seq(1,15),sens_reinf_all[,,1,i],60,2,ARI_plot[i],pars[i,2],pars[i,3]),
                  cbind(seq(1,15),sens_reinf_all[,,2,i],26,2,ARI_plot[i],pars[i,2],pars[i,3]),
                  cbind(seq(1,15),sens_reinf_n[,,1,i],60,3,ARI_plot[i],pars[i,2],pars[i,3]),
                  cbind(seq(1,15),sens_reinf_n[,,2,i],26,3,ARI_plot[i],pars[i,2],pars[i,3]))
  
}

colnames(sens_t) <- c("Month",seq(1,n_samp),"Threshold","Type","ARI","Protection","Sample")
sens_t <- as.data.frame((sens_t))
sens_t <- melt(sens_t,id.vars=c("Month","Threshold","Type","ARI","Protection","Sample"))

sens_t$Threshold_f = factor(sens_t$Threshold, levels=c('60','26'))
levels(sens_t$Threshold_f) <- c("RISK11 threshold = 60", "RISK11 threshold = 26")

sens_t$Type_f = factor(sens_t$Type, levels=c('1','2','3'))
levels(sens_t$Type_f) <- c("Observed", "Adjusted","Adjusted (RISK11- only)")

sens_t$Sample_f = factor(sens_t$Sample, levels=c('1','2'))
levels(sens_t$Sample_f) <- c("Two sample positive", "One or two sample postive")

# Plots for paper- figures 1 and 2

# Set ARI to use for fig 1
ARI_plot <- pars[,1]#*ARI_base

# fig 1 combines plot of cases over time and plot of sens over time
cases_plot <- ggplot(cases_to_plot)+
  geom_step(aes(x=as.numeric(as.character(Month)),y=as.numeric(as.character(Median)),colour=RISK11_status,lty=Type_f))+
  #geom_ribbon(data=cases_to_plot[cases_to_plot$Type=="Observed",],aes(x=as.numeric(as.character(Month)),min=as.numeric(as.character(Low)),max=as.numeric(as.character(High)),fill=RISK11_status),alpha=0.1)+
  #geom_ribbon(data=cases_to_plot[cases_to_plot$Type=="New exposure",],aes(x=as.numeric(as.character(Month)),min=as.numeric(as.character(Low)),max=as.numeric(as.character(High)),fill=RISK11_status),alpha=0.1)+
  #geom_step(aes(x=as.numeric(as.character(Month)),y=as.numeric(as.character(Low)),colour=RISK11_status))+
  #geom_step(aes(x=as.numeric(as.character(Month)),y=as.numeric(as.character(High)),colour=RISK11_status))+
  #geom_ribbon(data=cases_to_plot,aes(x=as.numeric(as.character(Month)),min=as.numeric(as.character(Low)),max=as.numeric(as.character(High)),fill=RISK11_status),alpha=0.3)+
  facet_wrap(~Threshold_f)+
  theme_bw()+
  scale_x_continuous("Month", seq(1,15), limits=c(0,15.1),expand = c(0, 0))+
  scale_y_continuous("Cumulative TB cases", seq(0,20,1), limits=c(0,20),expand = c(0, 0))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.title=element_blank())+
  theme(legend.position="bottom")+
  theme(strip.background =element_rect(fill="white"))

# Plot sens at 9,12,15 months with baseline ARI and protection, 60 and 26 thresholds - figure 2 in paper
sens_plot <- ggplot(sens_t[sens_t$Type_f%in%c("Observed","Adjusted")&
                            sens_t$Sample_f=="Two sample positive"&
                            sens_t$ARI==1&
                            sens_t$Protection==0.41&
                            sens_t$Month%in%c(9,12,15),],
                   aes(x=as.factor(Month),y=value,fill=Type_f))+
  geom_boxplot()+
  facet_grid(~Threshold_f)+
  theme_bw()+
  scale_x_discrete("Month")+
  scale_y_continuous("Sensitivty", seq(0,1,0.2), limits=c(0,1.02),expand = c(0, 0))+
  scale_fill_manual(values=c("white","grey"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.title=element_blank())+
  theme(legend.position="bottom")+
  theme(strip.background =element_rect(fill="white"))

# combine 
grid.arrange(cases_plot, sens_plot, nrow = 2)
