###########################################################################################################
RiskReductions <- function(scenario){
  #Sets the type-specific and age-dependent relative risk reductions 
  #according to a given scenario and for the different vaccines
  #For each vaccine and gender (women, men and MSM), 
  #RRR is a matrix of size 10x100
  #Rows correspond to HPV types (hr) and columns to ages 1:100
  
  #scenario = (scenario_HR, scenario_LR,AW/RRP)
  #scenario_HR =  1: direct protection only
  #               2: + herd immunity (base-case)
  #               3: + waning for all non16/18
  #               4: + waning for all non targeted types (cross-protection)
  #scenario_LR =  1: Base-case: protection equal to mean effect of type 18 = 0.84
  #               2: elimination
  #               3: direct protection only
  #warts/RRP =    0: neither
  #               1: include RRP only
  #               2: include AW only
  #               3: include both RRP and AW (base-case)
  
  sc_HR <- scenario[1]
  sc_LR <- scenario[2]
  
  #Gender neutral vaccination (2 dose) at age 10
  
  #Vaccine efficacy types for ongogenic types 16,18,31,33,35,39,45,51,52,58
  VE_2vCons <- c(0.98,0.98,0,0,0,0,0,0,0,0) #2v conservative
  VE_2vMid <- c(0.98,0.98,0.75,0.5,0,0,0.8,0,0,0) #2v middle
  VE_2vLib1 <- c(0.98,0.98,0.88,0.68,0,0.75,0.82,0.54,0,0) #2v very high cross-protection
  VE_2vLib2 <- c(0.98,0.98,0.66,0.41,0.40,0,0.81,0,0.36,0.30) #2v high cross-protection
  VE_9v <- c(0.98,0.98,0.97,0.97,0,0,0.97,0,0.97,0.97) #9v
  
  #Compute RRRs for the HR types
  if (sc_HR == 1){ 
    #scenario 1: All-or-nothing VE, life-long, no herd immunity
    
    RRR_2vCons_w <- matrix(vacc.up[1]*VE_2vCons,types,100)
    RRR_2vCons_m <- matrix(vacc.up[2]*VE_2vCons,types,100)
    RRR_2vCons_w[,1:a0] <- 0 #no reduction before vaccination
    RRR_2vCons_m[,1:a0] <- 0
    
    RRR_2vMid_w <- matrix(vacc.up[1]*VE_2vMid,types,100)
    RRR_2vMid_m <- matrix(vacc.up[2]*VE_2vMid,types,100)
    RRR_2vMid_w[,1:a0] <- 0 #no reduction before vaccination
    RRR_2vMid_m[,1:a0] <- 0
    
    RRR_2vLib1_w <- matrix(vacc.up[1]*VE_2vLib1,types,100)
    RRR_2vLib1_m <- matrix(vacc.up[2]*VE_2vLib1,types,100)
    RRR_2vLib1_w[,1:a0] <- 0 #no reduction before vaccination
    RRR_2vLib1_m[,1:a0] <- 0
    
    RRR_2vLib2_w <- matrix(vacc.up[1]*VE_2vLib2,types,100)
    RRR_2vLib2_m <- matrix(vacc.up[2]*VE_2vLib2,types,100)
    RRR_2vLib2_w[,1:a0] <- 0 #no reduction before vaccination
    RRR_2vLib2_m[,1:a0] <- 0
    
    RRR_4v_w <- RRR_2vCons_w
    RRR_4v_m <- RRR_2vCons_m
    
    RRR_9v_w <- matrix(vacc.up[1]*VE_9v,types,100)
    RRR_9v_m <- matrix(vacc.up[2]*VE_9v,types,100)
    RRR_9v_w[,1:a0] <- 0 #no reduction before vaccination
    RRR_9v_m[,1:a0] <- 0
    
    #For now
    RRR_2vCons_MSM <- RRR_2vCons_m
    RRR_2vMid_MSM <- RRR_2vMid_m
    RRR_2vLib1_MSM <- RRR_2vLib1_m
    RRR_2vLib2_MSM <- RRR_2vLib2_m
    RRR_4v_MSM <- RRR_4v_m
    RRR_9v_MSM <- RRR_9v_m
    
  } else if (sc_HR == 2){ 
    #scenario 2: All-or-nothing VE, life-long, + herd immunity
    
    #Relative risk reductions
    RRRs_w <- (HPVinc_pre_f_set - HPVinc_post_f)/HPVinc_pre_f_set
    RRRs_w[is.na(RRRs_w)] <- 0
    RRRs_m <- (HPVinc_pre_m_set - HPVinc_post_m)/HPVinc_pre_m_set
    RRRs_m[is.na(RRRs_m)] <- 0
    
    #matrices for the different vaccines
    RRR_2vCons_w <- matrix(0,types,100)
    RRR_2vCons_m <- matrix(0,types,100)
    RRR_2vMid_w <- matrix(0,types,100)
    RRR_2vMid_m <- matrix(0,types,100)
    RRR_2vLib1_w <- matrix(0,types,100)
    RRR_2vLib1_m <- matrix(0,types,100)
    RRR_2vLib2_w <- matrix(0,types,100)
    RRR_2vLib2_m <- matrix(0,types,100)
    RRR_9v_w <- matrix(0,types,100)
    RRR_9v_m <- matrix(0,types,100)
    
    #plugging in the correct RRRs
    RRR_2vCons_w[c(1,2),] <- RRRs_w[c(1,2),]
    RRR_2vCons_m[c(1,2),] <- RRRs_m[c(1,2),]
    RRR_2vMid_w[c(1,2,3,4,7),] <- RRRs_w[c(1,2,4,8,13),]
    RRR_2vMid_m[c(1,2,3,4,7),] <- RRRs_m[c(1,2,4,8,13),]
    RRR_2vLib1_w[c(1,2,3,4,6,7,8),] <- RRRs_w[c(1,2,5,9,12,15,17),]
    RRR_2vLib1_m[c(1,2,3,4,6,7,8),] <- RRRs_m[c(1,2,5,9,12,15,17),]
    RRR_2vLib2_w[c(1,2,3,4,5,7,9,10),] <- RRRs_w[c(1,2,3,7,11,14,18,20),]
    RRR_2vLib2_m[c(1,2,3,4,5,7,9,10),] <- RRRs_m[c(1,2,3,7,11,14,18,20),]
    RRR_4v_w <- RRR_2vCons_w
    RRR_4v_m <- RRR_2vCons_m
    RRR_9v_w[c(1,2,3,4,7,9,10),] <- RRRs_w[c(1,2,6,10,16,19,21),]
    RRR_9v_m[c(1,2,3,4,7,9,10),] <- RRRs_m[c(1,2,6,10,16,19,21),]
    
    #For now
    RRR_2vCons_MSM <- RRR_2vCons_m
    RRR_2vMid_MSM <- RRR_2vMid_m
    RRR_2vLib1_MSM <- RRR_2vLib1_m
    RRR_2vLib2_MSM <- RRR_2vLib2_m
    RRR_4v_MSM <- RRR_4v_m
    RRR_9v_MSM <- RRR_9v_m

  } else if (sc_HR == 3){ 
    #scenario 3: All-or-nothing VE, waning, + herd immunity
    
    #Extreme waning scenario: waning starts at age 20, protection very small at age 40
    #exp(-20*0.15) = 0.05
    #waning only for the non-16/18 HR types
    
    #Relative risk reductions
    RRRs_w <- (HPVinc_pre_f_set - HPVinc_post_f_waning015)/HPVinc_pre_f_set
    RRRs_w[is.na(RRRs_w)] <- 0
    RRRs_m <- (HPVinc_pre_m_set - HPVinc_post_m_waning015)/HPVinc_pre_m_set
    RRRs_m[is.na(RRRs_m)] <- 0
    
    #matrices for the different vaccines
    RRR_2vCons_w <- matrix(0,types,100)
    RRR_2vCons_m <- matrix(0,types,100)
    RRR_2vMid_w <- matrix(0,types,100)
    RRR_2vMid_m <- matrix(0,types,100)
    RRR_2vLib1_w <- matrix(0,types,100)
    RRR_2vLib1_m <- matrix(0,types,100)
    RRR_2vLib2_w <- matrix(0,types,100)
    RRR_2vLib2_m <- matrix(0,types,100)
    RRR_9v_w <- matrix(0,types,100)
    RRR_9v_m <- matrix(0,types,100)
    
    #plugging in the correct RRRs
    RRR_2vCons_w[c(1,2),] <- RRRs_w[c(1,2),]
    RRR_2vCons_m[c(1,2),] <- RRRs_m[c(1,2),]
    RRR_2vMid_w[c(1,2,3,4,7),] <- RRRs_w[c(1,2,4,8,13),]
    RRR_2vMid_m[c(1,2,3,4,7),] <- RRRs_m[c(1,2,4,8,13),]
    RRR_2vLib1_w[c(1,2,3,4,6,7,8),] <- RRRs_w[c(1,2,5,9,12,15,17),]
    RRR_2vLib1_m[c(1,2,3,4,6,7,8),] <- RRRs_m[c(1,2,5,9,12,15,17),]
    RRR_2vLib2_w[c(1,2,3,4,5,7,9,10),] <- RRRs_w[c(1,2,3,7,11,14,18,20),]
    RRR_2vLib2_m[c(1,2,3,4,5,7,9,10),] <- RRRs_m[c(1,2,3,7,11,14,18,20),]
    RRR_4v_w <- RRR_2vCons_w
    RRR_4v_m <- RRR_2vCons_m
    RRR_9v_w[c(1,2,3,4,7,9,10),] <- RRRs_w[c(1,2,6,10,16,19,21),]
    RRR_9v_m[c(1,2,3,4,7,9,10),] <- RRRs_m[c(1,2,6,10,16,19,21),]
    
    #For now
    RRR_2vCons_MSM <- RRR_2vCons_m
    RRR_2vMid_MSM <- RRR_2vMid_m
    RRR_2vLib1_MSM <- RRR_2vLib1_m
    RRR_2vLib2_MSM <- RRR_2vLib2_m
    RRR_4v_MSM <- RRR_4v_m
    RRR_9v_MSM <- RRR_9v_m
    
  } else if (sc_HR == 4){ 
    #scenario 4: All-or-nothing VE, waning (non-targeted types only), + herd immunity
    
    #Extreme waning scenario: waning starts at age 20, protection very small at age 40
    #exp(-20*0.15) = 0.05
    
    #waning only for the non-targeted types (cross protection)
    
    #matrices for the different vaccines
    RRR_2vCons_w <- matrix(0,types,100)
    RRR_2vCons_m <- matrix(0,types,100)
    RRR_2vMid_w <- matrix(0,types,100)
    RRR_2vMid_m <- matrix(0,types,100)
    RRR_2vLib1_w <- matrix(0,types,100)
    RRR_2vLib1_m <- matrix(0,types,100)
    RRR_2vLib2_w <- matrix(0,types,100)
    RRR_2vLib2_m <- matrix(0,types,100)
    RRR_9v_w <- matrix(0,types,100)
    RRR_9v_m <- matrix(0,types,100)
    
    
    #Relative risk reductions (without waning)
    RRRs_w <- (HPVinc_pre_f_set - HPVinc_post_f)/HPVinc_pre_f_set
    RRRs_w[is.na(RRRs_w)] <- 0
    RRRs_m <- (HPVinc_pre_m_set - HPVinc_post_m)/HPVinc_pre_m_set
    RRRs_m[is.na(RRRs_m)] <- 0
    
    #plugging in the correct RRRs
    RRR_2vCons_w[c(1,2),] <- RRRs_w[c(1,2),]
    RRR_2vCons_m[c(1,2),] <- RRRs_m[c(1,2),]
    RRR_4v_w <- RRR_2vCons_w
    RRR_4v_m <- RRR_2vCons_m
    RRR_9v_w[c(1,2,3,4,7,9,10),] <- RRRs_w[c(1,2,6,10,16,19,21),]
    RRR_9v_m[c(1,2,3,4,7,9,10),] <- RRRs_m[c(1,2,6,10,16,19,21),]
    
    #Relative risk reductions (with waning)
    RRRs_w <- (HPVinc_pre_f_set - HPVinc_post_f_waning015)/HPVinc_pre_f_set
    RRRs_w[is.na(RRRs_w)] <- 0
    RRRs_m <- (HPVinc_pre_m_set - HPVinc_post_m_waning015)/HPVinc_pre_m_set
    RRRs_m[is.na(RRRs_m)] <- 0
    
    #plugging in the correct RRRs
    RRR_2vMid_w[c(1,2,3,4,7),] <- RRRs_w[c(1,2,4,8,13),]
    RRR_2vMid_m[c(1,2,3,4,7),] <- RRRs_m[c(1,2,4,8,13),]
    RRR_2vLib1_w[c(1,2,3,4,6,7,8),] <- RRRs_w[c(1,2,5,9,12,15,17),]
    RRR_2vLib1_m[c(1,2,3,4,6,7,8),] <- RRRs_m[c(1,2,5,9,12,15,17),]
    RRR_2vLib2_w[c(1,2,3,4,5,7,9,10),] <- RRRs_w[c(1,2,3,7,11,14,18,20),]
    RRR_2vLib2_m[c(1,2,3,4,5,7,9,10),] <- RRRs_m[c(1,2,3,7,11,14,18,20),]

    #For now
    RRR_2vCons_MSM <- RRR_2vCons_m
    RRR_2vMid_MSM <- RRR_2vMid_m
    RRR_2vLib1_MSM <- RRR_2vLib1_m
    RRR_2vLib2_MSM <- RRR_2vLib2_m
    RRR_4v_MSM <- RRR_4v_m
    RRR_9v_MSM <- RRR_9v_m
  }
  
  #Compute RRRs for the LR types
  if (sc_LR == 1){ #base-case scenario: effect equal to mean effect type 18
    RRR_611_w <- 0.84
    RRR_611_m <- 0.84
  } else if (sc_LR == 2){ #elimination for types 6 and 11
    RRR_611_w <- 1
    RRR_611_m <- 1 
  } else if (sc_LR == 3){ #direct protection only for types 6 and 11
    VE_611 <- 0.98
    RRR_611_w <- vacc.up[1]*VE_611
    RRR_611_m <- vacc.up[2]*VE_611
  }
  
  RRR_2vCons <- list(w=RRR_2vCons_w, m=RRR_2vCons_m, MSM=RRR_2vCons_MSM)
  RRR_2vMid <- list(w=RRR_2vMid_w, m=RRR_2vMid_m, MSM=RRR_2vMid_MSM)
  RRR_2vLib1 <- list(w=RRR_2vLib1_w, m=RRR_2vLib1_m, MSM=RRR_2vLib1_MSM)
  RRR_2vLib2 <- list(w=RRR_2vLib2_w, m=RRR_2vLib2_m, MSM=RRR_2vLib2_MSM)
  RRR_4v <- list(w=RRR_4v_w, m=RRR_4v_m, MSM=RRR_4v_MSM)
  RRR_9v <- list(w=RRR_9v_w, m=RRR_9v_m, MSM=RRR_9v_MSM)
  RRR_611 <- list(w=RRR_611_w, m=RRR_611_m)
  
  return(list(RRR_2vCons=RRR_2vCons,RRR_2vMid=RRR_2vMid,RRR_2vLib1=RRR_2vLib1,RRR_2vLib2=RRR_2vLib2,
              RRR_4v=RRR_4v,RRR_9v=RRR_9v,
              RRR_611=RRR_611,vacc.up=vacc.up))
}
###########################################################################################################

RiskReductionScreening <- function(scenario,agegrps.scr){
  #Risk reductions for screening based on risk reductions in HPV prevalence
  sc_HR <- scenario[1]
  
  if (sc_HR == 1){ 
    #scenario 1: All-or-nothing VE, life-long, no herd immunity
    #Risk reductions are the same
    
    RRscreening_2v.cons <- RRR_2v.cons_w[,agegrps.scr]
    RRscreening_2v.mid <- RRR_2v.mid_w[,agegrps.scr]
    RRscreening_2v.lib1 <- RRR_2v.lib1_w[,agegrps.scr]
    RRscreening_2v.lib2 <- RRR_2v.lib2_w[,agegrps.scr]
    RRscreening_4v <- RRR_4v_w[,agegrps.scr]
    RRscreening_9v <- RRR_9v_w[,agegrps.scr]
    
  } else if (sc_HR == 2){ 
    #scenario 2: All-or-nothing VE, life-long, + herd immunity
    
    #Relative risk reductions
    RRRprev_w <- (HPVprev_pre_f_set - HPVprev_post_f)/HPVprev_pre_f_set
    RRRprev_w[is.na(RRRprev_w)] <- 0
    
    #matrices for the different vaccines
    RRscreening_2v.cons <- matrix(0,types,7)
    RRscreening_2v.mid <- matrix(0,types,7)
    RRscreening_2v.lib1 <- matrix(0,types,7)
    RRscreening_2v.lib2 <- matrix(0,types,7)
    RRscreening_9v <- matrix(0,types,7)
    
    #plugging in the correct RRRs
    RRscreening_2v.cons[c(1,2),] <- RRRprev_w[c(1,2),agegrps.scr]
    RRscreening_2v.mid[c(1,2,3,4,7),] <- RRRprev_w[c(1,2,4,8,13),agegrps.scr]
    RRscreening_2v.lib1[c(1,2,3,4,6,7,8),] <- RRRprev_w[c(1,2,5,9,12,15,17),agegrps.scr]
    RRscreening_2v.lib2[c(1,2,3,4,5,7,9,10),] <- RRRprev_w[c(1,2,3,7,11,14,18,20),agegrps.scr]
    RRscreening_4v <- RRscreening_2v.cons
    RRscreening_9v[c(1,2,3,4,7,9,10),] <- RRRprev_w[c(1,2,6,10,16,19,21),agegrps.scr]
    
  } else if (sc_HR == 3){ 
    #scenario 3: All-or-nothing VE, waning, + herd immunity
    
    #Relative risk reductions
    RRRprev_w <- (HPVprev_pre_f_set - HPVprev_post_f_waning015)/HPVprev_pre_f_set
    RRRprev_w[is.na(RRRprev_w)] <- 0
    
    #matrices for the different vaccines
    RRscreening_2v.cons <- matrix(0,types,7)
    RRscreening_2v.mid <- matrix(0,types,7)
    RRscreening_2v.lib1 <- matrix(0,types,7)
    RRscreening_2v.lib2 <- matrix(0,types,7)
    RRscreening_9v <- matrix(0,types,7)
    
    #plugging in the correct RRRs
    RRscreening_2v.cons[c(1,2),] <- RRRprev_w[c(1,2),agegrps.scr]
    RRscreening_2v.mid[c(1,2,3,4,7),] <- RRRprev_w[c(1,2,4,8,13),agegrps.scr]
    RRscreening_2v.lib1[c(1,2,3,4,6,7,8),] <- RRRprev_w[c(1,2,5,9,12,15,17),agegrps.scr]
    RRscreening_2v.lib2[c(1,2,3,4,5,7,9,10),] <- RRRprev_w[c(1,2,3,7,11,14,18,20),agegrps.scr]
    RRscreening_4v <- RRscreening_2v.cons
    RRscreening_9v[c(1,2,3,4,7,9,10),] <- RRRprev_w[c(1,2,6,10,16,19,21),agegrps.scr]
    
  } else if (sc_HR == 4){ 
    #scenario 4: All-or-nothing VE, waning, + herd immunity
    
    #Relative risk reductions
    RRRprev_w <- (HPVprev_pre_f_set - HPVprev_post_f)/HPVprev_pre_f_set
    RRRprev_w[is.na(RRRprev_w)] <- 0
    
    #matrices for the different vaccines
    RRscreening_2v.cons <- matrix(0,types,7)
    RRscreening_9v <- matrix(0,types,7)
    
    #plugging in the correct RRRs
    RRscreening_2v.cons[c(1,2),] <- RRRprev_w[c(1,2),agegrps.scr]
    RRscreening_4v <- RRscreening_2v.cons
    RRscreening_9v[c(1,2,3,4,7,9,10),] <- RRRprev_w[c(1,2,6,10,16,19,21),agegrps.scr]
    
    #Relative risk reductions
    RRRprev_w <- (HPVprev_pre_f_set - HPVprev_post_f_waning015)/HPVprev_pre_f_set
    RRRprev_w[is.na(RRRprev_w)] <- 0
    
    #matrices for the different vaccines
    RRscreening_2v.mid <- matrix(0,types,7)
    RRscreening_2v.lib1 <- matrix(0,types,7)
    RRscreening_2v.lib2 <- matrix(0,types,7)
    
    #plugging in the correct RRRs
    RRscreening_2v.mid[c(1,2,3,4,7),] <- RRRprev_w[c(1,2,4,8,13),agegrps.scr]
    RRscreening_2v.lib1[c(1,2,3,4,6,7,8),] <- RRRprev_w[c(1,2,5,9,12,15,17),agegrps.scr]
    RRscreening_2v.lib2[c(1,2,3,4,5,7,9,10),] <- RRRprev_w[c(1,2,3,7,11,14,18,20),agegrps.scr]
  }
  
  return(list(RR_2v.cons=RRscreening_2v.cons,RR_2v.mid=RRscreening_2v.mid,
              RR_2v.lib1=RRscreening_2v.lib1,RR_2v.lib2=RRscreening_2v.lib2,
              RR_4v=RRscreening_4v,RR_9v=RRscreening_9v))
}

###########################################################################################################
#Two functions that we need
gamma.distr <- function(a,b,sh,sc){
  b1 <- b[1]
  b2 <- b[2]
  pgamma(b2-a,shape=sh,scale=sc) - pgamma(b1-a,shape=sh,scale=sc)
}

RRR.ca <- function(b,age_vec,prob_vec,risk.reduction_type,gamma.shape,gamma.scale,col.nr){
  gvec <- gamma.distr(age_vec,b,gamma.shape,gamma.scale)
  cond.prob <- gvec*prob_vec/(sum(gvec*prob_vec))
  
  #lines(age_vec,cond.prob,col=col.nr)
  
  RRR.ca_value <- sum(cond.prob*risk.reduction_type[age_vec])
  return(RRR.ca_value)      
}

RiskReductionCancer <- function(risk.reduction,gamma.pars,HPV.inc){
  #Computation of probability distribution of age of infection given age of cancer diagnosis
  #and Risk reductions on cancer level
  
  age_vec <- 1:100 #doesn't matter that we include ages after 89
  B_ca <- cbind(as.vector(ca.ages),ca.ages+5)
  agegroups <- c('15-19','20-24','25-29','30-34','35-39','40-44','45-49','50-54','55-59','60-64','65-69','70-74','75-79','80-84','85-89')
  RRRmatrix.ca <- matrix(0,types,length(ca.ages))
  
  shape.par <- gamma.pars[1]
  scale.par <- gamma.pars[2]
  
  if (scenario[1] == 1){ #no herd immunity, no waning
    RRRmatrix.ca <- risk.reduction[,ca.ages]
  } else {
    #We compute the new RRRs (per cancer age-group) for each HPV type
    for (j in 1:types){
      
      #Distribution (discrete) of T1 (age of HPV infection, given you have infection before 90)
      prob_vec <- HPV.inc[j,age_vec]/sum(HPV.inc[j,age_vec])
      
      #plot(age_vec,prob_vec,ylim=c(0,0.3),xlim=c((age_vec[1]),tail(age_vec,1)),type='l',lwd=2,
           #xlab="age of infection",ylab="probability",
           #main='Age of infection given age of cancer')
      
      for (i in 1:length(ca.ages)){
        b <- B_ca[i,]
        RRRmatrix.ca[j,i] <- RRR.ca(b,age_vec,prob_vec,risk.reduction[j,],
                                    gamma.shape=shape.par,gamma.scale=scale.par,col.nr=i)
      }
      
      #legend("topright", legend = agegroups, col=1:length(ca.ages), pch=1, title='cancer age')
    }
  }
  
  return(RRRmatrix.ca)
}
