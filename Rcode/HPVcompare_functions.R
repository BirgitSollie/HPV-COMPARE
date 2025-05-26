#Functions needed for the computation

###########################################################################################################

WartsExpected <- function(AW.inc,cost_AW,S){
  
  #Warts diagnoses expected for the cohort
  nr.AW_exp <- cohort*S[1,]*AW.inc[a0:amax]*PAF.HPV_AW
  costs.AW_exp <- c(1,dvec_c)*cost_AW*nr.AW_exp
  QALYs.AW_exp <- c(1,dvec_h)*QALY_AW*nr.AW_exp
  
  list(nr_exp=nr.AW_exp,costs_exp=costs.AW_exp,QALYs_exp=QALYs.AW_exp)
}

###########################################################################################################

RRPExpected <- function(IncUnder18_RRP,IncAbove18_RRP,cost_RRP,QALY_RRP,S,gender){
  
  pop.size_age <- cohort*S[1,(18-a0+1):(amax-a0+1)]
  RRPdiags_adults <- sum(pop.size_age)*IncAbove18_RRP/100000
  nr.RRPpatients_exp <- PAF.HPV_RRP*RRPdiags_adults*(age.distr_RRP[18:amax]/sum(age.distr_RRP[18:amax]))

  #Make big Exp_prob_matrix
  x <- 0:(amax-1)
  Exp_probs <- exp(-x/nrYearsAverage_RRP)
  
  Exp_prob_mat <- matrix(rep(c(Exp_probs,0),(amax-1+1)),amax-1+1,amax-1+2)
  Exp_prob_mat[upper.tri(Exp_prob_mat,diag=FALSE)]=0 #set upper triangle to 0
  
  Exp_prob_mat_adults <- Exp_prob_mat[1:(amax-18+1),1:(amax-18+1)] #for ages 18-amax
  
  nr.RRP_exp <- apply(t(matrix(rep(nr.RRPpatients_exp,amax-18+1),amax-18+1,amax-18+1))*Exp_prob_mat_adults,1,sum)
  
  costs.RRP_exp <- cost_RRP*dvec_c[(18-a0):(amax-a0)]*nr.RRP_exp
  QALYs.RRP_exp <- QALY_RRP*dvec_h[(18-a0):(amax-a0)]*nr.RRP_exp
  
  nr.exp_total <- RRPdiags_adults
  
  #incidence children (assumption: via birth) #assume PAF.HPV = 1
  if (gender=="w"){
    nr.child_exp <- cohort*S[1,(16-a0+1):(49-a0+1)]*as.numeric(births_age_w[2:35,2]/1000)
    nr.RRPchildpatients_exp <- sum(sum(nr.child_exp)*(0.5*surv_m[1:18,]/100000 + 0.5*surv_w[1:18,]/100000))*IncUnder18_RRP/100000
    nr.RRPchildpatients_exp_agemother <- nr.RRPchildpatients_exp*nr.child_exp/sum(nr.child_exp)
    
    Exp_prob_mat_children <- Exp_prob_mat[,1:17] #for ages 1-17
    nr.RRPchildpatients_exp_age <- nr.RRPchildpatients_exp*(age.distr_RRP[1:17]/sum(age.distr_RRP[1:17]))
    nr.RRPchild_exp <- apply(t(matrix(rep(nr.RRPchildpatients_exp_age,amax),17,amax))*Exp_prob_mat_children,1,sum)
                  
    #discounting from age of mother
    dvec_c_long <- 1/(d_costs^(1:150))
    dvec_h_long <- 1/(d_health^(1:150))
    cost.discount_matrix <- matrix(0,99,34)
    health.discount_matrix <- matrix(0,99,34)
    for (i in 1:34){
      cost.discount_matrix[,i] <- dvec_c_long[(17+i-1-a0):(115+i-1-a0)]
      health.discount_matrix[,i] <- dvec_h_long[(17+i-1-a0):(115+i-1-a0)]
    }
    
    age.dist_mother <- nr.child_exp/sum(nr.child_exp)
    discount.age_matrix <- nr.RRPchild_exp %*% t(age.dist_mother)
    
    #discounted effects per age of mother
    costs.RRPchild_exp <- apply(discount.age_matrix*cost.discount_matrix*cost_RRP,2,sum)
    QALYs.RRPchild_exp <- apply(discount.age_matrix*cost.discount_matrix*QALY_RRP,2,sum)
    
    #age 16-99
    nr.RRPpatients_exp <- c(0,0,nr.RRPpatients_exp) + c(nr.RRPchildpatients_exp_agemother,numeric(50))
    costs.RRP_exp <- c(0,0,costs.RRP_exp) + c(costs.RRPchild_exp,numeric(50))
    QALYs.RRP_exp <- c(0,0,QALYs.RRP_exp) + c(QALYs.RRPchild_exp,numeric(50))
    
    nr.exp_total <- c(RRPdiags_adults,nr.RRPchildpatients_exp)
  }
  
  list(nr_exp=nr.RRPpatients_exp,costs_exp=costs.RRP_exp,QALYs_exp=QALYs.RRP_exp,exp_total=nr.exp_total)
}
###########################################################################################################

ScreeningExpected <- function(scr.up,fr.lesions,fr.cancer,AFs_CIN2plus_age1,AFs_CIN2plus_age2,AFs_cervix){
  #Computed the expected number of CIN2+ diagnosed in screening in a pre-vaccination setting.
  
  nr.agegrps.scr <- length(agegrps.scr)
  
  #AFs
  AFs_agegrps <- matrix(0,types,nr.agegrps.scr)
  AFs_agegrps[,1] <-  AFs_CIN2plus_age1
  AFs_agegrps[,2:nr.agegrps.scr] <- rep(AFs_CIN2plus_age2,(nr.agegrps.scr -1))
  AFs_agegrps_cascr <- matrix(rep(AFs_cervix,nr.agegrps.scr),types,nr.agegrps.scr)
  
  #percentage of a0 year-old cohort that live up to certain age group
  fr.surv_scr <- S_w[1,(agegrps.scr-a0+1)]
  
  #number of women in screening per age grp
  fr.screening <- fr.surv_scr*scr.up

  #Expected number of lesions pre vaccination, per age group and type
  nr.CIN2plus_exp <- cohort*(AFs_agegrps%*%diag(fr.lesions*fr.screening))
  nr.cancer_exp <- cohort*(AFs_agegrps_cascr%*%diag(fr.cancer*fr.screening))
  
  cost_CIN2plus_age <- dvec_c[(agegrps.scr-a0)]*cost_CIN2plus
  QALY_CIN2plus_age <- dvec_h[(agegrps.scr-a0)]*QALY_CIN2plus
  
  costs_CIN2plus_exp <- nr.CIN2plus_exp%*%diag(cost_CIN2plus_age)
  QALYs_CIN2plus_exp <- nr.CIN2plus_exp%*%diag(QALY_CIN2plus_age)
  
  #unnecessary colposcopies (leading to CIN0/1 diagnoses)
  nr.colpo_exp <- apply(nr.CIN2plus_exp,2,sum)/PPV_colpo_no
  cost_colpo_age <- dvec_c[(agegrps.scr-a0)]*cost_colpo
  costs_colpo_exp <- nr.colpo_exp*cost_colpo_age
  
  list(nr.les_exp=nr.CIN2plus_exp,nr.ca_exp=nr.cancer_exp,
       costs_exp=costs_CIN2plus_exp,QALYs_exp=QALYs_CIN2plus_exp,
       nr.colpo_exp=nr.colpo_exp,costs.colpo_exp=costs_colpo_exp,cost.colpo_vec=cost_colpo_age)
}
###########################################################################################################

CancerExpected <- function(S,ca.risk,ca_survdata,AFs_types,costs_ca,QALY_ca,PAF,HR){
  
  #Risk for time spend in the age group (approximately)
  ca.risk_grp <- 5*ca.risk
  
  #Probability of surviving to certain age group
  pr.surv_ca <- S[1,(agemids.ca-a0+1)]
  
  AFs_ca <- PAF*AFs_types
  
  nr.cancers_exp <- cohort*(AFs_ca %*% t(ca.risk_grp*pr.surv_ca))
  
  #Life-years lost (from age a0 till max(agemids))
  Lifeyearslost <- LifeyearsLost(max(agemids.ca),S,ca_survdata,PAF,HR)
  LYs.loss_ca <- Lifeyearslost$lifey.loss_ca
  ca_survdata_new <- Lifeyearslost$ca_survdata_new
  
  LYs.loss_age <- LYs.loss_ca[(agemids.ca-a0+1)]*dvec_h[(agemids.ca-a0)]
  LYs.loss_exp <- nr.cancers_exp%*%diag(LYs.loss_age)
  
  surv_grps <- c(rep(1,9),rep(2,2),rep(3,2),rep(4,2),rep(5,3))
  prob.surv_ca <- ca_survdata_new[surv_grps[agegrps.ca],10]
  
  #Costs treatment and palliative care
  costs1_ca_age <- dvec_c[(agemids.ca-a0)]*costs_ca[1]
  costs2_ca_age <- (1-prob.surv_ca)*dvec_c[(agemids.ca-a0)]*costs_ca[2]

  costs_ca1_exp <- nr.cancers_exp%*%diag(costs1_ca_age)
  costs_ca2_exp <- nr.cancers_exp%*%diag(costs2_ca_age)
  
  #Additional QALY loss per cancer diagnosis (irrespective of survival)
  QALYloss <- sum(nr.cancers_exp)*QALY_ca
  LYs.loss_exp <- LYs.loss_exp + QALYloss
  
  list(nr.ca_exp=nr.cancers_exp,costs1_exp=costs_ca1_exp,costs2_exp=costs_ca2_exp,LYs.loss_exp=LYs.loss_exp)
}

###########################################################################################################

LifeyearsLost <- function(amax_ca,S,ca_survdata,PAF,HR){
  #Computing for each age a (from a0 to amax_ca) the expected loss in LY if cancer diagnosed at age a 
  # = Expected (mean) lifetime (normal) - Expected (mean) lifetime after cancer diagnosis
  
  surv.age <- S[1:(amax_ca-a0+1),1:(amax-a0+1)] #normal survival
  surv.age_ca <- matrix(0,amax_ca-a0+1,amax-a0+1)
  
  #Increased survival for HPV+ cancer
  ca_survdata1 <- cbind(rep(1,5),ca_survdata)
  ca_survdata_new <- matrix(sapply(ca_survdata1,Increased_Surv_HPV,alpha=PAF,HR=HR),5,11)
  
  for (a in a0:amax_ca){
    surv.vec <- surv.age[a-a0+1,(a-a0+1):(amax-a0+1)] 
    
    #survival after cancer diagnosis (need to survive cancer AND other causes)
    if (a <= 44){
      surv.vec_ca <- c(ca_survdata_new[1,],rep(ca_survdata_new[1,11],amax-a-10))
      surv.age_ca[a-a0+1,(a-a0+1):(amax-a0+1)] <- surv.vec*surv.vec_ca
    } else if (a >= 45 & a <= 54){
      surv.vec_ca <- c(ca_survdata_new[2,],rep(ca_survdata_new[2,11],amax-a-10))
      surv.age_ca[a-a0+1,(a-a0+1):(amax-a0+1)] <- surv.vec*surv.vec_ca
    } else if (a >= 55 & a <= 64){
      surv.vec_ca <- c(ca_survdata_new[3,],rep(ca_survdata_new[3,11],amax-a-10))
      surv.age_ca[a-a0+1,(a-a0+1):(amax-a0+1)] <- surv.vec*surv.vec_ca
    } else if (a >= 65 & a<= 74){
      surv.vec_ca <- c(ca_survdata_new[4,],rep(ca_survdata_new[4,11],amax-a-10))
      surv.age_ca[a-a0+1,(a-a0+1):(amax-a0+1)] <- surv.vec*surv.vec_ca
    } else if (a >= 75 & a<= amax-10){
      surv.vec_ca <- c(ca_survdata_new[5,],rep(ca_survdata_new[5,11],amax-a-10))
      surv.age_ca[a-a0+1,(a-a0+1):(amax-a0+1)] <- surv.vec*surv.vec_ca
    } else if (a > amax-10){
      surv.vec_ca <- ca_survdata_new[5,1:(amax-a+1)]
      surv.age_ca[a-a0+1,(a-a0+1):(amax-a0+1)] <- surv.vec*surv.vec_ca
    }
  }
  
  #Expected number of LY
  #Note that we count the 1s on the diagonal, but we do that in both Exp.surv and Exp.surv_ca
  Exp.surv <- apply(surv.age,1,sum)
  Exp.surv_ca <- apply(surv.age_ca,1,sum)
  
  #Expected loss in life-years after cancer diagnosis, irrespective of surviving cancer or not
  lifey.loss_ca <- Exp.surv - Exp.surv_ca
  
  list(lifey.loss_ca=lifey.loss_ca,ca_survdata_new=ca_survdata_new)
}
###########################################################################################################

# Functions to compute the improved survival for HPV+ cancers
Surv_mixture <- function(s_HPV,s_i,alpha,HR){
  alpha*s_HPV + (1-alpha)*(s_HPV^(1/HR)) - s_i
}
###########################################################################################################

Increased_Surv_HPV <- function(s_i,alpha,HR){
  uniroot(Surv_mixture,c(0,1),tol=0.0001,s_i=s_i,alpha=alpha,HR=HR)[[1]]
}
###########################################################################################################
