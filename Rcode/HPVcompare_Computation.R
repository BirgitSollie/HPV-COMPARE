#Computation of costs, effects and ICERs

#Lists with empty matrices for data storage
scr.names <- c("colpos","lesions","cancer")
scr.prev_2v.cons <- lapply(1:3, matrix, data= 0, nrow=nr.sim, ncol=7)
scr.prev_2v.mid <- lapply(1:3, matrix, data= 0, nrow=nr.sim, ncol=7)
scr.prev_2v.lib1 <- lapply(1:3, matrix, data= 0, nrow=nr.sim, ncol=7)
scr.prev_2v.lib2 <- lapply(1:3, matrix, data= 0, nrow=nr.sim, ncol=7)
scr.prev_4v <- lapply(1:3, matrix, data= 0, nrow=nr.sim, ncol=7)
scr.prev_9v <- lapply(1:3, matrix, data= 0, nrow=nr.sim, ncol=7)
names(scr.prev_2v.cons) <- scr.names
names(scr.prev_2v.mid) <- scr.names
names(scr.prev_2v.lib1) <- scr.names
names(scr.prev_2v.lib2) <- scr.names
names(scr.prev_4v) <- scr.names
names(scr.prev_9v) <- scr.names

ca.names <- c("cervix","anus_w","anus_m","oroph_w","oroph_m","vulva","vagina","penis")
ca.prev_2v.cons <- lapply(1:8, matrix, data= 0, nrow=nr.sim, ncol=15)
ca.prev_2v.mid <- lapply(1:8, matrix, data= 0, nrow=nr.sim, ncol=15)
ca.prev_2v.lib1 <- lapply(1:8, matrix, data= 0, nrow=nr.sim, ncol=15)
ca.prev_2v.lib2 <- lapply(1:8, matrix, data= 0, nrow=nr.sim, ncol=15)
ca.prev_4v <- lapply(1:8, matrix, data= 0, nrow=nr.sim, ncol=15)
ca.prev_9v <- lapply(1:8, matrix, data= 0, nrow=nr.sim, ncol=15)
names(ca.prev_2v.cons) <- ca.names
names(ca.prev_2v.mid) <- ca.names
names(ca.prev_2v.lib1) <- ca.names
names(ca.prev_2v.lib2) <- ca.names
names(ca.prev_4v) <- ca.names
names(ca.prev_9v) <- ca.names

#total nr prevented, QALYs gained, savings treatment, savings pall.care
Effect.names <- c("warts_w","warts_m","RRP_w","RRP_m","colpos","lesions","cervix",
                  "anus_w","anus_m","oroph_w","oroph_m","vulva","vagina","penis")
Effects_2v.cons <- lapply(1:14, matrix, data= 0, nrow=nr.sim, ncol=4)
Effects_2v.mid <- lapply(1:14, matrix, data= 0, nrow=nr.sim, ncol=4)
Effects_2v.lib1 <- lapply(1:14, matrix, data= 0, nrow=nr.sim, ncol=4)
Effects_2v.lib2 <- lapply(1:14, matrix, data= 0, nrow=nr.sim, ncol=4)
Effects_4v <- lapply(1:14, matrix, data= 0, nrow=nr.sim, ncol=4)
Effects_9v <- lapply(1:14, matrix, data= 0, nrow=nr.sim, ncol=4)
names(Effects_2v.cons) <- Effect.names
names(Effects_2v.mid) <- Effect.names
names(Effects_2v.lib1) <- Effect.names
names(Effects_2v.lib2) <- Effect.names
names(Effects_4v) <- Effect.names
names(Effects_9v) <- Effect.names

Expected.total <- matrix(0,nr.sim,14)
colnames(Expected.total) <- Effect.names

#Vector just to save the total number of juvenile onset RRP
RRP_juv <- rep(0,nr.sim)

#matrix to save number of colpos expected in scenario with no vaccination
colpos_exp_no <- matrix(0,nr.sim,7)

for (i in 1:nr.sim){

  #############
  ### WARTS ###
  #############
  
  AW_inc_w <- AW.inc_w_sim[i,]
  AW_inc_m <- AW.inc_m_sim[i,]
  Expected_warts_w <- WartsExpected(AW_inc_w,cost_AW_w,S_w)
  Expected_warts_m <- WartsExpected(AW_inc_w,cost_AW_m,S_m)

  Expected.total[i,"warts_w"] <- sum(Expected_warts_w$nr_exp)
  Expected.total[i,"warts_m"] <- sum(Expected_warts_m$nr_exp)
  
  ###########
  ### RRP ###
  ###########
  
  PrevUnder18_RRP <- PrevUnder18_RRP_sim[i]
  PrevAbove18_RRP <- PrevAbove18_RRP_sim[i]
  IncUnder18_RRP <- PrevUnder18_RRP/6.94 #=125 cohorts/18 cohorts
  IncAbove18_RRP <- PrevAbove18_RRP/nrYearsAverage_RRP
  Expected_RRP_w <- RRPExpected(IncUnder18_RRP,IncAbove18_RRP,cost_RRP,QALY_RRP,S_w,"w")
  RRP_juv[i] <- Expected_RRP_w$exp_total[2]
  Expected_RRP_m <- RRPExpected(IncUnder18_RRP,IncAbove18_RRP,cost_RRP,QALY_RRP,S_m,"m")
  
  Expected.total[i,"RRP_w"] <- sum(Expected_RRP_w$nr_exp)
  Expected.total[i,"RRP_m"] <- sum(Expected_RRP_m$nr_exp)
    
  #################
  ### SCREENING ###
  #################

  scr.up <- scr.up_sim[i,]
  fr.lesions <- fr.lesions_sim[i,]
  fr.cancer <- fr.cancer_sim[i,]

  AFs_CIN2plus_age1 <- AFs_CIN2plus_age1_sim[i,]
  AFs_CIN2plus_age2 <- AFs_CIN2plus_age2_sim[i,]
  AFs_cervix <- AFs_cervix_sim[i,]

  Expected_screening <- ScreeningExpected(scr.up,fr.lesions,fr.cancer,AFs_CIN2plus_age1,AFs_CIN2plus_age2,AFs_cervix)
  colpos_exp_no[i,] <- Expected_screening$nr.colpo_exp
  Expected.total[i,"lesions"] <- sum(Expected_screening$nr.les_exp)
  Expected.total[i,"colpos"] <- sum(Expected_screening$nr.colpo_exp)

  ##############
  ### CANCER ###
  ##############

  #Cervix
  ca.risk_cervix <- ca.risk_cervix_sim[i,]
  surv_cervix <- t(apply(matrix(surv_probs_cervix_sim[i,],5,10),1,cumprod))
  Expected_cervix <- CancerExpected(S_w,ca.risk_cervix,surv_cervix,AFs_cervix,costs_cervix,QALY_cervix,PAF.HPV_cervix,HR_cervix)
  Expected.total[i,"cervix"] <- sum(Expected_cervix$nr.ca_exp)
  
  #Anus
  ca.risk_anus_w <- ca.risk_anus_w_sim[i,]
  ca.risk_anus_m <- ca.risk_anus_m_sim[i,]
  surv_anus <- t(apply(matrix(surv_probs_anus_sim[i,],5,10),1,cumprod))
  PAF.HPV_anus <- PAF.HPV_anus_sim[i]
  AFs_anus <- AFs_anus_sim[i,]
  HR_anus <- HR_anus_sim[i]
  Expected_anus_w <- CancerExpected(S_w,ca.risk_anus_w,surv_anus,AFs_anus,costs_anus_w,QALY_anus,PAF.HPV_anus,HR_anus)
  Expected_anus_m <- CancerExpected(S_m,ca.risk_anus_m,surv_anus,AFs_anus,costs_anus_m,QALY_anus,PAF.HPV_anus,HR_anus)
  Expected.total[i,"anus_w"] <- sum(Expected_anus_w$nr.ca_exp)
  Expected.total[i,"anus_m"] <- sum(Expected_anus_m$nr.ca_exp)
  
  #Oropharynx
  ca.risk_oroph_w <- ca.risk_oroph_w_sim[i,]
  ca.risk_oroph_m <- ca.risk_oroph_m_sim[i,]
  surv_oroph <- t(apply(matrix(surv_probs_oroph_sim[i,],5,10),1,cumprod))
  PAF.HPV_oroph <- PAF.HPV_oroph_sim[i]
  AFs_oroph <- AFs_oroph_sim[i,]
  HR_oroph <- HR_oroph_sim[i]
  Expected_oroph_w <- CancerExpected(S_w,ca.risk_oroph_w,surv_oroph,AFs_oroph,costs_oroph_w,QALY_oroph,PAF.HPV_oroph,HR_oroph)
  Expected_oroph_m <- CancerExpected(S_m,ca.risk_oroph_m,surv_oroph,AFs_oroph,costs_oroph_m,QALY_oroph,PAF.HPV_oroph,HR_oroph)
  Expected.total[i,"oroph_w"] <- sum(Expected_oroph_w$nr.ca_exp)
  Expected.total[i,"oroph_m"] <- sum(Expected_oroph_m$nr.ca_exp)
  
  #Vulva
  ca.risk_vulva <- ca.risk_vulva_sim[i,]
  surv_vulva <- t(apply(matrix(surv_probs_vulva_sim[i,],5,10),1,cumprod))
  PAF.HPV_vulva <- PAF.HPV_vulva_sim[i]
  AFs_vulva <- AFs_vulva_sim[i,]
  HR_vulva <- HR_vulva_sim[i]
  Expected_vulva <- CancerExpected(S_w,ca.risk_vulva,surv_vulva,AFs_vulva,costs_vulva,QALY_vulva,PAF.HPV_vulva,HR_vulva)
  Expected.total[i,"vulva"] <- sum(Expected_vulva$nr.ca_exp)
  
  #Vagina
  ca.risk_vagina <- ca.risk_vagina_sim[i,]
  surv_vagina <- t(apply(matrix(surv_probs_vagina_sim[i,],5,10),1,cumprod))
  PAF.HPV_vagina <- PAF.HPV_vagina_sim[i]
  AFs_vagina <- AFs_vagina_sim[i,]
  Expected_vagina <- CancerExpected(S_w,ca.risk_vagina,surv_vagina,AFs_vagina,costs_vagina,QALY_vagina,PAF.HPV_vagina,HR_vagina)
  Expected.total[i,"vagina"] <- sum(Expected_vagina$nr.ca_exp)
  
  #Penis
  ca.risk_penis <- ca.risk_penis_sim[i,]
  surv_penis <- t(apply(matrix(surv_probs_penis_sim[i,],5,10),1,cumprod))
  PAF.HPV_penis <- PAF.HPV_penis_sim[i]
  AFs_penis <- AFs_penis_sim[i,]
  HR_penis <- HR_penis_sim[i]
  Expected_penis <- CancerExpected(S_w,ca.risk_penis,surv_penis,AFs_penis,costs_penis,QALY_penis,PAF.HPV_penis,HR_penis)
  Expected.total[i,"penis"] <- sum(Expected_penis$nr.ca_exp)
  
  ###########################
  ### EFFECTS PER VACCINE ###
  ###########################
  
  ### WARTS ###
  
  #total nr prevented, QALYs gained, savings treatment, savings pall.care
  Effects_2v.cons$warts_w[i,] <- c(0,0,0,0)
  Effects_2v.mid$warts_w[i,] <- c(0,0,0,0)
  Effects_2v.lib1$warts_w[i,] <- c(0,0,0,0)
  Effects_2v.lib2$warts_w[i,] <- c(0,0,0,0)
  Effects_4v$warts_w[i,] <- c(RRR_611_w*sum(Expected_warts_w$nr_exp),
                              RRR_611_w*sum(Expected_warts_w$QALYs_exp),
                              RRR_611_w*sum(Expected_warts_w$costs_exp),
                              0)
  Effects_2v.cons$warts_m[i,] <- c(0,0,0,0)
  Effects_2v.mid$warts_m[i,] <- c(0,0,0,0)
  Effects_2v.lib1$warts_m[i,] <- c(0,0,0,0)
  Effects_2v.lib2$warts_m[i,] <- c(0,0,0,0)
  Effects_4v$warts_m[i,] <- c(RRR_611_m*sum(Expected_warts_m$nr_exp),
                              RRR_611_m*sum(Expected_warts_m$QALYs_exp),
                              RRR_611_m*sum(Expected_warts_m$costs_exp),
                              0)
  if (scenario[3] == 0 | scenario[3] == 1){
    Effects_4v$warts_w[i,] <- c(0,0,0,0)
    Effects_4v$warts_m[i,] <- c(0,0,0,0)
  }
  Effects_9v$warts_w[i,] <- Effects_4v$warts_w[i,]
  Effects_9v$warts_m[i,] <- Effects_4v$warts_m[i,]
  
  ### RRP ###
  
  #total nr prevented, QALYs gained, savings treatment, savings pall.care
  Effects_2v.cons$RRP_w[i,] <- c(0,0,0,0)
  Effects_2v.mid$RRP_w[i,] <- c(0,0,0,0)
  Effects_2v.lib1$RRP_w[i,] <- c(0,0,0,0)
  Effects_2v.lib2$RRP_w[i,] <- c(0,0,0,0)
  Effects_4v$RRP_w[i,] <- c(RRR_611_w*sum(Expected_RRP_w$nr_exp),
                              RRR_611_w*sum(Expected_RRP_w$QALYs_exp),
                              RRR_611_w*sum(Expected_RRP_w$costs_exp),
                              0)
  Effects_2v.cons$RRP_m[i,] <- c(0,0,0,0)
  Effects_2v.mid$RRP_m[i,] <- c(0,0,0,0)
  Effects_2v.lib1$RRP_m[i,] <- c(0,0,0,0)
  Effects_2v.lib2$RRP_m[i,] <- c(0,0,0,0)
  Effects_4v$RRP_m[i,] <- c(RRR_611_m*sum(Expected_RRP_m$nr_exp),
                            RRR_611_m*sum(Expected_RRP_m$QALYs_exp),
                            RRR_611_m*sum(Expected_RRP_m$costs_exp),
                            0)
  if (scenario[3] == 0 | scenario[3] == 2){
    Effects_4v$RRP_w[i,] <- c(0,0,0,0)
    Effects_4v$RRP_m[i,] <- c(0,0,0,0)
  }
  Effects_9v$RRP_w[i,] <- Effects_4v$RRP_w[i,]
  Effects_9v$RRP_m[i,] <- Effects_4v$RRP_m[i,]
  
  ### SCREENING ###
  
  #nr lesions prevented
  scr.prev_2v.cons$lesions[i,] <- apply(Expected_screening$nr.les_exp*RRscr_2v.cons,2,sum)
  scr.prev_2v.mid$lesions[i,] <- apply(Expected_screening$nr.les_exp*RRscr_2v.mid,2,sum)
  scr.prev_2v.lib1$lesions[i,] <- apply(Expected_screening$nr.les_exp*RRscr_2v.lib1,2,sum)
  scr.prev_2v.lib2$lesions[i,] <- apply(Expected_screening$nr.les_exp*RRscr_2v.lib2,2,sum)
  scr.prev_4v$lesions[i,] <- apply(Expected_screening$nr.les_exp*RRscr_4v,2,sum)
  scr.prev_9v$lesions[i,] <- apply(Expected_screening$nr.les_exp*RRscr_9v,2,sum)
  
  #nr colposcopies prevented
  scr.prev_2v.cons$colpos[i,] <- Expected_screening$nr.colpo_exp - 
  ((apply(Expected_screening$nr.les_exp,2,sum) - scr.prev_2v.cons$lesions[i,])/
    PPV_colpo_2v)
  scr.prev_2v.mid$colpos[i,] <- Expected_screening$nr.colpo_exp - 
    ((apply(Expected_screening$nr.les_exp,2,sum) - scr.prev_2v.mid$lesions[i,])/
       PPV_colpo_2v.cp)
  scr.prev_2v.lib1$colpos[i,] <- Expected_screening$nr.colpo_exp - 
    ((apply(Expected_screening$nr.les_exp,2,sum) - scr.prev_2v.lib1$lesions[i,])/
       PPV_colpo_2v.cp)
  scr.prev_2v.lib2$colpos[i,] <- Expected_screening$nr.colpo_exp - 
    ((apply(Expected_screening$nr.les_exp,2,sum) - scr.prev_2v.lib2$lesions[i,])/
       PPV_colpo_2v.cp)
  scr.prev_4v$colpos[i,] <- Expected_screening$nr.colpo_exp - 
    ((apply(Expected_screening$nr.les_exp,2,sum) - scr.prev_4v$lesions[i,])/
       PPV_colpo_2v)
  scr.prev_9v$colpos[i,] <- Expected_screening$nr.colpo_exp - 
    ((apply(Expected_screening$nr.les_exp,2,sum) - scr.prev_9v$lesions[i,])/
       PPV_colpo_9v)
  
  #cancer prevented in screening
  scr.prev_2v.cons$cancer[i,] <- apply(Expected_screening$nr.ca_exp*RRscr_2v.cons,2,sum)
  scr.prev_2v.mid$cancer[i,] <- apply(Expected_screening$nr.ca_exp*RRscr_2v.mid,2,sum)
  scr.prev_2v.lib1$cancer[i,] <- apply(Expected_screening$nr.ca_exp*RRscr_2v.lib1,2,sum)
  scr.prev_2v.lib2$cancer[i,] <- apply(Expected_screening$nr.ca_exp*RRscr_2v.lib2,2,sum)
  scr.prev_4v$cancer[i,] <- apply(Expected_screening$nr.ca_exp*RRscr_4v,2,sum)
  scr.prev_9v$cancer[i,] <- apply(Expected_screening$nr.ca_exp*RRscr_9v,2,sum)
  
  #total nr prevented, QALYs gained, savings treatment, savings pall.care =0 
  #(we do not include the cancer cases here)
  Effects_2v.cons$colpos[i,] <- c(sum(scr.prev_2v.cons$colpos[i,]),0,
                                  sum(scr.prev_2v.cons$colpos[i,]*Expected_screening$cost.colpo_vec),0)
  Effects_2v.mid$colpos[i,] <- c(sum(scr.prev_2v.mid$colpos[i,]),0,
                                  sum(scr.prev_2v.mid$colpos[i,]*Expected_screening$cost.colpo_vec),0)
  Effects_2v.lib1$colpos[i,] <- c(sum(scr.prev_2v.lib1$colpos[i,]),0,
                                 sum(scr.prev_2v.lib1$colpos[i,]*Expected_screening$cost.colpo_vec),0)
  Effects_2v.lib2$colpos[i,] <- c(sum(scr.prev_2v.lib2$colpos[i,]),0,
                                 sum(scr.prev_2v.lib2$colpos[i,]*Expected_screening$cost.colpo_vec),0)
  Effects_4v$colpos[i,] <- c(sum(scr.prev_4v$colpos[i,]),0,
                                 sum(scr.prev_4v$colpos[i,]*Expected_screening$cost.colpo_vec),0)
  Effects_9v$colpos[i,] <- c(sum(scr.prev_9v$colpos[i,]),0,
                                 sum(scr.prev_9v$colpos[i,]*Expected_screening$cost.colpo_vec),0)
  
  Effects_2v.cons$lesions[i,] <- c(sum(scr.prev_2v.cons$lesions[i,]),
                                   sum(Expected_screening$QALYs_exp*RRscr_2v.cons,na.rm=TRUE),
                                   sum(Expected_screening$costs_exp*RRscr_2v.cons,na.rm=TRUE),
                                   0)
  Effects_2v.mid$lesions[i,] <- c(sum(scr.prev_2v.mid$lesions[i,]),
                                  sum(Expected_screening$QALYs_exp*RRscr_2v.mid,na.rm=TRUE),
                                  sum(Expected_screening$costs_exp*RRscr_2v.mid,na.rm=TRUE),
                                  0)
  Effects_2v.lib1$lesions[i,] <- c(sum(scr.prev_2v.lib1$lesions[i,]),
                                  sum(Expected_screening$QALYs_exp*RRscr_2v.lib1,na.rm=TRUE),
                                  sum(Expected_screening$costs_exp*RRscr_2v.lib1,na.rm=TRUE),
                                  0)
  Effects_2v.lib2$lesions[i,] <- c(sum(scr.prev_2v.lib2$lesions[i,]),
                                  sum(Expected_screening$QALYs_exp*RRscr_2v.lib2,na.rm=TRUE),
                                  sum(Expected_screening$costs_exp*RRscr_2v.lib2,na.rm=TRUE),
                                  0)
  Effects_4v$lesions[i,] <- c(sum(scr.prev_4v$lesions[i,]),
                              sum(Expected_screening$QALYs_exp*RRscr_4v,na.rm=TRUE),
                              sum(Expected_screening$costs_exp*RRscr_4v,na.rm=TRUE),
                              0)
  Effects_9v$lesions[i,] <- c(sum(scr.prev_9v$lesions[i,]),
                              sum(Expected_screening$QALYs_exp*RRscr_9v,na.rm=TRUE),
                              sum(Expected_screening$costs_exp*RRscr_9v,na.rm=TRUE),
                              0)
  
  ### CANCER ###
  
  #Cervix
  ca.prev_2v.cons$cervix[i,] <- apply(Expected_cervix$nr.ca_exp*RRRca_2v.cons$cervix,2,sum,na.rm=TRUE)
  ca.prev_2v.mid$cervix[i,] <- apply(Expected_cervix$nr.ca_exp*RRRca_2v.mid$cervix,2,sum,na.rm=TRUE)
  ca.prev_2v.lib1$cervix[i,] <- apply(Expected_cervix$nr.ca_exp*RRRca_2v.lib1$cervix,2,sum,na.rm=TRUE)
  ca.prev_2v.lib2$cervix[i,] <- apply(Expected_cervix$nr.ca_exp*RRRca_2v.lib2$cervix,2,sum,na.rm=TRUE)
  ca.prev_4v$cervix[i,] <- apply(Expected_cervix$nr.ca_exp*RRRca_4v$cervix,2,sum,na.rm=TRUE)
  ca.prev_9v$cervix[i,] <- apply(Expected_cervix$nr.ca_exp*RRRca_9v$cervix,2,sum,na.rm=TRUE)
  
  #total nr prevented, QALYs gained, savings treatment, savings pall.care
  Effects_2v.cons$cervix[i,] <- c(sum(ca.prev_2v.cons$cervix[i,]),
                                  sum(Expected_cervix$LYs.loss_exp*RRRca_2v.cons$cervix,na.rm=TRUE),
                                  sum(Expected_cervix$costs1_exp*RRRca_2v.cons$cervix,na.rm=TRUE),
                                  sum(Expected_cervix$costs2_exp*RRRca_2v.cons$cervix,na.rm=TRUE))
  Effects_2v.mid$cervix[i,] <- c(sum(ca.prev_2v.mid$cervix[i,]),
                                 sum(Expected_cervix$LYs.loss_exp*RRRca_2v.mid$cervix,na.rm=TRUE),
                                 sum(Expected_cervix$costs1_exp*RRRca_2v.mid$cervix,na.rm=TRUE),
                                 sum(Expected_cervix$costs2_exp*RRRca_2v.mid$cervix,na.rm=TRUE))
  Effects_2v.lib1$cervix[i,] <- c(sum(ca.prev_2v.lib1$cervix[i,]),
                                 sum(Expected_cervix$LYs.loss_exp*RRRca_2v.lib1$cervix,na.rm=TRUE),
                                 sum(Expected_cervix$costs1_exp*RRRca_2v.lib1$cervix,na.rm=TRUE),
                                 sum(Expected_cervix$costs2_exp*RRRca_2v.lib1$cervix,na.rm=TRUE))
  Effects_2v.lib2$cervix[i,] <- c(sum(ca.prev_2v.lib2$cervix[i,]),
                                 sum(Expected_cervix$LYs.loss_exp*RRRca_2v.lib2$cervix,na.rm=TRUE),
                                 sum(Expected_cervix$costs1_exp*RRRca_2v.lib2$cervix,na.rm=TRUE),
                                 sum(Expected_cervix$costs2_exp*RRRca_2v.lib2$cervix,na.rm=TRUE))
  Effects_4v$cervix[i,] <- c(sum(ca.prev_4v$cervix[i,]),
                             sum(Expected_cervix$LYs.loss_exp*RRRca_4v$cervix,na.rm=TRUE),
                             sum(Expected_cervix$costs1_exp*RRRca_4v$cervix,na.rm=TRUE),
                             sum(Expected_cervix$costs2_exp*RRRca_4v$cervix,na.rm=TRUE))
  Effects_9v$cervix[i,] <- c(sum(ca.prev_9v$cervix[i,]),
                             sum(Expected_cervix$LYs.loss_exp*RRRca_9v$cervix,na.rm=TRUE),
                             sum(Expected_cervix$costs1_exp*RRRca_9v$cervix,na.rm=TRUE),
                             sum(Expected_cervix$costs2_exp*RRRca_9v$cervix,na.rm=TRUE))
  
  #Anus _w
  ca.prev_2v.cons$anus_w[i,] <- apply(Expected_anus_w$nr.ca_exp*RRRca_2v.cons$anus_w,2,sum,na.rm=TRUE)
  ca.prev_2v.mid$anus_w[i,] <- apply(Expected_anus_w$nr.ca_exp*RRRca_2v.mid$anus_w,2,sum,na.rm=TRUE)
  ca.prev_2v.lib1$anus_w[i,] <- apply(Expected_anus_w$nr.ca_exp*RRRca_2v.lib1$anus_w,2,sum,na.rm=TRUE)
  ca.prev_2v.lib2$anus_w[i,] <- apply(Expected_anus_w$nr.ca_exp*RRRca_2v.lib2$anus_w,2,sum,na.rm=TRUE)
  ca.prev_4v$anus_w[i,] <- apply(Expected_anus_w$nr.ca_exp*RRRca_4v$anus_w,2,sum,na.rm=TRUE)
  ca.prev_9v$anus_w[i,] <- apply(Expected_anus_w$nr.ca_exp*RRRca_9v$anus_w,2,sum,na.rm=TRUE)
  
  #total nr prevented, QALYs gained, savings treatment, savings pall.care
  Effects_2v.cons$anus_w[i,] <- c(sum(ca.prev_2v.cons$anus_w[i,]),
                                  sum(Expected_anus_w$LYs.loss_exp*RRRca_2v.cons$anus_w,na.rm=TRUE),
                                  sum(Expected_anus_w$costs1_exp*RRRca_2v.cons$anus_w,na.rm=TRUE),
                                  sum(Expected_anus_w$costs2_exp*RRRca_2v.cons$anus_w,na.rm=TRUE))
  Effects_2v.mid$anus_w[i,] <- c(sum(ca.prev_2v.mid$anus_w[i,]),
                                 sum(Expected_anus_w$LYs.loss_exp*RRRca_2v.mid$anus_w,na.rm=TRUE),
                                 sum(Expected_anus_w$costs1_exp*RRRca_2v.mid$anus_w,na.rm=TRUE),
                                 sum(Expected_anus_w$costs2_exp*RRRca_2v.mid$anus_w,na.rm=TRUE))
  Effects_2v.lib1$anus_w[i,] <- c(sum(ca.prev_2v.lib1$anus_w[i,]),
                                 sum(Expected_anus_w$LYs.loss_exp*RRRca_2v.lib1$anus_w,na.rm=TRUE),
                                 sum(Expected_anus_w$costs1_exp*RRRca_2v.lib1$anus_w,na.rm=TRUE),
                                 sum(Expected_anus_w$costs2_exp*RRRca_2v.lib1$anus_w,na.rm=TRUE))
  Effects_2v.lib2$anus_w[i,] <- c(sum(ca.prev_2v.lib2$anus_w[i,]),
                                 sum(Expected_anus_w$LYs.loss_exp*RRRca_2v.lib2$anus_w,na.rm=TRUE),
                                 sum(Expected_anus_w$costs1_exp*RRRca_2v.lib2$anus_w,na.rm=TRUE),
                                 sum(Expected_anus_w$costs2_exp*RRRca_2v.lib2$anus_w,na.rm=TRUE))
  Effects_4v$anus_w[i,] <- c(sum(ca.prev_4v$anus_w[i,]),
                             sum(Expected_anus_w$LYs.loss_exp*RRRca_4v$anus_w,na.rm=TRUE),
                             sum(Expected_anus_w$costs1_exp*RRRca_4v$anus_w,na.rm=TRUE),
                             sum(Expected_anus_w$costs2_exp*RRRca_4v$anus_w,na.rm=TRUE))
  Effects_9v$anus_w[i,] <- c(sum(ca.prev_9v$anus_w[i,]),
                             sum(Expected_anus_w$LYs.loss_exp*RRRca_9v$anus_w,na.rm=TRUE),
                             sum(Expected_anus_w$costs1_exp*RRRca_9v$anus_w,na.rm=TRUE),
                             sum(Expected_anus_w$costs2_exp*RRRca_9v$anus_w,na.rm=TRUE))
  
  #Anus _m #All MSM
  ca.prev_2v.cons$anus_m[i,] <- apply(Expected_anus_m$nr.ca_exp*RRRca_2v.cons$anus_m,2,sum,na.rm=TRUE)
  ca.prev_2v.mid$anus_m[i,] <- apply(Expected_anus_m$nr.ca_exp*RRRca_2v.mid$anus_m,2,sum,na.rm=TRUE)
  ca.prev_2v.lib1$anus_m[i,] <- apply(Expected_anus_m$nr.ca_exp*RRRca_2v.lib1$anus_m,2,sum,na.rm=TRUE)
  ca.prev_2v.lib2$anus_m[i,] <- apply(Expected_anus_m$nr.ca_exp*RRRca_2v.lib2$anus_m,2,sum,na.rm=TRUE)
  ca.prev_4v$anus_m[i,] <- apply(Expected_anus_m$nr.ca_exp*RRRca_4v$anus_m,2,sum,na.rm=TRUE)
  ca.prev_9v$anus_m[i,] <- apply(Expected_anus_m$nr.ca_exp*RRRca_9v$anus_m,2,sum,na.rm=TRUE)
  
  #total nr prevented, QALYs gained, savings treatment, savings pall.care
  Effects_2v.cons$anus_m[i,] <- c(sum(ca.prev_2v.cons$anus_m[i,]),
                                  sum(Expected_anus_m$LYs.loss_exp*RRRca_2v.cons$anus_m,na.rm=TRUE),
                                  sum(Expected_anus_m$costs1_exp*RRRca_2v.cons$anus_m,na.rm=TRUE),
                                  sum(Expected_anus_m$costs2_exp*RRRca_2v.cons$anus_m,na.rm=TRUE))
  Effects_2v.mid$anus_m[i,] <- c(sum(ca.prev_2v.mid$anus_m[i,]),
                                 sum(Expected_anus_m$LYs.loss_exp*RRRca_2v.mid$anus_m,na.rm=TRUE),
                                 sum(Expected_anus_m$costs1_exp*RRRca_2v.mid$anus_m,na.rm=TRUE),
                                 sum(Expected_anus_m$costs2_exp*RRRca_2v.mid$anus_m,na.rm=TRUE))
  Effects_2v.lib1$anus_m[i,] <- c(sum(ca.prev_2v.lib1$anus_m[i,]),
                                 sum(Expected_anus_m$LYs.loss_exp*RRRca_2v.lib1$anus_m,na.rm=TRUE),
                                 sum(Expected_anus_m$costs1_exp*RRRca_2v.lib1$anus_m,na.rm=TRUE),
                                 sum(Expected_anus_m$costs2_exp*RRRca_2v.lib1$anus_m,na.rm=TRUE))
  Effects_2v.lib2$anus_m[i,] <- c(sum(ca.prev_2v.lib2$anus_m[i,]),
                                 sum(Expected_anus_m$LYs.loss_exp*RRRca_2v.lib2$anus_m,na.rm=TRUE),
                                 sum(Expected_anus_m$costs1_exp*RRRca_2v.lib2$anus_m,na.rm=TRUE),
                                 sum(Expected_anus_m$costs2_exp*RRRca_2v.lib2$anus_m,na.rm=TRUE))
  Effects_4v$anus_m[i,] <- c(sum(ca.prev_4v$anus_m[i,]),
                             sum(Expected_anus_m$LYs.loss_exp*RRRca_4v$anus_m,na.rm=TRUE),
                             sum(Expected_anus_m$costs1_exp*RRRca_4v$anus_m,na.rm=TRUE),
                             sum(Expected_anus_m$costs2_exp*RRRca_4v$anus_m,na.rm=TRUE))
  Effects_9v$anus_m[i,] <- c(sum(ca.prev_9v$anus_m[i,]),
                             sum(Expected_anus_m$LYs.loss_exp*RRRca_9v$anus_m,na.rm=TRUE),
                             sum(Expected_anus_m$costs1_exp*RRRca_9v$anus_m,na.rm=TRUE),
                             sum(Expected_anus_m$costs2_exp*RRRca_9v$anus_m,na.rm=TRUE))
  
  #Oroph _w
  ca.prev_2v.cons$oroph_w[i,] <- apply(Expected_oroph_w$nr.ca_exp*RRRca_2v.cons$oroph_w,2,sum,na.rm=TRUE)
  ca.prev_2v.mid$oroph_w[i,] <- apply(Expected_oroph_w$nr.ca_exp*RRRca_2v.mid$oroph_w,2,sum,na.rm=TRUE)
  ca.prev_2v.lib1$oroph_w[i,] <- apply(Expected_oroph_w$nr.ca_exp*RRRca_2v.lib1$oroph_w,2,sum,na.rm=TRUE)
  ca.prev_2v.lib2$oroph_w[i,] <- apply(Expected_oroph_w$nr.ca_exp*RRRca_2v.lib2$oroph_w,2,sum,na.rm=TRUE)
  ca.prev_4v$oroph_w[i,] <- apply(Expected_oroph_w$nr.ca_exp*RRRca_4v$oroph_w,2,sum,na.rm=TRUE)
  ca.prev_9v$oroph_w[i,] <- apply(Expected_oroph_w$nr.ca_exp*RRRca_9v$oroph_w,2,sum,na.rm=TRUE)
  
  #total nr prevented, QALYs gained, savings treatment, savings pall.care
  Effects_2v.cons$oroph_w[i,] <- c(sum(ca.prev_2v.cons$oroph_w[i,]),
                                   sum(Expected_oroph_w$LYs.loss_exp*RRRca_2v.cons$oroph_w,na.rm=TRUE),
                                   sum(Expected_oroph_w$costs1_exp*RRRca_2v.cons$oroph_w,na.rm=TRUE),
                                   sum(Expected_oroph_w$costs2_exp*RRRca_2v.cons$oroph_w,na.rm=TRUE))
  Effects_2v.mid$oroph_w[i,] <- c(sum(ca.prev_2v.mid$oroph_w[i,]),
                                  sum(Expected_oroph_w$LYs.loss_exp*RRRca_2v.mid$oroph_w,na.rm=TRUE),
                                  sum(Expected_oroph_w$costs1_exp*RRRca_2v.mid$oroph_w,na.rm=TRUE),
                                  sum(Expected_oroph_w$costs2_exp*RRRca_2v.mid$oroph_w,na.rm=TRUE))
  Effects_2v.lib1$oroph_w[i,] <- c(sum(ca.prev_2v.lib1$oroph_w[i,]),
                                  sum(Expected_oroph_w$LYs.loss_exp*RRRca_2v.lib1$oroph_w,na.rm=TRUE),
                                  sum(Expected_oroph_w$costs1_exp*RRRca_2v.lib1$oroph_w,na.rm=TRUE),
                                  sum(Expected_oroph_w$costs2_exp*RRRca_2v.lib1$oroph_w,na.rm=TRUE))
  Effects_2v.lib2$oroph_w[i,] <- c(sum(ca.prev_2v.lib2$oroph_w[i,]),
                                  sum(Expected_oroph_w$LYs.loss_exp*RRRca_2v.lib2$oroph_w,na.rm=TRUE),
                                  sum(Expected_oroph_w$costs1_exp*RRRca_2v.lib2$oroph_w,na.rm=TRUE),
                                  sum(Expected_oroph_w$costs2_exp*RRRca_2v.lib2$oroph_w,na.rm=TRUE))
  Effects_4v$oroph_w[i,] <- c(sum(ca.prev_4v$oroph_w[i,]),
                              sum(Expected_oroph_w$LYs.loss_exp*RRRca_4v$oroph_w,na.rm=TRUE),
                              sum(Expected_oroph_w$costs1_exp*RRRca_4v$oroph_w,na.rm=TRUE),
                              sum(Expected_oroph_w$costs2_exp*RRRca_4v$oroph_w,na.rm=TRUE))
  Effects_9v$oroph_w[i,] <- c(sum(ca.prev_9v$oroph_w[i,]),
                              sum(Expected_oroph_w$LYs.loss_exp*RRRca_9v$oroph_w,na.rm=TRUE),
                              sum(Expected_oroph_w$costs1_exp*RRRca_9v$oroph_w,na.rm=TRUE),
                              sum(Expected_oroph_w$costs2_exp*RRRca_9v$oroph_w,na.rm=TRUE))
  
  #Oroph _m
  ca.prev_2v.cons$oroph_m[i,] <- apply(Expected_oroph_m$nr.ca_exp*RRRca_2v.cons$oroph_m,2,sum,na.rm=TRUE)
  ca.prev_2v.mid$oroph_m[i,] <- apply(Expected_oroph_m$nr.ca_exp*RRRca_2v.mid$oroph_m,2,sum,na.rm=TRUE)
  ca.prev_2v.lib1$oroph_m[i,] <- apply(Expected_oroph_m$nr.ca_exp*RRRca_2v.lib1$oroph_m,2,sum,na.rm=TRUE)
  ca.prev_2v.lib2$oroph_m[i,] <- apply(Expected_oroph_m$nr.ca_exp*RRRca_2v.lib2$oroph_m,2,sum,na.rm=TRUE)
  ca.prev_4v$oroph_m[i,] <- apply(Expected_oroph_m$nr.ca_exp*RRRca_4v$oroph_m,2,sum,na.rm=TRUE)
  ca.prev_9v$oroph_m[i,] <- apply(Expected_oroph_m$nr.ca_exp*RRRca_9v$oroph_m,2,sum,na.rm=TRUE)
  
  #total nr prevented, QALYs gained, savings treatment, savings pall.care
  Effects_2v.cons$oroph_m[i,] <- c(sum(ca.prev_2v.cons$oroph_m[i,]),
                                   sum(Expected_oroph_m$LYs.loss_exp*RRRca_2v.cons$oroph_m,na.rm=TRUE),
                                   sum(Expected_oroph_m$costs1_exp*RRRca_2v.cons$oroph_m,na.rm=TRUE),
                                   sum(Expected_oroph_m$costs2_exp*RRRca_2v.cons$oroph_m,na.rm=TRUE))
  Effects_2v.mid$oroph_m[i,] <- c(sum(ca.prev_2v.mid$oroph_m[i,]),
                                  sum(Expected_oroph_m$LYs.loss_exp*RRRca_2v.mid$oroph_m,na.rm=TRUE),
                                  sum(Expected_oroph_m$costs1_exp*RRRca_2v.mid$oroph_m,na.rm=TRUE),
                                  sum(Expected_oroph_m$costs2_exp*RRRca_2v.mid$oroph_m,na.rm=TRUE))
  Effects_2v.lib1$oroph_m[i,] <- c(sum(ca.prev_2v.lib1$oroph_m[i,]),
                                  sum(Expected_oroph_m$LYs.loss_exp*RRRca_2v.lib1$oroph_m,na.rm=TRUE),
                                  sum(Expected_oroph_m$costs1_exp*RRRca_2v.lib1$oroph_m,na.rm=TRUE),
                                  sum(Expected_oroph_m$costs2_exp*RRRca_2v.lib1$oroph_m,na.rm=TRUE))
  Effects_2v.lib2$oroph_m[i,] <- c(sum(ca.prev_2v.lib2$oroph_m[i,]),
                                  sum(Expected_oroph_m$LYs.loss_exp*RRRca_2v.lib2$oroph_m,na.rm=TRUE),
                                  sum(Expected_oroph_m$costs1_exp*RRRca_2v.lib2$oroph_m,na.rm=TRUE),
                                  sum(Expected_oroph_m$costs2_exp*RRRca_2v.lib2$oroph_m,na.rm=TRUE))
  Effects_4v$oroph_m[i,] <- c(sum(ca.prev_4v$oroph_m[i,]),
                              sum(Expected_oroph_m$LYs.loss_exp*RRRca_4v$oroph_m,na.rm=TRUE),
                              sum(Expected_oroph_m$costs1_exp*RRRca_4v$oroph_m,na.rm=TRUE),
                              sum(Expected_oroph_m$costs2_exp*RRRca_4v$oroph_m,na.rm=TRUE))
  Effects_9v$oroph_m[i,] <- c(sum(ca.prev_9v$oroph_m[i,]),
                              sum(Expected_oroph_m$LYs.loss_exp*RRRca_9v$oroph_m,na.rm=TRUE),
                              sum(Expected_oroph_m$costs1_exp*RRRca_9v$oroph_m,na.rm=TRUE),
                              sum(Expected_oroph_m$costs2_exp*RRRca_9v$oroph_m,na.rm=TRUE))
  
  #Vulva
  ca.prev_2v.cons$vulva[i,] <- apply(Expected_vulva$nr.ca_exp*RRRca_2v.cons$vulva,2,sum,na.rm=TRUE)
  ca.prev_2v.mid$vulva[i,] <- apply(Expected_vulva$nr.ca_exp*RRRca_2v.mid$vulva,2,sum,na.rm=TRUE)
  ca.prev_2v.lib1$vulva[i,] <- apply(Expected_vulva$nr.ca_exp*RRRca_2v.lib1$vulva,2,sum,na.rm=TRUE)
  ca.prev_2v.lib2$vulva[i,] <- apply(Expected_vulva$nr.ca_exp*RRRca_2v.lib2$vulva,2,sum,na.rm=TRUE)
  ca.prev_4v$vulva[i,] <- apply(Expected_vulva$nr.ca_exp*RRRca_4v$vulva,2,sum,na.rm=TRUE)
  ca.prev_9v$vulva[i,] <- apply(Expected_vulva$nr.ca_exp*RRRca_9v$vulva,2,sum,na.rm=TRUE)
  
  #total nr prevented, QALYs gained, savings treatment, savings pall.care
  Effects_2v.cons$vulva[i,] <- c(sum(ca.prev_2v.cons$vulva[i,]),
                                 sum(Expected_vulva$LYs.loss_exp*RRRca_2v.cons$vulva,na.rm=TRUE),
                                 sum(Expected_vulva$costs1_exp*RRRca_2v.cons$vulva,na.rm=TRUE),
                                 sum(Expected_vulva$costs2_exp*RRRca_2v.cons$vulva,na.rm=TRUE))
  Effects_2v.mid$vulva[i,] <- c(sum(ca.prev_2v.mid$vulva[i,]),
                                sum(Expected_vulva$LYs.loss_exp*RRRca_2v.mid$vulva,na.rm=TRUE),
                                sum(Expected_vulva$costs1_exp*RRRca_2v.mid$vulva,na.rm=TRUE),
                                sum(Expected_vulva$costs2_exp*RRRca_2v.mid$vulva,na.rm=TRUE))
  Effects_2v.lib1$vulva[i,] <- c(sum(ca.prev_2v.lib1$vulva[i,]),
                                sum(Expected_vulva$LYs.loss_exp*RRRca_2v.lib1$vulva,na.rm=TRUE),
                                sum(Expected_vulva$costs1_exp*RRRca_2v.lib1$vulva,na.rm=TRUE),
                                sum(Expected_vulva$costs2_exp*RRRca_2v.lib1$vulva,na.rm=TRUE))
  Effects_2v.lib2$vulva[i,] <- c(sum(ca.prev_2v.lib2$vulva[i,]),
                                sum(Expected_vulva$LYs.loss_exp*RRRca_2v.lib2$vulva,na.rm=TRUE),
                                sum(Expected_vulva$costs1_exp*RRRca_2v.lib2$vulva,na.rm=TRUE),
                                sum(Expected_vulva$costs2_exp*RRRca_2v.lib2$vulva,na.rm=TRUE))
  Effects_4v$vulva[i,] <- c(sum(ca.prev_4v$vulva[i,]),
                            sum(Expected_vulva$LYs.loss_exp*RRRca_4v$vulva,na.rm=TRUE),
                            sum(Expected_vulva$costs1_exp*RRRca_4v$vulva,na.rm=TRUE),
                            sum(Expected_vulva$costs2_exp*RRRca_4v$vulva,na.rm=TRUE))
  Effects_9v$vulva[i,] <- c(sum(ca.prev_9v$vulva[i,]),
                            sum(Expected_vulva$LYs.loss_exp*RRRca_9v$vulva,na.rm=TRUE),
                            sum(Expected_vulva$costs1_exp*RRRca_9v$vulva,na.rm=TRUE),
                            sum(Expected_vulva$costs2_exp*RRRca_9v$vulva,na.rm=TRUE))
  
  #Vagina
  ca.prev_2v.cons$vagina[i,] <- apply(Expected_vagina$nr.ca_exp*RRRca_2v.cons$vagina,2,sum,na.rm=TRUE)
  ca.prev_2v.mid$vagina[i,] <- apply(Expected_vagina$nr.ca_exp*RRRca_2v.mid$vagina,2,sum,na.rm=TRUE)
  ca.prev_2v.lib1$vagina[i,] <- apply(Expected_vagina$nr.ca_exp*RRRca_2v.lib1$vagina,2,sum,na.rm=TRUE)
  ca.prev_2v.lib2$vagina[i,] <- apply(Expected_vagina$nr.ca_exp*RRRca_2v.lib2$vagina,2,sum,na.rm=TRUE)
  ca.prev_4v$vagina[i,] <- apply(Expected_vagina$nr.ca_exp*RRRca_4v$vagina,2,sum,na.rm=TRUE)
  ca.prev_9v$vagina[i,] <- apply(Expected_vagina$nr.ca_exp*RRRca_9v$vagina,2,sum,na.rm=TRUE)
  
  #total nr prevented, QALYs gained, savings treatment, savings pall.care
  Effects_2v.cons$vagina[i,] <- c(sum(ca.prev_2v.cons$vagina[i,]),
                                  sum(Expected_vagina$LYs.loss_exp*RRRca_2v.cons$vagina,na.rm=TRUE),
                                  sum(Expected_vagina$costs1_exp*RRRca_2v.cons$vagina,na.rm=TRUE),
                                  sum(Expected_vagina$costs2_exp*RRRca_2v.cons$vagina,na.rm=TRUE))
  Effects_2v.mid$vagina[i,] <- c(sum(ca.prev_2v.mid$vagina[i,]),
                                 sum(Expected_vagina$LYs.loss_exp*RRRca_2v.mid$vagina,na.rm=TRUE),
                                 sum(Expected_vagina$costs1_exp*RRRca_2v.mid$vagina,na.rm=TRUE),
                                 sum(Expected_vagina$costs2_exp*RRRca_2v.mid$vagina,na.rm=TRUE))
  Effects_2v.lib1$vagina[i,] <- c(sum(ca.prev_2v.lib1$vagina[i,]),
                                 sum(Expected_vagina$LYs.loss_exp*RRRca_2v.lib1$vagina,na.rm=TRUE),
                                 sum(Expected_vagina$costs1_exp*RRRca_2v.lib1$vagina,na.rm=TRUE),
                                 sum(Expected_vagina$costs2_exp*RRRca_2v.lib1$vagina,na.rm=TRUE))
  Effects_2v.lib2$vagina[i,] <- c(sum(ca.prev_2v.lib2$vagina[i,]),
                                 sum(Expected_vagina$LYs.loss_exp*RRRca_2v.lib2$vagina,na.rm=TRUE),
                                 sum(Expected_vagina$costs1_exp*RRRca_2v.lib2$vagina,na.rm=TRUE),
                                 sum(Expected_vagina$costs2_exp*RRRca_2v.lib2$vagina,na.rm=TRUE))
  Effects_4v$vagina[i,] <- c(sum(ca.prev_4v$vagina[i,]),
                             sum(Expected_vagina$LYs.loss_exp*RRRca_4v$vagina,na.rm=TRUE),
                             sum(Expected_vagina$costs1_exp*RRRca_4v$vagina,na.rm=TRUE),
                             sum(Expected_vagina$costs2_exp*RRRca_4v$vagina,na.rm=TRUE))
  Effects_9v$vagina[i,] <- c(sum(ca.prev_9v$vagina[i,]),
                             sum(Expected_vagina$LYs.loss_exp*RRRca_9v$vagina,na.rm=TRUE),
                             sum(Expected_vagina$costs1_exp*RRRca_9v$vagina,na.rm=TRUE),
                             sum(Expected_vagina$costs2_exp*RRRca_9v$vagina,na.rm=TRUE))
  
  #Penis
  ca.prev_2v.cons$penis[i,] <- apply(Expected_penis$nr.ca_exp*RRRca_2v.cons$penis,2,sum,na.rm=TRUE)
  ca.prev_2v.mid$penis[i,] <- apply(Expected_penis$nr.ca_exp*RRRca_2v.mid$penis,2,sum,na.rm=TRUE)
  ca.prev_2v.lib1$penis[i,] <- apply(Expected_penis$nr.ca_exp*RRRca_2v.lib1$penis,2,sum,na.rm=TRUE)
  ca.prev_2v.lib2$penis[i,] <- apply(Expected_penis$nr.ca_exp*RRRca_2v.lib2$penis,2,sum,na.rm=TRUE)
  ca.prev_4v$penis[i,] <- apply(Expected_penis$nr.ca_exp*RRRca_4v$penis,2,sum,na.rm=TRUE)
  ca.prev_9v$penis[i,] <- apply(Expected_penis$nr.ca_exp*RRRca_9v$penis,2,sum,na.rm=TRUE)
  
  #total nr prevented, QALYs gained, savings treatment, savings pall.care
  Effects_2v.cons$penis[i,] <- c(sum(ca.prev_2v.cons$penis[i,]),
                                 sum(Expected_penis$LYs.loss_exp*RRRca_2v.cons$penis,na.rm=TRUE),
                                 sum(Expected_penis$costs1_exp*RRRca_2v.cons$penis,na.rm=TRUE),
                                 sum(Expected_penis$costs2_exp*RRRca_2v.cons$penis,na.rm=TRUE))
  Effects_2v.mid$penis[i,] <- c(sum(ca.prev_2v.mid$penis[i,]),
                                sum(Expected_penis$LYs.loss_exp*RRRca_2v.mid$penis,na.rm=TRUE),
                                sum(Expected_penis$costs1_exp*RRRca_2v.mid$penis,na.rm=TRUE),
                                sum(Expected_penis$costs2_exp*RRRca_2v.mid$penis,na.rm=TRUE))
  Effects_2v.lib1$penis[i,] <- c(sum(ca.prev_2v.lib1$penis[i,]),
                                sum(Expected_penis$LYs.loss_exp*RRRca_2v.lib1$penis,na.rm=TRUE),
                                sum(Expected_penis$costs1_exp*RRRca_2v.lib1$penis,na.rm=TRUE),
                                sum(Expected_penis$costs2_exp*RRRca_2v.lib1$penis,na.rm=TRUE))
  Effects_2v.lib2$penis[i,] <- c(sum(ca.prev_2v.lib2$penis[i,]),
                                sum(Expected_penis$LYs.loss_exp*RRRca_2v.lib2$penis,na.rm=TRUE),
                                sum(Expected_penis$costs1_exp*RRRca_2v.lib2$penis,na.rm=TRUE),
                                sum(Expected_penis$costs2_exp*RRRca_2v.lib2$penis,na.rm=TRUE))
  Effects_4v$penis[i,] <- c(sum(ca.prev_4v$penis[i,]),
                            sum(Expected_penis$LYs.loss_exp*RRRca_4v$penis,na.rm=TRUE),
                            sum(Expected_penis$costs1_exp*RRRca_4v$penis,na.rm=TRUE),
                            sum(Expected_penis$costs2_exp*RRRca_4v$penis,na.rm=TRUE))
  Effects_9v$penis[i,] <- c(sum(ca.prev_9v$penis[i,]),
                            sum(Expected_penis$LYs.loss_exp*RRRca_9v$penis,na.rm=TRUE),
                            sum(Expected_penis$costs1_exp*RRRca_9v$penis,na.rm=TRUE),
                            sum(Expected_penis$costs2_exp*RRRca_9v$penis,na.rm=TRUE))

}

#ICER computation
TotalEffects_2v.cons <- Reduce("+",Effects_2v.cons)
TotalEffects_2v.mid <- Reduce("+",Effects_2v.mid)
TotalEffects_2v.lib1 <- Reduce("+",Effects_2v.lib1)
TotalEffects_2v.lib2 <- Reduce("+",Effects_2v.lib2)
TotalEffects_4v <- Reduce("+",Effects_4v)
TotalEffects_9v <- Reduce("+",Effects_9v)
colnames(TotalEffects_2v.cons) <- c("nr.prev","QALYs","Savings1","Savings2")
colnames(TotalEffects_2v.mid) <- c("nr.prev","QALYs","Savings1","Savings2")
colnames(TotalEffects_2v.lib1) <- c("nr.prev","QALYs","Savings1","Savings2")
colnames(TotalEffects_2v.lib2) <- c("nr.prev","QALYs","Savings1","Savings2")
colnames(TotalEffects_4v) <- c("nr.prev","QALYs","Savings1","Savings2")
colnames(TotalEffects_9v) <- c("nr.prev","QALYs","Savings1","Savings2")

CostDiff_2v.cons <- sum(cohort*vacc.up*vacc.price_2v) - TotalEffects_2v.cons[,"Savings1"] - TotalEffects_2v.cons[,"Savings2"]
ICER_2v.cons <- CostDiff_2v.cons/TotalEffects_2v.cons[,"QALYs"]
summary(ICER_2v.cons)

CostDiff_2v.mid <- sum(cohort*vacc.up*vacc.price_2v) - TotalEffects_2v.mid[,"Savings1"] - TotalEffects_2v.mid[,"Savings2"]
ICER_2v.mid <- CostDiff_2v.mid/TotalEffects_2v.mid[,"QALYs"]
summary(ICER_2v.mid)

CostDiff_2v.lib1 <- sum(cohort*vacc.up*vacc.price_2v) - TotalEffects_2v.lib1[,"Savings1"] - TotalEffects_2v.lib1[,"Savings2"]
ICER_2v.lib1 <- CostDiff_2v.lib1/TotalEffects_2v.lib1[,"QALYs"]
summary(ICER_2v.lib1)

CostDiff_2v.lib2 <- sum(cohort*vacc.up*vacc.price_2v) - TotalEffects_2v.lib2[,"Savings1"] - TotalEffects_2v.lib2[,"Savings2"]
ICER_2v.lib2 <- CostDiff_2v.lib2/TotalEffects_2v.lib2[,"QALYs"]
summary(ICER_2v.lib2)

CostDiff_4v <- sum(cohort*vacc.up*vacc.price_4v) - TotalEffects_4v[,"Savings1"] - TotalEffects_4v[,"Savings2"]
ICER_4v <- CostDiff_4v/TotalEffects_4v[,"QALYs"]
summary(ICER_4v) 

CostDiff_9v <- sum(cohort*vacc.up*vacc.price_9v) - TotalEffects_9v[,"Savings1"] - TotalEffects_9v[,"Savings2"]
ICER_9v <- CostDiff_9v/TotalEffects_9v[,"QALYs"]
summary(ICER_9v)  

#Incremental
ICER_9v_vs_2v.cons <- (CostDiff_9v - CostDiff_2v.cons)/(TotalEffects_9v[,"QALYs"] - 
                                                          TotalEffects_2v.cons[,"QALYs"])
ICER_9v_vs_2v.mid <- (CostDiff_9v - CostDiff_2v.mid)/(TotalEffects_9v[,"QALYs"] - 
                                                        TotalEffects_2v.mid[,"QALYs"])
ICER_9v_vs_2v.lib1 <- (CostDiff_9v - CostDiff_2v.lib1)/(TotalEffects_9v[,"QALYs"] - 
                                                        TotalEffects_2v.lib1[,"QALYs"])
ICER_9v_vs_2v.lib2 <- (CostDiff_9v - CostDiff_2v.lib2)/(TotalEffects_9v[,"QALYs"] - 
                                                        TotalEffects_2v.lib2[,"QALYs"])
ICER_4v_vs_2v.cons <- (CostDiff_4v - CostDiff_2v.cons)/(TotalEffects_4v[,"QALYs"] - 
                                                          TotalEffects_2v.cons[,"QALYs"])
ICER_4v_vs_2v.mid <- (CostDiff_4v - CostDiff_2v.mid)/(TotalEffects_4v[,"QALYs"] - 
                                                        TotalEffects_2v.mid[,"QALYs"])
ICER_4v_vs_2v.lib1 <- (CostDiff_4v - CostDiff_2v.lib1)/(TotalEffects_4v[,"QALYs"] - 
                                                        TotalEffects_2v.lib1[,"QALYs"])
ICER_4v_vs_2v.lib2 <- (CostDiff_4v - CostDiff_2v.lib2)/(TotalEffects_4v[,"QALYs"] - 
                                                        TotalEffects_2v.lib2[,"QALYs"])
ICER_9v_vs_4v <- (CostDiff_9v - CostDiff_4v)/(TotalEffects_9v[,"QALYs"] - 
                                                TotalEffects_4v[,"QALYs"])

#Net Monetary Benefit
WTPthreshold <- 20000

NMB_9v_vs_2v.cons <- (TotalEffects_9v[,"QALYs"] - TotalEffects_2v.cons[,"QALYs"])*WTPthreshold -
  (CostDiff_9v - CostDiff_2v.cons)
NMB_9v_vs_2v.mid <- (TotalEffects_9v[,"QALYs"] - TotalEffects_2v.mid[,"QALYs"])*WTPthreshold -
  (CostDiff_9v - CostDiff_2v.mid)
NMB_9v_vs_2v.lib1 <- (TotalEffects_9v[,"QALYs"] - TotalEffects_2v.lib1[,"QALYs"])*WTPthreshold -
  (CostDiff_9v - CostDiff_2v.lib1)
NMB_9v_vs_2v.lib2 <- (TotalEffects_9v[,"QALYs"] - TotalEffects_2v.lib2[,"QALYs"])*WTPthreshold -
  (CostDiff_9v - CostDiff_2v.lib2)
NMB_4v_vs_2v.cons <- (TotalEffects_4v[,"QALYs"] - TotalEffects_2v.cons[,"QALYs"])*WTPthreshold -
  (CostDiff_4v - CostDiff_2v.cons)
NMB_4v_vs_2v.mid <- (TotalEffects_4v[,"QALYs"] - TotalEffects_2v.mid[,"QALYs"])*WTPthreshold -
  (CostDiff_4v - CostDiff_2v.mid)
NMB_4v_vs_2v.lib1 <- (TotalEffects_4v[,"QALYs"] - TotalEffects_2v.lib1[,"QALYs"])*WTPthreshold -
  (CostDiff_4v - CostDiff_2v.lib1)
NMB_4v_vs_2v.lib2 <- (TotalEffects_4v[,"QALYs"] - TotalEffects_2v.lib2[,"QALYs"])*WTPthreshold -
  (CostDiff_4v - CostDiff_2v.lib2)
NMB_9v_vs_4v <- (TotalEffects_9v[,"QALYs"] - TotalEffects_4v[,"QALYs"])*WTPthreshold -
  (CostDiff_9v - CostDiff_4v)


#save results in list
Effects <- list("2v.cons" = TotalEffects_2v.cons[,"QALYs"], "2v.mid" = TotalEffects_2v.mid[,"QALYs"],
                "2v.lib1" = TotalEffects_2v.lib1[,"QALYs"],"2v.lib2" = TotalEffects_2v.lib2[,"QALYs"], 
                "4v" = TotalEffects_4v[,"QALYs"],"9v" = TotalEffects_9v[,"QALYs"])
Costs <- list("2v.cons" = CostDiff_2v.cons,"2v.mid" = CostDiff_2v.mid,
              "2v.lib1" = CostDiff_2v.lib1,"2v.lib2" = CostDiff_2v.lib2,
              "4v" = CostDiff_4v,"9v" = CostDiff_9v)
ICERs <- list("2v.cons" = ICER_2v.cons, "2v.mid" = ICER_2v.mid,
              "2v.lib1" = ICER_2v.lib1, "2v.lib2" = ICER_2v.lib2,
              "4v" = ICER_4v, "9v" = ICER_9v,
              "9v_vs_2v.cons" = ICER_9v_vs_2v.cons, "9v_vs_2v.mid" = ICER_9v_vs_2v.mid,
              "9v_vs_2v.lib1" = ICER_9v_vs_2v.lib1,"9v_vs_2v.lib2" = ICER_9v_vs_2v.lib2,
              "4v_vs_2v.cons" = ICER_4v_vs_2v.cons,"4v_vs_2v.mid" = ICER_4v_vs_2v.mid,
              "4v_vs_2v.lib1" = ICER_4v_vs_2v.lib1,"4v_vs_2v.lib2" = ICER_4v_vs_2v.lib2,
              "9v_vs_4v" = ICER_9v_vs_4v)
NMB <- list("9v_vs_2v.cons" = NMB_9v_vs_2v.cons, "9v_vs_2v.mid" = NMB_9v_vs_2v.mid,
            "9v_vs_2v.lib1" = NMB_9v_vs_2v.lib1,"9v_vs_2v.lib2" = NMB_9v_vs_2v.lib2,
            "4v_vs_2v.cons" = NMB_4v_vs_2v.cons,"4v_vs_2v.mid" = NMB_4v_vs_2v.mid,
            "4v_vs_2v.lib1" = NMB_4v_vs_2v.lib1,"4v_vs_2v.lib2" = NMB_4v_vs_2v.lib2,
            "9v_vs_4v" = NMB_9v_vs_4v)
Results <- list(Effects_2v.cons=Effects_2v.cons,Effects_2v.mid=Effects_2v.mid,
                Effects_2v.lib1=Effects_2v.lib1,Effects_2v.lib2=Effects_2v.lib2,
                Effects_4v=Effects_4v,Effects_9v=Effects_9v,Expected.total=Expected.total,
                TotalEffects=Effects,Costs=Costs,ICERs=ICERs,NMB=NMB)
