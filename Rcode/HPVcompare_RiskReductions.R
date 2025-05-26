#Relative risk reductions (type by age specific) per vaccine and gender
RRRs <- RiskReductions(scenario)

RRR_2v.cons_w <- RRRs$RRR_2vCons$w
RRR_2v.cons_m <- RRRs$RRR_2vCons$m
RRR_2v.cons_MSM <- RRRs$RRR_2vCons$MSM
RRR_2v.mid_w <- RRRs$RRR_2vMid$w
RRR_2v.mid_m <- RRRs$RRR_2vMid$m
RRR_2v.mid_MSM <- RRRs$RRR_2vMid$MSM
RRR_2v.lib1_w <- RRRs$RRR_2vLib1$w
RRR_2v.lib1_m <- RRRs$RRR_2vLib1$m
RRR_2v.lib1_MSM <- RRRs$RRR_2vLib1$MSM
RRR_2v.lib2_w <- RRRs$RRR_2vLib2$w
RRR_2v.lib2_m <- RRRs$RRR_2vLib2$m
RRR_2v.lib2_MSM <- RRRs$RRR_2vLib2$MSM
RRR_4v_w <- RRRs$RRR_4v$w
RRR_4v_m <- RRRs$RRR_4v$m
RRR_4v_MSM <- RRRs$RRR_4v$MSM
RRR_9v_w <- RRRs$RRR_9v$w
RRR_9v_m <- RRRs$RRR_9v$m
RRR_9v_MSM <- RRRs$RRR_9v$MSM

RRR_611_w <- RRRs$RRR_611$w
RRR_611_m <- RRRs$RRR_611$m

## Risk reductions on screening
RRscreening <- RiskReductionScreening(scenario,agegrps.scr)

RRscr_2v.cons <- RRscreening$RR_2v.cons
RRscr_2v.mid <- RRscreening$RR_2v.mid
RRscr_2v.lib1 <- RRscreening$RR_2v.lib1
RRscr_2v.lib2 <- RRscreening$RR_2v.lib2
RRscr_4v <- RRscreening$RR_4v
RRscr_9v <- RRscreening$RR_9v

## Risk reductions on cancer
#Per vaccin and cancer type
ca.names <- c("cervix","anus_w","anus_m","oroph_w","oroph_m","vulva","vagina","penis")
RRRca_2v.cons <- lapply(1:8, matrix, data= 0, nrow=types, ncol=15)
RRRca_2v.mid <- lapply(1:8, matrix, data= 0, nrow=types, ncol=15)
RRRca_2v.lib1 <- lapply(1:8, matrix, data= 0, nrow=types, ncol=15)
RRRca_2v.lib2 <- lapply(1:8, matrix, data= 0, nrow=types, ncol=15)
RRRca_4v <- lapply(1:8, matrix, data= 0, nrow=types, ncol=15)
RRRca_9v <- lapply(1:8, matrix, data= 0, nrow=types, ncol=15)
names(RRRca_2v.cons) <- ca.names
names(RRRca_2v.mid) <- ca.names
names(RRRca_2v.lib1) <- ca.names
names(RRRca_2v.lib2) <- ca.names
names(RRRca_4v) <- ca.names
names(RRRca_9v) <- ca.names

#Part of male cancers that are MSM (which is not the same as PAF MSM)
prop.MSM_anus <- 1 #we assume all anal cancers are MSM
prop.MSM_oroph <- p.MSM*RR.MSM_oroph/(p.MSM*(RR.MSM_oroph - 1) +1)
prop.MSM_penis <- p.MSM*RR.MSM_penis/(p.MSM*(RR.MSM_penis - 1) +1)
prop.MSM <- c(0,0,prop.MSM_anus,0,prop.MSM_oroph,0,0,prop.MSM_penis)

#cancer for women
for (j in c(1,2,4,6,7)){
  RRRca_2v.cons[[j]] <- RiskReductionCancer(RRR_2v.cons_w,gamma.pars[j,],HPV.inc_w)
  RRRca_2v.mid[[j]] <- RiskReductionCancer(RRR_2v.mid_w,gamma.pars[j,],HPV.inc_w)
  RRRca_2v.lib1[[j]] <- RiskReductionCancer(RRR_2v.lib1_w,gamma.pars[j,],HPV.inc_w)
  RRRca_2v.lib2[[j]] <- RiskReductionCancer(RRR_2v.lib2_w,gamma.pars[j,],HPV.inc_w)
  RRRca_4v[[j]] <- RiskReductionCancer(RRR_4v_w,gamma.pars[j,],HPV.inc_w)
  RRRca_9v[[j]] <- RiskReductionCancer(RRR_9v_w,gamma.pars[j,],HPV.inc_w)
}
#cancer for men
for (j in c(3,5,8)){
  RRRca_2v.cons[[j]] <- prop.MSM[j]*RiskReductionCancer(RRR_2v.cons_MSM,gamma.pars[j,],HPV.inc_MSM) +
    (1-prop.MSM[j])*RiskReductionCancer(RRR_2v.cons_m,gamma.pars[j,],HPV.inc_m)
  RRRca_2v.mid[[j]] <- prop.MSM[j]*RiskReductionCancer(RRR_2v.mid_MSM,gamma.pars[j,],HPV.inc_MSM) +
    (1-prop.MSM[j])*RiskReductionCancer(RRR_2v.mid_m,gamma.pars[j,],HPV.inc_m)
  RRRca_2v.lib1[[j]] <- prop.MSM[j]*RiskReductionCancer(RRR_2v.lib1_MSM,gamma.pars[j,],HPV.inc_MSM) +
    (1-prop.MSM[j])*RiskReductionCancer(RRR_2v.lib1_m,gamma.pars[j,],HPV.inc_m)
  RRRca_2v.lib2[[j]] <- prop.MSM[j]*RiskReductionCancer(RRR_2v.lib2_MSM,gamma.pars[j,],HPV.inc_MSM) +
    (1-prop.MSM[j])*RiskReductionCancer(RRR_2v.lib2_m,gamma.pars[j,],HPV.inc_m)
  RRRca_4v[[j]] <- prop.MSM[j]*RiskReductionCancer(RRR_4v_MSM,gamma.pars[j,],HPV.inc_MSM) +
    (1-prop.MSM[j])*RiskReductionCancer(RRR_4v_m,gamma.pars[j,],HPV.inc_m)
  RRRca_9v[[j]] <- prop.MSM[j]*RiskReductionCancer(RRR_9v_MSM,gamma.pars[j,],HPV.inc_MSM) +
    (1-prop.MSM[j])*RiskReductionCancer(RRR_9v_m,gamma.pars[j,],HPV.inc_m)
}
