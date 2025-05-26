#Parameters and data needed for the computation

##################################
### GENERAL (FIXED) PARAMETERS ###
##################################

cohort <- 100000 #cohort size per gender (at age a0)
a0 <- 10 #age of vaccination
amax <- 99
adiff <- amax - a0
types <- 10 #16,18,31,33,35,39,45,51,52,58 (oncogenic types only)

#7 age groups in screening: 30-34, 35-39, 40-44, 45-49, 50-54, 55-59, 60+ (60-64)
agegrps.scr_start <- 30
agegrps.scr_end <- 60
agegrps.scr <- seq(agegrps.scr_start,agegrps.scr_end,by=5)

#15 age groups in cancer diagnoses: 15-19, 20-24, 25-29,..., 85-89
agegrps.ca_start <- floor(a0/5)+2
agegrps.ca_end <- floor(85/5)+1
agegrps.ca <- c(agegrps.ca_start:agegrps.ca_end)
ca.ages <- seq((floor(a0/5)+1)*5,85,5)
agemids.ca <- agegrps.ca*5 - 2

#Discount rates
#d_costs <- 1.04
d_costs <- 1.03 #new Dutch guidelines
d_health <- 1.015
dvec_c <- 1/(d_costs^(1:adiff))
dvec_h <- 1/(d_health^(1:adiff))

CEthreshold <- 20000

#Population/survival
surv_w <- read.table("Data/TablePop_Women2020.txt",header=T)
S_w <- diag(adiff + 1)
for (i in 1:adiff){
  for (j in (i+1):(adiff+1)){
    S_w[i,j] <- surv_w[a0+j,]/surv_w[a0+i,]
  }
}

surv_m <- read.table("Data/TablePop_Men2020.txt",header=T)
S_m <- diag(adiff + 1)
for (i in 1:adiff){
  for (j in (i+1):(adiff+1)){
    S_m[i,j] <- surv_m[a0+j,]/surv_m[a0+i,]
  }
}

#vaccine prices (2-dose)
vacc.price_2v <- 70 #(admin. costs now indexed to 2024) 
vacc.price_4v <- 80
vacc.price_9v <- 120 #price difference €50

##########################
### TRANSMISSION MODEL ###
##########################

load(paste("TransmissionModel/VaccineUptake",vacc.up[1],",",vacc.up[2],"/HPVinc_pre.RData", sep = ""))
load(paste("TransmissionModel/VaccineUptake",vacc.up[1],",",vacc.up[2],"/HPVinc_pre_set.RData", sep = "")) #pre vacc
load(paste("TransmissionModel/VaccineUptake",vacc.up[1],",",vacc.up[2],"/HPVinc_post.RData", sep = "")) #post vacc
load(paste("TransmissionModel/VaccineUptake",vacc.up[1],",",vacc.up[2],"/HPVinc_post_waning.RData", sep = "")) #post vacc with extreme waning
load(paste("TransmissionModel/VaccineUptake",vacc.up[1],",",vacc.up[2],"/HPVprev_pre_set.RData", sep = "")) #pre vacc
load(paste("TransmissionModel/VaccineUptake",vacc.up[1],",",vacc.up[2],"/HPVprev_post.RData", sep = "")) #post vacc
load(paste("TransmissionModel/VaccineUptake",vacc.up[1],",",vacc.up[2],"/HPVprev_post_waning.RData", sep = "")) #post vacc with extreme waning

HPV.inc_w <- HPVinc_pre_f[c(1,2,3,4,5,6,7,8,9,11),]
HPV.inc_m <- HPVinc_pre_m[c(1,2,3,4,5,6,7,8,9,11),]
HPV.inc_MSM <- HPV.inc_m

#############
### WARTS ###
#############

#File AnogenitalWarts_Incidence.R makes AW.inc_w and AW.inc_m; vectors of length 100 (ages 1:100) 
#Incidence per 1000 person-years
source('Rcode/AnogenitalWarts_Incidence.R')

cost_AW_m <- 138.6 #now indexed to 2024 #research by UMCG 2021
cost_AW_w <- 138.6
QALY_AW <- 0 #base-case LYs only
#QALY_AW <- 0.018 #Westra et al. 2013
PAF.HPV_AW <- 0.9

#Uncertainty on warts (data based on approx 10% of population)
#inc below age 12 and above age 84 is set to 0 (deterministic)
born2019 <- 0.1*169680/2
pop_w_age <- c(born2019*surv_w[13:85,]/surv_w[1,])
pop_m_age <- c(born2019*surv_m[13:85,]/surv_m[1,])
nr.warts_w <- AW.inc_w[12:84]*pop_w_age/1000
nr.warts_m <- AW.inc_m[12:84]*pop_m_age/1000
AW.inc_w_sim <- mapply(rbeta,shape1=nr.warts_w+0.5,shape2=pop_w_age-nr.warts_w+0.5,
                       SIMPLIFY = TRUE,MoreArgs=list(n=nr.sim))
AW.inc_m_sim <- mapply(rbeta,shape1=nr.warts_m+0.5,shape2=pop_m_age-nr.warts_m+0.5,
                       SIMPLIFY = TRUE,MoreArgs=list(n=nr.sim))
AW.inc_w_sim <- cbind(matrix(0,nr.sim,11),AW.inc_w_sim,matrix(0,nr.sim,15))
AW.inc_m_sim <- cbind(matrix(0,nr.sim,11),AW.inc_m_sim,matrix(0,nr.sim,15))

###########
### RRP ###
###########

#Prevalence per 100.000
PrevUnder18_RRP_sim <- rpois(nr.sim,lambda=2.5)
PrevAbove18_RRP_sim <- rpois(nr.sim,lambda=2)
nrYearsAverage_RRP <- 10

#costs indexed to 2024
cost_RRP <- 2038 #Yearly cost per RRP patient
QALY_RRPlifetime <- 0 #life time QALY loss. Base-case LYs only
#QALY_RRPlifetime <- 1.05
QALY_RRP <- QALY_RRPlifetime/nrYearsAverage_RRP #Yearly QALY loss
PAF.HPV_RRP <- 1

#Age distribution (based on San Giorgi et al 2015)
x <- 1:100
m2 <- dlnorm(x,meanlog = 3.57, sdlog = sqrt(0.09))
m3 <- dlnorm(x,meanlog = 4.15, sdlog = sqrt(0.02))
m1 <- dlnorm(x,meanlog = 1.94, sdlog = sqrt(1.39))
#plot(x,m1,type="l",col="black")
#lines(x,m2,col="red")
#lines(x,m3,col="green")
mixture.probs <- c(0.19,0.61,0.2)
age.distr_RRP <- mixture.probs[1]*m1 + mixture.probs[2]*m2 + mixture.probs[3]*m3
age.distr_RRP <- age.distr_RRP/(sum(age.distr_RRP))
#plot(x,age.distr_RRP,type="l",
#     main="Age distribution RRP",xlab="age",ylab="")

#Births per 1000 women
births_age_w <- read.table("Data/BirthsPer1000_AgeWomen2020.txt",header=T)

#################
### SCREENING ###
#################

#Colposcopies
#costs indexed to 2024
cost_colpo <- 609 #costs CIN0 Jansen2020
#PPV per age group
PPV_colpo_no <- c(0.659,rep(0.402,6)) #From Inturrisi IJC2020 table S3
PPV_colpo_2v <- c(0.548,rep(0.403,6))
PPV_colpo_2v.cp <- c(0.538,rep(0.384,6))
PPV_colpo_9v <- c(0.356,rep(0.342,6))

cost_CIN2 <- 1855 # costs from own adjusted computation
cost_CIN3 <- 2226
cost_CIN2plus <- 0.68*cost_CIN3 + 0.32*cost_CIN2
QALY_CIN2plus <- 0 #Base-case LYs only
#QALY_CIN2plus <- 0.03/12

#screening uptake per age group (beschermingsgraad from monitor bevolkinsonderzoek 2019)
scr.up_sim <- matrix(rep(c(0.648,0.707,0.733,0.729,0.769,0.771,0.762),nr.sim),nr.sim,7,byrow=TRUE)

#Fraction of lesions and cancer per age grp (from RIVM data, see RIVMdata_statistics.R)
# Lesions (CIN2 and CIN3) without cancer
nr.lesions <- c(2924, 1637, 1086, 1012,  740,  420,  248)
nr.cancer <- c(80, 50, 47, 49, 24,  9, 13)
nr.women <- c(61799, 61348, 63462, 75097, 79572, 79224, 66513)
fr.lesions_sim <- mapply(rbeta,shape1=nr.lesions+0.5,shape2=nr.women-nr.lesions+0.5,SIMPLIFY = TRUE,MoreArgs=list(n=nr.sim))
fr.cancer_sim <- mapply(rbeta,shape1=nr.cancer+0.5,shape2=nr.women-nr.cancer+0.5,SIMPLIFY = TRUE,MoreArgs=list(n=nr.sim))

#Type attribution as fraction of all CIN2+ lesions (from POBASCAM and VUSA data, see Computation_AFsCIN2.R)
#This is including cancers, because VUSA data only considers CIN3+ (not cancers separately)
n_age1 <- 1019 #nr of women in data <34
n_age2 <- 1990 #nr of women in data >33

#Last (10th) entry contains all other types
AFs_CIN2plus_age1_est <- c(0.5628, 0.0619, 0.1144, 0.0691, 0.0058, 0.0101, 0.0274, 0.0364, 0.0402, 0.0360, 0.0359)
AFs_CIN2plus_age2_est <- c(0.4698, 0.0807, 0.1317, 0.0665, 0.0313, 0.0253, 0.0233, 0.0482, 0.0569, 0.0427, 0.0236)
n_CIN2plus_age1 <- AFs_CIN2plus_age1_est*n_age1 #non integers, but not a problem
n_CIN2plus_age2 <- AFs_CIN2plus_age2_est*n_age2

alpha_dirichlet <- rep(0.5,types+1) #prior dirichlet distribution AFs
AFs_CIN2plus_age1_sim <- rdirichlet(nr.sim, (alpha_dirichlet+n_CIN2plus_age1))[,1:types]
AFs_CIN2plus_age2_sim <- rdirichlet(nr.sim, (alpha_dirichlet+n_CIN2plus_age2))[,1:types]

##############
### CANCER ###
##############

#General population size
pop_w <- read.table("Data/PopWomen_Agegrp.txt",header=T)
pop_m <- read.table("Data/PopMen_Agegrp.txt",header=T)

#Mean population over years 2015-2019
pop2015.2019_w <- as.numeric(apply(pop_w[6:10,agegrps.ca + 2],2,mean)) #last group only age 85-89
pop2015.2019_m <- as.numeric(apply(pop_m[6:10,agegrps.ca + 2],2,mean)) #last group only age 85-89

#Cancer incidence
inc_cervix <- read.table("Data/Inc_Cervix.txt",header=T)
inc_anus_w <- read.table("Data/Inc_Anus_Women.txt",header=T)
inc_anus_m <- read.table("Data/Inc_Anus_Men.txt",header=T)
inc_oroph_w <- read.table("Data/Inc_Oropharynx_Women.txt",header=T)
inc_oroph_m <- read.table("Data/Inc_Oropharynx_Men.txt",header=T)
inc_vulva <- read.table("Data/Inc_Vulva.txt",header=T)
inc_vagina <- read.table("Data/Inc_Vagina.txt",header=T)
inc_penis <- read.table("Data/Inc_Penis.txt",header=T)

#Mean incidence over years 2015-2019
inc2015.2019_cervix <- as.numeric(apply(inc_cervix[6:10,agegrps.ca + 2],2,mean))
inc2015.2019_anus_w <- as.numeric(apply(inc_anus_w[6:10,agegrps.ca + 2],2,mean))
inc2015.2019_anus_m <- as.numeric(apply(inc_anus_m[6:10,agegrps.ca + 2],2,mean))
inc2015.2019_oroph_w <- as.numeric(apply(inc_oroph_w[6:10,agegrps.ca + 2],2,mean))
inc2015.2019_oroph_m <- as.numeric(apply(inc_oroph_m[6:10,agegrps.ca + 2],2,mean))
inc2015.2019_vulva <- as.numeric(apply(inc_vulva[6:10,agegrps.ca + 2],2,mean))
inc2015.2019_vagina <- as.numeric(apply(inc_vagina[6:10,agegrps.ca + 2],2,mean))
inc2015.2019_penis <- as.numeric(apply(inc_penis[6:10,agegrps.ca + 2],2,mean))

#Yearly risk of cancer diagnosis
ca.risk_cervix_sim <- mapply(rbeta,shape1=inc2015.2019_cervix+0.5,shape2=pop2015.2019_w-inc2015.2019_cervix+0.5,SIMPLIFY = TRUE,MoreArgs=list(n=nr.sim))
ca.risk_anus_w_sim <- mapply(rbeta,shape1=inc2015.2019_anus_w+0.5,shape2=pop2015.2019_w-inc2015.2019_anus_w+0.5,SIMPLIFY = TRUE,MoreArgs=list(n=nr.sim))
ca.risk_anus_m_sim <- mapply(rbeta,shape1=inc2015.2019_anus_m+0.5,shape2=pop2015.2019_m-inc2015.2019_anus_m+0.5,SIMPLIFY = TRUE,MoreArgs=list(n=nr.sim))
ca.risk_oroph_w_sim <- mapply(rbeta,shape1=inc2015.2019_oroph_w+0.5,shape2=pop2015.2019_w-inc2015.2019_oroph_w+0.5,SIMPLIFY = TRUE,MoreArgs=list(n=nr.sim))
ca.risk_oroph_m_sim <- mapply(rbeta,shape1=inc2015.2019_oroph_m+0.5,shape2=pop2015.2019_m-inc2015.2019_oroph_m+0.5,SIMPLIFY = TRUE,MoreArgs=list(n=nr.sim))
ca.risk_vulva_sim <- mapply(rbeta,shape1=inc2015.2019_vulva+0.5,shape2=pop2015.2019_w-inc2015.2019_vulva+0.5,SIMPLIFY = TRUE,MoreArgs=list(n=nr.sim))
ca.risk_vagina_sim <- mapply(rbeta,shape1=inc2015.2019_vagina+0.5,shape2=pop2015.2019_w-inc2015.2019_vagina+0.5,SIMPLIFY = TRUE,MoreArgs=list(n=nr.sim))
ca.risk_penis_sim <- mapply(rbeta,shape1=inc2015.2019_penis+0.5,shape2=pop2015.2019_m-inc2015.2019_penis+0.5,SIMPLIFY = TRUE,MoreArgs=list(n=nr.sim))

#survival data
survdata_cervix <- read.table("Data/Surv_Cervix.txt",header=T)
survdata_anus <- read.table("Data/Surv_Anus.txt",header=T)
survdata_anus[2,13] <- 71 #must be non-increasing
survdata_anus[4,11:13] <- c(60,58,56)
survdata_anus[5,10] <- 46
survdata_oroph <- read.table("Data/Surv_Oropharynx.txt",header=T)
survdata_vulva <- read.table("Data/Surv_Vulva.txt",header=T)
survdata_vulva[1,11] <- 89
survdata_vulva[3,11:13] <- c(77,76,75)
survdata_vagina <- read.table("Data/Surv_Vagina.txt",header=T)
survdata_vagina[1:2,2] <- 50
survdata_vagina[1:3,9:13] <- rbind(c(68,68,68,67,67),c(68,68,68,67,67),c(68,68,68,67,67))
survdata_vagina[4,11:13] <- c(37,36,35)
survdata_vagina[5,9:13] <- c(29,27,25,24,23)
survdata_penis <- read.table("Data/Surv_Penis.txt",header=T)
survdata_penis[1,2] <- 50
survdata_penis[3,10:12] <- c(78,78,78)

#Simulate survival probabilities
#Every row in _sim contains one simulation of the 50 probabilities in the survival matrix
nr.risk_cervix <- as.vector(as.matrix(survdata_cervix[,3:12])*matrix(rep(survdata_cervix[,2],10),5,10)/100)
nr.surv_cervix <- as.vector(as.matrix(survdata_cervix[,4:13])*matrix(rep(survdata_cervix[,2],10),5,10)/100)
surv_probs_cervix_sim <- mapply(rbeta,shape1=nr.surv_cervix+0.5,shape2=nr.risk_cervix-nr.surv_cervix+0.5,SIMPLIFY = TRUE,MoreArgs=list(n=nr.sim))
nr.risk_anus <- as.vector(as.matrix(survdata_anus[,3:12])*matrix(rep(survdata_anus[,2],10),5,10)/100)
nr.surv_anus <- as.vector(as.matrix(survdata_anus[,4:13])*matrix(rep(survdata_anus[,2],10),5,10)/100)
surv_probs_anus_sim <- mapply(rbeta,shape1=nr.surv_anus+0.5,shape2=nr.risk_anus-nr.surv_anus+0.5,SIMPLIFY = TRUE,MoreArgs=list(n=nr.sim))
nr.risk_oroph <- as.vector(as.matrix(survdata_oroph[,3:12])*matrix(rep(survdata_oroph[,2],10),5,10)/100)
nr.surv_oroph <- as.vector(as.matrix(survdata_oroph[,4:13])*matrix(rep(survdata_oroph[,2],10),5,10)/100)
surv_probs_oroph_sim <- mapply(rbeta,shape1=nr.surv_oroph+0.5,shape2=nr.risk_oroph-nr.surv_oroph+0.5,SIMPLIFY = TRUE,MoreArgs=list(n=nr.sim))
nr.risk_vulva <- as.vector(as.matrix(survdata_vulva[,3:12])*matrix(rep(survdata_vulva[,2],10),5,10)/100)
nr.surv_vulva <- as.vector(as.matrix(survdata_vulva[,4:13])*matrix(rep(survdata_vulva[,2],10),5,10)/100)
surv_probs_vulva_sim <- mapply(rbeta,shape1=nr.surv_vulva+0.5,shape2=nr.risk_vulva-nr.surv_vulva+0.5,SIMPLIFY = TRUE,MoreArgs=list(n=nr.sim))
nr.risk_vagina <- as.vector(as.matrix(survdata_vagina[,3:12])*matrix(rep(as.numeric(survdata_vagina[,2]),10),5,10)/100)
nr.surv_vagina <- as.vector(as.matrix(survdata_vagina[,4:13])*matrix(rep(as.numeric(survdata_vagina[,2]),10),5,10)/100)
surv_probs_vagina_sim <- mapply(rbeta,shape1=nr.surv_vagina+0.5,shape2=nr.risk_vagina-nr.surv_vagina+0.5,SIMPLIFY = TRUE,MoreArgs=list(n=nr.sim))
nr.risk_penis <- as.vector(as.matrix(survdata_penis[,3:12])*matrix(rep(as.numeric(survdata_penis[,2]),10),5,10)/100)
nr.surv_penis <- as.vector(as.matrix(survdata_penis[,4:13])*matrix(rep(as.numeric(survdata_penis[,2]),10),5,10)/100)
surv_probs_penis_sim <- mapply(rbeta,shape1=nr.surv_penis+0.5,shape2=nr.risk_penis-nr.surv_penis+0.5,SIMPLIFY = TRUE,MoreArgs=list(n=nr.sim))

#Costs associated with treatment and palliative care respectively. Kok et al. 2011
#Indexing from 2011 -> 2024 (CPI 93.73% to 130.31% = *1.39)
costs_cervix <- c(8000,19600)*1.39
costs_anus_w <- c(5000,18800)*1.39
costs_anus_m <- c(5000,19500)*1.39
costs_oroph_w <- c(6000,19500)*1.39
costs_oroph_m <- c(6000,19600)*1.39
costs_vulva <- c(8000,16600)*1.39
costs_vagina <- c(8000,16600)*1.39
costs_penis <- c(4000,19500)*1.39

#QALY loss per cancer diagnosis: in base-case =0
QALY_cervix <- 0
QALY_anus <- 0
QALY_oroph <- 0
QALY_vulva <- 0
QALY_vagina <- 0
QALY_penis <- 0

#AF HPV to cancer and type-specific AFs
PAF.HPV_cervix <- 1
n_cervix <- c(1348,150,69,117,46,27,80,28,40,27,126) #From Sanjose et. al. 2010 (Europe)
AFs_cervix_sim <- rdirichlet(nr.sim, (alpha_dirichlet+n_cervix))[,1:types]

PAF.HPV_anus_sim <- rbeta(nr.sim,148+0.5,169-148+0.5)
n_anus <- c(354,16,8,12,7,2,4,0,3,8,24)
AFs_anus_sim <- rdirichlet(nr.sim, (alpha_dirichlet+n_anus))[,1:types]

PAF.HPV_oroph_sim <- rbeta(nr.sim,281+0.5,926-281+0.5)
n_oroph <- c(1650,33,6,44,24,3,7,0,3,11,90)
AFs_oroph_sim <- rdirichlet(nr.sim,(alpha_dirichlet+n_oroph))[,1:types]

PAF.HPV_vulva_sim <- rbeta(nr.sim,165+0.5,903-165+0.5)
n_vulva <- c(311,20,4,28,0,3,14,0,8,4,35)
AFs_vulva_sim <- rdirichlet(nr.sim,(alpha_dirichlet+n_vulva))[,1:types]

PAF.HPV_vagina_sim <- rbeta(nr.sim,108+0.5,152-108+0.5)
n_vagina <- c(178,15,16,15,3,6,11,7,9,11,32)
AFs_vagina_sim <- rdirichlet(nr.sim,(alpha_dirichlet+n_vagina))[,1:types]

PAF.HPV_penis_sim <- rbeta(nr.sim,135+0.5,419-135+0.5)
n_penis <- c(229,5,3,10,9,2,9,3,5,4,54)
AFs_penis_sim <- rdirichlet(nr.sim,(alpha_dirichlet+n_penis))[,1:types]

#HRs survival HPV+ vs HPV- cancer
HR_cervix <- 1
HR_vagina <- 1
#HR_anus <- 0.36 #Urbute et al. 2020 (95% CI = [0.22,0.58])
#HR_oroph <- 0.46 #O’Rorke et al. 2012 (95% CI = [0.37,0.57])
#HR_vulva <- 0.64 #Zhang et al. 2018 (95% CI = [0.47,0.87])
#HR_penis <- 0.2 #Djajadiningrat et al. 2015 (95% CI = [0.1,0.9])

compute.sd_HR <- function(HR,CI1,CI2){
  sd_HR <- (((log(CI2)-log(HR))/1.96) + ((log(HR)-log(CI1))/1.96))/2
  return(sd_HR)
}
#Sample logHR from normal distributions
HR_anus_sim <- exp(rnorm(nr.sim,mean=log(0.36),sd=compute.sd_HR(0.36,0.22,0.58)))
HR_oroph_sim <- exp(rnorm(nr.sim,mean=log(0.46),sd=compute.sd_HR(0.46,0.37,0.57)))
HR_vulva_sim <- exp(rnorm(nr.sim,mean=log(0.64),sd=compute.sd_HR(0.64,0.47,0.87)))
HR_penis_sim <- exp(rnorm(nr.sim,mean=log(0.2),sd=compute.sd_HR(0.2,0.1,0.9)))

#Distribution (gamma) of T2 (time from HPV to cancer)
#Estimated with Estimation_TimetoCancer.R
gamma.par_cervix <- c(5.39,3.62)
gamma.par_anus_w <- c(20.91,1.91)
gamma.par_anus_m <- c(13.09,2.92)
gamma.par_oroph_w <- c(24.78,1.58)
gamma.par_oroph_m <- c(25.73,1.47)
gamma.par_vulva <- c(8.19,5.00)
gamma.par_vagina <- c(10.33,3.90)
gamma.par_penis <- c(22.27,1.96)

gamma.pars <- rbind(gamma.par_cervix,gamma.par_anus_w,gamma.par_anus_m,
                    gamma.par_oroph_w,gamma.par_oroph_m,
                    gamma.par_vulva,gamma.par_vagina,gamma.par_penis)

#Parameters w.r.t MSM group
p.MSM <- 0.06 #proportion of MSM in male population

#Increased risk on cancers (see Frisch et al. 2003)
RR.MSM_oroph <- 2.5
RR.MSM_penis <- 1
#we assume all anal cancers are MSM

