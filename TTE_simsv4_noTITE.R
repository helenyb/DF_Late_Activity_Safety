##combo trials

##comparison simulations
source("CRM_2_combo_v12.R")


library(rjags)
library(doParallel)
library(mvtnorm)
registerDoParallel(cores=64)


##parameters for the scenarios

eff_scen1<-c(0.2,0.3,0.4,0.5,0.6,0.7)
eff_scen2<-c(0.3,0.4,0.5,0.5,0.5,0.5)
eff_scen3<-c(0.10,0.12,0.14,0.16,0.18,0.20)
eff_scen4<-c(0.10,0.15,0.20,0.30,0.50,0.70)

tox_scen1<-c(0.1,0.2,0.3,0.4,0.5,0.6)
tox_scen2<-c(0.1,0.13,0.16,0.2,0.25,0.4)
tox_scen3<-c(0.4,0.45,0.5,0.55,0.6,0.65)
tox_scen4<-c(0.3,0.4,0.45,0.5,0.55,0.6)
tox_scen5<-c(0.10,0.12,0.14,0.16,0.18,0.20)

eff.pattern1<-c(1,0,0)
# eff.pattern2<-c(1:3)/6
# eff.pattern3<-c(3:1)/6
ndoses<-6
ncycles<-3 

for(pattern in 1){
  for(eff_scen in 1:4){
    parMat<-matrix(nrow=2,ncol=ndoses)
    cycleMat<-     cyc_func_eff(cyc1_vec = get(paste(c("eff_scen",eff_scen),collapse="")),
                                split_vec = get(paste(c("eff.pattern",pattern),collapse="")))
    for(j in 1:ncol(cycleMat)){
      parMat[,j]<- find_lognormal_parms3(p1=cycleMat[1,j],p3=min(2*cycleMat[1,j],(1+cycleMat[1,j])/2),int2=seq(0.01,10,0.01))
      
    }
    assign(paste(c("eff_scen",eff_scen,".",pattern,"pars"),collapse=""),
           parMat)
  }
}

for(tox in 1:5){
  parMat<-matrix(nrow=2,ncol=ndoses)
  cycleMat<-cyc_func_tox(cyc1_vec=get(paste(c("tox_scen",tox),collapse="")),cyc_dec = 1/3,ncycles=ncycles)
  for(j in 1:ncol(cycleMat)){
    parMat[,j]<-  find_lognormal_parms3(p1=sum(cycleMat[,j]),p3=min(2*cycleMat[1,j],(1+cycleMat[1,j])/2),int2=seq(0.01,10,0.01))
    
    
  }
  
  assign(paste(c("tox_scen",tox,".pars"),collapse=""),
         parMat)
}

nsims<-1000
doses<-c(1.5,2.5,3.5,4.5,6.0,7.0)
#MTD.options<-c("MODEL_ONLY","CHOICE_ESCALTION")
stop.num<-29
BmeamT<-log(1/4)
BprecT<-0.5
ncycles<-1 # code as one cycle to use no TITE version

for(eff.pattern in 1){
  for(eff in 1:4){
    for(tox in 1:5){
      assign(paste(c("comparison_combo2v4noTITE_cyc3_eff",eff,".",eff.pattern,"_tox",tox),collapse=""),
             foreach(i=1:nsims, combine = list) %dopar% {
               ##function
               
               
               CRM_2_combo_v12(seed=i,tru.E.pars = get(paste(c("eff_scen",eff,".",pattern,"pars"),collapse="")),
                               tru.T.pars=get(paste(c("tox_scen",tox,".pars"),collapse="")),
                               tru.corET=-0.5,co_size=3,ncohorts=20 ,target=0.391,targetE=0.3,doses=doses,p_a_meanT=log(1/16),p_b_meanT=BmeamT,
                               p_a_precT=1,p_b_precT=BprecT,ncycles=1,dose.skipping=T,
                               sufficient.information=T,sufficient.information.lim=stop.num,hard.safety.rule=95,safety.stopping.low.unsafe=T,
                               safety.stopping.high.toosafe=T,kfold.skipping=T,kfold=2,precision.stopping="SAFETY_AND_EFFICACY",precision.stopping.level=stop.num,initial.one.cycle=T,
                               p_a_meanE=-3,p_b_meanE=-0.2,effbound=0.2,
                               p_a_precE=1,p_b_precE=1,p_phi_mean=0,p_phi_prec=0.01, w1=0.33, w2=1.09, upper.tox=0.392,TITE=F)
               
               
               
             }
             
             
      )
      
      save.image(paste(c("Comparison_combo_cycles_TTE_v4_noTITE.RData"),collapse=""))
    }
  }
}

