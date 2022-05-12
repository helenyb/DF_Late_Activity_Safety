#Liu & Johnson
library(rjags)
source("data_generation_TTE.R")
##pre-amble functions
#for two parameter logistic


current_patient_data_frame<-function(current_time,patient.dataframe,follow_up){
  #input the whole data and convert to format for TITECRM
  #each row is a patient, with entries c("entry.time","dose.level","DLT","DLT.time","Eff","Eff.time") 
  num_patients<-max(patient.dataframe$patient_ID)
  patient.dataframe<-patient.dataframe[patient.dataframe$time_of<=current_time,]
  patient.dataframe<-patient.dataframe[patient.dataframe$cycle_num<=follow_up,]
  
  
  patient_dataframe<-matrix(NA,nrow=num_patients,ncol=6)
  for (i in 1:num_patients){
    patient_data_ind<-patient.dataframe[patient.dataframe$patient_ID==i,]
    patient_DLT<-max(patient_data_ind$DLT)
    patient_Eff<-max(patient_data_ind$Eff)
    patient_entry.time<-patient_data_ind$entry_time[1]
    patient_dose.level<-patient_data_ind$dose_level[1]
    if(patient_DLT==1){
      patient_DLT.time<-max(patient_data_ind$DLT.time,na.rm=T)
    }else{
      patient_DLT.time<-NA
    }
    
    if(patient_Eff==1){
      patient_Eff.time<-max(patient_data_ind$Eff.time,na.rm=T)
    }else{
      patient_Eff.time<-NA
    }
    patient_dataframe[i,]<-c(patient_entry.time,patient_dose.level,patient_DLT,patient_DLT.time,patient_Eff,patient_Eff.time)
  }
  
  patient_dataframe<-data.frame(patient_dataframe)
  names(patient_dataframe)<-c("entry.time","dose.level","DLT","DLT.time","Eff","Eff.time") 
  follow_up_time<-function(x) min(follow_up,x)
  patient_follow_up<- sapply(current_time-patient_dataframe$entry.time,follow_up_time)
  patient_weights_tox<-patient_follow_up/follow_up
  patient_weights_eff<-patient_follow_up/follow_up
  
  #only DLT if we have seen it at the current time
  current.DLT<-patient_dataframe$DLT
  current.DLT[which(current.DLT==1)]<-patient_dataframe$DLT.time[which(current.DLT==1)]<=current_time
  
  #only DLT time observed if current DLT is true
  current.DLT.time<-patient_dataframe$DLT.time
  current.DLT.time[which(current.DLT==0)]<-NA
  
  #only Eff if we have seen it at the current time
  current.Eff<-patient_dataframe$Eff
  current.Eff[which(current.Eff==1)]<-patient_dataframe$Eff.time[which(current.Eff==1)]<=current_time
  
  #only Eff time observed if current Eff is true
  current.Eff.time<-patient_dataframe$Eff.time
  current.Eff.time[which(current.Eff==0)]<-NA
  
  entry.time<-patient_dataframe$entry.time
  dose.level<-patient_dataframe$dose.level
  
  patient_weights_tox[which(current.DLT==1)]<-1
  patient_weights_eff[which(current.DLT==1)]<-(current.DLT.time[which(current.DLT==1)]-entry.time[which(current.DLT==1)])/follow_up
  patient_weights_eff[which(current.Eff==1)]<-1
  
  
  
  return(data.frame(entry.time,dose.level,current.DLT,current.DLT.time,current.Eff,current.Eff.time,patient_weights_tox,patient_weights_eff,patient_follow_up))
  
}

gibbs_sampler.LJ.combo<-function(data_list,iter,model_string){
  model1.spec<-textConnection(model_string)
  jags <- jags.model(model1.spec,data =data_list,n.chains=1,n.adapt=1000,quiet=T)
  update(jags, 1000,progress.bar="none")
  tt<-jags.samples(jags,c('betaT','betaE'),iter,progress.bar="none")
  return(tt)
}


#create hard safety matrix

hard.safety.mat.function<-function(perc,co.size,max.cohorts){
  hard.safety.mat<-matrix(nrow=2,ncol=max.cohorts)
  hard.safety.mat[2,]<-seq(from=co.size,to=max.cohorts*co.size, by=co.size)
  for(i in 1:max.cohorts){
    numi<- hard.safety.mat[2,i]
    try_vec<-c(1:numi)
    probs_vec<-1-pbeta(0.3,1+try_vec,1+numi-try_vec)
    hard.safety.mat[1,i]<- min( which(100*probs_vec>perc))
  }                     
  return(hard.safety.mat)
}


##INPUT:
#seed=seed for reproducability
#tru.T.pars=parameters for data generation of DLT times 
#tru.E.pars=parameters for data generation of efficacy times 
#tru.corET=correlation between DLT times and efficacy times in data generation
#co_size=cohort size
#ncohorts=maximum number of cohorts in study
#target=target toxicity (all cycles)
#targetE=target efficacy (all cycles)
#doses=real doses
#pa_T= vector of prior a's for toxicity (length=ndoses)
#pb_T= vector of prior b's for toxicity (length=ndoses)
#pa_E= vector of prior a's for efficacy (length=ndoses)
#pb_E= vector of prior b's for efficacy (length=ndoses)
#ncycles=number of cycles
#dose.skipping=enforce dose skipping rule? T=no dose skipping (default)
#sufficient.information==enforce stopping for sufficient information? no more than 9 patients per dose. T=enforce stopping for sufficient information (default)
#sufficient.information.lim= number of patients needed before stopping when the next assignment is the same. 
#hard.safety.rule=percentage for hard safety rule based on Beta(1.1)? 
#e.g. : 85= 2/3,3/6,4/9. 90=2/3,4/6,5/9, 95=3/3,4/6,5/9, <50 means no hard safety enforced
#safety.stopping.low.unsafe= (T=stop when P(p1>0.3)>0.8 (cycle 1))
#safety.stopping.high.toosafe= (T= stop when P(pJ>0.3)>0.8 (cycle 1))
#kfold.skipping: Is a kfold dose skipping rule implemented?
#kfold: the k-fold rise thsat is allowable in dose skipping
#initial.one.cycle: Is the initial period based on one cycle at a time? (F=wait until all cycles completed before next dose in initial period)
#C_eff: Dose is "admissible in efficacy if P(P(efficacy)>effbound)>C_eff
#C_tox: Dose is "admissible in safety if P(P(DLT)<toxbound)>C_tox
#effbound: Dose is "admissible in efficacy if P(P(efficacy)>effbound)>C_eff
#toxbound: Dose is "admissible in safety if P(P(DLT)<toxbound)>C_tox
#w1 & w2 : weights for utility function
#uppertox: value above which we penalise toxicity in utility
#backfill: should doses deemed safe be backfilled? (Default FALSE)
#backfill.num: How many cohorts to add as backfilling?


##OUTPUT:
# dose.rec: dose recommendation
# num.pat: number of patients
# dose.ass: number of cohorts per dose (vector)
# stop.code: stopping reason 
#1= No admissible doses
#2= Precision
#3= Max patients
#4= Sufficient information
#5= Lowest dose fails hard safety
#6= Lowest dose unsafe
#7= Highest Dose too safe
# num.DLT: number of DLTS (total)
# DLT.mat: 2xnumdoses matrix. rows for grades 1=no DLT, 2=DLT. cols for doses.
# Eff.mat: 2xnumdoses matrix. rows for grades 1=noEff, 2=Eff. cols for doses.
# duration: total trial duration (until all recruited patients are fully observed)
# max_admissable max admissable dose at the end of the trial according to hard safety only (has hard safety eliminated any?)
# dosevec: Dose assignment of cohorts in sequence order
# DLT.vec: Binary sequence of DLT outcomes (all cycles)
# EFF.vec: Binary sequence of efficacy outcomes (all cycles)
# MTD.rec: MTD recommendation (note this is not subject to safety rules, will always output a dose)


LiuJohnson_v2<-function(seed,tru.E.pars,tru.T.pars,tru.corET,co_size,ncohorts ,target,targetE,doses,
                     pa_T,pb_T,pa_E,pb_E,ncycles,dose.skipping=T,
                     sufficient.information=T,sufficient.information.lim=0,hard.safety.rule=95,safety.stopping.low.unsafe=T,
                     safety.stopping.high.toosafe=T,kfold.skipping=T,kfold=2,initial.one.cycle=T,
                     w1, w2, upper.tox,C_eff=0.2,C_tox=0.2,effbound=0.2,toxbound=0.391,
                     backfill=F,backfill.num=2){
  set.seed(seed)
  stop_vec<-rep(0,7)
  ndoses<-length(doses)
  dosevec<-c()
  patient_ID1<-1
  current.time<-0
  
  #check that the data generation inputs tally with the number of doses and cycles
  if(!all(dim(tru.T.pars)==dim(tru.E.pars))){
    stop("Incompatible truth matrices for toxicity and efficacy")
  }
  
  
  if(ndoses!=ncol(tru.T.pars)){
    stop("Incompatable truth matrix with number of doses")
  }
  
  #define all doses as admissable before any are dropped for safety
  max_admissable<-ndoses
  
  if(hard.safety.rule>50){
    hard.safety<-T
    hard.safety.mat<-hard.safety.mat.function(perc=hard.safety.rule,co.size = co_size,max.cohorts = ncohorts)
  }else{
    hard.safety<-F
  }
  
  
  expand_vec<-rep(0, ndoses)
  dose_rec<-NA
  nextdose<-1
  stop<-0
  
  
  initial<-1
  model.LJ.string <-"
model {


piT[1] = betaT[1]
piE[1] = betaE[1]

for(j in 2:ndoses){

piT[j] = piT[j-1] + (1-piT[j-1])*betaT[j]

piE[j] = piE[j-1] + (1-piE[j-1])*betaE[j]
}

for(i in 1:length(yE)){

G_piE[i]=weightE[i]*piE[patdoses[i]]
G_piT[i]=weightT[i]*piT[patdoses[i]]


yE[i] ~ dbinom(G_piE[i],1)
yT[i] ~ dbinom(G_piT[i],1)

}





for(j in 1:ndoses){
betaT[j] ~ dbeta(pa_T[j],pb_T[j])
betaE[j] ~ dbeta(pa_E[j],pb_E[j])
}

}
"


while(stop==0){
  
  
  #DATA GENERATION
  
  if(current.time==0){
    #first cohort
    all.data<-multiple_patient_generation_TTE(patient_ID1=patient_ID1,efficacy_pars =tru.E.pars,tox_pars = tru.T.pars,corET=tru.corET,ncycles=ncycles ,thenextdose=nextdose,
                                              entry_time=current.time,num_patients = co_size)
    current.time<-current.time+1
    patient_ID1<-max(all.data$patient_ID)+1
    current.data<-all.data[all.data$time_of<=current.time,]
    dosevec[current.time]<-nextdose
  }else{
    
    #subsequent cohorts
    #escalation
    all.data<-rbind(all.data,multiple_patient_generation_TTE(patient_ID1=patient_ID1,efficacy_pars =tru.E.pars,tox_pars = tru.T.pars,corET=tru.corET,ncycles=ncycles ,thenextdose=nextdose,
                                                             entry_time=current.time,num_patients = co_size)
    )
    patient_ID1<-max(all.data$patient_ID)+1
    if(backfill==T){
      #expansion
      if(nextdose>1){ #if there are doses below the next dose
        below.doses<-c(1:(nextdose-1)) #which doses are below? 
        expand.below<- expand_vec[below.doses] #have the below doses been expanded?
        # which.expand<-which[expand.below==0] #which doses should now be expanded?
        for(expand.doses in below.doses){ #for each dose that needs expanding
          if(expand.below[expand.doses]==0){
            #generate data for 2 cohorts on the dose
            all.data<-rbind(all.data,multiple_patient_generation_TTE(patient_ID1=patient_ID1,efficacy_pars =tru.E.pars,tox_pars = tru.T.pars,corET=tru.corET,ncycles=ncycles ,thenextdose=expand.doses,
                                                                     entry_time=current.time,num_patients = backfill.num*co_size))
            #update patient ID for next assignment
            patient_ID1<-max(all.data$patient_ID)+1
            #update the expansion vector to say this dose has been expanded
            expand_vec[expand.doses]<-1
          }
        }
      }
      
    }
    
    current.time<-current.time+1
    
    current.data<-all.data[all.data$time_of<=current.time,]
    dosevec[current.time]<-nextdose
    
  }
  #  patient.data<-patient_data_frame(all.data)
  
  current.patient.data<-current_patient_data_frame(current_time=current.time,patient.dataframe=all.data,follow_up=ncycles)
  
  
  ##initial period where we keep escalating until we see a DLT
  
  if(initial==1){
    
    nextdose<-min(max(current.patient.data$dose.level)+1,length(doses))
    
    ##sufficient information 
    if(sufficient.information==T){
      npats_doses<-sapply(c(1:ndoses), function(x) sum(current.patient.data$dose.level==x))
      num_ass<-tabulate(dosevec,nbins = ndoses) #number of assignments
      
      if(npats_doses[nextdose]>=sufficient.information.lim){
        stop<-4
        stop_vec[4]<-1
        dose_rec<-nextdose
        break
      }
    }
    
    if(initial.one.cycle==T){
      if(sum(current.patient.data$current.DLT)>0){
        initial<-0
      }
      
    }else{
      
      for (cyc in 1:(ncycles-1)){
        current.patient.data<-current_patient_data_frame(current_time=current.time,patient.dataframe=all.data,follow_up=ncycles)
        current.time<-current.time+1
        
        if(sum(current.patient.data$current.DLT)>0){
          initial<-0
          current.time<-current.time-1
        }
        
      }
    }
    
  }
  if(initial==0){
    
    
    
    
    #posterior 
    
    patdoses<-current.patient.data$dose.level
    weightE<-current.patient.data$patient_weights_eff
    weightT<-current.patient.data$patient_weights_tox
    
    
    
    current_data_list<-list(yE=current.patient.data$current.Eff,
                            yT=current.patient.data$current.DLT,patdoses=patdoses,
                            weightT=weightT,weightE=weightE,pa_E=pa_E,pb_E=pb_E,
                            pa_T=pa_T,pb_T=pb_T,ndoses=ndoses)
    
    
    gibbs_out<-gibbs_sampler.LJ.combo(data_list=current_data_list,iter=5000,model_string=model.LJ.string)
    betaE_mean<-rowMeans(gibbs_out$betaE)
    betaT_mean<-rowMeans(gibbs_out$betaT)
    #distribution of pT and pE
    pT_out<-matrix(nrow=ndoses,ncol=5000)
    pE_out<-matrix(nrow=ndoses,ncol=5000)
    for(dose.level in 1:ndoses){
      for(itr in 1:5000){
        pT_out[dose.level,itr]<-1-prod(1-gibbs_out$betaT[c(1:dose.level),itr,1])
        pE_out[dose.level,itr]<-1-prod(1-gibbs_out$betaE[c(1:dose.level),itr,1])
      }
    }
    
    

    max.dose.tried<-max(dosevec,na.rm=T)
    #admissable doses
    admiss<-probeff<-probtox<-c()
    for(dose.level in 1:length(doses)){
      probeff[dose.level]<-mean(pE_out[dose.level,]>effbound)
      probtox[dose.level]<-mean(pT_out[dose.level,]<toxbound)
      
      admiss[dose.level]<-  (probeff[dose.level]>C_eff)&(probtox[dose.level]>C_tox)
    }
    admiss.doses<-which(admiss)
    
    utility<-c()
    for(dose.level in admiss.doses){
      utility[dose.level]<-(1-prod(1-betaE_mean[c(1:dose.level)])) - 
        w1*(1-prod(1-betaT_mean[c(1:dose.level)])) -
        w2*(1-prod(1-betaT_mean[c(1:dose.level)]))*as.numeric((1-prod(1-betaT_mean[c(1:dose.level)]))>upper.tox)
      
    }
    

    ##choosing the next dose:
    
    nextdose<-which.max(utility)
    

    ##if no admissable doses
    if(length(nextdose)==0){
      stop<-1
      stop_vec[1]<-1
      dose_rec<-NA
      break
    }    
    
    #k fold dose increase on EXPERIMENTED DOSES
    if(kfold.skipping==T){
      
      prevdoseval<-doses[max(dosevec,na.rm=T)]
      nextdoseval<-doses[nextdose]
      
      if(nextdoseval>(kfold*prevdoseval)){
        maxdoseval<-kfold*prevdoseval
        max_dose<-max(which(doses<=maxdoseval))
        #ensures we don't escalate to a less "optimal dose" by the dose skipping restriction
        nextdose<-which.max(utility[c(1:max_dose)])
      }
    }
    
    #stopping based on beta-binomial
    npats_doses<-sapply(c(1:ndoses), function(x) sum(current.patient.data$dose.level==x))
    
    if(((safety.stopping.low.unsafe==T)&(npats_doses[1]>0))|((safety.stopping.high.toosafe==T)&(npats_doses[ndoses]>0))){
      
      
      current.patient.data_cyc1<-current_patient_data_frame(current_time=current.time,patient.dataframe=all.data,follow_up=1)
      

      weightE1<-weightT1<-c()
      patdoses1<-current.patient.data_cyc1$dose.level
      
      
      #formatting for gibbs sampler
      for(obser in 1:nrow(current.patient.data_cyc1)){
        weightE1[obser]<-current.patient.data_cyc1$patient_weights_eff[obser]
        weightT1[obser]<-current.patient.data_cyc1$patient_weights_tox[obser]
        
      }
      
      
      current_data_list1<-list(yE=current.patient.data_cyc1$current.Eff,
                               yT=current.patient.data_cyc1$current.DLT,patdoses=patdoses1,
                               weightT=weightT1,weightE=weightE1,pa_E=pa_E,pb_E=pb_E,
                               pa_T=pa_T,pb_T=pb_T,ndoses=ndoses)
      
      
      gibbs_out1<-gibbs_sampler.LJ.combo(data_list=current_data_list1,iter=5000,model_string=model.LJ.string)
      betaE_mean1<-rowMeans(gibbs_out1$betaE)
      betaT_mean1<-rowMeans(gibbs_out1$betaT)
      pT_out1<-matrix(nrow=ndoses,ncol=5000)
      pE_out1<-matrix(nrow=ndoses,ncol=5000)
      for(dose.level in 1:ndoses){
        for(itr in 1:5000){
          pT_out1[dose.level,itr]<-1-prod(1-gibbs_out1$betaT[c(1:dose.level),itr,1])
          pE_out1[dose.level,itr]<-1-prod(1-gibbs_out1$betaE[c(1:dose.level),itr,1])
        }
      }
      
      
      
      #only cycle 1 
   
      if((safety.stopping.low.unsafe==T)&(npats_doses[1]>0)){
        posterior1<- pT_out1[1,]
        
        posterior1<-posterior1[1:length(posterior1)]
        cyc1_0.3g<-mean(posterior1>0.3,na.rm=T)
        
        if(cyc1_0.3g>0.8){
          stop<-6
          stop_vec[6]<-1
          nextdose<-NA
          # break
        }
      }
      
      if((safety.stopping.high.toosafe==T)&(npats_doses[ndoses]>0)){
        posteriorJ<-pT_out1[ndoses,]
        posteriorJ<-posteriorJ[1:length(posteriorJ)]
        cycJ_0.3l<-mean(posteriorJ<0.3,na.rm=T)
        if(cycJ_0.3l>0.8){

          stop<-7
          stop_vec[7]<-1
          nextdose<-NA
          #  break
        }
      }
      
    }
    
    
    
    
    #number of DLTs per dose level
    nDLTs_doses<-sapply(c(1:ndoses),function(x) sum(current.data$DLT[(current.data$dose_level==x)&(current.data$patient>0)]))
    nDLTs_doses_c1<-sapply(c(1:ndoses),function(x) sum(current.data$DLT[(current.data$dose_level==x)&(current.data$patient>0)&(current.data$cycle==1)]))
    
    
    ##hard safety
    if(hard.safety==T){
      
      explored<-c(1:ndoses)[npats_doses>0]

      for(do in explored){
        
        if(nDLTs_doses_c1[do]>=hard.safety.mat[1, which(hard.safety.mat[2,]==npats_doses[do])]){
          
          max_admissable<-min(c(max_admissable,do-1))
        }
      }
      if(!is.na(nextdose)){
        if(nextdose>max_admissable){
          nextdose<-max_admissable
        }
      }
      if(max_admissable==0){
        stop<-5
        stop_vec[5]<-1
        dose_rec<-NA
        nextdose<-NA

      }
      
      
    }
    
    
    
    
    ##sufficient information 
    if((sufficient.information==T)&(!is.na(nextdose))){
      if(npats_doses[nextdose]>=sufficient.information.lim){
        
        stop<-4
        stop_vec[4]<-1
        dose_rec<-nextdose

      }
    }

    ##precision stopping not applicable
    
  }
  if(nrow(current.patient.data)==(ncohorts*co_size)){ #max patients reached
    
    
    current.patient.data<-current_patient_data_frame(current_time=(current.time+ncycles),patient.dataframe=all.data,follow_up = ncycles)
    patdoses<-current.patient.data$dose.level
    
    patdoses<-current.patient.data$dose.level
    weightE<-current.patient.data$patient_weights_eff
    weightT<-current.patient.data$patient_weights_tox
    
    
    
    current_data_list<-list(yE=current.patient.data$current.Eff,
                            yT=current.patient.data$current.DLT,patdoses=patdoses,
                            weightT=weightT,weightE=weightE,pa_E=pa_E,pb_E=pb_E,
                            pa_T=pa_T,pb_T=pb_T,ndoses=ndoses)
    
    
    gibbs_out<-gibbs_sampler.LJ.combo(data_list=current_data_list,iter=5000,model_string=model.LJ.string)
    betaE_mean<-rowMeans(gibbs_out$betaE)
    betaT_mean<-rowMeans(gibbs_out$betaT)
    pT_out<-matrix(nrow=ndoses,ncol=5000)
    pE_out<-matrix(nrow=ndoses,ncol=5000)
    for(dose.level in 1:ndoses){
      for(itr in 1:5000){
        pT_out[dose.level,itr]<-1-prod(1-gibbs_out$betaT[c(1:dose.level),itr,1])
        pE_out[dose.level,itr]<-1-prod(1-gibbs_out$betaE[c(1:dose.level),itr,1])
      }
    }
    
    
    max.dose.tried<-max(dosevec,na.rm=T)
    #admissable doses
    admiss<-probeff<-probtox<-c()
    for(dose.level in 1:length(doses)){
      probeff[dose.level]<-mean(pE_out[dose.level,]>effbound)
      probtox[dose.level]<-mean(pT_out[dose.level,]<toxbound)
      
      admiss[dose.level]<-  (probeff[dose.level]>C_eff)&(probtox[dose.level]>C_tox)
    }
    admiss.doses<-which(admiss)
    
    utility<-c()
    for(dose.level in admiss.doses){
      utility[dose.level]<-(1-prod(1-betaE_mean[c(1:dose.level)])) - 
        w1*(1-prod(1-betaT_mean[c(1:dose.level)])) -
        w2*(1-prod(1-betaT_mean[c(1:dose.level)]))*as.numeric((1-prod(1-betaT_mean[c(1:dose.level)]))>upper.tox)
      
    }
    
    
    ##choosing the dose:

    dose_rec<-which.max(utility)
    stop<-3
    stop_vec[3]<-1
  }
  
}


#follow up for all patients
current.data<-all.data
current.time<-max(current.data$time_of)
current.patient.data<-current_patient_data_frame(current_time = current.time,patient.dataframe = all.data,follow_up = ncycles)



#grade matrix out
DLT.matrix.out<-matrix(0,ncol=ndoses,nrow=2)
for (pat in 1:(max(current.data$patient))){
  pat_max_grade<-max(current.data[current.data$patient_ID==pat,]$DLT)
  pat_dose<-current.data[current.data$patient_ID==pat,]$dose_level[1]
  DLT.matrix.out[pat_max_grade+1,pat_dose]<- DLT.matrix.out[pat_max_grade+1,pat_dose]+1
}

#grade matrix out
Eff.matrix.out<-matrix(0,ncol=ndoses,nrow=2)
for (pat in 1:(max(current.data$patient))){
  pat_max_grade<-max(current.data[current.data$patient_ID==pat,]$Eff)
  pat_dose<-current.data[current.data$patient_ID==pat,]$dose_level[1]
  Eff.matrix.out[pat_max_grade+1,pat_dose]<- Eff.matrix.out[pat_max_grade+1,pat_dose]+1
}

MTD.rec<-min(c(max(dosevec),max_admissable,which.min(abs(target-sapply(c(1:ndoses),function(x) (1-prod(1-betaT_mean[c(1:x)])))))))


#output list
output<-list(
  dose.rec=dose_rec ,#dose recommendation
  num.pat=max(current.data$patient) ,#: number of patients
  dose.ass=tabulate(dosevec,nbins = ndoses) ,# : number of cohorts per dose (vector)
  stop.code=stop_vec, #: stopping reason 
  num.DLT= sum(current.data$DLT[current.data$patient_ID>0]),#: number of DLTS (total)
  DLT.mat=DLT.matrix.out ,#: 2xnumdoses matrix. rows for grades 1=no DLT, 2=DLT. cols for doses.
  Eff.mat=Eff.matrix.out ,#: 2xnumdoses matrix. rows for grades 1=noEff, 2=Eff. cols for doses.
  duration= current.time,#: total trial duration (until all recruited patients are fully observed)
  max_admissable=max_admissable, # max admissable dose at the end of the trial (has hard safety eliminated any?)
  dosevec=dosevec,
  DLT.vec=current.patient.data$current.DLT,
  EFF.vec=current.patient.data$current.Eff,
  MTD.rec=MTD.rec
)


return(output)

}

