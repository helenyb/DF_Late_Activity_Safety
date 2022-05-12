##a function to do simulations for Yuan & Yin (2008) method in comparison setting

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



full2current<-function(fulldata,current.time,follow.up){
  patient.number<-fulldata$patient.number
  entry.time<-fulldata$entry.time
  dose.level<-fulldata$dose.level
  cycle.obs<-10*sapply(current.time-entry.time, function(x) min(x,follow.up))
  y_E<-sapply(patient.number, function(x) min(c(fulldata$t_E[x],0.1*cycle.obs[x])))
  y_T<-sapply(patient.number, function(x) min(c(fulldata$t_T[x],0.1*cycle.obs[x])))
  delta_E<-fulldata$t_E<=0.1*cycle.obs
  delta_T<-fulldata$t_T<=0.1*cycle.obs
  for(pats in 1:nrow(fulldata)){
    if((delta_T[pats]==T)&(y_T[pats]<y_E[pats])){
      y_E[pats]<-y_T[pats]
      delta_E[pats]<-F
    }
  }
  
  return(data.frame(patient.number, y_T, y_E,  delta_T,delta_E, entry.time, cycle.obs, dose.level))
}




##preamble function
gibbs_sampler.yuan2<-function(data_list,model_string){
  
  model1.spec<-textConnection(model_string)
  jags <- jags.model(model1.spec,data =data_list,n.chains=1,quiet=T, inits= list(alpha_T=1, beta_T=0.5,lambda_T=0.01, alpha_E=1, beta_E=0.5,lambda_E=0.01,pi=0.8,phi=1.5,".RNG.name" = "base::Wichmann-Hill"))
  update(jags, 500,progress.bar="none")
  tt<-jags.samples(jags,c('beta_T','alpha_T','lambda_T','beta_E','alpha_E','lambda_E','pi','phi'),n.iter=2000,progress.bar="none")
  return(tt)
}


#generic function for S_E and S_T

S_function<-function(x, alpha, beta, lambda, Z){
  exp(-lambda*(x^alpha)*exp(beta*Z))
}

S_star<-function(x, alpha, beta, lambda, Z, pi){
  1- pi + pi*S_function(x, alpha, beta, lambda, Z)
}


#formatting for gibbs sampler
#yOUT 
#1= no E no T
#2= no E yes T
#3= yes E no T
#4= tes E yes T
category_function<-function(E_vec,T_vec){
  length_out<-length(E_vec)
  length_out2<-length(T_vec)
  if(length_out!=length_out2){
    # browser()
    stop("incompatible lengths")
  }
  # browser()
  out_vec<-c()
  for(i in 1:length_out){
    
    if(E_vec[i]+T_vec[i]==2){
      out_vec[i]<-4
    }else if(E_vec[i]+T_vec[i]==0){
      out_vec[i]<-1
    }else if((E_vec[i]==1)&(T_vec[i]==0)){
      out_vec[i]<-3
    }else if((E_vec[i]==0)&(T_vec[i]==1)){
      out_vec[i]<-2
    }
    
  }
  
  return(out_vec)
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



#convert our dataframe into format for Yuan (current time)

fulldata2YuanFormat_current<-function(full_data,current.time){
  dose_level<-patient.number<-t_T<-t_E<-y_T<-y_E<-delta_T<-delta_E<-entry.time<-cycle_obs<-c()
  full_data<-full_data[full_data$time_of<=current.time,]
  #input is dataframe with names c("patient_ID","cycle_num","dose_level","entry_time","time_of","DLT","Eff")
  # browser()
  for(patient in 1:max(full_data$patient_ID)){
    patient_data<-full_data[full_data$patient_ID==patient,]
    patient.number[patient]<-patient
    cycle_obs[patient]<-max(patient_data$cycle_num)
    delta_T[patient]<-as.numeric(sum(patient_data$DLT)>0)
    delta_E[patient]<-as.numeric(sum(patient_data$Eff)>0)
    dose_level[patient]<-patient_data$dose_level[1]
    if(sum(patient_data$DLT)>0){
      y_T[patient]<-max(patient_data$DLT.time,na.rm=T)-patient_data$entry_time[1]
    }else{
      y_T[patient]<-length(patient_data$DLT)
    }
    if(sum(patient_data$Eff)>0){
      y_E[patient]<- max(patient_data$Eff.time,na.rm=T)-patient_data$entry_time[1]
    }else{
      y_E[patient]<-y_T[patient]
    }
    entry.time[patient]<-patient_data$entry_time[1]
  }
  return(data.frame(patient.number=patient.number,y_T=0.1*y_T,y_E=0.1*y_E,entry.time=0.1*entry.time,delta_T=delta_T,delta_E=delta_E,dose.level=dose_level))
}
# 

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
#precision.stopping: Is a precision stopping rule implemented? (Default=TRUE)
#precision.stopping.level:  number of assigned cohorts needed before precision stopping can be triggered
#initial.one.cycle: Is the initial period based on one cycle at a time? (F=wait until all cycles completed before next dose in initial period)
#C_eff: Dose is "admissible in efficacy if P(P(efficacy)>effbound)>C_eff
#C_tox: Dose is "admissible in safety if P(P(DLT)<toxbound)>C_tox
#effbound: Dose is "admissible in efficacy if P(P(efficacy)>effbound)>C_eff
#toxbound: Dose is "admissible in safety if P(P(DLT)<toxbound)>C_tox
#p_cross: Dose is escalated in initial period if P(P(DLT)<toxbound)>p_cross
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

Yuan2008_comp_v4<-function(seed,tru.E.pars,tru.T.pars,tru.corET,co_size,ncohorts ,target,targetE,doses,ncycles,dose.skipping=T,
                           sufficient.information=T,sufficient.information.lim=29,hard.safety.rule=95,safety.stopping.low.unsafe=T,
                           safety.stopping.high.toosafe=T,kfold.skipping=T,kfold=2,precision.stopping=T,precision.stopping.level=29,initial.one.cycle=T,
                           C_eff=0.2,C_tox=0.2,effbound=0.2,toxbound=0.391,p_cross=0.8,
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
  #model for jags
  model.string <-"
model {


C<-100

for(i in 1:length(y_E)){


S_y_T_val[i]<-   exp(-lambda_T*(y_T[i]^alpha_T)*exp(beta_T*Z[i]))
S_y_E_val[i]<-   exp(-lambda_E*(y_E[i]^alpha_E)*exp(beta_E*Z[i]))

S_joint_val[i]<-  ( (S_y_T_val[i])^(-1/phi) + (S_y_E_val[i])^(-1/phi) - 1 )^(-phi)




dS_dy_T_val[i]<- (-lambda_T*exp(beta_T*Z[i]))*alpha_T*(y_T[i]^(alpha_T-1))*exp((y_T[i]^alpha_T)*(-lambda_T*exp(beta_T*Z[i])))
dS_dy_E_val[i]<- (-lambda_E*exp(beta_E*Z[i]))*alpha_E*(y_E[i]^(alpha_E-1))*exp((y_E[i]^alpha_E)*(-lambda_E*exp(beta_E*Z[i])))


dSJ_dy_E_val[i]<- dS_dy_E_val[i]* (S_y_E_val[i]^(-(1+phi)/phi))*(S_y_E_val[i]^(-1/phi) + S_y_T_val[i]^(-1/phi) - 1 )^(-(phi+1))
dSJ_dy_T_val[i]<- dS_dy_T_val[i]* (S_y_T_val[i]^(-(1+phi)/phi))*(S_y_T_val[i]^(-1/phi) + S_y_E_val[i]^(-1/phi) - 1 )^(-(phi+1))


d2SJ_yTyE_val[i]<- dS_dy_T_val[i]*dS_dy_E_val[i]*((phi+1)/phi)*(S_y_T_val[i]^(-(1+phi)/phi))*(S_y_E_val[i]^(-(1+phi)/phi))*((S_y_E_val[i]^(-1/phi) + S_y_T_val[i]^(-1/phi) -1 )^(-(phi+2)))



L1[i]<- pi*d2SJ_yTyE_val[i]

L2[i]<- -(1-pi)*dS_dy_T_val[i] - pi*dSJ_dy_T_val[i]

L3[i]<- -pi*dSJ_dy_E_val[i]

L4[i]<- (1-pi)*S_y_T_val[i] + pi*S_joint_val[i]


Lik[i]<-(L1[i]^(delta_T[i]*delta_E[i]))*(L2[i]^(delta_T[i]*(1-delta_E[i])))*(L3[i]^((1-delta_T[i])*delta_E[i]))*(L4[i]^((1-delta_T[i])*(1-delta_E[i])))

p[i]<-Lik[i]/C

ones[i] ~ dbern(p[i])

}






lambda_T ~ dgamma(0.1,0.1) T(0.001,2)
lambda_E ~ dgamma(0.1,0.1) T(0.001,3)

alpha_T ~ dgamma(0.1,0.1) T(0.1,3)
alpha_E ~ dgamma(0.1,0.1) T(0.1,4)

beta_T ~ dgamma(0.1,0.1) T(0.01,1) 
beta_E ~ dgamma(0.1,0.1) T(0.01,1)

phi ~ dunif(0.1,5) 
pi ~ dunif(0.6,0.99)


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
     #format into dataset for current time (Yuan format)
    current_data<-fulldata2YuanFormat_current(full_data=all.data,current.time=current.time)
    #format into dataset for current time (Standard format)
    current.data<-all.data[all.data$time_of<=current.time,]
    dosevec[current.time]<-nextdose
  }else{
    
    #subsequent cohorts
    #escalation
    all.data<-rbind(all.data,multiple_patient_generation_TTE(patient_ID1=patient_ID1,efficacy_pars =tru.E.pars,tox_pars = tru.T.pars,corET=tru.corET,ncycles=ncycles ,thenextdose=nextdose,
                                                             entry_time=current.time,num_patients = co_size)
    )
    current_data<-fulldata2YuanFormat_current(full_data=all.data,current.time=current.time)
    
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
            # all.data.Yuan<- fulldata2YuanFormat(full_data=all.data,followup=ncycles)
            
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
    #format into dataset for current time (Yuan format)
    current_data<-fulldata2YuanFormat_current(full_data=all.data,current.time=current.time)
    
    dosevec[current.time]<-nextdose
    
  }

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
    
    dose_obs_lev<-current_data$dose.level
    dose_obs<-sapply(dose_obs_lev,function(x) doses[x]) 
    #format current data for gibbs sampler
    currentdata_list<-list(y_E=current_data$y_E , y_T= current_data$y_T, delta_T= current_data$delta_T, delta_E=current_data$delta_E , Z= dose_obs, ones=rep(1,length(dose_obs)))

    gibbs.out<-gibbs_sampler.yuan2(currentdata_list,model_string=model.string)
    max.dose.tried<-max(current_data$dose.level)
    
    post.tox.prob<-mean((1-S_function(0.1*ncycles,gibbs.out$alpha_T, gibbs.out$beta_T, gibbs.out$lambda_T,doses[max.dose.tried]))<toxbound)
    #should we escalate?
    if(post.tox.prob>p_cross){
      nextdose<-min(c(length(doses),max.dose.tried+1))

    }else{
      #otherwise choose using ATAE
      
      admiss<-probeff<-probtox<-c()
      for(dose.level in 1:length(doses)){
        probeff[dose.level]<-mean((1-S_star(0.1*ncycles,gibbs.out$alpha_E, gibbs.out$beta_E, gibbs.out$lambda_E,doses[dose.level],gibbs.out$pi))>effbound)
        probtox[dose.level]<-mean((1-S_function(0.1*ncycles,gibbs.out$alpha_T, gibbs.out$beta_T, gibbs.out$lambda_T,doses[dose.level]))<toxbound)
        
        admiss[dose.level]<-  (probeff[dose.level]>C_eff)&(probtox[dose.level]>C_tox)
      }
      admiss.doses<-which(admiss)
      
      if(length(admiss.doses)==0){
        nextdose<-numeric(0)
      }else{
        ATAE<-((1/mean(gibbs.out$alpha_T))*((mean(gibbs.out$lambda_T)*exp(doses[admiss.doses]*mean(gibbs.out$beta_T)))^(-1/mean(gibbs.out$alpha_T)))*pgamma((mean(gibbs.out$lambda_T)*exp(doses[admiss.doses]*mean(gibbs.out$beta_T))*0.1*ncycles^mean(gibbs.out$alpha_T)),1/mean(gibbs.out$alpha_T))*gamma(1/mean(gibbs.out$alpha_T)))/
          ((1-mean(gibbs.out$pi))*0.1*ncycles + mean(gibbs.out$pi)*((1/mean(gibbs.out$alpha_E))*((mean(gibbs.out$lambda_E)*exp(doses[admiss.doses]*mean(gibbs.out$beta_E)))^(-1/mean(gibbs.out$alpha_E)))*pgamma((mean(gibbs.out$lambda_E)*exp(doses[admiss.doses]*mean(gibbs.out$beta_E))*0.1*ncycles^mean(gibbs.out$alpha_E)),1/mean(gibbs.out$alpha_E))*gamma(1/mean(gibbs.out$alpha_E))))
        
        ##choosing the next dose:
        
        nextdose<-admiss.doses[which(ATAE==max(ATAE,na.rm = T))]
      }
    }
    
    

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
        admiss.doses.maxed<-admiss.doses[admiss.doses<=max_dose]
        ATAE.maxed<-ATAE[admiss.doses<=max_dose]
        nextdose<-admiss.doses.maxed[which.max(ATAE.maxed)]
        if(length(nextdose)==0){ #if only admissable doses are too high
          
          if(probtox[max_dose]>C_tox){
            nextdose<-max_dose #can escalate to a non-admissable dose if safe (in order to get to the higher doses eventually)
          }else{ #otherwise stop for no admissible doses
          stop<-1
          stop_vec[1]<-1
          dose_rec<-NA
          break
          }
        }
      }
    }
    
    #number of patients on doses
    npats_doses<-sapply(c(1:ndoses), function(x) sum(current.patient.data$dose.level==x))
    
    
    #safety stopping for lowest dose unsafe and highest dose too safe
    if(((safety.stopping.low.unsafe==T)&(npats_doses[1]>0))|((safety.stopping.high.toosafe==T)&(npats_doses[ndoses]>0))){
      
      
      #only uses cycle 1 data
      current.patient.data_cyc1<-all.data[(all.data$cycle_num==1)&(all.data$time_of<=current.time),]
      current_data.cyc1Yuan<-fulldata2YuanFormat_current(full_data=current.patient.data_cyc1,current.time=current.time)
      
      
      
      currentdata_list1<-list(y_E=current_data$y_E , y_T= current_data$y_T, delta_T= current_data$delta_T, delta_E=current_data$delta_E , Z= dose_obs, ones=rep(1,length(dose_obs)))
      
      
      #only cycle 1 
      
      currentdata_list1<-list(y_E=current_data.cyc1Yuan$y_E , y_T= current_data.cyc1Yuan$y_T, delta_T= current_data.cyc1Yuan$delta_T, delta_E=current_data.cyc1Yuan$delta_E , Z= dose_obs, ones=rep(1,length(dose_obs)))

      gibbs_out1<-gibbs_sampler.yuan2(currentdata_list1,model_string=model.string)
      #if lowest unsafe rule is being used and patients on lowest dose
      if((safety.stopping.low.unsafe==T)&(npats_doses[1]>0)){
        posterior1<- (1-S_function(0.1,gibbs_out1$alpha_T, gibbs_out1$beta_T, gibbs_out1$lambda_T,doses[1]))
        
        
        #expit((gibbs_out1$alphaT)+doses[1]*(gibbs_out1$betaT))
        
        posterior1<-posterior1[1:length(posterior1)]
        cyc1_0.3g<-mean(posterior1>0.3,na.rm=T)
        
        if(cyc1_0.3g>0.8){
          stop<-6
          stop_vec[6]<-1
          nextdose<-NA
          # break
        }
      }
      #if highest too safe rule is being used and patients on highest dose
      
      if((safety.stopping.high.toosafe==T)&(npats_doses[ndoses]>0)){
        posteriorJ<-(1-S_function(0.1,gibbs_out1$alpha_T, gibbs_out1$beta_T, gibbs_out1$lambda_T,doses[ndoses]))
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
    
    
    
    
    #number of DLTs per dose level (all and c1)
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
        # break
      }
      
      
    }
    
    
    
    
    ##sufficient information 
    if((sufficient.information==T)&(!is.na(nextdose))){
      if(npats_doses[nextdose]>=sufficient.information.lim){
        
        stop<-4
        stop_vec[4]<-1
        dose_rec<-nextdose
        #   break
      }
    }
    
    ##precision stopping
    if(precision.stopping==T){
    if((sum(npats_doses)>=precision.stopping.level)&(!is.na(nextdose))){
      if(npats_doses[nextdose]>0){
        
        MTD_vectorT<-(1/gibbs.out$beta_T)*log(log(1-target)/(-gibbs.out$lambda_T*(0.1*ncycles)^gibbs.out$alpha_T))
        CV_T<-mad(MTD_vectorT, center = median(MTD_vectorT), constant = 1.4826, 
                  na.rm = FALSE, low = FALSE, high = FALSE)/median(MTD_vectorT)
        
        eff_vectorE<-(1/gibbs.out$beta_E)*log(log((gibbs.out$pi-targetE)/gibbs.out$pi)/(-gibbs.out$lambda_E*(0.1*ncycles)^gibbs.out$alpha_E))
        
        CV_E<-mad(eff_vectorE, center = median(eff_vectorE), constant = 1.4826, 
                  na.rm = FALSE, low = FALSE, high = FALSE)/median(eff_vectorE)
        if((CV_E<0.3)&(CV_T<0.3)){
          dose_rec<-nextdose
          stop<-2
          stop_vec[2]<-1
        }
        
        
        
        
      }
      
    }
    }
  }
  if(nrow(current.patient.data)==(ncohorts*co_size)){ #max patients reached
    
    
    current.patient.data<-current_patient_data_frame(current_time=(current.time+ncycles),patient.dataframe=all.data,follow_up = ncycles)
    dose_obs_lev<-current_data$dose.level
    dose_obs<-sapply(dose_obs_lev,function(x) doses[x]) 
    
    currentdata_list<-list(y_E=current_data$y_E , y_T= current_data$y_T, delta_T= current_data$delta_T, delta_E=current_data$delta_E , Z= dose_obs, ones=rep(1,length(dose_obs)))
    

    gibbs.out<-gibbs_sampler.yuan2(currentdata_list,model_string=model.string)
    max.dose.tried<-max(current_data$dose.level)
    
    
    
    
    
    admiss<-probeff<-probtox<-c()
    for(dose.level in 1:length(doses)){
      probeff[dose.level]<-mean((1-S_star(0.1*ncycles,gibbs.out$alpha_E, gibbs.out$beta_E, gibbs.out$lambda_E,doses[dose.level],gibbs.out$pi))>effbound)
      probtox[dose.level]<-mean((1-S_function(0.1*ncycles,gibbs.out$alpha_T, gibbs.out$beta_T, gibbs.out$lambda_T,doses[dose.level]))<toxbound)
      admiss[dose.level]<-  (probeff[dose.level]>C_eff)&(probtox[dose.level]>C_tox)
      
    }
    admiss.doses<-which(admiss)
    if(sum(admiss.doses)==0){
      stop<-1

      break
    }
    ATAE<-((1/mean(gibbs.out$alpha_T))*((mean(gibbs.out$lambda_T)*exp(doses[admiss.doses]*mean(gibbs.out$beta_T)))^(-1/mean(gibbs.out$alpha_T)))*pgamma((mean(gibbs.out$lambda_T)*exp(doses[admiss.doses]*mean(gibbs.out$beta_T))*0.1*ncycles^mean(gibbs.out$alpha_T)),1/mean(gibbs.out$alpha_T))*gamma(1/mean(gibbs.out$alpha_T)))/
      ((1-mean(gibbs.out$pi))*0.1*ncycles + mean(gibbs.out$pi)*((1/mean(gibbs.out$alpha_E))*((mean(gibbs.out$lambda_E)*exp(doses[admiss.doses]*mean(gibbs.out$beta_E)))^(-1/mean(gibbs.out$alpha_E)))*pgamma((mean(gibbs.out$lambda_E)*exp(doses[admiss.doses]*mean(gibbs.out$beta_E))*0.1*ncycles^mean(gibbs.out$alpha_E)),1/mean(gibbs.out$alpha_E))*gamma(1/mean(gibbs.out$alpha_E))))
    
    ##choosing the recommended dose:
    dose_rec<-admiss.doses[which(ATAE==max(ATAE,na.rm = T))]
    
    stop<-3
    stop_vec[3]<-1
  }
  
}


#follow up for all patients
current.data<-all.data
current.time<-max(current.data$time_of)
current.patient.data<-current_patient_data_frame(current_time = current.time,patient.dataframe = all.data,follow_up = ncycles)



#DLT matrix out
DLT.matrix.out<-matrix(0,ncol=ndoses,nrow=2)
for (pat in 1:(max(current.data$patient))){
  pat_max_grade<-max(current.data[current.data$patient_ID==pat,]$DLT)
  pat_dose<-current.data[current.data$patient_ID==pat,]$dose_level[1]
  DLT.matrix.out[pat_max_grade+1,pat_dose]<- DLT.matrix.out[pat_max_grade+1,pat_dose]+1
}

#Efficacy matrix out
Eff.matrix.out<-matrix(0,ncol=ndoses,nrow=2)
for (pat in 1:(max(current.data$patient))){
  pat_max_grade<-max(current.data[current.data$patient_ID==pat,]$Eff)
  pat_dose<-current.data[current.data$patient_ID==pat,]$dose_level[1]
  Eff.matrix.out[pat_max_grade+1,pat_dose]<- Eff.matrix.out[pat_max_grade+1,pat_dose]+1
}
#calculate MTD for additional output
MTD.rec<-min(c(max(dosevec),max_admissable,which.min(abs(target-(1-exp(-mean(gibbs.out$lambda_T)*((0.1*ncycles)^(mean(gibbs.out$alpha_T)))*exp(mean(gibbs.out$beta_T)*doses)))))))


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

