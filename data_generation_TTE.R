#function to give cumulative sums of rows of matrix
matrix_cumulative_row_sums_function<-function(matrix_in){
  nrows<-nrow(matrix_in)
  ncols<-ncol(matrix_in)
  
  matrix_out<-matrix(nrow=nrows,ncol=ncols)
  
  for(i in 1:ncols){
    matrix_out[,i]<-cumsum(matrix_in[,i])
  }
  
  return(matrix_out)
}


##find the parameters

find_lognormal_parms3<-function(p1,p3,int1=seq(-20,20,0.05),int2=seq(0.05,20,0.05)){
  #p1=P(X<=1), p3=P(X<=3)
  
  lognormal_cdf<-function(parms,quantile,prob){
    mu<-parms[1]
    sig<-parms[2]
    prob_q<-pnorm(q=log(quantile),mean=mu,sd=sig)
    return(abs(prob_q-prob))
  }
  mean_i<-0
  matrix1<-matrix(0,nrow=length(int1),ncol=length(int2))
  for(mean_val in int1){
    mean_i<-mean_i+1
    sd_i<-0
    for(sd_val in int2){
      sd_i<-sd_i+1
      if(lognormal_cdf(parms=c(mean_val,sd_val),quantile = 1,prob=p1)<0.01){
        matrix1[mean_i,sd_i]<-1
        if(lognormal_cdf(parms=c(mean_val,sd_val),quantile = 3,prob=p3)<0.01){
          matrix1[mean_i,sd_i]<-2
        }
      }
    }
  }
  
  
  
  closes<-which(matrix1==2,arr.ind = T)
  close_vals1<-close_vals3<-c()
  for(closes_i in 1:nrow(closes)){
    close_vals1[closes_i]<-lognormal_cdf(parms=c(int1[closes[closes_i,1]],int2[closes[closes_i,2]]),quantile = 1,prob=p1)
    close_vals3[closes_i]<-lognormal_cdf(parms=c(int1[closes[closes_i,1]],int2[closes[closes_i,2]]),quantile = 3,prob=p3)
  }
  final_choice<-closes[which.min(close_vals1+close_vals3),]
  return(c(int1[final_choice[1]],
           int2[final_choice[2]]))
  
}

#generates a single patient's COMPLETE outcome as a data-frame
#patient: patient ID for index of rands
#means and sd for Efficacy and Toxicity
#ncycles: number of cycles
#thenext dose: dose level
#entry_time: time the patient enters the trial


single_patient_generation_TTE<-function(patient,efficacy_pars,tox_pars,corET,ncycles,thenextdose,entry_time){
  meanE<-efficacy_pars[1,thenextdose]
    meanT<-tox_pars[1,thenextdose]
    sdE<-efficacy_pars[2,thenextdose]
    sdT<-tox_pars[2,thenextdose]
   #first element is Efficacy, second is Toxicity
  mean_vec<-c(meanE,meanT)
  var_mat<-matrix(c(sdE^2,sdE*sdT*corET,sdE*sdT*corET,sdT^2),nrow=2)
  timesET<-exp(rmvnorm(1,mean=mean_vec,sigma=var_mat))
  if(timesET[2]<timesET[1]){#if toxicity occurs before efficacy, we don't see the efficacy
    timesET[1]<-ncycles+1
  }
  
  max_cycle_obs<-min(c(ncycles,ceiling(timesET[2])))
  
  
  Eff.time<-DLT.time<-rep(NA,max_cycle_obs)
  
  
  Eff.obs<-as.numeric(timesET[1]<c(1:max_cycle_obs))
  DLT.obs<-as.numeric(timesET[2]<c(1:max_cycle_obs))
  
  
  Eff.time[Eff.obs==1]<-timesET[1]+entry_time
  DLT.time[DLT.obs==1]<-timesET[2]+entry_time
  
  #browser()
  out_data<-data.frame(patient_ID=rep(patient,max_cycle_obs),cycle_num=c(1:max_cycle_obs),dose_level=rep(thenextdose,max_cycle_obs),
                       entry_time=rep(entry_time,max_cycle_obs),time_of=entry_time+c(1:max_cycle_obs),
                       DLT=DLT.obs,Eff=Eff.obs,DLT.time=DLT.time,Eff.time=Eff.time)
  
  
  return(out_data)
}



#generates multiple patients' COMPLETE outcome as a data-frame
multiple_patient_generation_TTE<-function(patient_ID1,efficacy_pars,tox_pars,corET,ncycles,thenextdose,entry_time,num_patients){
  meanE<-efficacy_pars[1,thenextdose]
  meanT<-tox_pars[1,thenextdose]
  sdE<-efficacy_pars[2,thenextdose]
  sdT<-tox_pars[2,thenextdose]
  
   out_data<-single_patient_generation_TTE(patient=patient_ID1,efficacy_pars=efficacy_pars,tox_pars=tox_pars,corET=corET,ncycles=ncycles,thenextdose=thenextdose,entry_time=entry_time)
  for(i in 2:num_patients){
    out_data<-rbind(out_data,single_patient_generation_TTE(patient=patient_ID1+i-1,efficacy_pars=efficacy_pars,tox_pars=tox_pars,corET=corET,ncycles=ncycles,thenextdose=thenextdose,
                                                           entry_time=entry_time)
    )
  }
  return(out_data)
}




## translate the 1st cycle true P(DLT) into a matrix
cyc_func_tox<-function(cyc1_vec,cyc_dec,ncycles){
  out_mat<-matrix(nrow=ncycles,ncol=length(cyc1_vec))
  for(j in 1:ncol(out_mat)){
    cyc_all_vec<-c()
    cyc_all_vec[1]<-cyc1_vec[j]
    for ( i in 2:ncycles){
      cyc_all_vec[i]<-prod(1-cyc_all_vec[1:(i-1)])*cyc_all_vec[1]*(cyc_dec)^(i-1)
    }
    out_mat[,j]<-cyc_all_vec
  }
  
  return(out_mat)
}

##translate efficacy vector into matrix - cycle per row
cyc_func_eff<-function(cyc1_vec,split_vec){
  if(sum(split_vec)!=1){
    stop("Split of Efficacy probabilities does not sum to 1")
  }
  
  out_mat<-matrix(nrow=length(split_vec),ncol=length(cyc1_vec))
  for(i in 1:length(cyc1_vec)){
    out_mat[,i]<-split_vec*cyc1_vec[i]
  }
  return(out_mat)
}