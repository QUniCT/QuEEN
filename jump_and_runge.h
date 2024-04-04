///////////////////////////////////////////////////////////////////////
int jump_and_runge()
{
 int i,j,ran_int=0;
 long midum=1000;
 float per_thread;
//,rand_num;
 int *ran_seed;
///////////////////////////////////////////////////////////////////////
//////  split the trajectories into seperate threads  /////////////////
///////////////////////////////////////////////////////////////////////
 #pragma omp parallel
 {
  int i;
  #pragma omp single
  {  
   i=omp_get_num_threads();
   if(i<num_threads||num_threads==0)num_threads=i;
   if(num_threads>Ntrajectories){
    num_threads=Ntrajectories;
    per_thread=1;
   }else{
    per_thread=Ntrajectories/num_threads;
    Ntrajectories=num_threads*per_thread;
   }
   printf("Number of trajectories is %d to split evenly over the %d threads: %f per thread\n",Ntrajectories,num_threads,per_thread); 
   Ntrajectories=per_thread;
  }
 }
 ran_seed=(int*)calloc(num_threads,sizeof(int));
 if (ran_seed==NULL){
  printf("\nAllocating error for *ran_seed\n");
  printf("\nRequired size is %d sizeof(double)",Dim);
  exit(0);
 } 
///////////////////////////////////////////////////////////////////////
//////  set seeds for each thread  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////
 for(i=0;i<num_threads;i++){
  ran_seed[i]=midum;
  for(j=0;j<(Ntrajectories*Nsteph*NstepDt);j++)seed2(&midum);
 }
///////////////////////////////////////////////////////////////////////
//////  start the parallel section  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////
  #pragma omp parallel\
   num_threads(num_threads)\
   default(none)\
   shared(Ntrajectories,Dim,x0,NstepDt,Nsteph,Dt,\
    add_white_noise_dephasing,dph_gamma,\
    add_ou_process_dephasing,oup_sigma,oup_tau,\
    add_static_noise_dephasing,snoise,\
    add_filtered_noise_dephasing,filtered_noise_width,filtered_noise_site,filtered_noise_freqres,filtered_noise_order,filtered_noise_processes,\
    add_pl_wave_dephasing,pwd_amp,pwd_freq,pwd_phase,\
    add_fm_signal_dephasing,fmn_amp,fmn_base,fmn_signal,fmn_dev,\
    add_sq_wave_dephasing,sqw_amp,sqw_freq,sqw_samp,\
    add_sitecorr_wn,add_sitecorr_sn,add_sitecorr_ou,site_corr,site_corr_rate,site_anticorr_rate,\
    use_pulses,pulses_count,pulses_site_i,pulses_site_j,pulses_amp,pulses_now,pulses_wth,\
    use_qa,qa_a_n,qa_a_bigA,qa_a_tp1,qa_a_bigT1,qa_a_a,qa_a_b,qa_a_c,qa_a_tp2,qa_a_bigT2,qa_MatrixQ,\
    decay,decay_1,decay_rate,h,receptive_state,donor_state,vectorgammasum,\
    count_jumps_intosink,count_jumps,count_dephase,population_sink,\
    per_thread,num_threads,\
    ran_seed,ran_int,\
    rho_sys_pointer,MatrixA,\
    calculate_sd,standard_deviation)
 {
  int jump,jlvl,jindex,sink_thread,count_ou=0,sqw_num,sqw_k; 
  long int count_jumps_thread=0,count_jumps_intosink_thread=0,count_dephase_thread=0;
  int seed=Ntrajectories*Nsteph*NstepDt*omp_get_thread_num()*time(NULL);
  long randk;
  int i,j,n,m,k,p;
  float rand_num;
  double **MatrixH,**MatrixHqa,*pb_thread,norm,aux,*population_sink_thread,**standevs,**standevm,*save_ou_number,*ou_const,*ou_corr,zero=10e-10,time_point;
  double *filtered_noise,**filtered_noise_memory;
 _Complex double *k1,*k2,*k3,*k4,*auxx,partial,*xp_thread,*xn_thread,*rho_sys_thread,***rho_sys_thread_pointer;
  static long idum2=123456789;
  static long iy=0;
  static long iv[32];
  float temp;
  long midum=1000;
  long IM1=2147483563,IM2=2147483399,IMM1=2147483562;
  int IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,NTAB=32;
  float EPS=1.2e-7,RNMX=1.0-1.2e-7;
  double AM=1.0/2147483563,NDIV=1+2147483562/32;
  double pulses_x;
  double qa_At;
///////////////////////////////////////////////////////////////////////
//////  set random seeds  /////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
//  #pragma omp critical
//  {
//   seed=ran_seed[ran_int];
//   ran_int++;
//  }
///////////////////////////////////////////////////////////////////////
//////  allocate thread arrays  ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////
//  hamiltonian  //////////////////////////////////////////////////////
//  seed=rand_r(&seed);
  MatrixH=(double**)calloc(Dim,sizeof(double));
  if(MatrixH==NULL){
   printf("\nAllocating error for *MatrixH\n");
   printf("\nRequired size is %d sizeof(double)",Dim);
   exit(0);
  }
  for(i=0;i<Dim;i++){
   *(MatrixH+i)=(double*)calloc(Dim,sizeof(double));
   if(*(MatrixH+i)==NULL){
    printf("\nAllocating error for *(MatrixH+i)\n");
    printf("\nRequired size is %d sizeof(double)",Dim);
    exit(0);
   }
  }
  MatrixHqa=(double**)calloc(Dim,sizeof(double));
  if(MatrixHqa==NULL){
   printf("\nAllocating error for *MatrixH\n");
   printf("\nRequired size is %d sizeof(double)",Dim);
   exit(0);
  }
  for(i=0;i<Dim;i++){
   *(MatrixHqa+i)=(double*)calloc(Dim,sizeof(double));
   if(*(MatrixHqa+i)==NULL){
    printf("\nAllocating error for *(MatrixH+i)\n");
    printf("\nRequired size is %d sizeof(double)",Dim);
    exit(0);
   }
  }
//  wavefunction  /////////////////////////////////////////////////////
  for(i=0;i<Dim;i++)for(j=0;j<Dim;j++)MatrixH[i][j]=MatrixA[i][j];
  xp_thread=(_Complex double*)calloc(Dim,sizeof(_Complex double));
  if (xp_thread==NULL){
   printf("\nAllocating error for *xp_thread\n");
   printf("\nRequired size is %d sizeof(double)",Dim);
   exit(0);
  } 
  xn_thread=(_Complex double*)calloc(Dim,sizeof(_Complex double));
  if (xn_thread==NULL){
   printf("\nAllocating error for *xn_thread\n");
   printf("\nRequired size is %d sizeof(double)",Dim);
   exit(0);
  } 
//  runge-kutta parameters  ///////////////////////////////////////////
  k1=(_Complex double*)calloc(Dim,sizeof(_Complex double));
  if (k1==NULL){
   printf("\nAllocating error for *k1\n");
   printf("\nRequired size is %d sizeof(double)",Dim);
   exit(0);
  } 
  k2=(_Complex double*)calloc(Dim,sizeof(_Complex double));
  if (k2==NULL){
   printf("\nAllocating error for *k2\n");
   printf("\nRequired size is %d sizeof(double)",Dim);
   exit(0);
  } 
  k3=(_Complex double*)calloc(Dim,sizeof(_Complex double));
  if (k3==NULL){
   printf("\nAllocating error for *Dim\n");
   printf("\nRequired size is %d sizeof(double)",Dim);
   exit(0);
  } 
  k4=(_Complex double*)calloc(Dim,sizeof(_Complex double));
  if (k4==NULL){
   printf("\nAllocating error for *k4\n");
   printf("\nRequired size is %d sizeof(double)",Dim);
   exit(0);
  } 
  auxx=(_Complex double*)calloc(Dim,sizeof(_Complex double));
  if (auxx==NULL){
   printf("\nAllocating error for *auxx\n");
   printf("\nRequired size is %d sizeof(double)",Dim);
   exit(0);
  } 
//  jump probabilities  ///////////////////////////////////////////////
  pb_thread=(double*)calloc(decay,sizeof(double));
  if(pb_thread==NULL){
   printf("\nAllocating error for *pb\n");
   printf("\nRequired size is %d sizeof(double)",Dim);
   exit(0);
  } 
//  density matrix  ///////////////////////////////////////////////////
 p=0.5*Dim*(Dim+1);
 rho_sys_thread=(_Complex double*)calloc((NstepDt+1)*p,sizeof(_Complex double));
 if(rho_sys_thread==NULL){
  printf("\nAllocating error for *rho_sys_thread(t)\n");
  printf("\nRequired size is %d sizeof(double)",(NstepDt+1)*p);
  exit(0);
 } 
 rho_sys_thread_pointer=(_Complex double***)calloc(NstepDt+1,sizeof(_Complex double**));
 for(i=0;i<=NstepDt;i++){
  rho_sys_thread_pointer[i]=(_Complex double**)calloc(Dim,sizeof(_Complex double*));
  if(rho_sys_thread_pointer==NULL){
   printf("\nAllocating error for *rho_sys_thread_pointer(t)\n");
   printf("\nRequired size is %d sizeof(double)",p);
   exit(0);
  }
 }
 for(i=0;i<=NstepDt;i++){
  for(j=0;j<Dim;j++){
   m=0; for(k=0;k<=j;k++)m+=k;
   rho_sys_thread_pointer[i][j]=&rho_sys_thread[i*p+m];
  }
 }
//  sink population  //////////////////////////////////////////////////
 population_sink_thread=(double*)calloc(NstepDt+1,sizeof(double));
 if(population_sink_thread==NULL){
  printf("\nAllocating error for population_sink_thread\n");
  printf("\nRequired size is %d sizeof(double)",NstepDt+1);
  exit(0);
 }
//  standard deviation  ///////////////////////////////////////////////
 if(calculate_sd){
  standevm=(double**)calloc(NstepDt+1,sizeof(double));
  for(i=0;i<NstepDt+1;i++){
   standevm[i]=(double*)calloc(Dim,sizeof(double));
   if(standevm==NULL){
    printf("\nAllocating error for standevm\n");
    printf("\nRequired size is %d sizeof(double)",Dim);
    exit(0);
   }
  }
  standevs=(double**)calloc(NstepDt+1,sizeof(double));
  for(i=0;i<NstepDt+1;i++){
   standevs[i]=(double*)calloc(Dim,sizeof(double));
   if(standevs==NULL){
    printf("\nAllocating error for standevs\n");
    printf("\nRequired size is %d sizeof(double)",Dim);
    exit(0);
   }
  }
  for(i=0;i<Dim;i++)standevm[0][i]=creal(x0[i]*~x0[i]);
 }
//  filtered noise arrays  ////////////////////////////////////////////
 if(add_filtered_noise_dephasing){
  filtered_noise_memory=(double**)calloc(filtered_noise_processes,sizeof(double));
  for(j=0;j<filtered_noise_processes;j++){
   filtered_noise_memory[j]=(double*)calloc(filtered_noise_order[filtered_noise_site[j]],sizeof(double));
   if(filtered_noise_memory==NULL){
    printf("\nAllocating error for filtered_noise_freqres\n");
    printf("\nRequired size is %d sizeof(double)",filtered_noise_order[filtered_noise_site[j]]);
    exit(0);
   }
  }
  for(j=0;j<filtered_noise_processes;j++)for(i=0;i<filtered_noise_order[filtered_noise_site[j]];i++){
   filtered_noise_memory[j][i]=getgaussrandom()*filtered_noise_width[filtered_noise_site[j]];
  }
  filtered_noise=(double*)calloc(filtered_noise_processes,sizeof(double));
  if(filtered_noise==NULL){
   printf("\nAllocating error for filtered_noise\n");
   printf("\nRequired size is %d sizeof(double)",filtered_noise_processes);
   exit(0);
  }
 }
//  ou arrays  ////////////////////////////////////////////////////////
 if(add_ou_process_dephasing){
  save_ou_number=(double*)calloc(Dim,sizeof(double));
  if(save_ou_number==NULL){
   printf("\nAllocating error for save_ou_number\n");
   printf("\nRequired size is %d sizeof(double)",Dim);
   exit(0);
  }
  ou_const=(double*)calloc(Dim,sizeof(double));
  if(ou_const==NULL){
   printf("\nAllocating error for ou_const\n");
   printf("\nRequired size is %d sizeof(double)",Dim);
   exit(0);
  }
  ou_corr=(double*)calloc(Dim,sizeof(double));
  if(ou_corr==NULL){
   printf("\nAllocating error for ou_corr\n");
   printf("\nRequired size is %d sizeof(double)",Dim);
   exit(0);
  }
  for(i=0;i<Dim;i++)ou_corr[i]=exp(-h/oup_tau[i]);
  for(i=0;i<Dim;i++)ou_const[i]=sqrt(1-ou_corr[i]*ou_corr[i])*oup_sigma[i];
 }
///////////////////////////////////////////////////////////////////////
//////  start the trajectories  ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////
  for(k=1;k<=Ntrajectories; k++){//-----------------------------//For1


///////////////////////////////////////////////////////////////////////
//----  nice noises  ----//////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
   if(add_static_noise_dephasing){
    if(add_sitecorr_sn){
     for(i=0;i<Dim;i++)for(j=0;j<Dim;j++)MatrixH[i][j]=MatrixA[i][j];
     for(p=0;p<Dim;p++)MatrixH[p][p]=(snoise[p][p]>zero||snoise[p][p]<-zero)?getgaussrandom():0.0;
     for(p=0;p<Dim;p++){
      if(snoise[p][p]>zero||snoise[p][p]<-zero){
       if(site_anticorr_rate[p]<zero){
        MatrixH[p][p]=MatrixH[site_corr[p]][site_corr[p]];         
       }else if(site_corr_rate[p]<zero){
        MatrixH[p][p]=MatrixH[p][p];
       }else{
        MatrixH[p][p]=site_anticorr_rate[p]*MatrixH[p][p]+site_corr_rate[p]*MatrixH[site_corr[p]][site_corr[p]];
       }
      }
     }
     for(p=0;p<Dim;p++)MatrixH[p][p]=(snoise[p][p]>zero||snoise[p][p]<-zero)?MatrixH[p][p]*snoise[p][p]:0.0;
     for(p=0;p<Dim;p++)MatrixH[p][p]+=MatrixA[p][p];
    }else{
     for(i=0;i<Dim;i++)for(j=0;j<Dim;j++)MatrixH[i][j]=(snoise[i][j]>zero||snoise[i][j]<-zero)?MatrixA[i][j]+getgaussrandom()*snoise[i][j]:MatrixA[i][j];
    }
   }else{
    for(i=0;i<Dim;i++)for(j=0;j<Dim;j++)MatrixH[i][j]=MatrixA[i][j];
   }
   if(add_ou_process_dephasing)for(p=0;p<Dim;p++)if(oup_sigma[p]>0.0)save_ou_number[p]=oup_sigma[p]*getgaussrandom();


   for(i=0;i<Dim;i++)xp_thread[i]=x0[i];//Refresh initial conditions
   sink_thread=0;//sink population, refreshed each trajectory

   if(use_qa)for(i=0;i<Dim;i++)for(j=0;j<Dim;j++)MatrixHqa[i][j]=MatrixH[i][j];

   for(m=1; m<=NstepDt; m++){//----------------------------------//For2
    if(!sink_thread){//-------------------------------------//if sink 1
     for(n=0; n<Nsteph; n++){//----------------------------------//For3


///////////////////////////////////////////////////////////////////////
//----  nasty noises  ----/////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
// ou process /////////////////////////////////////////////////////////
      if(add_ou_process_dephasing){//-------------------------//if ou 1
       for(p=0;p<Dim;p++)if(oup_sigma[p]>0.0)save_ou_number[p]=ou_corr[p]*save_ou_number[p]+ou_const[p]*getgaussrandom();
       if(add_sitecorr_ou){
        for(p=0;p<Dim;p++){
         if(site_anticorr_rate[p]<zero){
          if(oup_sigma[site_corr[p]]>0.0)MatrixH[p][p]=MatrixA[p][p]+save_ou_number[site_corr[p]];         
         }else if(site_corr_rate[p]<zero){
          if(oup_sigma[p]>0.0)MatrixH[p][p]=MatrixA[p][p]+save_ou_number[p];
         }else{
          if(oup_sigma[site_corr[p]]>0.0||oup_sigma[p]>0.0)MatrixH[p][p]=MatrixA[p][p]+site_anticorr_rate[p]*save_ou_number[p]+site_corr_rate[p]*save_ou_number[site_corr[p]];
         }
        }
       }else{
        for(p=0;p<Dim;p++)if(oup_sigma[p]>0.0)MatrixH[p][p]=MatrixA[p][p]+save_ou_number[p];
       }
      }//-----------------------------------------------------//if ou 1
// filtered noise /////////////////////////////////////////////////////
      if(add_filtered_noise_dephasing){//-------------------------//if fn
       for(p=0;p<filtered_noise_processes;p++){
        for(j=filtered_noise_order[filtered_noise_site[p]]-1;j>0;j--)filtered_noise_memory[p][j]=filtered_noise_memory[p][j-1];
        filtered_noise_memory[p][0]=getgaussrandom()*filtered_noise_width[filtered_noise_site[p]];
       }
       for(p=0;p<filtered_noise_processes;p++)filtered_noise[p]=0;
       for(p=0;p<filtered_noise_processes;p++)for(j=0;j<filtered_noise_order[filtered_noise_site[p]];j++)filtered_noise[p]+=filtered_noise_memory[p][j]*filtered_noise_freqres[p][j];
       for(p=0;p<filtered_noise_processes;p++)MatrixH[filtered_noise_site[p]][filtered_noise_site[p]]=MatrixA[filtered_noise_site[p]][filtered_noise_site[p]]+filtered_noise[p];
      }//-------------------------------------------------------//if fn
// time_point dependent things //////////////////////////////////////////////
      if(add_pl_wave_dephasing||add_fm_signal_dephasing||add_sq_wave_dephasing||use_pulses||use_qa)time_point=(m-1)*Dt+n*h;
// plane wave /////////////////////////////////////////////////////////
      if(add_pl_wave_dephasing){
       for(p=0;p<Dim;p++)if(pwd_amp[p]>0.0)MatrixH[p][p]=pwd_amp[p]*sin((pwd_freq[p]+pwd_phase[p])*time_point)+MatrixA[p][p];
      }
// frequency modulated ////////////////////////////////////////////////
      if(add_fm_signal_dephasing){
       for(p=0;p<Dim;p++)if(fmn_amp[p]>0.0)MatrixH[p][p]=MatrixA[p][p]+fmn_amp[p]*cos(fmn_base[p]*time_point+fmn_dev[p]*cos(fmn_signal[p]*time_point));
      }
// square wave ////////////////////////////////////////////////////////
      if(add_sq_wave_dephasing){//-----------------------------//if sqw
       for(p=0;p<Dim;p++)if(sqw_amp[p]>0.0){
        MatrixH[p][p]=0;
        for(sqw_num=1;sqw_num<sqw_samp[p];sqw_num++){
         sqw_k=2*sqw_num-1;
         MatrixH[p][p]+=sin(sqw_freq[p]*sqw_k*time_point)/sqw_k;
        }
        MatrixH[p][p]*=sqw_amp[p];
        MatrixH[p][p]+=MatrixA[p][p];
       }
      }//------------------------------------------------------//if sqw
// pulses /////////////////////////////////////////////////////////////
      if(use_pulses){//----------------------------------------//if pls
       for(p=0;p<pulses_count;p++){
        pulses_x=(time_point-pulses_now[p])/pulses_wth[p];
        MatrixH[pulses_site_i[p]][pulses_site_j[p]]=MatrixA[pulses_site_i[p]][pulses_site_j[p]]+pulses_amp[p]*exp(-(pulses_x*pulses_x));
        MatrixH[pulses_site_j[p]][pulses_site_i[p]]=MatrixH[pulses_site_i[p]][pulses_site_j[p]];
       }
      }//------------------------------------------------------//if pls
// QA algorithm ///////////////////////////////////////////////////////
      if(use_qa){//---------------------------------------------//if qa
       qa_At=0;
       for(p=0;p<qa_a_n+1;p++)qa_At+=qa_a_bigA[p]*cos(qa_a_a[p]*time_point+qa_a_b[p]*erf((time_point-qa_a_tp2[p])/qa_a_bigT2[p])+qa_a_c[p])*exp(-((time_point-qa_a_tp1[p])/qa_a_bigT1[p])*((time_point-qa_a_tp1[p])/qa_a_bigT1[p]));
       for(p=0;p<Dim;p++)for(j=0;j<Dim;j++)MatrixH[p][j]=0;
       for(p=0;p<Dim;p++)for(j=p;j<Dim;j++)if(qa_MatrixQ[p][j]!=0||MatrixHqa[p][j]!=0)MatrixH[p][j]=MatrixHqa[p][j]+qa_At*qa_MatrixQ[p][j];
       for(p=0;p<Dim;p++)for(j=p+1;j<Dim;j++)MatrixH[j][p]=MatrixH[p][j];
      }//-------------------------------------------------------//if qa
///////////////////////////////////////////////////////////////////////
//////  jumps  ////////////////////////////////////////////////////////
////  modified version of method from  Giuliano, Giulio, and       ////
////  Strini's book                                                ////
///////////////////////////////////////////////////////////////////////
      jump=0;
      if(decay&&!sink_thread){//run jumps
       if(add_sitecorr_wn){
        for(p=0;p<decay_1;p++)pb_thread[p]=decay_rate[p]*h*xp_thread[donor_state[p]]*~xp_thread[donor_state[p]];
        for(p=decay_1;p<decay;p++){
         pb_thread[p]=decay_rate[p]*h*creal(xp_thread[donor_state[p]]*~xp_thread[site_corr[donor_state[p]]]);
        }
       }else{
        for(p=0;p<decay;p++)pb_thread[p]=decay_rate[p]*h*xp_thread[donor_state[p]]*~xp_thread[donor_state[p]];
       }
       for(p=1;p<decay;p++)pb_thread[p]+=pb_thread[p-1];
       rand_num=((double)rand_r(&seed)/RAND_MAX);
///////////////////////////////////////////////////////////////////////
//////  run jump test  ////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
       if(pb_thread[decay-1]>rand_num){//jump if
        jlvl=-1; do{jlvl++;}while(rand_num>pb_thread[jlvl]);
        if(receptive_state[jlvl]==donor_state[jlvl]){// sort jumps(sink dispersion)
         for(p=0;p<Dim;p++)xp_thread[p]=0; 
         count_jumps_intosink_thread++;
         sink_thread=1;
        }else if(receptive_state[jlvl]==Dim){// sort jumps(white noise dephasing)
         for(p=0;p<Dim;p++)xp_thread[p]=0;
         xp_thread[donor_state[jlvl]]=1;
         count_dephase_thread++;
         jump=1;
        }else if(receptive_state[jlvl]>Dim){// sort jumps(white noise dephasing correlation)// not working - please ignore for now
         xp_thread[donor_state[jlvl]]*=I;
         xp_thread[site_corr[donor_state[jlvl]]]*=I;
         count_dephase_thread++;
         jump=1;
        }else{// sort jumps(in-system dispersion)
         for(p=0;p<Dim;p++)xp_thread[p]=0; 
         xp_thread[receptive_state[jlvl]]=1;
         count_jumps_thread++;
         jump=1;
        }// sort jumps
       }//jump if
      }//run jumps


///////////////////////////////////////////////////////////////////////
//////  no jumps  /////////////////////////////////////////////////////
////  Runge-Kutta algorithm (fourth order)                         ////
///////////////////////////////////////////////////////////////////////
      if(!sink_thread&&!jump){//------------------------//if not-sink 3
       for(i=0;i<Dim;i++){
        partial=0;
        for(j=0;j<Dim;j++)if(MatrixH[i][j]!=0)partial-=I*MatrixH[i][j]*xp_thread[j];
        if(vectorgammasum[i]!=0)partial-=vectorgammasum[i]*xp_thread[i]; 
        k1[i]=h*partial;
       }
       for(i=0;i<Dim;i++)auxx[i]=xp_thread[i]+0.5*k1[i];
       for(i=0;i<Dim;i++){
        partial=0;
        for(j=0;j<Dim;j++)if(MatrixH[i][j]!=0)partial-=I*MatrixH[i][j]*auxx[j];
        if(vectorgammasum[i]!=0)partial-=vectorgammasum[i]*auxx[i]; 
        k2[i]=h*partial;
       }
       for(i=0;i<Dim;i++)auxx[i]=xp_thread[i]+0.5*k2[i];
       for(i=0;i<Dim;i++){
        partial=0;
        for(j=0;j<Dim;j++)if(MatrixH[i][j]!=0)partial-=I*MatrixH[i][j]*auxx[j];
        if(vectorgammasum[i]!=0)partial-=vectorgammasum[i]*auxx[i]; 
        k3[i]=h*partial;
       }
       for(i=0;i<Dim;i++)auxx[i]=xp_thread[i]+k3[i];  
       for(i=0;i<Dim;i++){
        partial=0;
        for(j=0;j<Dim;j++)if(MatrixH[i][j]!=0)partial-=I*MatrixH[i][j]*auxx[j];
        if(vectorgammasum[i]!=0)partial-=vectorgammasum[i]*auxx[i]; 
        k4[i]=h*partial;
       }
       for(i=0;i<Dim;i++)xn_thread[i]=xp_thread[i]+(k1[i]+2*(k2[i]+k3[i])+k4[i])/6;
       if(decay||add_white_noise_dephasing){
        aux=0;
        for(i=0;i<Dim;i++) aux+=xn_thread[i]*~xn_thread[i];
        norm=sqrt(aux);
        for(i=0;i<Dim;i++)xn_thread[i]/=norm; 
       }
       for(i=0;i<Dim;i++)xp_thread[i]=xn_thread[i];
      }
     }//---------------------------------------------------//Close For3
    }else{//-----------------------------------------//Chiude if sink 1
     for(i=0;i<Dim;i++){xn_thread[i]=0; xp_thread[i]=xn_thread[i];}
    }
    for(i=0;i<Dim;i++){
     for(j=0;j<=i;j++)*(rho_sys_thread_pointer[m][i]+j)+=xp_thread[i]*~xp_thread[j];
     if(calculate_sd){
      standevs[m][i]+=(((*(rho_sys_thread_pointer[m][i]+i)/k)-standevm[m][i])*((*(rho_sys_thread_pointer[m][i]+i)/k)-standevm[m][i])*k)/(k+1);
      standevm[m][i]+=((*(rho_sys_thread_pointer[m][i]+i)/k)-standevm[m][i])/(k+1);
     }
    }
    population_sink_thread[m]+=sink_thread;
   }//-----------------------------------------------------//Close For2
  }//------------------------------------------------------//Close For1


///////////////////////////////////////////////////////////////////////
//////  collection  ///////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
  for(m=1;m<=NstepDt;m++)for(i=0;i<Dim;i++)for(j=0;j<=i;j++){
   #pragma omp critical(densitymatrix)
   {
    *(rho_sys_pointer[m][i]+j)+=*(rho_sys_thread_pointer[m][i]+j);
   }
  }
  #pragma omp critical(sinkpop)
  {
   for(m=1;m<=NstepDt;m++)population_sink[m]+=population_sink_thread[m];
  }
  #pragma omp critical(systemjumps)
  {
   count_jumps+=count_jumps_thread;
  }
  #pragma omp critical(sinkjumps)
  {
   count_jumps_intosink+=count_jumps_intosink_thread;
  }
  #pragma omp critical(dephasejumps)
  {
   count_dephase+=count_dephase_thread;
  }
  if(calculate_sd){
  #pragma omp critical(standarddev)
   {
    for(m=1;m<=NstepDt;m++)for(i=0;i<Dim;i++)standard_deviation[i][m]+=standevs[m][i]/(Ntrajectories-1);
   }
  }
 }
 printf("sink jump rate %f\n",(double)count_jumps_intosink/((double)Ntrajectories*(double)num_threads));
 printf("level jump rate %f\n",(double)count_jumps/((double)Ntrajectories*(double)num_threads));
 if(add_white_noise_dephasing)printf("dephase jump rate %f\n",(double)count_dephase/((double)Ntrajectories*(double)num_threads));
}
