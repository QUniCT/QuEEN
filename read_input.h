//Dichiarazioni funzioni
int read_input();

int read_input()
{
 int livelli; //number of levels
 int i,j,k,l,err;
 double component,component2,component3,component4;
 double pi=3.141592;
/////////////////////////////////////////////////////////////////////// 
//////  read system dimensions  ///////////////////////////////////////
////  Hilbert space contains the excited states and a ground       ////
////  state. All input files must start with the number of         ////
////  excited states                                               ////
/////////////////////////////////////////////////////////////////////// 
 err=fscanf(hamiltonian_file,"%d",&livelli);
 if(err!=1){printf("problem reading hamiltonian file\n");exit(0);}
 err=fscanf(dissipator_file,"%d",&i);
 if(err!=1){printf("problem reading dissipator file\n");exit(0);}
 if(i!=livelli){
  printf("\nThe hamiltonian and the dissipator have disparate dimensions\n N(d)=%d and N(h)=%d\n",i,livelli);
  exit(0);
 }
 err=fscanf(initcon_file,"%d",&i);
 if(err!=1){printf("problem reading initcon file\n");exit(0);}
 if(i!=livelli){
  printf("\nThe hamiltonian and the initial conditions have disparate dimensions\n N(ic)=%d and N(h)=%d\n",i,livelli);
  exit(0);
 }
 if(add_white_noise_dephasing){
  err=fscanf(wnd_file,"%d",&i);
  if(err!=1){printf("problem reading .wnd file\n");exit(0);}
  if(i!=livelli){
   printf("\nThe hamiltonian and the wnd have disparate dimensions\n N(a)=%d and N(h)=%d\n",i,livelli);
   exit(0);
  }
 }
 if(add_ou_process_dephasing){
  err=fscanf(oud_file,"%d",&i);
  if(err!=1){printf("problem reading .oup file\n");exit(0);}
  if(i!=livelli){
   printf("\nThe hamiltonian and the oup have disparate dimensions\n N(a)=%d and N(h)=%d\n",i,livelli);
   exit(0);
  }
 }
 if(add_filtered_noise_dephasing){
  err=fscanf(fnd_file,"%d",&i);
  if(err!=1){printf("problem reading .fnd file\n");exit(0);}
  if(i!=livelli){
   printf("\nThe hamiltonian and the fnd have disparate dimensions\n N(a)=%d and N(h)=%d\n",i,livelli);
   exit(0);
  }
 }
 if(add_pl_wave_dephasing){
  err=fscanf(pwd_file,"%d",&i);
  if(err!=1){printf("problem reading .pwd file\n");exit(0);}
  if(i!=livelli){
   printf("\nThe hamiltonian and the pwd have disparate dimensions\n N(a)=%d and N(h)=%d\n",i,livelli);
   exit(0);
  }
 }
 if(add_fm_signal_dephasing){
  err=fscanf(fmd_file,"%d",&i);
  if(err!=1){printf("problem reading .fmn file\n");exit(0);}
  if(i!=livelli){
   printf("\nThe hamiltonian and the fmn have disparate dimensions\n N(a)=%d and N(h)=%d\n",i,livelli);
   exit(0);
  }
 }
 if(add_sq_wave_dephasing){
  err=fscanf(sqd_file,"%d",&i);
  if(err!=1){printf("problem reading .sqw file\n");exit(0);}
  if(i!=livelli){
   printf("\nThe hamiltonian and the sqw have disparate dimensions\n N(a)=%d and N(h)=%d\n",i,livelli);
   exit(0);
  }
 }
 if(add_sitecorr_wn||add_sitecorr_sn||add_sitecorr_ou){
  err=fscanf(sitecorr_file,"%d",&i);
  if(err!=1){printf("problem reading sitecorr file\n");exit(0);}
  if(i!=livelli){
   printf("\nThe hamiltonian and the sitecorr have disparate dimensions\n N(a)=%d and N(h)=%d\n",i,livelli);
   exit(0);
  }
 }
 if(use_pulses){
  err=fscanf(pulses_file,"%d",&i);
  if(err!=1){printf("problem reading sitecorr file\n");exit(0);}
  if(i!=livelli){
   printf("\nThe hamiltonian and the pulses have disparate dimensions\n N(a)=%d and N(h)=%d\n",i,livelli);
   exit(0);
  }
 }
 if(use_qa){
  err=fscanf(qa_q_file,"%d",&i);
  if(err!=1){printf("problem reading sitecorr file\n");exit(0);}
  if(i!=livelli){
   printf("\nThe hamiltonian and the Q Matrix have disparate dimensions\n N(a)=%d and N(h)=%d\n",i,livelli);
   exit(0);
  }
 }
 if(add_static_noise_dephasing){
  err=fscanf(snd_file,"%d",&i);
  if(err!=1){printf("problem reading snoise file\n");exit(0);}
  if(i!=livelli){
   printf("\nThe hamiltonian and the snoise have disparate dimensions\n N(a)=%d and N(h)=%d\n",i,livelli);
   exit(0);
  }
 }
 Dim=livelli+1; //Hilbert space dimension (doesn't include the sink)
 jump_dim=Dim*Dim; //Maximum number of jump operators
 printf("\n Excited states=%d; --- All states=%d\n",livelli,Dim);
///////////////////////////////////////////////////////////////////////
//////  read hamiltonian  /////////////////////////////////////////////
////  fill the upper triangle of the matrix (i<j), and the code    ////
////  will symmeterize the hamiltonian automatically               ////
////  All components are zero by default                           ////
///////////////////////////////////////////////////////////////////////
////  File format:                                                 ////
////  i j H{i,j}                                                   ////
///////////////////////////////////////////////////////////////////////
 genera_MatrixA();
 while(fscanf(hamiltonian_file,"%d %d %lf",&i,&j,&component)!=EOF){
  if(i<Dim&&j<Dim){
   MatrixA[i][j]=component;
  }else{
   printf("the hamiltonian index has exceeded the dimensions\n Dimension=%d  i=%d  j=%d",livelli,i,j);
   exit(0);
  }
 }
 fclose(hamiltonian_file);
 simmetrizza_MatrixA();
 if(Dim<25){
  printf("\nhamiltonian:");
  for(i=0;i<Dim;i++){
   printf("\n");
   for(j=0;j<Dim;j++)printf("%.2f ",MatrixA[i][j]);
  }
 }
///////////////////////////////////////////////////////////////////////
//////  read dissipator  //////////////////////////////////////////////
////  fill the full matrix with relaxation rates for jumps from    ////
////  site j to site i. If i=j, site j jumps into a sink - which   ////
////  is not part of the hamiltonian                               ////
////  All rates are zero by default                                ////
///////////////////////////////////////////////////////////////////////
////  File format:                                                 ////
////  i j Gamma{i,j}                                               ////
///////////////////////////////////////////////////////////////////////
 genera_VectorGamma();
 while(fscanf(dissipator_file,"%d %d %lf",&i,&j,&component)!=EOF){
  if(i<Dim&&j<Dim){
   *(vectorgamma_pointer[i]+j)=component;
  }else{
   printf("\n\nError: the dissipator index has exceeded the dimensions\n Dimension=%d  i=%d  j=%d",livelli,i,j);
   exit(0);
  }
 }
 fclose(dissipator_file);
 if(Dim<25){
  printf("\n\ndissipator:");
  for(i=0;i<Dim;i++){
   printf("\n");
   for(j=0;j<Dim;j++)printf("%.2f ",*(vectorgamma_pointer[i]+j));
  }
 }
///////////////////////////////////////////////////////////////////////
//////  read initial conditions  //////////////////////////////////////
////  fill the initial values of the wavefunction coeffiecients    ////
////  if the command line flag -nic is used, these values are      ////
////  normalized automatically by the code                         ////
////  All coeffiecients are zero by default                        ////
///////////////////////////////////////////////////////////////////////
////  File format:                                                 ////
////  j C{j}                                                       ////
///////////////////////////////////////////////////////////////////////
 genera_initcon();
 while(fscanf(initcon_file,"%d %lf",&j,&component)!=EOF){
  if(j<Dim){
   x0[j]=component;
  }else{
   printf("\n\nError: the initial condition index has exceeded the dimensions\n Dimension=%d  j=%d",livelli,j);
   exit(0);
  }
 }
 fclose(initcon_file);
 component=0; for(j=0;j<Dim;j++)component+=creal(x0[j])*creal(x0[j]);
 if(component<0.00001){
  printf("\n\nError: the initial state is empty, sum(rho[j,j])=%f",component);
  exit(0);
 }
 component=1/sqrt(component);
 for(j=0;j<Dim;j++)x0[j]=component*x0[j];
 if(Dim<25){
  printf("\n\ninitial state:");
  for(j=0;j<Dim;j++)printf("\nsite %d, %.2f ",j,creal(x0[j]));
 }
///////////////////////////////////////////////////////////////////////
//////  maybe read wnd  ///////////////////////////////////////////////
////  fill a vector with pure dephasing parameter (gamma) for      ////
////  each level. It is only read if the command line flag         ////
////  -dephase is used                                             ////
////  All components are zero by default                           ////
///////////////////////////////////////////////////////////////////////
////  File format:                                                 ////
////  j gamma{j}                                                   ////
///////////////////////////////////////////////////////////////////////
 if(add_white_noise_dephasing){
  dph_gamma=(double *)calloc(Dim,sizeof(double));
   if(dph_gamma==NULL){
    printf("\nAllocating error for dph_gamma\n");
    printf("\nRequired size is %d sizeof(double)\n",Dim);
    exit(0);
   } 
  while(fscanf(wnd_file,"%d %lf",&j,&component)!=EOF){
   if(j<Dim){
    dph_gamma[j]=component;
   }else{
    printf("\n\nError: the dephase index has exceeded the dimensions\n Dimension=%d j=%d",livelli,j);
    exit(0);
   }
  }
 }
///////////////////////////////////////////////////////////////////////
//////  maybe read oup  ///////////////////////////////////////////////
////  fill a vector with the ou parameters - noise width, sigma    ////
////  and correlation time, tau. Only called in the flag -ou is    ////
////  used                                                         ////
////  All components are zero by default                           ////
///////////////////////////////////////////////////////////////////////
////  File format:                                                 ////
////  j sigma{j} tau{j}                                            ////
///////////////////////////////////////////////////////////////////////
 if(add_ou_process_dephasing){
  oup_sigma=(double *)calloc(Dim,sizeof(double));
   if(oup_sigma==NULL){
    printf("\nAllocating error for dephase\n");
    printf("\nRequired size is %d sizeof(double)\n",Dim);
    exit(0);
   } 
  oup_tau=(double *)calloc(Dim,sizeof(double));
   if(oup_tau==NULL){
    printf("\nAllocating error for dph_gamma\n");
    printf("\nRequired size is %d sizeof(double)\n",Dim);
    exit(0);
   } 
  while(fscanf(oud_file,"%d %lf %lf",&j,&component,&component2)!=EOF){
   if(j<Dim){
    oup_sigma[j]=sqrt(component);
    oup_tau[j]=component2;
   }else{
    printf("\n\nError: the dephase index has exceeded the dimensions\n Dimension=%d j=%d",livelli,j);
    exit(0);
   }
  }
 }
///////////////////////////////////////////////////////////////////////
//////  maybe read fnd  ///////////////////////////////////////////////
////  fill a vector with filtered noise parameters (order, width,  ////
////  and frequncy boundaries) for each level. It is only read if  ////
////  the command line flag -fn is used                            ////
////  All components are zero by default                           ////
///////////////////////////////////////////////////////////////////////
////  File format:                                                 ////
////  j order{j} frequency_low{j} frequency_high{j}                ////
///////////////////////////////////////////////////////////////////////
 if(add_filtered_noise_dephasing){
  filtered_noise_width=(double *)calloc(Dim,sizeof(double));
   if(filtered_noise_width==NULL){
    printf("\nAllocating error for filtered_noise_width\n");
    printf("\nRequired size is %d sizeof(double)\n",Dim);
    exit(0);
   } 
  filtered_noise_order=(int *)calloc(Dim,sizeof(int));
   if(filtered_noise_order==NULL){
    printf("\nAllocating error for filtered_noise_order\n");
    printf("\nRequired size is %d sizeof(int)",Dim);
    exit(0);
   } 
  filtered_noise_lfreq=(double *)calloc(Dim,sizeof(double));
   if(filtered_noise_lfreq==NULL){
    printf("\nAllocating error for filtered_noise_lfreq\n");
    printf("\nRequired size is %d sizeof(double)\n",Dim);
    exit(0);
   } 
  filtered_noise_hfreq=(double *)calloc(Dim,sizeof(double));
   if(filtered_noise_hfreq==NULL){
    printf("\nAllocating error for filtered_noise_lfreq\n");
    printf("\nRequired size is %d sizeof(double)\n",Dim);
    exit(0);
   }
  while(fscanf(fnd_file,"%d %lf %d %lf %lf",&j,&component,&i,&component2,&component3)!=EOF){
   if(j<Dim){
    filtered_noise_width[j]=component;
    filtered_noise_order[j]=i;
    filtered_noise_lfreq[j]=component2*h;
    filtered_noise_hfreq[j]=component3*h;
   }else{
    printf("\n\nError: the dephase index has exceeded the dimensions\n Dimension=%d j=%d",livelli,j);
    exit(0);
   }
  }
///////////////////////////////////////////////////////////////////////
////  now we set up some other arrays used by the filter  /////////////
///////////////////////////////////////////////////////////////////////
  filtered_noise_processes=0;
  for(j=0;j<Dim;j++)if(filtered_noise_order[j]>0)filtered_noise_processes+=1;
  if(filtered_noise_processes<1){
   printf("\nError: If you flag filtered noise, you need to add filter process into the file *.fnd");
   exit(0);
  }
  filtered_noise_site=(int *)calloc(filtered_noise_processes,sizeof(int));
  if(filtered_noise_site==NULL){
   printf("\nAllocating error for filtered_noise_site\n");
   printf("\nRequired size is %d sizeof(int)",Dim);
   exit(0);
  }
  i=0;
  for(j=0;j<Dim;j++)if(filtered_noise_order[j]>0){
   filtered_noise_site[i]=j;
   i+=1;
  }
  filtered_noise_freqres=(double**)calloc(filtered_noise_processes,sizeof(double));
  for(j=0;j<filtered_noise_processes;j++){
   filtered_noise_freqres[j]=(double*)calloc(filtered_noise_order[filtered_noise_site[j]],sizeof(double));
   if(filtered_noise_freqres==NULL){
    printf("\nAllocating error for filtered_noise_freqres\n");
    printf("\nRequired size is %d sizeof(double)\n",filtered_noise_order[filtered_noise_site[j]]);
    exit(0);
   }
  }
///////////////////////////////////////////////////////////////////////
////  set up the frequency responses  /////////////////////////////////
///////////////////////////////////////////////////////////////////////
  for(j=0;j<filtered_noise_processes;j++){//-------------------openfor1
   if(filtered_noise_hfreq[filtered_noise_site[j]]<0.00001){//--openif1
    if(filtered_noise_lfreq[filtered_noise_site[j]]<0.00001){//-openif2
//-- no filter ------------------------------------------------------//
     filtered_noise_freqres[j][0]=1;
    }else{//----------------------------------------------------elseif2
//-- lowpass --------------------------------------------------------//
     for(i=0;i<filtered_noise_order[filtered_noise_site[j]];i++){//for2
      k=(i-filtered_noise_order[filtered_noise_site[j]]*0.5)*pi; 
      if(abs(k)>0.00001){//-------------------------------------openif3
       filtered_noise_freqres[j][i]=sin(2*k*filtered_noise_lfreq[filtered_noise_site[j]])/k;
      }else{//--------------------------------------------------elseif3
       filtered_noise_freqres[j][i]=2*k*filtered_noise_lfreq[filtered_noise_site[j]];
      }//------------------------------------------------------closeif3
     }//------------------------------------------------------closefor2
    }//--------------------------------------------------------closeif2
   }else{//-----------------------------------------------------elseif1
    if(filtered_noise_lfreq[filtered_noise_site[j]]<0.00001){//-openif2
//-- highpass -------------------------------------------------------//
     for(i=0;i<filtered_noise_order[filtered_noise_site[j]];i++){//for2
      k=(i-filtered_noise_order[filtered_noise_site[j]]*0.5)*pi; 
      if(abs(k)>0.00001){//-------------------------------------openif3
       filtered_noise_freqres[j][i]=(-sin(2*k*filtered_noise_hfreq[filtered_noise_site[j]]))/k;
//       filtered_noise_freqres[j][i]=(sin(k)-sin(2*k*filtered_noise_hfreq[filtered_noise_site[j]]))/k;
      }else{//--------------------------------------------------elseif3
       filtered_noise_freqres[j][i]=1-2*k*filtered_noise_hfreq[filtered_noise_site[j]];
      }//------------------------------------------------------closeif3
     }//------------------------------------------------------closefor2
    }else{//----------------------------------------------------elseif2
     if(filtered_noise_lfreq[filtered_noise_site[j]]>filtered_noise_hfreq[filtered_noise_site[j]]){
//--------------------------------------------------------------openif3
//-- bandpass -------------------------------------------------------//
      for(i=0;i<filtered_noise_order[filtered_noise_site[j]];i++){//for2
       k=(i-filtered_noise_order[filtered_noise_site[j]]*0.5)*pi;
       if(abs(k)>0.00001){//------------------------------------openif4
        filtered_noise_freqres[j][i]=(sin(2*k*filtered_noise_hfreq[filtered_noise_site[j]])-sin(2*k*filtered_noise_lfreq[filtered_noise_site[j]]))/k;
       }else{//-------------------------------------------------elseif4
        filtered_noise_freqres[j][i]=2*(filtered_noise_hfreq[filtered_noise_site[j]]-filtered_noise_lfreq[filtered_noise_site[j]]);
       }//-----------------------------------------------------closeif4
      }//-----------------------------------------------------closefor2
     }else{//---------------------------------------------------elseif3
//-- bandstop -------------------------------------------------------//
      for(i=0;i<filtered_noise_order[filtered_noise_site[j]];i++){//for2
       k=(i-filtered_noise_order[filtered_noise_site[j]]*0.5)*pi;
       if(abs(k)>0.00001){//------------------------------------openif4
//        filtered_noise_freqres[j][i]=(sin(k)-sin(2*k*filtered_noise_hfreq[filtered_noise_site[j]])+sin(2*k*filtered_noise_lfreq[filtered_noise_site[j]]))/k;
        filtered_noise_freqres[j][i]=(-sin(2*k*filtered_noise_hfreq[filtered_noise_site[j]])+sin(2*k*filtered_noise_lfreq[filtered_noise_site[j]]))/k;
       }else{//-------------------------------------------------elseif4
        filtered_noise_freqres[j][i]=1-2*(filtered_noise_hfreq[filtered_noise_site[j]]-filtered_noise_lfreq[filtered_noise_site[j]]);
       }//-----------------------------------------------------closeif4
      }//-----------------------------------------------------closefor2
     }//-------------------------------------------------------closeif3
    }//--------------------------------------------------------closeif2
   }//---------------------------------------------------------closeif1
  }//---------------------------------------------------------closefor1
 }
///////////////////////////////////////////////////////////////////////
//////  maybe read pwd  ///////////////////////////////////////////////
////  fill a vector with the plane-wave parameters - amplitude,    ////
////  frequency, and phase - A, f, and p respectively - for each   ////
////  site, j                                                      ////
////  All components are zero by default                           ////
///////////////////////////////////////////////////////////////////////
////  File format:                                                 ////
////  j A{j} f{j} p{j}                                             ////
///////////////////////////////////////////////////////////////////////
 if(add_pl_wave_dephasing){
  pwd_amp=(double *)calloc(Dim,sizeof(double));
   if(pwd_amp==NULL){
    printf("\nAllocating error for pwd_amp\n");
    printf("\nRequired size is %d sizeof(double)\n",Dim);
    exit(0);
   } 
  pwd_freq=(double *)calloc(Dim,sizeof(double));
   if(pwd_freq==NULL){
    printf("\nAllocating error for pwd_freq\n");
    printf("\nRequired size is %d sizeof(double)\n",Dim);
    exit(0);
   } 
  pwd_phase=(double *)calloc(Dim,sizeof(double));
   if(pwd_phase==NULL){
    printf("\nAllocating error for pwd_phase\n");
    printf("\nRequired size is %d sizeof(int)",Dim);
    exit(0);
   } 
  while(fscanf(pwd_file,"%d %lf %lf %lf",&j,&component,&component2,&component3)!=EOF){
   if(j<Dim){
    pwd_amp[j]=4*component/pi;
    pwd_freq[j]=2*pi*component2;
    pwd_phase[j]=2*pi*component3;
   }else{
    printf("\n\nError: the pwd index has exceeded the dimensions\n Dimension=%d j=%d",livelli,j);
    exit(0);
   }
  }
 }
///////////////////////////////////////////////////////////////////////
//////  maybe read fmn  ///////////////////////////////////////////////
////  fill a vector with the fm parameters - amplitude, base       ////
////  frequency, signal frequency, and peak deviation:             ////
////          A, fc, fm, and fd respectively                       ////
////  The process is generated using:                              ////
////          x(t)=A*cos(2*pi*fc*t+(fd/fm)*cos(2*pi*fm*t))         ////
////  for each site, j                                             ////
////  All components are zero by default                           ////
///////////////////////////////////////////////////////////////////////
////  File format:                                                 ////
////  j A{j} fc{j} fm{j} fd{d}                                     ////
///////////////////////////////////////////////////////////////////////
 if(add_fm_signal_dephasing){
  fmn_amp=(double *)calloc(Dim,sizeof(double));
   if(fmn_amp==NULL){
    printf("\nAllocating error for fmn_amp\n");
    printf("\nRequired size is %d sizeof(double)\n",Dim);
    exit(0);
   } 
  fmn_base=(double *)calloc(Dim,sizeof(double));
   if(fmn_base==NULL){
    printf("\nAllocating error for fmn_base\n");
    printf("\nRequired size is %d sizeof(double)\n",Dim);
    exit(0);
   } 
  fmn_signal=(double *)calloc(Dim,sizeof(double));
   if(fmn_signal==NULL){
    printf("\nAllocating error for fmn_signal\n");
    printf("\nRequired size is %d sizeof(double)\n",Dim);
    exit(0);
   } 
  fmn_dev=(double *)calloc(Dim,sizeof(double));
   if(fmn_dev==NULL){
    printf("\nAllocating error for fmn_dev\n");
    printf("\nRequired size is %d sizeof(double)\n",Dim);
    exit(0);
   } 
  while(fscanf(fmd_file,"%d %lf %lf %lf %lf",&j,&component,&component2,&component3,&component4)!=EOF){
   if(j<Dim){
    fmn_amp[j]=component;
    fmn_base[j]=2*pi*component2;
    fmn_signal[j]=2*pi*component3;
    fmn_dev[j]=component4/component3;
   }else{
    printf("\n\nError: the fmn index has exceeded the dimensions\n Dimension=%d j=%d",livelli,j);
    exit(0);
   }
  }
 }
///////////////////////////////////////////////////////////////////////
//////  maybe read sqw  ///////////////////////////////////////////////
////  fill a vector with the square-wave parameters - amplitude,   ////
////  frequency, and fourier sample - A, f, and s respectively -   ////
////  for each site, j                                             ////
////  All components are zero by default                           ////
///////////////////////////////////////////////////////////////////////
////  File format:                                                 ////
////  j A{j} f{j} s{j}                                             ////
///////////////////////////////////////////////////////////////////////
 if(add_sq_wave_dephasing){
  sqw_amp=(double *)calloc(Dim,sizeof(double));
   if(sqw_amp==NULL){
    printf("\nAllocating error for sqw_amp\n");
    printf("\nRequired size is %d sizeof(double)\n",Dim);
    exit(0);
   } 
  sqw_freq=(double *)calloc(Dim,sizeof(double));
   if(sqw_freq==NULL){
    printf("\nAllocating error for sqw_freq\n");
    printf("\nRequired size is %d sizeof(double)\n",Dim);
    exit(0);
   } 
  sqw_samp=(int *)calloc(Dim,sizeof(int));
   if(sqw_samp==NULL){
    printf("\nAllocating error for sqw_samp\n");
    printf("\nRequired size is %d sizeof(int)",Dim);
    exit(0);
   } 
  while(fscanf(sqd_file,"%d %lf %lf %d",&j,&component,&component2,&i)!=EOF){
   if(j<Dim){
    sqw_amp[j]=4*component/pi;
    sqw_freq[j]=2*pi*component2;
    sqw_samp[j]=i+1;
   }else{
    printf("\n\nError: the sqw index has exceeded the dimensions\n Dimension=%d j=%d",livelli,j);
    exit(0);
   }
  }
 }
///////////////////////////////////////////////////////////////////////
//////  maybe read sitecorr  //////////////////////////////////////////
////  fill the upper triangle of a matrix with correlation         ////
////  parameters for the site noise dephasing: Corr_{i,j}          ////
////  The correlation uses the formula:                            ////
////    x_j=(1-corr_{i,j})*x_j+corr_{i,j}*x_i                      ////
////  All components are zero by default                           ////
///////////////////////////////////////////////////////////////////////
////  file format:                                                 ////
////  i j site_corr[i][j]                                          ////
////  NB: only fill for j>i, everything else will be ignored       ////
////  Also, corr[i][j] must be between 0 and 1                     ////
////  Also, you can only correlate to one other site (i'll fix     ////
////  it later)                                                    ////
///////////////////////////////////////////////////////////////////////
 if(add_sitecorr_wn||add_sitecorr_sn||add_sitecorr_ou){
  site_corr_rate=(double *)calloc(Dim,sizeof(double));
   if(site_corr_rate==NULL){
    printf("\nAllocating error for site_corr_rate\n");
    printf("\nRequired size is %d sizeof(double)\n",Dim);
    exit(0);
   } 
  site_anticorr_rate=(double *)calloc(Dim,sizeof(double));
   if(site_anticorr_rate==NULL){
    printf("\nAllocating error for site_anticorr_rate\n");
    printf("\nRequired size is %d sizeof(double)\n",Dim);
    exit(0);
   }
  site_corr=(int *)calloc(Dim,sizeof(int));
   if(site_corr==NULL){
    printf("\nAllocating error for site_corr\n");
    printf("\nRequired size is %d sizeof(int)",Dim);
    exit(0);
   } 
  while(fscanf(sitecorr_file,"%d %d %lf",&i,&j,&component)!=EOF){
   if(j>i){
    if(component<0||component>1){
     printf("\n\nError: your site correlation parameter isn't between 0 and 1. My algorithm does wierd and irrelevant things outside that range. Try again.\n\n");
     exit(0);
    }
    site_corr[j]=i;
    site_corr_rate[j]=component;
   }else{
    printf("\n\nError: you should only fill the upper triangle (j>i) for the sitecorr file");
    exit(0);
   }
  }  
 for(j=0;j<Dim;j++)site_anticorr_rate[j]=1-site_corr_rate[j];
 }
///////////////////////////////////////////////////////////////////////
//////  maybe read staticnoise  ///////////////////////////////////////
////  fill a matrix with static gaussian noise width (sigma)       ////
////  for each hamiltonian components                              ////
////  It is only read if the command line flag -snoise is used     ////
////  All components are zero by default                           ////
///////////////////////////////////////////////////////////////////////
////  File format:                                                 ////
////  i j noise[i][j]                                              ////
////  NB: this is necessarily symetric, so only fill for j>i       ////
///////////////////////////////////////////////////////////////////////
 if(add_static_noise_dephasing){
  snoise=(double**)calloc(Dim,sizeof(double));
  if(snoise==NULL){
   printf("\nAllocating error for *MatrixA\n");
   printf("\nRequired size is %d sizeof(double)\n",Dim);
   exit(0);
  } 
  for(i=0;i<Dim;i++){
   *(snoise+i)=(double *) calloc(Dim, sizeof(double));
   if(*(snoise+i)==NULL){
    printf("\nAllocating error for *(snoise+i)\n");
    printf("\nRequired size is %d sizeof(double)\n",Dim);
    exit(0);
   }
  } 
  while(fscanf(snd_file,"%d %d %lf",&i,&j,&component)!=EOF){
   if(j>=i){
    if(component<0){
     snoise[i][j]=-sqrt(-component);
    }else{
     snoise[i][j]=sqrt(component);
    }
    if(j>i)snoise[j][i]=snoise[i][j];
   }else{
    printf("\n\nError: you should only fill the upper triangle (j>i) for the snoise file");
    exit(0);
   }
  }
 }
///////////////////////////////////////////////////////////////////////
//////  maybe read pulses  ////////////////////////////////////////////
////  fill the upper triangle of a matrix with pulse variables     ////
////  All components are zero by default                           ////
///////////////////////////////////////////////////////////////////////
////  file format:                                                 ////
////  i j W_{ij} tau_{ij} T_{ij}                                   ////
////  NB: only fill for j>i, everything else will be ignored       ////
////  Also: time is alway positive, so don't expect a negative tau ////
////        to produce anything but gibberish                      ////
////        T and W can be negative.                               ////
///////////////////////////////////////////////////////////////////////
 if(use_pulses){
  pulses_omega=(double**)calloc(Dim,sizeof(double));
  if(pulses_omega==NULL){
   printf("\nAllocating error for *pulses_omega\n");
   printf("\nRequired size is %d sizeof(double)\n",Dim);
   exit(0);
  }
  for(i=0;i<Dim;i++){
   *(pulses_omega+i)=(double*)calloc(Dim,sizeof(double));
   if(*(pulses_omega+i)==NULL){
    printf("\nAllocating error for *(pulses_omega+i)\n");
    printf("\nRequired size is %d sizeof(double)\n",Dim);
    exit(0);
   }
  }
  pulses_tau=(double**)calloc(Dim,sizeof(double));
  if(pulses_tau==NULL){
   printf("\nAllocating error for *pulses_tau\n");
   printf("\nRequired size is %d sizeof(double)\n",Dim);
   exit(0);
  }
  for(i=0;i<Dim;i++){
   *(pulses_tau+i)=(double*)calloc(Dim,sizeof(double));
   if(*(pulses_tau+i)==NULL){
    printf("\nAllocating error for *(pulses_tau+i)\n");
    printf("\nRequired size is %d sizeof(double)\n",Dim);
    exit(0);
   }
  }
  pulses_width=(double**)calloc(Dim,sizeof(double));
  if(pulses_width==NULL){
   printf("\nAllocating error for *pulses_width\n");
   printf("\nRequired size is %d sizeof(double)\n",Dim);
   exit(0);
  }
  for(i=0;i<Dim;i++){
   *(pulses_width+i)=(double*)calloc(Dim,sizeof(double));
   if(*(pulses_width+i)==NULL){
    printf("\nAllocating error for *(pulses_width+i)\n");
    printf("\nRequired size is %d sizeof(double)\n",Dim);
    exit(0);
   }
  }
  while(fscanf(pulses_file,"%d %d %lf %lf %lf",&i,&j,&component,&component2,&component3)!=EOF){
   if(j<Dim){
    if(j>i){
     pulses_omega[i][j]=component;
     pulses_tau[i][j]=component2;
     pulses_width[i][j]=component3;
    }else{
     printf("\n\nError: you should only fill the upper triangle (j>i) for the pulses file");
     exit(0);
    }
   }else{
    printf("\n\nError: pulses file has exceded the system dimensions");
    exit(0);
   }
  }
  for(i=0;i<Dim;i++)for(j=0;j<Dim;j++)if(pulses_omega[i][j]!=0&&pulses_width[i][j]!=0)pulses_count++;
  if(pulses_count==0){
   printf("\n\nError: no pulse found");
   exit(0);
  }
  pulses_site_i=(int*)calloc(pulses_count,sizeof(int));
  if(pulses_site_i==NULL){
   printf("\nAllocating error for pulses_site_i\n");
   printf("\nRequired size is %d sizeof(double)\n",pulses_count);
   exit(0);
  } 
  pulses_site_j=(int*)calloc(pulses_count,sizeof(int));
  if(pulses_site_j==NULL){
   printf("\nAllocating error for pulses_site_j\n");
   printf("\nRequired size is %d sizeof(double)\n",pulses_count);
   exit(0);
  } 
  pulses_amp=(double*)calloc(pulses_count,sizeof(double));
  if(pulses_amp==NULL){
   printf("\nAllocating error for pulses_amp\n");
   printf("\nRequired size is %d sizeof(double)\n",pulses_count);
   exit(0);
  } 
  pulses_now=(double*)calloc(pulses_count,sizeof(double));
  if(pulses_now==NULL){
   printf("\nAllocating error for pulses_now\n");
   printf("\nRequired size is %d sizeof(double)\n",pulses_count);
   exit(0);
  } 
  pulses_wth=(double*)calloc(pulses_count,sizeof(double));
  if(pulses_wth==NULL){
   printf("\nAllocating error for pulses_wth\n");
   printf("\nRequired size is %d sizeof(double)\n",pulses_count);
   exit(0);
  } 
  pulses_count=0;
  for(i=0;i<Dim;i++)for(j=0;j<Dim;j++)if(pulses_omega[i][j]!=0&&pulses_width[i][j]!=0){
   pulses_site_i[pulses_count]=i;
   pulses_site_j[pulses_count]=j;
   pulses_amp[pulses_count]=pulses_omega[i][j];
   pulses_now[pulses_count]=pulses_tau[i][j];
   pulses_wth[pulses_count]=pulses_width[i][j];
   pulses_count++;
  }
 }
///////////////////////////////////////////////////////////////////////
//////  maybe read pulses  ////////////////////////////////////////////
////  fill the upper triangle of a matrix with Q matrix variables  ////
////  also define the parameters of the A(t) term                  ////
////  All components are zero by default                           ////
///////////////////////////////////////////////////////////////////////
////  file format (*.qa_q):                                        ////
////  i j Q_{ij}                                                   ////
////  NB: only fill for j>i, everything else will be ignored - the ////
////      Q matrix is expected to be hermitian                     ////
////                                                               ////
////  file format (*.qa_a):                                        ////
////  N_{A}                                                        ////
////  A_{i} tp1_{i} T_{i} a_{i} b_{i} c_{i} tp2_{i}                ////
////  Also: time is alway positive, so don't expect a negative tp* ////
////        to produce anything but gibberish                      ////
///////////////////////////////////////////////////////////////////////
 if(use_qa){
  err=fscanf(qa_a_file,"%d",&qa_a_n); 
  qa_a_n--;
  if(err!=1){printf("problem reading .oup file\n");exit(0);}
  qa_a_bigA=(double*)calloc(qa_a_n,sizeof(double));
  if(qa_a_bigA==NULL){
   printf("\nAllocating error for pulses_wth\n");
   printf("\nRequired size is %d sizeof(double)\n",qa_a_n);
   exit(0);
  } 
  qa_a_tp1=(double*)calloc(qa_a_n,sizeof(double));
  if(qa_a_tp1==NULL){
   printf("\nAllocating error for pulses_wth\n");
   printf("\nRequired size is %d sizeof(double)\n",qa_a_n);
   exit(0);
  } 
  qa_a_bigT1=(double*)calloc(qa_a_n,sizeof(double));
  if(qa_a_bigT1==NULL){
   printf("\nAllocating error for pulses_wth\n");
   printf("\nRequired size is %d sizeof(double)\n",qa_a_n);
   exit(0);
  } 
  qa_a_a=(double*)calloc(qa_a_n,sizeof(double));
  if(qa_a_a==NULL){
   printf("\nAllocating error for pulses_wth\n");
   printf("\nRequired size is %d sizeof(double)\n",qa_a_n);
   exit(0);
  } 
  qa_a_b=(double*)calloc(qa_a_n,sizeof(double));
  if(qa_a_b==NULL){
   printf("\nAllocating error for pulses_wth\n");
   printf("\nRequired size is %d sizeof(double)\n",qa_a_n);
   exit(0);
  } 
  qa_a_c=(double*)calloc(qa_a_n,sizeof(double));
  if(qa_a_c==NULL){
   printf("\nAllocating error for pulses_wth\n");
   printf("\nRequired size is %d sizeof(double)\n",qa_a_n);
   exit(0);
  } 
  qa_a_tp2=(double*)calloc(qa_a_n,sizeof(double));
  if(qa_a_tp2==NULL){
   printf("\nAllocating error for pulses_wth\n");
   printf("\nRequired size is %d sizeof(double)\n",qa_a_n);
   exit(0);
  } 
  qa_a_bigT2=(double*)calloc(qa_a_n,sizeof(double));
  if(qa_a_bigT2==NULL){
   printf("\nAllocating error for pulses_wth\n");
   printf("\nRequired size is %d sizeof(double)\n",qa_a_n);
   exit(0);
  } 
  i=0;
  while(fscanf(qa_a_file,"%lf %lf %lf %lf %lf %lf %lf %lf",&qa_a_bigA[i],&qa_a_tp1[i],&qa_a_bigT1[i],&qa_a_a[i],&qa_a_b[i],&qa_a_c[i],&qa_a_tp2[i],&qa_a_bigT2[i])!=EOF){
   i++;
   if(i>qa_a_n+1){
    printf("\n\nError: you have to check the *.qa_a file!");
    exit(0);
   }
  }
  qa_MatrixQ=(double **) calloc(Dim, sizeof(double));
  if (qa_MatrixQ==NULL){
   printf("\nAllocating error for *MatrixA\n");
   printf("\nRequired size is %d sizeof(double)",Dim);
   exit(0);
  } 
  for(i=0;i<Dim;i++){
   *(qa_MatrixQ+i)=(double *) calloc(Dim, sizeof(double));
   if(*(qa_MatrixQ+i)==NULL){
    printf("\nAllocating error for *(MatrixA+i)\n");
    printf("\nRequired size is %d sizeof(double)",Dim);
    exit(0);
   }
  } 
  while(fscanf(qa_q_file,"%d %d %lf",&i,&j,&component)!=EOF){
   if(i<Dim&&j<Dim){
    qa_MatrixQ[i][j]=component;
   }else{
    printf("the Q-Matrix index has exceeded the dimensions\n Dimension=%d  i=%d  j=%d",livelli,i,j);
    exit(0);
   }
  }
  for(i=0;i<Dim;i++)for(j=i+1;j<Dim;j++)qa_MatrixQ[j][i]=qa_MatrixQ[i][j];
 }
///////////////////////////////////////////////////////////////////////
}
///////////////////////////////////////////////////////////////////////
