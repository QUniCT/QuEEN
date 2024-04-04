//Declare functions
int apertura_file_scrittura();
int genera_MatrixA();
int genera_VectorGamma();
int alloca_variabili();
int alloca_operatore_densita();
int ricava_op_densita(int passo);
int normalizza_op_densita(int passo);
int scrivi_op_densita(int passo);
int scrivi_population(int passo);
int normalizza_ket();
double machine_epsilon();
//Define functions
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
int apertura_file_scrittura()
{
 file_name=malloc(snprintf(NULL,0,"%s%s",file_label,".hamiltonian")+1);
//////  open label.hamiltonian  ///////////////////////////////////////
 sprintf(file_name,"%s%s",file_label,".hamiltonian");
 if((hamiltonian_file=fopen(file_name,"r"))==NULL){
  printf("\nError::%s is not there\n",file_name);
  exit(0);
 }
//////  open label.dissipator  ////////////////////////////////////////
 sprintf(file_name,"%s%s",file_label,".dissipator");
 if((dissipator_file=fopen(file_name,"r"))==NULL){
  printf("\nError::%s is not there\n",file_name);
  exit(0);
 }
//////  open label.initcon  ///////////////////////////////////////////
 sprintf(file_name,"%s%s",file_label,".initcon");
 if((initcon_file=fopen(file_name,"r"))==NULL){
  printf("\nError::%s is not there\n",file_name);
  exit(0);
 }
//////  maybe open label.wnd  /////////////////////////////////////////
 if(add_white_noise_dephasing){
  sprintf(file_name,"%s%s",file_label,".wnd");
  if((wnd_file=fopen(file_name,"r"))==NULL){
   printf("\nError::%s is not there\n",file_name);
   exit(0);
  }
 }
//////  maybe open label.oud  /////////////////////////////////////////
 if(add_ou_process_dephasing){
  sprintf(file_name,"%s%s",file_label,".oud");
  if((oud_file=fopen(file_name,"r"))==NULL){
   printf("\nError::%s is not there\n",file_name);
   exit(0);
  }
 }
//////  maybe open label.snd  /////////////////////////////////////////
 if(add_static_noise_dephasing){
  sprintf(file_name,"%s%s",file_label,".snd");
  if((snd_file=fopen(file_name,"r"))==NULL){
   printf("\nError::%s is not there\n",file_name);
   exit(0);
  }
 }
//////  maybe open label.fnd  /////////////////////////////////////////
 if(add_filtered_noise_dephasing){
  sprintf(file_name,"%s%s",file_label,".fnd");
  if((fnd_file=fopen(file_name,"r"))==NULL){
   printf("\nError::%s is not there\n",file_name);
   exit(0);
  }
 }
//////  maybe open label.pwd  /////////////////////////////////////////
 if(add_pl_wave_dephasing){
  sprintf(file_name,"%s%s",file_label,".pwd");
  if((pwd_file=fopen(file_name,"r"))==NULL){
   printf("\nError::%s is not there\n",file_name);
   exit(0);
  }
 }
//////  maybe open label.fmd  /////////////////////////////////////////
 if(add_fm_signal_dephasing){
  sprintf(file_name,"%s%s",file_label,".fmd");
  if((fmd_file=fopen(file_name,"r"))==NULL){
   printf("\nError::%s is not there\n",file_name);
   exit(0);
  }
 }
//////  maybe open label.sqd  /////////////////////////////////////////
 if(add_sq_wave_dephasing){
  sprintf(file_name,"%s%s",file_label,".sqd");
  if((sqd_file=fopen(file_name,"r"))==NULL){
   printf("\nError::%s is not there\n",file_name);
   exit(0);
  }
 }
//////  maybe open label.sitecorr  ////////////////////////////////////
 if(add_sitecorr_wn||add_sitecorr_sn||add_sitecorr_ou){
  sprintf(file_name,"%s%s",file_label,".sitecorr");
  if((sitecorr_file=fopen(file_name,"r"))==NULL){
   printf("\nError::%s is not there\n",file_name);
   exit(0);
  }
 }
//////  maybe open label.pulses  //////////////////////////////////////
 if(use_pulses){
  sprintf(file_name,"%s%s",file_label,".pulses");
  if((pulses_file=fopen(file_name,"r"))==NULL){
   printf("\nError::%s is not there\n",file_name);
   exit(0);
  }
 }
//////  maybe open label.qa  //////////////////////////////////////
 if(use_qa){
  sprintf(file_name,"%s%s",file_label,".qa_q");
  if((qa_q_file=fopen(file_name,"r"))==NULL){
   printf("\nError::%s is not there\n",file_name);
   exit(0);
  }
  sprintf(file_name,"%s%s",file_label,".qa_a");
  if((qa_a_file=fopen(file_name,"r"))==NULL){
   printf("\nError::%s is not there\n",file_name);
   exit(0);
  }
 }
//////  open label.dm  ////////////////////////////////////////////////
 sprintf(file_name,"%s%s",file_label,".dm");
 if((density_file=fopen(file_name,"w"))==NULL){
  printf("\nError::problem opening %s\n",file_name);
  exit(0);
 }
//////  open label.purity  ////////////////////////////////////////////
 if(calc_purity){
  sprintf(file_name,"%s%s",file_label,".purity");
  if((purity_file=fopen(file_name,"w"))==NULL){
   printf("\nError::problem opening %s\n",file_name);
   exit(0);
  }
 }
//////  open label.population  ////////////////////////////////////////
 sprintf(file_name,"%s%s",file_label,".population");
 if((population_file=fopen(file_name,"w"))==NULL){
  printf("\nError::problem opening %s\n",file_name);
  exit(0);
 }
//////  open label.sd  ////////////////////////////////////////////////
 if(calculate_sd){
  sprintf(file_name,"%s%s",file_label,".convergence");
  if((sd_file=fopen(file_name,"w"))==NULL){
   printf("\nError::%s is not there\n",file_name);
   exit(0);
  }
 }
return 1;
}
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
int genera_MatrixA()
{
 int i;
 MatrixA=(double **) calloc(Dim, sizeof(double));
 if (MatrixA==NULL){
  printf("\nAllocating error for *MatrixA\n");
  printf("\nRequired size is %d sizeof(double)",Dim);
  exit(0);
 } 
 for(i=0;i<Dim;i++){
  *(MatrixA+i)=(double *) calloc(Dim, sizeof(double));
  if(*(MatrixA+i)==NULL){
   printf("\nAllocating error for *(MatrixA+i)\n");
   printf("\nRequired size is %d sizeof(double)",Dim);
   exit(0);
  }
 } 
 return 0; 
}
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
int simmetrizza_MatrixA()
{
 int i,j;
 for(i=0;i<Dim;i++)for(j=i;j<Dim;j++)MatrixA[j][i]=MatrixA[i][j];
 return 0; 
}
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
int genera_VectorGamma()
{
 int i,j;
 VectorGamma=(double *) calloc(jump_dim, sizeof(double));
 vectorgamma_pointer=(double**) calloc(Dim, sizeof(double*));
 vectorgammasum=(double *) calloc(Dim, sizeof(double));
 if (VectorGamma==NULL){
  printf("\nAllocating error for *VectorGamma\n");
  printf("\nRequired size is %d sizeof(double)",Dim);
  exit(0);
 } 
 if(vectorgamma_pointer==NULL){
  printf("\nAllocating error for vectorgamma_pointer\n");
  printf("\nRequired size is %d sizeof(double)",NstepDt+1);
  exit(0);
 } 
 if(vectorgammasum==NULL){
  printf("\nAllocating error for vectorgammasum\n");
  printf("\nRequired size is %d sizeof(double)",NstepDt+1);
  exit(0);
 } 
 for(i=0;i<Dim;i++)vectorgamma_pointer[i]=&VectorGamma[i*Dim];
 return 0; 
}
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
int genera_initcon()
{
 x0=(_Complex double *) calloc(Dim, sizeof(_Complex double));
 if (x0==NULL){
  printf("\nAllocating error for *x0\n");
  printf("\nRequired size is %d sizeof(double)",Dim);
  exit(0);
 } 
}
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
int alloca_operatore_densita()
{
 int i,j,k,m,p;
///////////////////////////////////////////////////////////////////////
////// allocate the density matrix  ///////////////////////////////////
////  the density matrix, rho_sys[i][j](t), has three dimensions:  ////
////  two dimensions for the density matrix and one for time, so   ////
////  that the dynamics of the density matrix are store.           ////
////  Since the density matrix is hermitian, there is no need to   ////
////  store the upper diagonal elements:                           ////
////              rho_sys[i][j](t)=rho_sys[j][i](t)*               ////
///////////////////////////////////////////////////////////////////////
 p=0.5*Dim*(Dim+1);
 rho_sys=(_Complex double*)calloc((NstepDt+1)*p,sizeof(_Complex double));
 if(rho_sys==NULL){
  printf("\nAllocating error for *rho_sys(t)\n");
  printf("\nRequired size is %d sizeof(double)",(NstepDt+1)*p);
  exit(0);
 } 
 rho_sys_pointer=(_Complex double***)calloc(NstepDt+1,sizeof(_Complex double**));
 for(i=0;i<=NstepDt;i++){
  rho_sys_pointer[i]=(_Complex double**)calloc(Dim,sizeof(_Complex double*));
  if(rho_sys_pointer==NULL){
   printf("\nAllocating error for *rho_sys_pointer(t)\n");
   printf("\nRequired size is %d sizeof(double)",p);
   exit(0);
  }
 }
 for(i=0;i<=NstepDt;i++){
  for(j=0;j<Dim;j++){
   m=0; for(k=0;k<=j;k++)m+=k;
   rho_sys_pointer[i][j]=&rho_sys[i*p+m];
  }
 }
 population_sink=(double*)calloc(NstepDt+1,sizeof(double));
 if(population_sink==NULL){
  printf("\nAllocating error for population_sink\n");
  printf("\nRequired size is %d sizeof(double)",Dim);
  exit(0);
 } 
 if(calc_purity){
  purity=(double *)calloc(NstepDt+1,sizeof(double));
  if(purity==NULL){
   printf("\nAllocating error for purity\n");
   printf("\nRequired size is %d sizeof(double)",NstepDt+1);
   exit(0);
  } 
 }
 return 0; 
}
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
int calculate_purity(int passo)
{
 static int i,j;
 for(i=0;i<Dim;i++)for(j=0;j<i;j++)purity[passo]+=2*(creal(*(rho_sys_pointer[passo][i]+j))*creal(*(rho_sys_pointer[passo][i]+j))-cimag(*(rho_sys_pointer[passo][i]+j))*cimag(*(rho_sys_pointer[passo][i]+j)));
 for(j=0;j<Dim;j++)purity[passo]+=(creal(*(rho_sys_pointer[passo][j]+j))*creal(*(rho_sys_pointer[passo][j]+j))+cimag(*(rho_sys_pointer[passo][j]+j))*cimag(*(rho_sys_pointer[passo][j]+j)));
 fprintf(purity_file,"%d %f\n",passo,purity[passo]);
 return 1;
}
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
int normalizza_op_densita(int passo)
{
 static int i, j;
 for(i=0;i<Dim;i++)for(j=0;j<=i;j++){
  *(rho_sys_pointer[passo][i]+j)/=(Ntrajectories*num_threads);
 }
 population_sink[passo]/=(Ntrajectories*num_threads);
return 1;
}
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
int scrivi_op_densita(int passo)
{
 int i,j;
 fprintf(density_file,"%.7e",tmin+passo*Dt);
 for(i=0;i<Dim;i++)for(j=0;j<=i;j++)fprintf(density_file," %.7f %.7f",creal(*(rho_sys_pointer[passo][i]+j)),cimag(*(rho_sys_pointer[passo][i]+j)));
 fprintf(density_file," %.7f\n",population_sink[passo]);
 return 1;
}
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
int scrivi_population(int passo)
{
 int i;
 fprintf(population_file,"%.7e",tmin+passo*Dt);
 for(i=0;i<Dim;i++)fprintf(population_file," %.7f",creal(*(rho_sys_pointer[passo][i]+i)));
 fprintf(population_file," %.7f\n",population_sink[passo]);
 return 1;
}
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
int normalizza_ket(_Complex double *x)
{ 
 static int i; 
 static double norm, aux;
 aux=0;
 for(i=0;i<Dim;i++) aux+=creal(x[i])*creal(x[i])+cimag(x[i])*cimag(x[i]);
 norm=sqrt(aux);
 for(i=0;i<Dim;i++)x[i]/=norm; 
 return 1;
}
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
int allocate_dm_save()
{
 int i,j;
 dm_save=(double***)calloc(Dim,sizeof(double**));
 if(dm_save==NULL){
  printf("\nAllocating error for *dm_save(t)\n");
  printf("\nRequired size is %d sizeof(double)",Dim);
  exit(0);
 }
 for(i=0;i<Dim;i++){
  *(dm_save+i)=(double**)calloc(NstepDt+1,sizeof(double *));
  if (*(dm_save+i)==NULL){
   printf("\nAllocating error for dm_save[%d]\n",i);
   printf("\nRequired size is %d sizeof(double)",NstepDt);
   exit(0);
  } 
  for(j=0;j<NstepDt+1;j++){
   *(*(dm_save+i)+j)=(double*)calloc(calculate_sd,sizeof(double));
   if (*(*(dm_save+i)+j)==NULL){
    printf("\nAllocating error for rho_sys[%d][%d]\n",i,j);
    printf("\nRequired size is %d sizeof(double)",NstepDt);
    exit(0);
   } 
  }
 }
 standard_deviation=(double**)calloc(Dim,sizeof(double*));
 if(standard_deviation==NULL){
  printf("\nAllocating error for standard_deviation\n");
  printf("\nRequired size is %d sizeof(double)",Dim);
  exit(0);
 } 
 for(i=0;i<Dim;i++){
  *(standard_deviation+i)=(double*)calloc(NstepDt+1,sizeof(double));
  if (*(standard_deviation+i)==NULL){
   printf("\nAllocating error for standard_deviation[%d]\n",i);
   printf("\nRequired size is %d sizeof(double)",NstepDt);
   exit(0);
  } 
 }
 return 0; 
}
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
 double getgaussrandom(){
 int i;
 double random[2],w;
 do{
  for(i=0;i<2;i++){random[i]=((double)rand()/RAND_MAX)*2-1;}
  w=random[0]*random[0]+random[1]*random[1];
 }while(w>=1.0||w==0);
 return random[0]*sqrt((-2.0*log(w))/w);
}
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
float seed2(long *idum)
{
 register int j;
 long k;
 float temp;
  static long idum2=123456789;
  static long iy=0;
  static long iv[32];
 long midum=1000;
 long IM1=2147483563,IM2=2147483399,IMM1=2147483562;
 int IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,NTAB=32;
 float EPS=1.2e-7,RNMX=1.0-1.2e-7;
 double AM=1.0/2147483563,NDIV=1+2147483562/32;

 if(*idum<=0){
  *idum=(-*idum<1)?1:-*idum;
  for(j=NTAB+7;j>=0;j--){
   k=*idum/IQ1;
   *idum=IA1*(*idum-k*IQ1)-k*IR1;
   if(*idum<0) *idum+=IM1;
  }
 }
 k=*idum/IQ1;
 *idum=IA1*(*idum-k*IQ1)-k*IR1;
 if(*idum<0) *idum+=IM1;
 return 0;
}

