#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <unistd.h>
#include <fcntl.h>
#include <complex.h>
#include <ctype.h>

#include "./Global-variables.h"
#include "./service_routines.h"
#include "./read_input.h"
#include "./jump_and_runge.h"

//////  gcc QuEEN.c -O3 -lm -fopenmp -o QuEEN
int main(int argc, char *argv[])
{
 register int i,m,p,sd,jlvl,Ntrajectories_save;
 long int count_jumps=0, count_jumps_intosink=0;
 double pop_tot=0;
// printf("\n Program is running!\n\n");
 printf("\n And we're off!\n\n");
///////////////////////////////////////////////////////////////////////
//////  set defaults  /////////////////////////////////////////////////
////  most of these can be customised via the command line         ////
///////////////////////////////////////////////////////////////////////
// midum=1000; //Set the seed of the pseudo number generator 
 Ntrajectories=10000;//-----------Set the number of trajectories
 num_threads=1;//-----------------use 1 thread unless asked to use more	
 calc_purity=0;//-----------------don't calculate purity unless asked
 add_white_noise_dephasing=0;//---don't add adiabatic noise unless asked
 add_static_noise_dephasing=0;//--don't add adiabatic noise unless asked
 add_ou_process_dephasing=0;//----don't generate an ou process unless asked
 add_pl_wave_dephasing=0;//-------don't generate a plase-wave process unless asked
 add_fm_signal_dephasing=0;//-----don't generate frequency modulated noise process unless asked
 add_sq_wave_dephasing=0;//-------don't generate a square-wave process unless asked
 add_filtered_noise_dephasing=0;//don't generate filtered noise unless asked
 add_sitecorr_wn=0;//-------------don't correlate the markovian site dephasing unless asked
 add_sitecorr_sn=0;//-------------don't correlate the static site dephasing unless asked
 add_sitecorr_ou=0;//-------------don't correlate the ou site dephasing unless asked
 use_pulses=0;//------------------don't use pulses unless asked
 use_pulses=0;//------------------don't use qa algorithm unless asked
 calculate_sd=0;//----------------don't calculate the standard deviation unless asked
 h=0.01;//------------------------Setting of Runge-Kutta step
 tmin=0;//------------------------Setting of temporal parameters
 tmax=10;//-----------------------Setting of temporal parameters
 Dt=1; //-------------------------coarse grid (it must be larger than h!) 
///////////////////////////////////////////////////////////////////////
//////  search the command line  //////////////////////////////////////
////  contains algorithm information: sample size, temporal        ////
////  parameters, flags to add noise or calculate purity etc.      ////
////  All system information is contained in the input files       ////
///////////////////////////////////////////////////////////////////////
 if(argc>1){
  for(i=1;i<argc;i++){
   if(strcmp(argv[i],"-label")==0)file_label=argv[i+1];
   if(strcmp(argv[i],"-sample")==0)Ntrajectories=atoi(argv[i+1]);
   if(strcmp(argv[i],"-tmax")==0)tmax=strtod(argv[i+1],NULL);
   if(strcmp(argv[i],"-coarse_grid")==0||strcmp(argv[i],"-cg")==0)Dt=strtod(argv[i+1],NULL);
   if(strcmp(argv[i],"-fine_grid")==0||strcmp(argv[i],"-fg")==0)h=strtod(argv[i+1],NULL);
   if(strcmp(argv[i],"-threads")==0)num_threads=atoi(argv[i+1]);
   if(strcmp(argv[i],"-white_noise")==0||strcmp(argv[i],"-wn")==0)add_white_noise_dephasing=1;
   if(strcmp(argv[i],"-filtered_noise")==0||strcmp(argv[i],"-fn")==0)add_filtered_noise_dephasing=1;
   if(strcmp(argv[i],"-static_noise")==0||strcmp(argv[i],"-sn")==0)add_static_noise_dephasing=1;
   if(strcmp(argv[i],"-ornstein_uhlenbeck")==0||strcmp(argv[i],"-ou")==0)add_ou_process_dephasing=1;
   if(strcmp(argv[i],"-plane_wave")==0||strcmp(argv[i],"-pw")==0)add_pl_wave_dephasing=1;
   if(strcmp(argv[i],"-frequency_modulated")==0||strcmp(argv[i],"-fm")==0)add_fm_signal_dephasing=1;
   if(strcmp(argv[i],"-square_wave")==0||strcmp(argv[i],"-sq")==0)add_sq_wave_dephasing=1;
   if(strcmp(argv[i],"-correlated_white_noise")==0||strcmp(argv[i],"-cwn")==0){
    add_white_noise_dephasing=1;
    add_sitecorr_wn=1;
   }
   if(strcmp(argv[i],"-correlated_static_noise")==0||strcmp(argv[i],"-csn")==0){
    add_static_noise_dephasing=1;
    add_sitecorr_sn=1;
   }
   if(strcmp(argv[i],"-correlated_ornstein_uhlenbeck")==0||strcmp(argv[i],"-cou")==0){
    add_ou_process_dephasing=1;
    add_sitecorr_ou=1;
   }
   if(strcmp(argv[i],"-pulses")==0||strcmp(argv[i],"-pu")==0)use_pulses=1;
   if(strcmp(argv[i],"-qa")==0)use_qa=1;
   if(strcmp(argv[i],"-purity")==0)calc_purity=1;
   if(strcmp(argv[i],"-convergence")==0||strcmp(argv[i],"-con")==0)calculate_sd=1;
  }
 }
 if(file_label==NULL){
  printf("Please select a system label with the command line flag -label. That's where the system info should be. Don't make me ask again\n");
  exit(0);
 }
///////////////////////////////////////////////////////////////////////
//////  open the files  ///////////////////////////////////////////////
////  *.hamiltonian, *.dissipator, *.initcon                       ////
////   are manditory input files                                   ////
////  *.dph, *.sgn, *,oup, and *.sitecorr                          ////
////   are optional input files                                    ////
////  *.dm, *.purity, *sd, and *.population                        ////
////   are output files                                            ////
///////////////////////////////////////////////////////////////////////
 apertura_file_scrittura();
 read_input();
///////////////////////////////////////////////////////////////////////
//////  setup jump information  ///////////////////////////////////////
////  this indexes all non-zero relaxation rates so the algorithm  ////
////  can loop over the non-zero components only when it           ////
////  calculates jumps                                             ////
///////////////////////////////////////////////////////////////////////
 for(i=0;i<Dim;i++)for(p=0;p<Dim;p++)vectorgammasum[p]+=*(vectorgamma_pointer[i]+p);
 if(add_white_noise_dephasing)for(p=0;p<Dim;p++)vectorgammasum[p]+=dph_gamma[p];
// if(add_sitecorr&&add_white_noise_dephasing){
//  for(i=0;i<Dim;i++)if(site_corr_rate[i]!=0){
 //  *(vectorgamma_pointer[site_corr[i]]+i)+=site_corr_rate[i];
 //  *(vectorgamma_pointer[i]+site_corr[i])+=site_corr_rate[i];
//  }
// }
 decay=0;
 for(i=0;i<jump_dim;i++)if(VectorGamma[i]!=0)decay++;
 if(add_white_noise_dephasing)for(i=0;i<Dim;i++)if(dph_gamma[i]!=0)decay++;
 if(add_sitecorr_wn)for(i=0;i<Dim;i++)if(site_corr_rate[i]!=0)decay++;
 if(decay){
  pb=(double*)calloc(decay,sizeof(double));
  if(pb==NULL){
   printf("\nAllocating error for pb\n");
   printf("\nRequired size is %d sizeof(double)",decay);
   exit(0);
  } 
  decay_rate=(double*)calloc(decay,sizeof(double));
  if(decay_rate==NULL){
   printf("\nAllocating error for decay_rate\n");
   printf("\nRequired size is %d sizeof(double)",decay);
   exit(0);
  } 
  receptive_state=(int*)calloc(decay,sizeof(int));
  if(receptive_state==NULL){
   printf("\nAllocating error for receptive_state\n");
   printf("\nRequired size is %d sizeof(double)",decay);
   exit(0);
  } 
  donor_state=(int*)calloc(decay,sizeof(int));
  if(donor_state==NULL){
   printf("\nAllocating error for donor_state\n");
   printf("\nRequired size is %d sizeof(double)",decay);
   exit(0);
  } 
  decay=0;
  for(jlvl=0;jlvl<Dim;jlvl++)for(p=0;p<Dim;p++)if(*(vectorgamma_pointer[jlvl]+p)!=0){
   receptive_state[decay]=jlvl;
   donor_state[decay]=p;
   decay_rate[decay]=*(vectorgamma_pointer[jlvl]+p);
   decay++;
  } 
  if(add_white_noise_dephasing)for(p=0;p<Dim;p++)if(dph_gamma[p]!=0){
   receptive_state[decay]=Dim;
   donor_state[decay]=p;
   decay_rate[decay]=dph_gamma[p];
   decay++;
  }
  if(add_sitecorr_wn){
   decay_1=decay;
   for(p=0;p<Dim;p++)if(site_corr_rate[p]!=0){
    receptive_state[decay]=Dim+1;
    donor_state[decay]=p;
    decay_rate[decay]=site_corr_rate[p];
    decay++;
   }
  }
//  for(p=0;p<decay;p++)vectorgammasum[donor_state[p]]+=decay_rate[p];
  for(i=0;i<Dim;i++)vectorgammasum[i]*=0.5;
 }
 if(!add_ou_process_dephasing&&!add_static_noise_dephasing&&!decay&&!add_filtered_noise_dephasing)Ntrajectories=1;
///////////////////////////////////////////////////////////////////////
//////  adjust temporal parameters  ///////////////////////////////////
////  create an integer number of steps for the coarse loop and    ////
////  the fine sub-loops                                           ////
///////////////////////////////////////////////////////////////////////
 Nsteph=Dt/h; 
 if((Dt-Nsteph*h)!=0)Dt=Nsteph*h;
 NstepDt=tmax/Dt;
 if((tmax-NstepDt*Dt)!=0)tmax=NstepDt*Dt;
 printf("\n\nDt: %.7e and tmax: %.7e\n",Dt,tmax);
///////////////////////////////////////////////////////////////////////
//////  Allocate variables  ///////////////////////////////////////////
////  allocate all the variables that need NstepDt components in   ////
////  one of its dimensions                                        ////
///////////////////////////////////////////////////////////////////////
 alloca_operatore_densita();
 if(calculate_sd)allocate_dm_save();
 for(i=0;i<Dim;i++)for(m=0;m<=i;m++)*(rho_sys_pointer[0][i]+m)+=x0[i]*~x0[m];
 scrivi_op_densita(0);
 Ntrajectories_save=Ntrajectories;
 jump_and_runge();
 for(m=1;m<=NstepDt;m++)normalizza_op_densita(m);
///////////////////////////////////////////////////////////////////////
//////  print the density matrix  /////////////////////////////////////
///////////////////////////////////////////////////////////////////////
 if(calc_purity)for(m=1;m<=NstepDt;m++)calculate_purity(m);
 if(calculate_sd)for(m=1;m<=NstepDt;m++){
  fprintf(sd_file,"%f",tmin+m*Dt);
  for(i=0;i<Dim;i++)fprintf(sd_file," %.7f",standard_deviation[i][m]);
  fprintf(sd_file,"\n");
 }
 for(m=1;m<=NstepDt;m++)scrivi_op_densita(m);
 for(m=1;m<=NstepDt;m++)scrivi_population(m);
///////////////////////////////////////////////////////////////////////
//////  print the final report  ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////
 for(i=0;i<Dim;i++)pop_tot+=creal(*(rho_sys_pointer[NstepDt][i]+i));
 pop_tot+=population_sink[NstepDt];
 printf("\nErrore sulle populazioni al tempo finale:%.2e\n",1-pop_tot);
 return 0;
}
