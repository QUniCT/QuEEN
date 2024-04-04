// system stuff //
double **MatrixA;
double *VectorGamma;
double **vectorgamma_pointer;
double *vectorgammasum;
double *pb;
double **pb_pointer;
double *decay_rate;
int *receptive_state;
int *donor_state;
_Complex double *x0;
double *population_sink;
int Ntrajectories;
int decay,decay_1;
int Dim;
int jump_dim;
int num_threads;
//  dephasing //
int add_white_noise_dephasing;
int add_static_noise_dephasing;
int add_ou_process_dephasing;
int add_pl_wave_dephasing;
int add_fm_signal_dephasing;
int add_sq_wave_dephasing, *sqw_samp;
int add_filtered_noise_dephasing, *filtered_noise_order, filtered_noise_processes, *filtered_noise_site;
double *dph_gamma;
double **snoise;
double *oup_sigma, *oup_tau;
double *pwd_amp, *pwd_freq, *pwd_phase;
double *fmn_amp, *fmn_base, *fmn_signal, *fmn_dev;
double *sqw_amp, *sqw_freq;
double *filtered_noise_width, *filtered_noise_lfreq, *filtered_noise_hfreq, **filtered_noise_freqres;
// correlate dephasing //
int add_sitecorr_wn;
int add_sitecorr_sn;
int add_sitecorr_ou;
double *site_corr_rate;
double *site_anticorr_rate;
int *site_corr;
// temporal variables //
double h;
double tmin;
double tmax;
double Dt;
int NstepDt; 
int Nsteph;
//  density and purity variables  //
long int count_jumps=0,count_jumps_intosink=0,count_dephase=0;
_Complex double *rho_sys,***rho_sys_pointer;
int calculate_sd;
double ***dm_save;
double **standard_deviation;
int calc_purity;
double *purity;
// pulse variables //
int use_pulses,pulses_count=0;
double **pulses_omega;
double **pulses_width;
double **pulses_tau;
int *pulses_site_i;
int *pulses_site_j;
double *pulses_amp;
double *pulses_now;
double *pulses_wth;
// QA variables //
int use_qa,qa_a_n;
double *qa_a_bigA,*qa_a_tp1,*qa_a_bigT1,*qa_a_a,*qa_a_b,*qa_a_c,*qa_a_tp2,*qa_a_bigT2;
double **qa_MatrixQ;


//  input/output files  //
char *file_label;
char *file_name;
// in //
FILE *hamiltonian_file;
FILE *dissipator_file;
FILE *wnd_file;
FILE *snd_file;
FILE *oud_file;
FILE *pwd_file;
FILE *fmd_file;
FILE *sqd_file;
FILE *fnd_file;
FILE *initcon_file;
FILE *sitecorr_file;
FILE *pulses_file;
FILE *qa_q_file;
FILE *qa_a_file;
// out //
FILE *density_file;
FILE *purity_file;
FILE *population_file;
FILE *sd_file;
FILE *noise_file;
