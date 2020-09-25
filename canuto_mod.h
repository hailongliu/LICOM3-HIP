#ifndef INCLUDED_CANUTO
#define INCLUDED_CANUTO

#include  "precision_mod.h"
#include <math.h>
#include <stdbool.h>

#define nmodel_ 1
#define ntbl 501
#define ifback 5
#define ifsali 1
#define ri0 (double)-4.0
//#define e (double)2.71828182845904509
#define pi_const double(3.14159265358979312)
//#define pi_const (double)3.14159265358979323846
#define s2min (double)1.0e-14
#define ifshearmin true
#define ifzeroshear true
#define icondear 0
#define ifepson2 2
#define epson2__ (double)0.288e0
#define ifbotenhance 0
#define eps_bot0 (double)2.0e-5
#define scale_bot (double)5.0e4
#define ifdeeplat 1
#define eplatidependmin (double)7.0e-2
#define ifcheckbottomeps 0
//extern double eps_bot__under;
#define ifrafgmax 1
#define ifsalback 5
#define nextrtbl0 62
#define ifexpabstable 1	
#define ifast 1
#define ifastexpabs ifast*ifexpabstable
#define nextrtbl1 1000
#define nextrtbl  1062
#define nposapprox 101
#define mt0 400
#define mt 1462
extern double canuto_mod_mp_and2on2a1_[2*mt+1];
extern double canuto_mod_mp_amtaun2a1_[2*mt+1];
extern double canuto_mod_mp_dand2on2_;
extern double canuto_mod_mp_sma1_[2*mt+1];
extern double canuto_mod_mp_sha1_[2*mt+1];
extern double canuto_mod_mp_ssa1_[2*mt+1];
extern double canuto_mod_mp_rri_,canuto_mod_mp_rnd2on2_,canuto_mod_mp_dri_, canuto_mod_mp_deltheta_r_,canuto_mod_mp_b1_,canuto_mod_mp_theta_rcrp_, canuto_mod_mp_theta_rcrn_,canuto_mod_mp_ako_,canuto_mod_mp_back_l_0_;

	
#define ifchengcon 0
#define idefmld 0
#define deltemld (double)0.5e0
#define delrhmld (double)0.125e-3
extern int canuto_mod_mp_ilomega_, canuto_mod_mp_ifirst_;
#define ilomega 0
#define amldminlom (double)5.0e4
#define nbig 100   

extern double canuto_mod_mp_rib_[2*mt+1];
extern double canuto_mod_mp_ridb_[2*mt+1];
extern double (canuto_mod_mp_slq2b_)[2*mt+1][2*mt+1];
extern double (canuto_mod_mp_smb_)[2*mt+1][2*mt+1];
extern double (canuto_mod_mp_shb_)[2*mt+1][2*mt+1];
extern double (canuto_mod_mp_ssb_)[2*mt+1][2*mt+1];
extern int canuto_mod_mp_irimax_[2*mt+1];
#define mt_ra_r (nposapprox-1)
#define n_theta_r_oct 75
extern double canuto_mod_mp_sisamax_[4*n_theta_r_oct+1];
extern double canuto_mod_mp_ra_rmax_[4*n_theta_r_oct+1];
extern double canuto_mod_mp_c_y_r0_[4*n_theta_r_oct+1];
extern double canuto_mod_mp_back_ra_r_[4*n_theta_r_oct+1];
extern double canuto_mod_mp_sm_r1_[4*n_theta_r_oct+1];
extern double canuto_mod_mp_sh_r1_[4*n_theta_r_oct+1];
extern double canuto_mod_mp_ss_r1_[4*n_theta_r_oct+1];
extern double canuto_mod_mp_slq2_r1_[4*n_theta_r_oct+1];
extern double canuto_mod_mp_rimax_;
#define ifpolartablewrite 0
#define ifbg_theta_interp 1
#define back_ph_0 (double)(6.0e-5)*(1.0e2/(2.0e0*pi_const))
#define adjust_gargett (double)1.0e0
#define back_k_0 (double)(0.1e0)*(2.0e0)*pi_const*(1.0e-2)*adjust_gargett
#define back_del_0 (double)pi_const/back_k_0
#define back_s2 (double)back_ph_0*back_k_0
#define back_sm2 (double)1.0e0/back_s2
#define v_back0 (double)0.01
#define t_back0 (double)0.01
#define s_back0 (double)0.01
#define ri_internal (double)1.0e0
#define backfrac (double)85.0e-2
#define backfact (double)exp(-1)
extern double canuto_mod_mp_ria_[ntbl],canuto_mod_mp_slq2a_[ntbl],canuto_mod_mp_sma_[ntbl],canuto_mod_mp_sha_[ntbl];
extern double canuto_mod_mp_and2on2a1_[2*mt+1],canuto_mod_mp_amtaun2a1_[2*mt+1];
#endif 
