#ifndef INCLUDE_PCONST
#define INCLUDE_PCONST

#include <stdbool.h>
#include <def-undef.h>
#include "precision_mod.h"
#include "param_mod.h"

//     -----------------------------------------------------------
//     Index Fields
//     -----------------------------------------------------------
//YU   real,dimension(:),allocatable:: dyr_global
    extern int pconst_mod_mp_ix_;
    extern int pconst_mod_mp_iy_;
    extern int pconst_mod_mp_i_global_[imt];
    extern int pconst_mod_mp_j_global_[jmt];
    extern double pconst_mod_mp_ahv_back_[jmt_global];
    extern double (*pconst_mod_mp_viv_)[km][jmt][imt];
    extern double (*pconst_mod_mp_vit_)[km][jmt][imt];
//     real(r8),dimension(imt_global,jmt,km):: vit_1d,viv_1d
//      integer,dimension(imt,jmt,max_blocks_clinic):: basin
//     integer,dimension(imt,jmt,max_blocks_clinic):: itnu
//lhl1204
    extern double pconst_mod_mp_dfricmx_;
    extern double pconst_mod_mp_dwndmix_;
    extern double pconst_mod_mp_pr_number_;
    extern double pconst_mod_mp_na_[max_blocks_clinic][jmt][imt];
#ifdef BCKMEX
    extern double pconst_mod_mp_diff_back_eq_;
    extern double pconst_mod_mp_diff_back_coef_;
    extern double pconst_mod_mp_diff_back_coef_max_;
#endif
// lihuimin, 2012.7.15
//     integer :: i_num   // actual grid number in i direction in this process
//     integer :: j_num   // actual grid number in j direction in this process
//     integer :: i_f_num // formal grid number in i dircetion 
//     integer :: j_f_num // formal grid number in j direction //lhl1204


//     -----------------------------------------------------------
//     Grids
//     -----------------------------------------------------------
#if (defined NETCDF) || (defined ALL)
    extern float pconst_mod_mp_lon_[imt_global];
    extern float pconst_mod_mp_lat_[jmt_global];
    extern float pconst_mod_mp_lon_o_[jmt_global][imt_global];
    extern float pconst_mod_mp_lat_o_[jmt_global][imt_global];
    extern float pconst_mod_mp_ulon_o_[jmt_global][imt_global];
    extern float pconst_mod_mp_ulat_o_[jmt_global][imt_global];
    extern float pconst_mod_mp_lev_[km];
    extern float pconst_mod_mp_lev1_[km+1];
#endif
//lhl090729
    extern double pconst_mod_mp_s_lon_[s_imt];
    extern double pconst_mod_mp_s_lat_[s_jmt];
//lhl090729
    extern double pconst_mod_mp_zkt_[km];
    extern double pconst_mod_mp_dzp_[km];
    extern double pconst_mod_mp_odzp_[km];
    extern double pconst_mod_mp_odzt_[km];
    extern double pconst_mod_mp_zkp_[kmp1];
    extern double (*pconst_mod_mp_ebea_)[jmt][imt];
    extern double (*pconst_mod_mp_ebeb_)[jmt][imt];
    extern double (*pconst_mod_mp_ebla_)[jmt][imt];
    extern double (*pconst_mod_mp_eblb_)[jmt][imt];
    extern double (*pconst_mod_mp_epea_)[jmt][imt];
    extern double (*pconst_mod_mp_epeb_)[jmt][imt];
    extern double (*pconst_mod_mp_epla_)[jmt][imt];
    extern double (*pconst_mod_mp_eplb_)[jmt][imt];
    extern double (*pconst_mod_mp_rrd1_)[jmt][imt];
    extern double (*pconst_mod_mp_rrd2_)[jmt][imt];
//lhl060506
// #if ( defined SMAG)
// #endif

//Yu
    extern double (*pconst_mod_mp_ohbt_)[jmt][imt];
    extern double (*pconst_mod_mp_ohbu_)[jmt][imt];
    extern double (*pconst_mod_mp_dzph_)[jmt][imt];
    extern double (*pconst_mod_mp_hbx_)[jmt][imt];
    extern double (*pconst_mod_mp_hby_)[jmt][imt];
    extern double (*pconst_mod_mp_snlat_)[jmt][imt];
//     real(r8),dimension(jmt):: COSU,COST
//     real(r8),dimension(:),allocatable::COSU_global,COST_global
    extern int *pconst_mod_mp_i_start_;
    extern int *pconst_mod_mp_j_start_;
//     real(r8),dimension(imt)::CF1,CF2,SF1,SF2


//     -----------------------------------------------------------
//     Reference T S & coefficients for calculation of d(density)
//     -----------------------------------------------------------
    extern double pconst_mod_mp_to_[km];
    extern double pconst_mod_mp_so_[km];
    extern double pconst_mod_mp_c_[9][km];
    extern double pconst_mod_mp_po_[km];


//lhl1204
//YU
//     -----------------------------------------------------------
//     Control Parameter
//     -----------------------------------------------------------
    extern int pconst_mod_mp_isop_;


//     -----------------------------------------------------------
//     Phycical Parameter
//     -----------------------------------------------------------
    extern double pconst_mod_mp_amv_;
    extern double pconst_mod_mp_ahv_;
    extern double pconst_mod_mp_ahice_;
//lhl1204
    extern double pconst_mod_mp_akmu_[max_blocks_clinic][km][jmt][imt];
    extern double pconst_mod_mp_akmt_[max_blocks_clinic][km][jmt][imt];
    extern double pconst_mod_mp_akt_[max_blocks_clinic][ntra][km][jmt][imt];
//lhl1204
    extern double pconst_mod_mp_am_[jmt];
    extern double pconst_mod_mp_ah_[jmt];
    extern double pconst_mod_mp_am3_[max_blocks_clinic][km][jmt][imt];
    extern double pconst_mod_mp_ah3_[max_blocks_clinic][km][jmt][imt];
//lhl
    extern double pconst_mod_mp_amx_[max_blocks_clinic][km][jmt][imt];
    extern double pconst_mod_mp_amy_[max_blocks_clinic][km][jmt][imt];
//lhl
    extern double pconst_mod_mp_gamma_;

#if ( defined SMAG)
    //   real(r8):: D0,CP,C0F,TBICE,OD0,SAG,CAG,OD0CP,ASEA, &
    //                 VSEA,AFB1,AFB2,AFC1,AFC2,AFT1,AFT2,KARMAN,RR
#else
    //   real(r8)::  D0,CP,C0F,TBICE,OD0,SAG,CAG,OD0CP,ASEA, &
    //                 VSEA,AFB1,AFB2,AFC1,AFC2,AFT1,AFT2
#endif
    extern double pconst_mod_mp_d0_;   
    extern double pconst_mod_mp_cp_;
    extern double pconst_mod_mp_c0f_;  
    extern double pconst_mod_mp_tbice_;
    extern double pconst_mod_mp_od0_;  
    extern double pconst_mod_mp_sag_;
    extern double pconst_mod_mp_cag_;
    extern double pconst_mod_mp_od0cp_; 
    extern double pconst_mod_mp_asea_;
    extern double pconst_mod_mp_vsea_;
    extern double pconst_mod_mp_afb1_;
    extern double pconst_mod_mp_afb2_;
    extern double pconst_mod_mp_afc1_; 
    extern double pconst_mod_mp_afc2_; 
    extern double pconst_mod_mp_aft1_; 
    extern double pconst_mod_mp_aft2_; 

    extern bool pconst_mod_mp_diag_msf_;
    extern bool pconst_mod_mp_diag_bsf_;
    extern bool pconst_mod_mp_diag_mth_;
    extern bool pconst_mod_mp_diag_budget_;
    extern bool pconst_mod_mp_test_input_;
//
    extern char pconst_mod_mp_abmon_[3][12];
    extern char pconst_mod_mp_abmon1_[3][12];
    extern char pconst_mod_mp_out_dir_[80];
//
    extern double pconst_mod_mp_dtb_;
    extern double pconst_mod_mp_dtc_;
    extern double pconst_mod_mp_dts_;
    extern double pconst_mod_mp_dtb2_;
    extern double pconst_mod_mp_dtc2_;
    extern double pconst_mod_mp_onbb_;
    extern double pconst_mod_mp_onbc_;
    extern double pconst_mod_mp_oncc_;
    //   REAL(r8):: DTB,DTC,DTS,DTB2,DTC2,ONBB,ONBC,ONCC
    extern int pconst_mod_mp_nbb_;
    extern int pconst_mod_mp_ncc_;
    extern int pconst_mod_mp_nss_;
    extern int pconst_mod_mp_isb_;
    extern int pconst_mod_mp_isc_;
    extern int pconst_mod_mp_ist_;
    extern int pconst_mod_mp_month_;
    extern int pconst_mod_mp_number_day_;
    extern int pconst_mod_mp_number_month_;
    extern int pconst_mod_mp_nmonth_[12];
    extern int pconst_mod_mp_nnmonth_[12];
//lihuimin 2012.6.18, add REFDATE
      extern int pconst_mod_mp_number_;
      extern int pconst_mod_mp_nstart_;
      extern int pconst_mod_mp_iyo_;
      extern int pconst_mod_mp_iyfm_;
      extern int pconst_mod_mp_mon0_;
      extern int pconst_mod_mp_mend_;
      extern int pconst_mod_mp_imd_;
      extern int pconst_mod_mp_iday_;
      extern int pconst_mod_mp_ii_;
      extern int pconst_mod_mp_jj_;
      extern int pconst_mod_mp_io_hist_;
      extern int pconst_mod_mp_io_rest_;
      extern int pconst_mod_mp_rest_freq_;
      extern int pconst_mod_mp_hist_freq_;
      extern int pconst_mod_mp_refdate_;
      extern int pconst_mod_mp_boundary_restore_;
      extern int pconst_mod_mp_klv_;
//      extern int pconst_mod_mp_kvt_; //use the layer for simple assimilation
//
    extern char pconst_mod_mp_adv_momentum_[80];
    extern char pconst_mod_mp_adv_tracer_[80];
    //   character (len=80) :: adv_momentum, adv_tracer
    extern int pconst_mod_mp_curr_ymd_licom_;  // Current date YYYYMMDD
    extern int pconst_mod_mp_curr_ymd_cpl_;
    extern int pconst_mod_mp_yy_licom_;            // year, month, day
    extern int pconst_mod_mp_mm_licom_;
    extern int pconst_mod_mp_dd_licom_;
    extern int pconst_mod_mp_tod_licom_;
    extern int pconst_mod_mp_yy_cpl_;
    extern int pconst_mod_mp_mm_cpl_;
    extern int pconst_mod_mp_dd_cpl_;
    extern int pconst_mod_mp_tod_cpl_;
    extern int pconst_mod_mp_imonth_;
    extern int pconst_mod_mp_ihour_;
    extern int pconst_mod_mp_iminute_;
    extern int pconst_mod_mp_isecond_;
    extern int pconst_mod_mp_tod_;
    extern int pconst_mod_mp_ymd_;
    extern int pconst_mod_mp_ymd_sync_;
    extern int pconst_mod_mp_tod_sync_;
//lhl20130419
    extern int pconst_mod_mp_ocn_cpl_dt_;
    extern int pconst_mod_mp_licom_cpl_dt_;
    extern int pconst_mod_mp_yearadd_;   //number year and day will be added to the restart file 
    extern int pconst_mod_mp_dayadd_;    //number year and day will be added to the restart file 
    extern int pconst_mod_mp_iyfmfnew_;  //number forcing yr will be added to the restart file 
    extern int pconst_mod_mp_first_step_; //test if first step 
    extern int pconst_mod_mp_num_step_per_day_;  //test if first step 
    extern int pconst_mod_mp_num_output_;         //number of accumulation for output 
    extern int pconst_mod_mp_num_outputacc_;    //number of accumulation for output 
    extern int pconst_mod_mp_num_dailyoutput_;  //number of accumulation for daily output 

    extern bool pconst_mod_mp_dts_accum_;     //test if daily accumulation or T/S step accumula
    extern bool pconst_mod_mp_daily_accum_;   //test if daily accumulation or T/S step accumula
    extern bool pconst_mod_mp_dailybudget_accum_;  //test if daily accumulation or T/S step accumula
    extern bool pconst_mod_mp_simple_assm_;  //using the simple assimilation(restoring) or not 

#if ( defined TIDEMIX )
    //yuzp-2016/12/8
    #define local_mixing_fraction (double)0.33
    #define mixing_ef (double)0.20
    #define shelf_cutoff (double)1500.
    #define decay_scale (double)1000.
    #define back_tidalmixing (double)1.0e-6
    #define max_tidalmixing (double)1.0e-2
    #define max_wavedis (double)0.1
    extern double (*pconst_mod_mp_fz_tide_)[km][jmt][imt];
    extern double (*pconst_mod_mp_ak_tide_)[km][jmt][imt];
    extern double pconst_mod_mp_fztidal_[max_blocks_clinic][km][jmt][imt];         //yuzp-2016/11/13
    extern double pconst_mod_mp_richardson_[max_blocks_clinic][km][jmt][imt];
    extern double pconst_mod_mp_wp3_tidal_[max_blocks_clinic][km][jmt][imt];
    extern double pconst_mod_mp_ak_tide1_[max_blocks_clinic][km][jmt][imt];    //yuzp-2016/11/19
#endif

#if ( defined CANUTOMIXOUT )
    extern double pconst_mod_mp_wp1_canuto_[max_blocks_clinic][km][jmt][imt];  //yuzp-2016/12/4
    extern double pconst_mod_mp_wp2_canuto_[max_blocks_clinic][km][jmt][imt];
    extern double pconst_mod_mp_wp3_canuto_[max_blocks_clinic][km][jmt][imt];
    extern double pconst_mod_mp_wp4_canuto_[max_blocks_clinic][km][jmt][imt];
    extern double pconst_mod_mp_wp5_canuto_[max_blocks_clinic][km][jmt][imt];
    extern double pconst_mod_mp_wp6_canuto_[max_blocks_clinic][km][jmt][imt];
    extern double pconst_mod_mp_wp7_canuto_[max_blocks_clinic][km][jmt][imt];
    extern double pconst_mod_mp_wp8_canuto_[max_blocks_clinic][km][jmt][imt];
    //    real(r8),dimension(imt,jmt,km,max_blocks_clinic) ::wp1_canuto,wp2_canuto,wp3_canuto,wp4_canuto,&
    //                                                       wp5_canuto,wp6_canuto,wp7_canuto,wp8_canuto     //yuzp-2016/12/4
    extern double pconst_mod_mp_wk1_canuto_[max_blocks_clinic][km][jmt][imt];
    extern double pconst_mod_mp_wk2_canuto_[max_blocks_clinic][km][jmt][imt];
    extern double pconst_mod_mp_wk3_canuto_[max_blocks_clinic][km][jmt][imt];
    extern double pconst_mod_mp_wk4_canuto_[max_blocks_clinic][km][jmt][imt];
    //    real(r8),dimension(imt,jmt,km,max_blocks_clinic) ::wk1_canuto,wk2_canuto,wk3_canuto,wk4_canuto     //yuzp-2016/12/4
    extern double pconst_mod_mp_fcort_canuto_[max_blocks_clinic][km][jmt][imt];//yuzp-2016/12/4
    extern double pconst_mod_mp_fcor_canuto_[max_blocks_clinic][km][jmt][imt];
    //    real(r8),dimension(imt,jmt,max_blocks_clinic) ::fcor_canuto,fcort_canuto     //yuzp-2016/12/4
    extern double pconst_mod_mp_wp12_canuto_[max_blocks_clinic][km][jmt][imt];
    extern double pconst_mod_mp_wp13_canuto_[max_blocks_clinic][km][jmt][imt];
    extern double pconst_mod_mp_wp10_canuto_[max_blocks_clinic][jmt][imt];
    extern double pconst_mod_mp_wp11_canuto_[max_blocks_clinic][jmt][imt];
    extern double pconst_mod_mp_alpha_canuto_[max_blocks_clinic][km][jmt][imt];
    extern double pconst_mod_mp_beta_canuto_[max_blocks_clinic][km][jmt][imt];
    //    real(r8),dimension(imt,jmt,km,max_blocks_clinic) ::wp12_canuto,wp13_canuto     //yuzp-2016/12/4
    //    real(r8),dimension(imt,jmt,max_blocks_clinic) ::wp10_canuto,wp11_canuto     //yuzp-2016/12/4
    //    real(r8),dimension(imt,jmt,km,max_blocks_clinic) ::alpha_canuto,beta_canuto     //yuzp-2016/12/4
#endif

#endif // //INCLUDE_PCONST
