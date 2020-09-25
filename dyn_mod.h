#ifndef INCLUDE_DYN
#define INCLUDE_DYN

#include "precision_mod.h"
#include "param_mod.h"

extern double(*buffer)[imt_global];
extern double(*dyn_mod_mp_h0_)[jmt][imt];
extern double(*dyn_mod_mp_dlub_)[jmt][imt];
extern double(*dyn_mod_mp_dlvb_)[jmt][imt];

extern double(*dyn_mod_mp_gg_)[km][jmt][imt];
extern double(*dyn_mod_mp_dlu_)[km][jmt][imt];
extern double(*dyn_mod_mp_dlv_)[km][jmt][imt];
extern double(*dyn_mod_mp_u_)[km][jmt][imt];
extern double(*dyn_mod_mp_v_)[km][jmt][imt];

extern double (*dyn_mod_mp_ub_)[jmt][imt];
extern double (*dyn_mod_mp_vb_)[jmt][imt];
extern double (*dyn_mod_mp_ubp_)[jmt][imt];
extern double (*dyn_mod_mp_vbp_)[jmt][imt];
extern double (*dyn_mod_mp_up_)[km][jmt][imt];
extern double (*dyn_mod_mp_vp_)[km][jmt][imt];
extern double (*dyn_mod_mp_utf_)[km][jmt][imt];
extern double (*dyn_mod_mp_vtf_)[km][jmt][imt];
extern double (*dyn_mod_mp_utl_)[km][jmt][imt];
extern double (*dyn_mod_mp_vtl_)[km][jmt][imt];

extern double (*dyn_mod_mp_ws_)[kmp1][jmt][imt];

extern double (*dyn_mod_mp_h0p_)[jmt][imt];
extern double (*dyn_mod_mp_h0f_)[jmt][imt];
extern double (*dyn_mod_mp_h0l_)[jmt][imt];
extern double (*dyn_mod_mp_h0bf_)[jmt][imt];
extern double (*dyn_mod_mp_h0bl_)[jmt][imt];

extern double (*dyn_mod_mp_sbcx_)[jmt][imt];
extern double (*dyn_mod_mp_bbcx_)[jmt][imt];
extern double (*dyn_mod_mp_sbcy_)[jmt][imt];
extern double (*dyn_mod_mp_bbcy_)[jmt][imt];

#endif // !INCLUDE_DYN
