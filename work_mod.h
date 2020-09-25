#ifndef INCLUDE_WORK
#define INCLUDE_WORK

#include "precision_mod.h"
#include "param_mod.h"

extern double work_mod_mp_wkk_[kmp1];

extern double(*work_mod_mp_work1_g_)[imt_global];
// extern double(*work_mod_mp_work2_g_)[imt_global];
// extern double(*work_mod_mp_work3_g_)[imt_global];
extern double(*work_mod_mp_stf_)[jmt][imt];
extern double(*work_mod_mp_tf_)[km][jmt][imt];
extern double(*work_mod_mp_wkb_)[km][jmt][imt];
extern double(*work_mod_mp_wkc_)[km][jmt][imt];
extern double(*work_mod_mp_wkd_)[km][jmt][imt];
// extern double(*work_mod_mp_uk_)[km][jmt][imt];
// extern double(*work_mod_mp_vk_)[km][jmt][imt];

extern double (*work_mod_mp_wka_)[km][jmt][imt];
extern double (*work_mod_mp_work_)[jmt][imt];
extern double work_mod_mp_wgp_[max_blocks_clinic][jmt][imt];
extern double work_mod_mp_pax_[max_blocks_clinic][jmt][imt];
extern double work_mod_mp_pay_[max_blocks_clinic][jmt][imt];
extern double work_mod_mp_pxb_[max_blocks_clinic][jmt][imt];
extern double work_mod_mp_pyb_[max_blocks_clinic][jmt][imt];
extern double work_mod_mp_whx_[max_blocks_clinic][jmt][imt];
extern double work_mod_mp_why_[max_blocks_clinic][jmt][imt];


#endif // !INCLUDE_WORK
