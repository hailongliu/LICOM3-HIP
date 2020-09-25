#ifndef INCLUDE_TRACER
#define INCLUDE_TRACER

#include "precision_mod.h"
#include "param_mod.h"

extern double tracer_mod_mp_fw_norm2_;

extern double (*tracer_mod_mp_atb_)[ntra][km+1][jmt][imt]; //注意，fortran是0:km

extern double (*tracer_mod_mp_at_)[ntra][km][jmt][imt];

extern double (*tracer_mod_mp_amld_)[jmt][imt];
extern double (*tracer_mod_mp_licomqice_)[jmt][imt];
extern double tracer_mod_mp_penetrate_[max_blocks_clinic][km][jmt][imt];
extern double tracer_mod_mp_pdensity_[max_blocks_clinic][km+1][jmt][imt];

extern double (*tracer_mod_mp_net_)[ntra][jmt][imt];
extern double tracer_mod_mp_ax_[max_blocks_clinic][ntra][km][jmt][imt];
extern double tracer_mod_mp_ay_[max_blocks_clinic][ntra][km][jmt][imt];
extern double tracer_mod_mp_az_[max_blocks_clinic][ntra][km][jmt][imt];
extern double tracer_mod_mp_dx_[max_blocks_clinic][ntra][km][jmt][imt];
extern double tracer_mod_mp_dy_[max_blocks_clinic][ntra][km][jmt][imt];
extern double tracer_mod_mp_dz_[max_blocks_clinic][ntra][km][jmt][imt];
extern double tracer_mod_mp_ddy_[max_blocks_clinic][ntra][km][jmt][imt];
extern double tracer_mod_mp_tend_[max_blocks_clinic][ntra][km][jmt][imt];
extern double tracer_mod_mp_dt_diff_[max_blocks_clinic][ntra][km][jmt][imt];
extern double tracer_mod_mp_dt_conv_[max_blocks_clinic][ntra][km][jmt][imt];
extern double tracer_mod_mp_dt_restore_[max_blocks_clinic][ntra][km][jmt][imt];

#ifdef ISO
    extern double tracer_mod_mp_aay_iso_[max_blocks_clinic][ntra][km][jmt][imt];
    extern double tracer_mod_mp_ddy_iso_[max_blocks_clinic][ntra][km][jmt][imt];
    extern double tracer_mod_mp_ax_iso_[max_blocks_clinic][ntra][km][jmt][imt];
    extern double tracer_mod_mp_ay_iso_[max_blocks_clinic][ntra][km][jmt][imt];
    extern double tracer_mod_mp_az_iso_[max_blocks_clinic][ntra][km][jmt][imt];
    extern double tracer_mod_mp_dx_iso_[max_blocks_clinic][ntra][km][jmt][imt];
    extern double tracer_mod_mp_dy_iso_[max_blocks_clinic][ntra][km][jmt][imt];
    extern double tracer_mod_mp_dz_iso_[max_blocks_clinic][ntra][km][jmt][imt];
#endif

#endif // !INCLUDE_TRACER
