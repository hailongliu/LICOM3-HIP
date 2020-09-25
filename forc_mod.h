#ifndef INCLUDE_FORC
#define INCLUDE_FORC

#include "precision_mod.h"
#include "param_mod.h"

extern double (*forc_mod_mp_ssf_)[jmt][imt];
extern double (*forc_mod_mp_sss_)[jmt][imt];
extern double (*forc_mod_mp_su_)[jmt][imt];
extern double (*forc_mod_mp_sv_)[jmt][imt];
extern double (*forc_mod_mp_ustar_)[jmt][imt];
extern double (*forc_mod_mp_psa_)[jmt][imt];
extern double forc_mod_mp_buoytur_[max_blocks_clinic][jmt][imt];
extern double forc_mod_mp_buoysol_[max_blocks_clinic][jmt][imt];
extern double forc_mod_mp_sst_[max_blocks_clinic][jmt][imt];
extern double forc_mod_mp_dqdt_[max_blocks_clinic][jmt][imt];
extern double forc_mod_mp_seaice_[max_blocks_clinic][jmt][imt];
extern double forc_mod_mp_fresh_[max_blocks_clinic][jmt][imt];
extern double (*forc_mod_mp_nswv_)[jmt][imt];
extern double (*forc_mod_mp_swv_)[jmt][imt];
extern double (*forc_mod_mp_wave_dis_)[jmt][imt];
extern double (*forc_mod_mp_tsf_)[jmt][imt];
extern double (*forc_mod_mp_restore_)[ntra][km][jmt][imt];
#endif // !INCLUDE_FORC
