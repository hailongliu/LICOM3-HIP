#ifndef IMCLUDE_PMIX
#define IMCLUDE_PMIX

#include "precision_mod.h"
#include "param_mod.h"

extern double pmix_mod_mp_pen_[kmm1];
extern double pmix_mod_mp_pen_chl_[max_blocks_clinic][km][jmt][imt];

extern double pmix_mod_mp_s2t_[max_blocks_clinic][kmm1][jmt][imt];
extern double pmix_mod_mp_s2u_[max_blocks_clinic][kmm1][jmt][imt];
extern double pmix_mod_mp_ricdttms_[max_blocks_clinic][kmm1][jmt][imt];
extern double pmix_mod_mp_ricdt_[max_blocks_clinic][kmm1][jmt][imt];
extern double pmix_mod_mp_ridt_[max_blocks_clinic][kmm1][jmt][imt];
extern double (*pmix_mod_mp_ric_)[kmm1][jmt][imt];
extern double (*pmix_mod_mp_rict_)[kmm1][jmt][imt];
extern double (*pmix_mod_mp_rict_ref_)[jmt][imt];
extern double (*pmix_mod_mp_rict_replace_)[kmm1][jmt][imt];
extern double (* pmix_mod_mp_rit_)[kmm1][jmt][imt];

extern double (* pmix_mod_mp_riu_)[km+1][jmt][imt];

#endif // !IMCLUDE_PMIX
