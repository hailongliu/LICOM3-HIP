#ifndef INCLUDE_HMIX_DEL2
#define INCLUDE_HMIX_DEL2

#include<string.h>
#include "pconst_mod.h"
#include "precision_mod.h"
#include "param_mod.h"
#include "constant_mod.h"
#include "grid.h"
#include "domain.h"

extern double hmix_del2_mp_am_;
extern double hmix_del2_mp_ah_;


extern double (*hmix_del2_mp_dtn_)[ny_block][nx_block];
extern double (*hmix_del2_mp_dts_)[ny_block][nx_block];
extern double (*hmix_del2_mp_dte_)[ny_block][nx_block];
extern double (*hmix_del2_mp_dtw_)[ny_block][nx_block];
extern double (*hmix_del2_mp_duc_)[ny_block][nx_block];
extern double (*hmix_del2_mp_dum_)[ny_block][nx_block];
extern double (*hmix_del2_mp_dun_)[ny_block][nx_block];
extern double (*hmix_del2_mp_dus_)[ny_block][nx_block];
extern double (*hmix_del2_mp_due_)[ny_block][nx_block];
extern double (*hmix_del2_mp_duw_)[ny_block][nx_block];
extern double (*hmix_del2_mp_dmc_)[ny_block][nx_block];
extern double (*hmix_del2_mp_dmn_)[ny_block][nx_block];
extern double (*hmix_del2_mp_dms_)[ny_block][nx_block];
extern double (*hmix_del2_mp_dme_)[ny_block][nx_block];
extern double (*hmix_del2_mp_dmw_)[ny_block][nx_block];

//extern void hmix_del2_mp_hdiffu_del2_(int *, double[][imt], double[][imt], double[][imt], double[][imt], struct block *);
//extern void hmix_del2_mp_hdifft_del2_(int *, double[][imt], double[][imt], struct block *);

//extern void hmix_del2_hdiffu_del2_(int *, double[ny_block][nx_block], double[ny_block][nx_block], double[ny_block][nx_block], double[ny_block][nx_block], int *, int *, int *, int *, int *);
//extern void hmix_del2_hdifft_del2_(int *, double[ny_block][nx_block], double[ny_block][nx_block], int *, int *, int *, int *, int *);

#endif // !INCLUDE_HMIX_DEL2
