#ifndef INCLUDE_HMIX_DEL4
#define INCLUDE_HMIX_DEL4

#include<string.h>
#include "pconst_mod.h"
#include "precision_mod.h"
#include "param_mod.h"
#include "constant_mod.h"
#include "grid.h"
#include "domain.h"


extern double hmix_del4_mp_am_;
extern double hmix_del4_mp_ah_;

extern double (*hmix_del4_mp_dtn_)[ny_block][nx_block];
extern double (*hmix_del4_mp_dts_)[ny_block][nx_block];
extern double (*hmix_del4_mp_dte_)[ny_block][nx_block];
extern double (*hmix_del4_mp_dtw_)[ny_block][nx_block];
extern double (*hmix_del4_mp_duc_)[ny_block][nx_block];
extern double (*hmix_del4_mp_dum_)[ny_block][nx_block];
extern double (*hmix_del4_mp_dun_)[ny_block][nx_block];
extern double (*hmix_del4_mp_dus_)[ny_block][nx_block];
extern double (*hmix_del4_mp_due_)[ny_block][nx_block];
extern double (*hmix_del4_mp_duw_)[ny_block][nx_block];
extern double (*hmix_del4_mp_dmc_)[ny_block][nx_block];
extern double (*hmix_del4_mp_dmn_)[ny_block][nx_block];
extern double (*hmix_del4_mp_dms_)[ny_block][nx_block];
extern double (*hmix_del4_mp_dme_)[ny_block][nx_block];
extern double (*hmix_del4_mp_dmw_)[ny_block][nx_block];
extern double (*hmix_del4_mp_amf_)[ny_block][nx_block];
extern double (*hmix_del4_mp_ahf_)[ny_block][nx_block];


//extern void hmix_del4_mp_hdiffu_del4_(int *, double[][imt], double[][imt], double[][imt], double[][imt], struct block *);
//extern void hmix_del4_mp_hdifft_del4_(int *, double[][imt], double[][imt], struct block *);

extern void hmix_del4_hdiffu_del4_(int *, double[][imt], double[][imt], double[][imt], double[][imt],int *, int *, int *, int *, int *);
extern void hmix_del4_hdifft_del4_(int *, double[][imt], double[][imt], int *, int *, int *, int *, int *);

#endif // !INCLUDE_HMIX_DEL4
