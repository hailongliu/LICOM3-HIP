#ifndef INCLUDE_DOMAIN
#define INCLUDE_DOMAIN

#include "precision_mod.h"
#include "param_mod.h"
#include "POP_HaloMod.h"
#include "POP_GridHorzMod.h"
#include "LICOM_Error_mod.h"
#include "blocks.h"
#include "distribution.h"

//#define domain_mp_nblocks_clinic_ 1
extern int domain_mp_nblocks_clinic_;
extern int domain_mp_errorcode_ ;
extern int *domain_mp_blocks_clinic_;

extern struct distrb domain_mp_distrb_clinic_;
//extern struct pop_halo domain_mp_pop_haloclinic_;


#endif // !INCLUDE_DOMAIN

