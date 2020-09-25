#ifndef INCLUDE_GATHER_SCATTER
#define INCLUDE_GATHER_SCATTER

#include "precision_mod.h"
#include "param_mod.h"
#include "blocks.h"
#include "distribution.h"
#include "constant_mod.h"
#include "LICOM_Error_mod.h"

extern void gather_scatter_mp_gather_global_dbl_(double[][imt_global], double[][jmt][imt], int *,struct distrb *);

#endif // !INCLUDE_GATHER_SCATTER
