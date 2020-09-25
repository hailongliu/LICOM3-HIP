#ifndef INCLUDE_ADVECTION
#define INCLUDE_ADVECTION

#include "precision_mod.h"
#include "param_mod.h"
#include "blocks.h"

extern void advection_mp_advection_momentum_ (double[][jmt][imt], double[][jmt][imt], double[][jmt][imt], double[][jmt][imt], double[][jmt][imt], int *) ;
extern void advection_mp_advection_tracer_ (double[][jmt][imt], double[][jmt][imt], double[][jmt][imt], double[][jmt][imt], double[][jmt][imt], int *, int*, struct block *) ;

extern void advection_advection_momentum_ (double[][jmt][imt], double[][jmt][imt], double[][jmt][imt], double[][jmt][imt], double[][jmt][imt], int *) ;
extern void advection_advection_tracer_ (double[][jmt][imt], double[][jmt][imt], double[][jmt][imt], double[][jmt][imt], double[][jmt][imt], int *, int*) ;


#endif // !INCLUDE_ADVECTION
