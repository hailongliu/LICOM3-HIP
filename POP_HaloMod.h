#ifndef INCLUDE_POP_HALOMOD
#define INCLUDE_POP_HALOMOD

#include "precision_mod.h"
#include "param_mod.h"

//extern int POP_HaloMod_mp_bufSizeSend_, POP_HaloMod_mp_bufSizeRecv_;

// struct pop_halo
// {
// 	int communicator, nummsgsend, nummsgrecv, numlocalcopies;
// 	int *recvtask, *sendtask, *sizesend, *sizerecv;
// 	int **srclocaladdr, **dstlocaladdr;
// 	int ***sendaddr, ***recvaddr;
// };

// extern void pop_halomod_mp_pop_haloupdate2dr4_(double[][jmt][imt], struct pop_halo *, char [], char [], int *, double *);
// extern void pop_halomod_mp_pop_haloupdate2dr8_(double[][jmt][imt], struct pop_halo *, char[], char[], int *, double *);
// extern void pop_halomod_mp_pop_haloupdate2di4_(double[][jmt][imt], struct pop_halo *, char[], char[], int *, int *);
// extern void pop_halomod_mp_pop_haloupdate3dr4_(double[][km][jmt][imt], struct pop_halo *, char[], char[], int *, double *);
// extern void pop_halomod_mp_pop_haloupdate3dr8_(double[][ntra][jmt][imt], struct pop_halo *, char[], char[], int *, double *);  //用fortran包装

// extern void pop_halomod_mp_pop_haloupdate3dr8_(double[][km][jmt][imt], struct pop_halo *, char[], char[], int *, double *);
// extern void pop_halomod_mp_pop_haloupdate3di4_(double[][km][jmt][imt], struct pop_halo *, char[], char[], int *, int *);
// extern void pop_halomod_mp_pop_haloupdate4dr4_(double[][ntra[[km][jmt][imt], struct pop_halo *, char[], char[], int *, double *);
// extern void pop_halomod_mp_pop_haloupdate4dr8_(double[][ntra][km][jmt][imt], struct pop_halo *, char[], char[], int *, double *);
// extern void pop_halomod_mp_pop_haloupdate4di4_(double[][ntra][km][jmt][imt], struct pop_halo *, char[], char[], int *, int *);

#endif // !INCLUDE_POP_HALOMOD
