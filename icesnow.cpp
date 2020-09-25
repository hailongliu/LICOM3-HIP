//  CVS: $Id: icesnow.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
//     ==================
#include "hip/hip_runtime.h"
#include "cuda_data.h"

#include <def-undef.h>
#include "param_mod.h"
#include "pconst_mod.h"
#include "shr_const_mod.h"
/*
__device__ double min(double x, double y) {
    return (x) < (y) ? (x) : (y);
}*/

__global__ void icesnow_cu(double *d_at, int *d_kmt, double *d_licomqice, double *d_dzp,
                           double d_tbice, double d_cp, double heat_ice_fusion) {

    double sal_ice, sal_ocn, tdiff;
    int i, j, k;
//    __shared__ double s_tbice;
//	if(threadIdx.x==0&&threadIdx.y==0){        
//		s_tbice= d_tbice;
  //  	}	
    i = (blockIdx.x) * blockDim.x + threadIdx.x;
    j = (blockIdx.y) * blockDim.y + threadIdx.y;

    sal_ocn = 34.7e0;                    //Reference salinity for ocean
    sal_ice = 4.0e0;                     //Reference salinity for sea ice
    tdiff = 0.0e0;

    if (j < jmt && i < imt) {
//wpf 0309, in sd do k=1,1; in old COUP version 1,15
        for (k = 0; k < 1; k++) {
            if (d_at[k * jmt * imt + j * imt + i] < d_tbice && d_kmt[j * imt + i] > 0) {
                tdiff = d_tbice - d_at[k * jmt * imt + j * imt + i];

                d_licomqice[j * imt + i] += tdiff * d_dzp[k] / d_dzp[0];

                d_at[1 * km * jmt * imt + k * jmt * imt + j * imt + i] +=
                        tdiff * (sal_ocn - sal_ice) * d_cp / heat_ice_fusion * 0.001e0;

                d_at[k * jmt * imt + j * imt + i] = d_tbice;
            }
        }
    }

    if (j < jmt && i < imt) {
        if (d_licomqice[j * imt + i] > 0 && d_at[j * imt + i] > d_tbice) {
            tdiff = fmin((d_at[j * imt + i] - d_tbice), d_licomqice[j * imt + i]);

            d_licomqice[j * imt + i] -= tdiff;

            d_at[j * imt + i] -= tdiff;

            d_at[1 * km * jmt * imt + j * imt + i] -=
                    tdiff * (sal_ocn - sal_ice) * d_cp / heat_ice_fusion *
                    0.001e0;
        }
    }

    for (k = 0; k < km; k++) {
        if (j < jmt && i < imt) {
            if (d_at[k * jmt * imt + j * imt + i] < d_tbice) {
                d_at[k * jmt * imt + j * imt + i] = d_tbice;
            }
        }
    }
}

extern "C" void icesnow_() {

    dim3 threadsPerBlock(8, 8, 1);
    dim3 numBlocks((imt + threadsPerBlock.x - 1) / threadsPerBlock.x,
                   (jmt + threadsPerBlock.y - 1) / threadsPerBlock.y);

#ifdef FRC_CORE
           CHECK(hipMemset(d_licomqice, 0.0e0, dataSize2d))
#endif
    hipLaunchKernelGGL(icesnow_cu, dim3(numBlocks), dim3(threadsPerBlock ), 0, 0, d_at, d_kmt, d_licomqice, d_dzp,
                                       pconst_mod_mp_tbice_, pconst_mod_mp_cp_, SHR_CONST_LATICE);


    return;

}
