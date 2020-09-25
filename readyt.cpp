// by Wpf, 2019.1.13
// Add icc compile option "-fp-model precise" to get same result 

#include "hip/hip_runtime.h"
#include <math.h>
#include <string.h>
#include <stdio.h>
#include<sys/timeb.h>
#include "cuda_data.h"
#include "precision_mod.h"
#include "param_mod.h"
#include "pconst_mod.h"
#include "dyn_mod.h"
#include "tracer_mod.h"
#include "work_mod.h"
#include "pmix_mod.h"
#include "forc_mod.h"
#include "domain.h"
#include "grid.h"
#include "blocks.h"
#include "constant_mod.h"
#include "common.h"

extern "C" void allocate_readyt_();
__device__ double dens_c_cu(double *tq, double *sq, int *kk, double *d_c) {
    double outcome;
    int k_tmpt = *kk;

    outcome = (d_c[k_tmpt] + (d_c[3 * km + k_tmpt] + d_c[6 * km + k_tmpt] * (*sq)) * (*sq) +
               (d_c[2 * km + k_tmpt] + d_c[7 * km + k_tmpt] * (*sq) + d_c[5 * km + k_tmpt] * (*tq)) *
               (*tq)) * (*tq) + (d_c[1 * km + k_tmpt] + (d_c[4 * km + k_tmpt] + d_c[8 * km + k_tmpt] * (*sq)) *
                                                        (*sq)) * (*sq);

    return outcome;
}

__global__ void dens_cu(double *d_dlu, double *d_dlv, double *d_to, double *d_so, double *d_ric, double *d_viv,
                        double *d_odzt, double *d_c, double d_od0, const double d_g, 
			double *d_rict, double *d_vit, double *d_atb, int d_km,
			int d_kmm1, int d_jmt, int d_imt) {
    int i, j, k;

    i = blockIdx.x * blockDim.x + threadIdx.x;
    j = blockIdx.y * blockDim.y + threadIdx.y;
    k = blockIdx.z * blockDim.z + threadIdx.z;

    if (k < d_kmm1 && j < d_jmt - 1 && 1 <= i && i < d_imt) {
        double tup, sup, tlo, slo;
        double rholo, rhoup;
        int kp1;

        tup = d_dlu[k * d_jmt * d_imt + j * d_imt + i] - d_to[k + 1];
        sup = d_dlv[k * d_jmt * d_imt + j * d_imt + i] - d_so[k + 1];
        tlo = d_dlu[(k + 1) * d_jmt * d_imt + j * d_imt + i] - d_to[k + 1];
        slo = d_dlv[(k + 1) * d_jmt * d_imt + j * d_imt + i] - d_so[k + 1];
        kp1 = k + 1;

        rhoup = dens_c_cu(&tup, &sup, &kp1, d_c); //---------------funcion 2
        rholo = dens_c_cu(&tlo, &slo, &kp1, d_c);//---------------funcion 2

        d_ric[k * d_jmt * d_imt + j * d_imt + i] = d_viv[(k + 1) * d_jmt * d_imt + j * d_imt + i] *
                                                   d_od0 * d_g * (rholo - rhoup) * d_odzt[k + 1];

        tup = d_atb[(k + 1) * d_jmt * d_imt + j * d_imt + i] - d_to[k + 1];
        sup = d_atb[1 * (d_km + 1) * d_jmt * d_imt + (k + 1) * d_jmt * d_imt + j * d_imt + i] - d_so[k + 1];
        tlo = d_atb[(k + 2) * d_jmt * d_imt + j * d_imt + i] - d_to[k + 1];
        slo = d_atb[1 * (d_km + 1) * d_jmt * d_imt + (k + 2) * d_jmt * d_imt + j * d_imt + i] - d_so[k + 1];

        kp1 = k + 1;
        rhoup = dens_c_cu(&tup, &sup, &kp1, d_c);//---------------funcion 2
        rholo = dens_c_cu(&tlo, &slo, &kp1, d_c);//---------------funcion 2

        d_rict[k * d_jmt * d_imt + j * d_imt + i] = d_vit[(k + 1) * d_jmt * d_imt + j * d_imt + i] *
                                                    d_od0 * d_g * (rholo - rhoup) * d_odzt[k + 1];
    }
}

__global__ void dens_cu2(double *d_to, double *d_so, double *d_rict, double *d_vit, double *d_odzt, double *d_c,
                         double *d_atb, double d_od0, const double d_g, int d_kmm1, int d_km, int d_jmt, int d_imt) {

    int i, j, k;

    i = blockIdx.x * blockDim.x + threadIdx.x;
    j = blockIdx.y * blockDim.y + threadIdx.y;
    k = blockIdx.z * blockDim.z + threadIdx.z;

    if (k < d_kmm1 && j < d_jmt && i < d_imt) {
        double tup, sup, tlo, slo;
        double rholo, rhoup;
        int kp1;

        tup = d_atb[(k + 1) * d_jmt * d_imt + j * d_imt + i] - d_to[k + 1];
        sup = d_atb[1 * (d_km + 1) * d_jmt * d_imt + (k + 1) * d_jmt * d_imt + j * d_imt + i] - d_so[k + 1];
        tlo = d_atb[(k + 2) * d_jmt * d_imt + j * d_imt + i] - d_to[k + 1];
        slo = d_atb[1 * (d_km + 1) * d_jmt * d_imt + (k + 2) * d_jmt * d_imt + j * d_imt + i] - d_so[k + 1];

        kp1 = k + 1;
        rhoup = dens_c_cu(&tup, &sup, &kp1, d_c);//---------------funcion 2
        rholo = dens_c_cu(&tlo, &slo, &kp1, d_c);//---------------funcion 2

        d_rict[k * d_jmt * d_imt + j * d_imt + i] = d_vit[(k + 1) * d_jmt * d_imt + j * d_imt + i] *
                                                    d_od0 * d_g * (rholo - rhoup) * d_odzt[k + 1];
    }
}

__global__ void readyt_cu_1(double *d_h0l, double *d_h0f, double *d_h0, double *d_utl, double *d_vtl, double *d_utf, double *d_vtf, double *d_u, double *d_v,
                            double *d_dlu, double *d_dlv, double *d_atb,
                            double *d_au0, double *d_aus, double *d_auw, double *d_ausw) {
    int i, j, k, n;

    i = blockIdx.x * blockDim.x + threadIdx.x;
    j = blockIdx.y * blockDim.y + threadIdx.y;
    k = blockIdx.z * blockDim.z + threadIdx.z;

        if(k==0 && jst-1<=j && j<jet && 0<= i && i<imt ) {
			
            d_h0l[j*imt +i]= d_h0f[j*imt +i];
            d_h0f[j*imt +i]= d_h0[j*imt +i];
        }
    

            if(0<=k && k<km && jst-1<=j && j<jet && 0<= i && i<imt) {
				
                     d_utl[k*jmt*imt +j*imt +i]= d_utf[k*jmt*imt +j*imt +i];
                     d_vtl[k*jmt*imt +j*imt +i]= d_vtf[k*jmt*imt +j*imt +i];
                     d_utf[k*jmt*imt +j*imt +i]= d_u[k*jmt*imt +j*imt +i];
                     d_vtf[k*jmt*imt +j*imt +i]= d_v[k*jmt*imt +j*imt +i];
            }      

    if (k < km && j < jmt && i < imt) {
        if (j < (jmt - 1) && i >= 1) {
            d_dlu[k * jmt * imt + j * imt + i] =
                    d_au0[j * imt + i] *
                    d_atb[(k + 1) * jmt * imt + j * imt + i] +
                    d_aus[j * imt + i] *
                    d_atb[(k + 1) * jmt * imt + (j + 1) * imt + i] +
                    d_auw[j * imt + i] *
                    d_atb[(k + 1) * jmt * imt + j * imt + i - 1] +
                    d_ausw[j * imt + i] *
                    d_atb[(k + 1) * jmt * imt + (j + 1) * imt + i - 1];

            d_dlv[k * jmt * imt + j * imt + i] =
                    d_au0[j * imt + i] *
                    d_atb[1 * (km + 1) * jmt * imt + (k + 1) * jmt * imt + j * imt + i] +
                    d_aus[j * imt + i] *
                    d_atb[1 * (km + 1) * jmt * imt + (k + 1) * jmt * imt + (j + 1) * imt + i] +
                    d_auw[j * imt + i] *
                    d_atb[1 * (km + 1) * jmt * imt + (k + 1) * jmt * imt + j * imt + i - 1] +
                    d_ausw[j * imt + i] *
                    d_atb[1 * (km + 1) * jmt * imt + (k + 1) * jmt * imt + (j + 1) * imt + i - 1];
        } else {
            d_dlu[k * jmt * imt + j * imt + i] = (double) 0.0e0;
            d_dlv[k * jmt * imt + j * imt + i] = (double) 0.0e0;
        }
    }

    return;
}

__global__ void readyt_cu_2(double *d_to, double *d_so, double *d_c, double *d_po,
                            double *d_atb, double *d_rict, double *d_rict_replace, double *d_gg, double *d_pdensity, double *d_vit,
                            double *d_pp, double *d_ppa, double *d_ppb, double *d_ppc, double *d_psa, double *d_at,
                            double *d_dzp, double d_od0, const double d_g, int d_kmm1, int d_km, int d_jmt, int d_imt) {
    int i, j, k;

    i = blockIdx.x * blockDim.x + threadIdx.x;
    j = blockIdx.y * blockDim.y + threadIdx.y;
    k = blockIdx.z * blockDim.z + threadIdx.z;

    if (k < d_km && j < d_jmt && i < d_imt) {
        if (d_vit[k * d_jmt * d_imt + j * d_imt + i] > 0.0) {
            double tq, sq;

//            tq = d_atb[(k + 1) * d_jmt * d_imt + j * d_imt + i] - d_to[k];
            tq = d_at[k  * d_jmt * d_imt + j * d_imt + i] - d_to[k];
            sq = d_at[1 * d_km  * d_jmt * d_imt + k * d_jmt * d_imt + j * d_imt + i] - d_so[k];
//            sq = d_atb[1 * (d_km + 1) * d_jmt * d_imt + (k + 1) * d_jmt * d_imt + j * d_imt + i] - d_so[k];

            d_pdensity[k * d_jmt * d_imt + j * d_imt + i] = 1.0e+3 + d_po[k] +
                                                            (d_c[0 * d_km + k] +
                                                             (d_c[3 * d_km + k] + d_c[6 * d_km + k] * sq) * sq +
                                                             (d_c[2 * d_km + k] + d_c[7 * d_km + k] * sq +
                                                              d_c[5 * d_km + k] * tq) * tq) * tq +
                                                            (d_c[1 * d_km + k] +
                                                             (d_c[4 * d_km + k] + d_c[8 * d_km + k] * sq) * sq) * sq;
        } else {
            d_pdensity[k * d_jmt * d_imt + j * d_imt + i] = 0.0;
        }

        d_gg[k * d_jmt * d_imt + j * d_imt + i] =
                -d_od0 * d_g *
                d_pdensity[k * d_jmt * d_imt + j * d_imt + i] *
                d_vit[k * d_jmt * d_imt + j * d_imt + i];

    if (k >= 1) {
        d_ppb[k * d_jmt * d_imt + j * d_imt + i] =
                d_vit[k * d_jmt * d_imt + j * d_imt + i] *
                (d_at[(k - 1) * d_jmt * d_imt + j * d_imt + i] -
                 (d_at[(k - 1) * d_jmt * d_imt + j * d_imt + i] -
                  d_at[k * d_jmt * d_imt + j * d_imt + i]) *
                 d_dzp[k - 1] / (d_dzp[k - 1] + d_dzp[k]));

        d_ppc[k * d_jmt * d_imt + j * d_imt + i] =
                d_vit[k * d_jmt * d_imt + j * d_imt + i] *
                (d_at[1 * d_km * d_jmt * d_imt + (k - 1) * d_jmt * d_imt + j * d_imt + i] -
                 (d_at[1 * d_km * d_jmt * d_imt + (k - 1) * d_jmt * d_imt + j * d_imt + i] -
                  d_at[1 * d_km * d_jmt * d_imt + k * d_jmt * d_imt + j * d_imt + i]) *
                 d_dzp[k - 1] / (d_dzp[k - 1] + d_dzp[k]));
    }

    }

    return;
}

__global__ void readyt_cu_2_2(double *d_psa, double *d_ppb, double *d_ppc, double *d_at, double *d_gg, double *d_vit, double *d_pp, double *d_ppa, double *d_dzp, int d_km,
                              int d_jmt, int d_imt) {
    int i, j, k;

    i = blockIdx.x * blockDim.x + threadIdx.x;
    j = blockIdx.y * blockDim.y + threadIdx.y;

    if (j < d_jmt && i < d_imt) {
            d_pp[j * d_imt + i] = 0.5e0 * d_dzp[0] *
                                  d_gg[j * d_imt + i] *
                                  d_vit[j * d_imt + i];

            d_ppa[j * d_imt + i] = d_psa[j * d_imt + i] *
                                   d_vit[j * d_imt + i];

            d_ppb[j * d_imt + i] = d_at[j * d_imt + i] *
                                   d_vit[j * d_imt + i];

            d_ppc[j * d_imt + i] = d_at[1 * d_km * d_jmt * d_imt + j * d_imt + i] *
                                   d_vit[j * d_imt + i];
//             d_h0[j * d_imt + i]=d_pp[j * d_imt + i];
        for (k = 1; k < d_km; k++) {
            d_pp[k * d_jmt * d_imt + j * d_imt + i] =
                    d_vit[k * d_jmt * d_imt + j * d_imt + i] *
                    (d_pp[(k - 1) * d_jmt * d_imt + j * d_imt + i] +
                     0.5e0 * (d_gg[k * d_jmt * d_imt + j * d_imt + i] *
                              d_dzp[k] + d_gg[(k - 1) * d_jmt * d_imt + j * d_imt + i] *
                                         d_dzp[k - 1]));

            d_ppa[k * d_jmt * d_imt + j * d_imt + i] =
                    d_vit[k * d_jmt * d_imt + j * d_imt + i] *
                    (d_ppa[(k - 1) * d_jmt * d_imt + j * d_imt + i] +
                     d_gg[(k - 1) * d_jmt * d_imt + j * d_imt + i] * d_dzp[k - 1]);
        }
    }

    return;
}

__global__ void thermal_cu(double *d_ppb, double *d_ppc, double *d_ppa, double *d_alpha,
                           double *d_beta, double *d_vit, double d_od0, int d_imt, int d_jmt, int d_km) {
    int i, j, k;

    i = (blockIdx.x) * blockDim.x + threadIdx.x;
    j = (blockIdx.y) * blockDim.y + threadIdx.y;
    k = blockIdx.z * blockDim.z + threadIdx.z;

    if (k < d_km && i < d_imt - 1 && i>0 && j < d_jmt) {
        double tt1, tt2, tt3, ss1, ss2, tm1, tm2, tm3;

        tt1 = d_ppb[k * d_imt * d_jmt + j * d_imt + i];
        tt2 = tt1 * tt1;
        tt3 = tt2 * tt1;

        ss1 = d_ppc[k * d_imt * d_jmt + j * d_imt + i] * 1000.0e0;
        ss2 = ss1 * ss1;

        tm1 = -d_ppa[k * d_imt * d_jmt + j * d_imt + i] / d_od0 / 10000.e0 *
              d_vit[k * d_imt * d_jmt + j * d_imt + i];
        tm2 = tm1 * tm1;
        tm3 = tm2 * tm1;

        d_beta[k * d_imt * d_jmt + j * d_imt + i] =
                (0.785567e-3 - 0.301985e-5 * tt1 + 0.555579e-7 * tt2 - 0.415613e-9 * tt3 +
                 ss1 * (-0.356603e-6 + 0.788212e-8 * tt1 + 0.408195e-10 * tm1 - 0.602281e-15 * tm2) +
                 ss2 * (0.515032e-8) + tm1 * (-0.121555e-7 + 0.192867e-9 * tt1 - 0.2131127e-11 * tt2) +
                 tm2 * (0.176621e-12 - 0.175379e-14 * tt1) + tm3 * (0.12155e-17)) *
                d_vit[k * d_imt * d_jmt + j * d_imt + i];

        d_alpha[k * d_imt * d_jmt + j * d_imt + i] =
                (0.665157e-1 + 0.170907e-1 * tt1 - 0.203814e-3 * tt2 + 0.298357e-5 * tt3 - 0.255019e-7 * tt3 * tt1 +
                 ss1 * (0.378110e-2 - 0.846960e-4 * tt1 - 0.164759e-6 * tm1 - 0.251520e-11 * tm2) +
                 ss2 * (-0.678662e-5) + tm1 * (0.380374e-4 - 0.933746e-6 * tt1 + 0.791325e-8 * tt2) +
                 0.512857e-12 * tm2 * tt2 - 0.302285e-13 * tm3) * d_beta[k * d_imt * d_jmt + j * d_imt + i] *
                d_vit[k * d_imt * d_jmt + j * d_imt + i];
    }/*else{
	d_beta[k * d_imt * d_jmt + j * d_imt + i] = 0.0;
	d_alpha[k * d_imt * d_jmt + j * d_imt + i] = 0.0;
    }*/

    return;
}

__global__ void readyt_cu_3(double *d_ricdttms, double *d_ricdt, double *d_gg, double *d_pdensity,
                            double *d_vit, double *d_pp, double *d_at, double *d_po, double *d_dlu, double *d_dlv,
                            double *d_dxur, double *d_dyur, int *d_kmu, double *d_buoytur,
                            double *d_buoysol, double *d_nswv, double *d_swv, double *d_alpha, double *d_beta,
                            double *d_odzt, double d_od0cp, double d_od0,
                            const double d_g, int d_kmm1, int d_km, int d_jmt, int d_imt) {

    int i, j, k;

    i = blockIdx.x * blockDim.x + threadIdx.x;
    j = blockIdx.y * blockDim.y + threadIdx.y;
    k = blockIdx.z * blockDim.z + threadIdx.z;

    if (k < d_km && j < d_jmt && i < d_imt) {
        double epsln = 1.0e-25;

        if (k == 0) {
            d_buoytur[j * d_imt + i] =
                    d_vit[j * d_imt + i] *
                    d_nswv[j * d_imt + i] *
                    d_alpha[j * d_imt + i] *
                    d_g * d_od0cp;

            d_buoysol[j * d_imt + i] =
                    d_vit[j * d_imt + i] *
                    d_swv[j * d_imt + i] *
                    d_alpha[j * d_imt + i] *
                    d_g * d_od0cp;
        }

        if (k < d_kmm1 && 1 <= i && i < d_imt - 1) {
            d_ricdttms[k * d_jmt * d_imt + j * d_imt + i] =
                    d_vit[(k + 1) * d_jmt * d_imt + j * d_imt + i] *
                    d_g *
                    ((d_at[k * d_jmt * d_imt + j * d_imt + i] -
                      d_at[(k + 1) * d_jmt * d_imt + j * d_imt + i]) *
                     d_alpha[(k + 1) * d_jmt * d_imt + j * d_imt + i] +
                     1000.0e0 * (d_at[1 * d_km * d_jmt * d_imt + k * d_jmt * d_imt + j * d_imt + i] -
                                 d_at[1 * d_km * d_jmt * d_imt + (k + 1) * d_jmt * d_imt + j * d_imt + i]) *
                     d_beta[(k + 1) * d_jmt * d_imt + j * d_imt + i]) * d_odzt[k + 1];

            d_ricdt[k * d_jmt * d_imt + j * d_imt + i] =
                    d_vit[(k + 1) * d_jmt * d_imt + j * d_imt + i] /
                    ((d_at[k * d_jmt * d_imt + j * d_imt + i] -
                      d_at[(k + 1) * d_jmt * d_imt + j * d_imt + i] + epsln) *
                     d_alpha[(k + 1) * d_jmt * d_imt + j * d_imt + i]) * 1000.0e0 *
                    (d_at[1 * d_km * d_jmt * d_imt + k * d_jmt * d_imt + j * d_imt + i] -
                     d_at[1 * d_km * d_jmt * d_imt + (k + 1) * d_jmt * d_imt + j * d_imt + i]) *
                    d_beta[(k + 1) * d_jmt * d_imt + j * d_imt + i];
        }
                #if ( defined CANUTOMIXOUT )
                        pconst_mod_mp_alpha_canuto_[iblock][k][j][i]=alpha[iblock][k][j][i];
                        pconst_mod_mp_beta_canuto_[iblock][k][j][i] =beta[iblock][k][j][i];
                 #endif


        d_gg[k * d_jmt * d_imt + j * d_imt + i] =
                -d_od0 * d_g *
                (d_pdensity[k * d_jmt * d_imt + j * d_imt + i] - d_po[k] - 1000.0e0) *
                d_vit[k * d_jmt * d_imt + j * d_imt + i];


        d_dlu[k * d_jmt * d_imt + j * d_imt + i] = 0.0;//c0
        d_dlv[k * d_jmt * d_imt + j * d_imt + i] = 0.0;//c0

        if ((1 <= i) && (j < (d_jmt - 1))) {
            if (k <= (d_kmu[j * d_imt + i]-1) ) { //jjr bug
                d_dlu[k * d_jmt * d_imt + j * d_imt + i] =
                        d_dxur[j * d_imt + i] * p5 *
                        (d_pp[k * d_jmt * d_imt + (j + 1) * d_imt + i] -
                         d_pp[k * d_jmt * d_imt + j * d_imt + i - 1] -
                         d_pp[k * d_jmt * d_imt + (j + 1) * d_imt + i - 1] +
                         d_pp[k * d_jmt * d_imt + j * d_imt + i]);

                d_dlv[k * d_jmt * d_imt + j * d_imt + i] =
                        d_dyur[j * d_imt + i] * p5 *
                        (d_pp[k * d_jmt * d_imt + (j + 1) * d_imt + i] -
                         d_pp[k * d_jmt * d_imt + j * d_imt + i - 1] +
                         d_pp[k * d_jmt * d_imt + (j + 1) * d_imt + i - 1] -
                         d_pp[k * d_jmt * d_imt + j * d_imt + i]);
            }
        }
//      if(k==0 )   d_h0[j * d_imt + i]=d_dlu[j * d_imt + i];
  }

    return;
}

__global__ void readyt_cu_3_2(double *d_gg, double *d_dzp, double *d_dlu, double *d_dlv,
                              double *d_ohbu, double *d_ohbt, double *d_zkt, double *d_viv, double *d_pxb,
                              double *d_pyb, int d_km, int d_jmt, int d_imt) {

    int i, j, k;
    double abcd;

    i = blockIdx.x * blockDim.x + threadIdx.x;
    j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i < d_imt && j < d_jmt) {
        d_pxb[j * d_imt + i] = 0.0;
        d_pyb[j * d_imt + i] = 0.0;

        for (k = 0; k < d_km; k++) {
            d_pxb[j * d_imt + i] += d_dzp[k] *
                                    d_ohbu[j * d_imt + i] *
                                    d_dlu[k * d_jmt * d_imt + j * d_imt + i] *
                                    d_viv[k * d_jmt * d_imt + j * d_imt + i];

            d_pyb[j * d_imt + i] += d_dzp[k] *
                                    d_ohbu[j * d_imt + i] *
                                    d_dlv[k * d_jmt * d_imt + j * d_imt + i] *
                                    d_viv[k * d_jmt * d_imt + j * d_imt + i];

//if ((1 <= i) && (j < (d_jmt - 1)))    d_h0[j * d_imt + i]+= d_dlu[k * d_jmt * d_imt+j * d_imt + i] ;
        }

        d_dlu[j * d_imt + i] = 0.0e0;
        d_dlu[1 * d_jmt * d_imt + j * d_imt + i] = 0.0e0;

        for (k = 0; k < d_km; k++) {
            abcd = d_gg[k * d_jmt * d_imt + j * d_imt + i] *
                   d_ohbt[j * d_imt + i] * d_dzp[k];
            d_dlu[j * d_imt + i] += abcd;
            d_dlu[1 * d_jmt * d_imt + j * d_imt + i] += abcd * d_zkt[k];
        }
    }

    return;
}

__global__ void readyt_cu_4(double *d_dlu, double *d_dlv, double *d_ohbt, const double d_g, int d_jmt, int d_imt) {

    int i, j;
    i = blockIdx.x * blockDim.x + threadIdx.x;
    j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i < d_imt && j < d_jmt) {
        d_dlv[j * d_imt + i] =
                (d_dlu[j * d_imt + i] +
                 d_dlu[1 * d_jmt * d_imt + j * d_imt + i] *
                 d_ohbt[j * d_imt + i]) / d_g;

        d_dlv[1 * d_jmt * d_imt + j * d_imt + i] =
                d_dlu[1 * d_jmt * d_imt + j * d_imt + i] *
                d_ohbt[j * d_imt + i] *
                d_ohbt[j * d_imt + i];
    }

    return;
}

__global__ void readyt_cu_5(double *d_dlv, double *d_au0, double *d_aus, double *d_auw, double *d_ausw,
                            double *d_wgp, double *d_work, double *d_hbx, double *d_hby, double *d_viv,
                            double *d_whx, double *d_why, double *d_work1, double *d_work2, double *d_psa,
                            double *d_dxur, double *d_dyur, int *d_kmu, double *d_pax, double *d_pay, double d_od0,
                            int d_jmt, int d_imt) {

    int i, j;
    i = blockIdx.x * blockDim.x + threadIdx.x;
    j = blockIdx.y * blockDim.y + threadIdx.y;

    if (j < d_jmt && i < d_imt) {
        if (j < (d_jmt - 1) && 1 <= i) {
            d_wgp[j * d_imt + i] = d_au0[j * d_imt + i] * d_dlv[j * d_imt + i] +
                                   d_aus[j * d_imt + i] * d_dlv[(j + 1) * d_imt + i] +
                                   d_auw[j * d_imt + i] * d_dlv[j * d_imt + i - 1] +
                                   d_ausw[j * d_imt + i] * d_dlv[(j + 1) * d_imt + i - 1];

            d_work[j * d_imt + i] = d_au0[j * d_imt + i] * d_dlv[1 * d_jmt * d_imt + j * d_imt + i] +
                                    d_aus[j * d_imt + i] * d_dlv[1 * d_jmt * d_imt + (j + 1) * d_imt + i] +
                                    d_auw[j * d_imt + i] * d_dlv[1 * d_jmt * d_imt + j * d_imt + i - 1] +
                                    d_ausw[j * d_imt + i] * d_dlv[1 * d_jmt * d_imt + (j + 1) * d_imt + i - 1];
        } else {
            d_wgp[j * d_imt + i] = (double) 0.0e0;
            d_work[j * d_imt + i] = (double) 0.0e0;
        }

        d_whx[j * d_imt + i] =
                d_hbx[j * d_imt + i] *
                d_work[j * d_imt + i] *
                d_viv[j * d_imt + i];
        d_why[j * d_imt + i] =
                d_hby[j * d_imt + i] *
                d_work[j * d_imt + i] *
                d_viv[j * d_imt + i];

        d_work1[j * d_imt + i] = 0.0;
        d_work2[j * d_imt + i] = 0.0;

        if ((1 <= i) && (j < (d_jmt - 1))) {
            if (0 <= d_kmu[j * d_imt + i] - 1) { //!jjr bug
                d_work1[j * d_imt + i] = d_dxur[j * d_imt + i] * p5 *
                                         (d_psa[(j + 1) * d_imt + i] -
                                          d_psa[j * d_imt + i - 1] -
                                          d_psa[(j + 1) * d_imt + i - 1] +
                                          d_psa[j * d_imt + i]);

                d_work2[j * d_imt + i] = d_dyur[j * d_imt + i] * p5 *
                                         (d_psa[(j + 1) * d_imt + i] -
                                          d_psa[j * d_imt + i - 1] +
                                          d_psa[(j + 1) * d_imt + i - 1] -
                                          d_psa[j * d_imt + i]);
            }
        }

        d_pay[j * d_imt + i] = -d_od0 * d_work2[j * d_imt + i];
        d_pax[j * d_imt + i] = -d_od0 * d_work1[j * d_imt + i];
    }

    return;
}

extern "C" void readyt_() {
    //allocate_readyt_();
    dim3
    threadsPerBlock(8, 8, 1);
    dim3
    numBlocks((imt + threadsPerBlock.x - 1) / threadsPerBlock.x,
              (jmt + threadsPerBlock.y - 1) / threadsPerBlock.y);

    dim3
    threadsPerBlock_km(8, 8, 4);
    dim3
    numBlocks_km((imt + threadsPerBlock_km.x - 1) / threadsPerBlock_km.x,
                 (jmt + threadsPerBlock_km.y - 1) / threadsPerBlock_km.y,
                 (km + threadsPerBlock_km.z - 1) / threadsPerBlock_km.z);

    dim3
    threadsPerBlock_kmm1(8, 8, 4);
    dim3
    numBlocks_kmm1((imt + threadsPerBlock_kmm1.x - 1) / threadsPerBlock_kmm1.x,
                   (jmt + threadsPerBlock_kmm1.y - 1) / threadsPerBlock_kmm1.y,
                   (kmm1 + threadsPerBlock_kmm1.z - 1) / threadsPerBlock_kmm1.z);

//    hipMemset(d_pp, c0, dataSize3d);
//    hipMemset(d_ppa, c0, dataSize3d);
//    hipMemset(d_ppb, c0, dataSize3d);
//    hipMemset(d_ppc, c0, dataSize3d);
    hipMemset(d_alpha, c0, dataSize3d);
    hipMemset(d_beta, c0, dataSize3d);

    hipMemset(d_rit, c0, dataSize3333d);
    hipMemset(d_ric, c0, dataSize3333d);
    hipMemset(d_rict, c0, dataSize3333d);
    hipMemset(d_ricdttms, c0, dataSize3333d);
    hipMemset(d_ricdt, c0, dataSize3333d);
    hipMemset(d_akt, c0, dataSize4d);


            iCheck();
    hipLaunchKernelGGL(readyt_cu_1, dim3(numBlocks_km), dim3(threadsPerBlock_km ), 0, 0, d_h0l, d_h0f, d_h0, d_utl, d_vtl, d_utf, d_vtf, d_u, d_v,
    d_dlu, d_dlv, d_atb, d_au0, d_aus, d_auw, d_ausw);
            iCheck();

    hipLaunchKernelGGL(dens_cu, dim3(numBlocks_kmm1), dim3(threadsPerBlock_kmm1 ), 0, 0, d_dlu, d_dlv, d_to, d_so, d_ric,
            d_viv, d_odzt, d_c, pconst_mod_mp_od0_, g, d_rict, d_vit, d_atb, km,
	    kmm1, jmt, imt);
            iCheck();

//    hipLaunchKernelGGL(dens_cu2, dim3(numBlocks_kmm1), dim3(threadsPerBlock_kmm1 ), 0, 0, d_to, d_so, d_rict, d_vit, d_odzt, d_c, d_atb,
//            pconst_mod_mp_od0_, g, kmm1, km, jmt, imt);

//    hipMemcpy(d_rict_replace, d_rict, dataSize3333d, hipMemcpyDeviceToDevice);
    hipLaunchKernelGGL(readyt_cu_2, dim3(numBlocks_km), dim3(threadsPerBlock_km ), 0, 0, d_to, d_so, d_c, d_po, d_atb, d_rict, d_rict_replace, d_gg, d_pdensity, d_vit,
            d_pp, d_ppa, d_ppb, d_ppc, d_psa, d_at, d_dzp, pconst_mod_mp_od0_, g, kmm1, km, jmt, imt);
            iCheck();

    hipLaunchKernelGGL(readyt_cu_2_2, dim3(numBlocks), dim3(threadsPerBlock ), 0, 0, d_psa, d_ppb, d_ppc, d_at, d_gg, d_vit, d_pp, d_ppa, d_dzp, km, jmt, imt);
            iCheck();

    hipLaunchKernelGGL(thermal_cu, dim3(numBlocks_km), dim3(threadsPerBlock_km ), 0, 0, d_ppb, d_ppc, d_ppa, d_alpha, d_beta, d_vit,
            pconst_mod_mp_od0_, imt, jmt, km);
            iCheck();

    hipLaunchKernelGGL(readyt_cu_3, dim3(numBlocks_km), dim3(threadsPerBlock_km ), 0, 0, d_ricdttms, d_ricdt, d_gg, d_pdensity, d_vit,
            d_pp, d_at, d_po, d_dlu, d_dlv, d_dxur, d_dyur, d_kmu, d_buoytur, d_buoysol, d_nswv, d_swv,
            d_alpha, d_beta, d_odzt, pconst_mod_mp_od0cp_, pconst_mod_mp_od0_, g, kmm1, km, jmt, imt);
            iCheck();

    hipLaunchKernelGGL(readyt_cu_3_2, dim3(numBlocks), dim3(threadsPerBlock ), 0, 0, d_gg, d_dzp,
            d_dlu, d_dlv, d_ohbu, d_ohbt, d_zkt, d_viv, d_pxb, d_pyb, km, jmt, imt);
            iCheck();

    hipLaunchKernelGGL(readyt_cu_4, dim3(numBlocks), dim3(threadsPerBlock ), 0, 0, d_dlu, d_dlv, d_ohbt, g, jmt, imt);
            iCheck();

    hipLaunchKernelGGL(readyt_cu_5, dim3(numBlocks), dim3(threadsPerBlock ), 0, 0, d_dlv, d_au0, d_aus, d_auw, d_ausw,
            d_wgp, d_work, d_hbx, d_hby, d_viv, d_whx, d_why, d_work1, d_work2, d_psa, d_dxur, d_dyur,
            d_kmu, d_pax, d_pay, pconst_mod_mp_od0_, jmt, imt);
            iCheck();

    if (pconst_mod_mp_ist_ == 0) {
        hipMemset(d_ax, 0, dataSize4d);
        hipMemset(d_ay, 0, dataSize4d);
        hipMemset(d_az, 0, dataSize4d);
    }

    return;
}

