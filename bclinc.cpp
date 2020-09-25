/*by zly 2018.11*/
// Fix  by Wpf, 2019.1.10
#include "hip/hip_runtime.h"
#include "cuda_data.h"
#include<stdio.h>

#include "precision_mod.h"
#include "param_mod.h"
#include "pconst_mod.h"
#include "constant_mod.h"
#include "dyn_mod.h"
#include "work_mod.h"
#include "domain.h"
#include "grid.h"
#include "blocks.h"
#include "operators.h"
#include "forc_mod.h"
#include "common.h"
#include <math.h>
#include <math_functions.h>

extern "C" void mpi2_();

extern "C" void pop_haloupdate_bclinc2_(int *);//d_U/d_V

extern "C" void pop_haloupdate_bclinc3_(int *);//d_wka

__global__ void bclinc12_cu(double d_od0, double d_c0f, double d_cag,
        double *d_snlat, double d_sag,  double *d_up, double *d_vp, double *d_bbcy,
        double *d_sbcy, double *d_sbcx, double *d_bbcx, double *d_su, double *d_sv,
        int *d_kmu, int d_imt, int d_jmt, int d_km) {
    int i, j;

    i = (blockIdx.x) * blockDim.x + threadIdx.x;
    j = (blockIdx.y) * blockDim.y + threadIdx.y;

    if (i >= 1 && i < d_imt - 1 && j >= 1 && j < d_jmt - 1) {
        int kmb = d_kmu[j * d_imt + i] - 1;

        if (kmb >= 0) {
            double temp = d_c0f * sqrt(d_up[kmb * d_jmt * d_imt + j * d_imt + i] *
                               d_up[kmb * d_jmt * d_imt + j * d_imt + i] +
                               d_vp[kmb * d_jmt * d_imt + j * d_imt + i] *
                               d_vp[kmb * d_jmt * d_imt + j * d_imt + i]);

            d_sbcx[j * d_imt + i] = d_su[j * d_imt + i] * d_od0;
            d_sbcy[j * d_imt + i] = d_sv[j * d_imt + i] * d_od0;

            d_bbcx[j * d_imt + i] = temp *
                    (d_up[kmb * d_jmt * d_imt + j * d_imt + i] * d_cag +
                     d_snlat[j * d_imt + i] *
                     d_vp[kmb * d_jmt * d_imt + j * d_imt + i] * d_sag);

            d_bbcy[j * d_imt + i] = temp *
                    (d_vp[kmb * d_jmt * d_imt + j * d_imt + i] * d_cag -
                     d_snlat[j * d_imt + i] *
                     d_up[kmb * d_jmt * d_imt + j * d_imt + i] * d_sag);
        }
    }

    return;
}

__global__ void bclinc11_cu(double d_od0, int d_isc, double d_onbb, double *d_ohbt, double *d_vit,
        double *d_dzp,double *d_h0bf, double *d_h0bl, double *d_gg,  double *d_psa,
        double *d_wka, double *d_wka1, double *d_wka2, double *d_zkt,
        int d_imt, int d_jmt, int d_km) {
    int i, j;

    i = (blockIdx.x) * blockDim.x + threadIdx.x;
    j = (blockIdx.y) * blockDim.y + threadIdx.y;

    if (i < d_imt && j < d_jmt) {
        double aa = 0, wkk[km + 1], work;
        int k;

        if (d_isc != 0) aa = 0.5;

        d_h0bf[j * d_imt + i] *= d_onbb;
        work = aa * d_h0bf[j * d_imt + i] + (1.0 - aa) * d_h0bl[j * d_imt + i];
        wkk[0] = (d_psa[j * d_imt + i] * d_od0 +
                  work * 9.806e0) *
                 d_vit[j * d_imt + i];

        for (k = 0; k < d_km; k++) {
            wkk[k + 1] = wkk[k] -
                         d_gg[k * d_jmt * d_imt + j * d_imt + i] *
                         d_dzp[k] *
                         d_vit[k * d_jmt * d_imt + j * d_imt + i];
        }

        for (k = 0; k < d_km; k++) {
            d_wka1[k * d_jmt * d_imt + j * d_imt + i] = 0.5 * (wkk[k] + wkk[k + 1]);
            d_wka2[k * d_jmt * d_imt + j * d_imt + i] =
                    (1.0 + d_ohbt[j * d_imt + i] * d_zkt[k]) *
                    work *
                    d_vit[k * d_jmt * d_imt + j * d_imt + i];

            d_wka[k * d_jmt * d_imt + j * d_imt + i] = d_wka2[k * d_jmt * d_imt + j * d_imt + i];
        }
    }

    return;
}

__global__ void bclinc2_cu(int d_isc,
        double *d_epea, double *d_epeb, double *d_epla, double *d_eplb, double *d_up,
        double *d_dlu, double *d_gg, double *d_fcor,
        double *d_wka1, double *d_wka2, double *d_dxur, double *d_dyur,
        double *d_viv, double *d_akmu, double d_dtc2, double *d_vp, double *d_bbcy,
        double *d_dlv, double *d_sbcy, int *d_kmu, double *d_wka, double *d_odzt, double *d_odzp,
        double *d_a8, double *d_b8, double *d_c8, double *d_d8, double *d_e8, double *d_f8,
        int d_imt, int d_jmt, int d_km) {
    int i, j, k;

    i = (blockIdx.x) * blockDim.x + threadIdx.x;
    j = (blockIdx.y) * blockDim.y + threadIdx.y;
    k = (blockIdx.z) * blockDim.z + threadIdx.z;

    if (i >= 1 && i < d_imt && j < d_jmt - 1  && k < d_km) {
        int kmb = d_kmu[j * d_imt + i] - 1;
        double gradx, grady, ggu, wk1, wk2;

        gradx = d_fcor[j * d_imt + i] * d_vp[k * d_jmt * d_imt + j * d_imt + i];
        grady = -d_fcor[j * d_imt + i] * d_up[k * d_jmt * d_imt + j * d_imt + i];

        if (kmb >= k) {
            ggu = 0.25 * (d_gg[k * d_jmt * d_imt + j * d_imt + i] +
                          d_gg[k * d_jmt * d_imt + j * d_imt + i - 1] +
                          d_gg[k * d_jmt * d_imt + (j + 1) * d_imt + i] +
                          d_gg[k * d_jmt * d_imt + (j + 1) * d_imt + i - 1]);

            gradx += d_dxur[j * d_imt + i] * 0.5 *
                     (d_wka1[k * d_jmt * d_imt + (j + 1) * d_imt + i] -
                      d_wka1[k * d_jmt * d_imt + j * d_imt + i - 1] -
                      d_wka1[k * d_jmt * d_imt + (j + 1) * d_imt + i - 1] +
                      d_wka1[k * d_jmt * d_imt + j * d_imt + i]);

            grady += d_dyur[j * d_imt + i] * 0.5 *
                     (d_wka1[k * d_jmt * d_imt + (j + 1) * d_imt + i] -
                      d_wka1[k * d_jmt * d_imt + j * d_imt + i - 1] +
                      d_wka1[k * d_jmt * d_imt + (j + 1) * d_imt + i - 1] -
                      d_wka1[k * d_jmt * d_imt + j * d_imt + i]);

            gradx -= ggu * d_dxur[j * d_imt + i] * 0.5 *
                     (d_wka2[k * d_jmt * d_imt + (j + 1) * d_imt + i] -
                      d_wka2[k * d_jmt * d_imt + j * d_imt + i - 1] -
                      d_wka2[k * d_jmt * d_imt + (j + 1) * d_imt + i - 1] +
                      d_wka2[k * d_jmt * d_imt + j * d_imt + i]);

            grady -= ggu * d_dyur[j * d_imt + i] * 0.5 *
                     (d_wka2[k * d_jmt * d_imt + (j + 1) * d_imt + i] -
                      d_wka2[k * d_jmt * d_imt + j * d_imt + i - 1] +
                      d_wka2[k * d_jmt * d_imt + (j + 1) * d_imt + i - 1] -
                      d_wka2[k * d_jmt * d_imt + j * d_imt + i]);
        }

        grady -= d_dlv[k * d_jmt * d_imt + j * d_imt + i];
        gradx -= d_dlu[k * d_jmt * d_imt + j * d_imt + i];

        if (d_isc == 0) {
            wk1 = -(d_epea[j * d_imt + i] * grady +
                    d_epeb[j * d_imt + i] * gradx);
            wk2 = -(d_epea[j * d_imt + i] * gradx -
                    d_epeb[j * d_imt + i] * grady);
        } else {
            wk1 = -(d_epla[j * d_imt + i] * grady +
                    d_eplb[j * d_imt + i] * gradx);
            wk2 = -(d_epla[j * d_imt + i] * gradx -
                    d_eplb[j * d_imt + i] * grady);
        }
        d_dlv[k * d_jmt * d_imt + j * d_imt + i] = wk1 * d_viv[k * d_jmt * d_imt + j * d_imt + i];
        d_dlu[k * d_jmt * d_imt + j * d_imt + i] = wk2 * d_viv[k * d_jmt * d_imt + j * d_imt + i];
    }

    if (i >= 2 && i < d_imt - 2 && j >= 2 && j < d_jmt - 2 && k < d_km) {

        d_wka[k * d_jmt * d_imt + j * d_imt + i] =
                d_vp[k * d_jmt * d_imt + j * d_imt + i] +
                d_dlv[k * d_jmt * d_imt + j * d_imt + i] * d_dtc2;

    }

    return;
}

__global__ void bclinc3_cu(double *d_viv, double *d_akmu, double d_dtc2, double *d_vp, double *d_bbcy,
        double *d_dlv, double *d_sbcy, int *d_kmu, double *d_wka, double *d_odzt, double *d_odzp,
        double *d_a8, double *d_b8, double *d_c8, double *d_d8, double *d_e8, double *d_f8,
        int d_imt, int d_jmt, int d_km) {
    int i, j;

    i = (blockIdx.x) * blockDim.x + threadIdx.x;
    j = (blockIdx.y) * blockDim.y + threadIdx.y;

    if (i >= 2 && i < d_imt - 2 && j >= 2 && j < d_jmt - 2) {
        int kmb = d_kmu[j * d_imt + i] - 1, k;

        if (kmb >= 0) {
            double g0;
                for (k = 1; k <= kmb; k++) {
                    d_a8[k * jmt * imt + j * imt + i] =
                            d_akmu[(k - 1) * jmt * imt + j * imt + i] *
                            d_odzt[k] * d_odzp[k] * d_dtc2 * 0.5;
                    d_d8[k * jmt * imt + j * imt + i] =
                            d_wka[k * jmt * imt + j * imt + i];
                }
                for (k = 1; k < kmb; k++) {
                    d_c8[k * jmt * imt + j * imt + i] =
                            d_akmu[k * jmt * imt + j * imt + i] *
                            d_odzt[k + 1] * d_odzp[k] * d_dtc2 * 0.5;
                    d_b8[k * jmt * imt + j * imt + i] =
                            1.0 + d_a8[k * jmt * imt + j * imt + i] +
                            d_c8[k * jmt * imt + j * imt + i];
                    d_e8[k * jmt * imt + j * imt + i] = 0.0;
                    d_f8[k * jmt * imt + j * imt + i] = 0.0;
                }

                //  B. C. AT TOP
            k = 0;
            d_a8[k * jmt * imt + j * imt + i] = d_odzp[k] * d_dtc2 * 0.5;
            d_c8[k * jmt * imt + j * imt + i] =
                    d_akmu[k * jmt * imt + j * imt + i] *
                    d_odzt[k + 1] * d_odzp[k] * d_dtc2 * 0.5;
            d_b8[k * jmt * imt + j * imt + i] = 1.0 + d_c8[k * jmt * imt + j * imt + i];
            d_d8[k * jmt * imt + j * imt + i] = d_wka[k * jmt * imt + j * imt + i];

            //B. C. AT BOTTOM
            d_b8[kmb * jmt * imt + j * imt + i] = 1.0 + d_a8[kmb * jmt * imt + j * imt + i];
            d_c8[kmb * jmt * imt + j * imt + i] = d_odzp[kmb] * d_dtc2 * 0.5;
            d_e8[(kmb + 1) * jmt * imt + j * imt + i] = 0.0;
            d_f8[(kmb + 1) * jmt * imt + j * imt + i] = 0.0;
            d_d8[kmb * jmt * imt + j * imt + i] =
                    d_wka[kmb * jmt * imt + j * imt + i] -
                    d_bbcy[j * imt + i] * d_odzp[kmb] * d_dtc2 * 0.5;

            //NOW INVERT
            for (k = kmb; k >= 0; k--) {
                g0 = 1.0 / (d_b8[k * jmt * imt + j * imt + i] -
                            d_c8[k * jmt * imt + j * imt + i] *
                            d_e8[(k + 1) * jmt * imt + j * imt + i]);
                d_e8[k * jmt * imt + j * imt + i] = d_a8[k * jmt * imt + j * imt + i] * g0;
                d_f8[k * jmt * imt + j * imt + i] =
                        (d_d8[k * jmt * imt + j * imt + i] +
                         d_c8[k * jmt * imt + j * imt + i] *
                         d_f8[(k + 1) * jmt * imt + j * imt + i]) * g0;
            }

            //B.C. AT SURFACE
            d_wka[j * imt + i] =
                    (d_e8[j * imt + i] * d_sbcy[j * imt + i] +
                     d_f8[j * imt + i]) * d_viv[j * imt + i];
            for (k = 1; k <= kmb; k++) {
                d_wka[k * jmt * imt + j * imt + i] =
                        (d_e8[k * jmt * imt + j * imt + i] *
                         d_wka[(k - 1) * jmt * imt + j * imt + i] +
                         d_f8[k * jmt * imt + j * imt + i]) *
                        d_viv[k * jmt * imt + j * imt + i];
            }
        }
    }

    return;
}

__global__ void bclinc4_cu(double *d_ohbu, double *d_viv, double *d_dzp,
        double *d_odzt, double *d_odzp, double *d_akmu, double d_afc1,
        double d_afc2, double d_dtc2, double *d_v, double *d_up, double *d_vp,
        double *d_vb, double *d_dlu, double *d_sbcx, double *d_bbcx,
        int *d_kmu, double *d_wka,
        double *d_a8, double *d_b8, double *d_c8, double *d_d8, double *d_e8, double *d_f8,
        int d_imt, int d_jmt, int d_km) {

    int i, j, k;

    i = (blockIdx.x) * blockDim.x + threadIdx.x;
    j = (blockIdx.y) * blockDim.y + threadIdx.y;

    if (i < d_imt && j < d_jmt) {
        double work = 0;

        for (k = 0; k < d_km; k++) {
            work += d_dzp[k] *
                    d_ohbu[j * d_imt + i] *
                    d_wka[k * d_jmt * d_imt + j * d_imt + i] *
                    d_viv[k * d_jmt * d_imt + j * d_imt + i];
        }

        for (k = 0; k < d_km; k++) {
            d_wka[k * d_jmt * d_imt + j * d_imt + i] =
                    (d_wka[k * d_jmt * d_imt + j * d_imt + i] -
                     work + d_vb[j * d_imt + i]) *
                    d_viv[k * d_jmt * d_imt + j * d_imt + i];
            d_vp[k * d_jmt * d_imt + j * d_imt + i] =
                    d_afc2 * d_v[k * d_jmt * d_imt + j * d_imt + i] +
                    d_afc1 * (d_vp[k * d_jmt * d_imt + j * d_imt + i] +
                                 d_wka[k * d_jmt * d_imt + j * d_imt + i]);
            d_v[k * d_jmt * d_imt + j * d_imt + i] =
                    d_wka[k * d_jmt * d_imt + j * d_imt + i];
        }

        if (i >= 2 && i < d_imt - 2 && j >= 2 && j < d_jmt - 2) {
            for (k = 0; k < d_km; k++) {
                d_wka[k * d_jmt * d_imt + j * d_imt + i] =
                        d_up[k * d_jmt * d_imt + j * d_imt + i] +
                        d_dlu[k * d_jmt * d_imt + j * d_imt + i] * d_dtc2;
            }

            double g0;
            int kz = d_kmu[j * imt + i];
            if (kz > 0) {
                for (k = 1; k < kz; k++) {
                    d_a8[k * jmt * imt + j * imt + i] =
                            d_akmu[(k - 1) * jmt * imt + j * imt + i] *
                            d_odzt[k] * d_odzp[k] * d_dtc2 * 0.5;
                    d_d8[k * jmt * imt + j * imt + i] =
                            d_wka[k * jmt * imt + j * imt + i];
                }
                for (k = 1; k < kz - 1; k++) {
                    d_c8[k * jmt * imt + j * imt + i] =
                            d_akmu[k * jmt * imt + j * imt + i] *
                            d_odzt[k + 1] * d_odzp[k] * d_dtc2 * 0.5;
                    d_b8[k * jmt * imt + j * imt + i] =
                            1.0 + d_a8[k * jmt * imt + j * imt + i] +
                            d_c8[k * jmt * imt + j * imt + i];
                    d_e8[k * jmt * imt + j * imt + i] = 0.0;
                    d_f8[k * jmt * imt + j * imt + i] = 0.0;
                }

                //  B. C. AT TOP
                k = 0;
                d_a8[k * jmt * imt + j * imt + i] = d_odzp[k] * d_dtc2 * 0.5;
                d_c8[k * jmt * imt + j * imt + i] =
                        d_akmu[k * jmt * imt + j * imt + i] *
                        d_odzt[k + 1] * d_odzp[k] * d_dtc2 * 0.5;
                d_b8[k * jmt * imt + j * imt + i] = 1.0 + d_c8[k * jmt * imt + j * imt + i];
                d_d8[k * jmt * imt + j * imt + i] = d_wka[k * jmt * imt + j * imt + i];
                d_e8[k * jmt * imt + j * imt + i] = 0.0;
                d_f8[k * jmt * imt + j * imt + i] = 0.0;

                //B. C. AT BOTTOM
                d_b8[(kz - 1) * jmt * imt + j * imt + i] = 1.0 + d_a8[(kz - 1) * jmt * imt + j * imt + i];
                d_c8[(kz - 1) * jmt * imt + j * imt + i] = d_odzp[kz - 1] * d_dtc2 * 0.5;
                d_e8[kz * jmt * imt + j * imt + i] = 0.0;
                d_f8[kz * jmt * imt + j * imt + i] = 0.0;
                d_d8[(kz - 1) * jmt * imt + j * imt + i] =
                        d_wka[(kz - 1) * jmt * imt + j * imt + i] -
                        d_bbcx[j * imt + i] * d_odzp[kz - 1] * d_dtc2 * 0.5;

                //NOW INVERT
                for (k = kz - 1; k >= 0; k--) {
                    g0 = 1.0 / (d_b8[k * jmt * imt + j * imt + i] -
                                d_c8[k * jmt * imt + j * imt + i] *
                                d_e8[(k + 1) * jmt * imt + j * imt + i]);
                    d_e8[k * jmt * imt + j * imt + i] =
                            d_a8[k * jmt * imt + j * imt + i] * g0;
                    d_f8[k * jmt * imt + j * imt + i] =
                            (d_d8[k * jmt * imt + j * imt + i] +
                             d_c8[k * jmt * imt + j * imt + i] *
                             d_f8[(k + 1) * jmt * imt + j * imt + i]) * g0;
                }

                //B.C. AT SURFACE
                d_wka[j * imt + i] =
                        (d_e8[j * imt + i] * d_sbcx[j * imt + i] +
                         d_f8[j * imt + i]) * d_viv[j * imt + i];
                for (k = 1; k < kz; k++) {
                    d_wka[k * jmt * imt + j * imt + i] =
                            (d_e8[k * jmt * imt + j * imt + i] *
                             d_wka[(k - 1) * jmt * imt + j * imt + i] +
                             d_f8[k * jmt * imt + j * imt + i]) *
                            d_viv[k * jmt * imt + j * imt + i];
                }
            }
        }
    }

    return;
}

__global__ void bclinc5_cu(double *d_viv, double *d_dzp, double d_afc1, double d_afc2, double *d_u, double *d_v,
        double *d_up, double *d_ub, double *d_wka, double *d_ohbu, double *d_vtf, double *d_utf,
        int d_imt, int d_jmt, int d_km) {

    int i, j, k;

    i = (blockIdx.x) * blockDim.x + threadIdx.x;
    j = (blockIdx.y) * blockDim.y + threadIdx.y;
    k = (blockIdx.z) * blockDim.z + threadIdx.z;

    if (i < d_imt && j < d_jmt && k < d_km) {
        double work = 0.0e0, wka;
        int t;

        for (t = 0; t < d_km; t++) {
            work += d_dzp[t] *
                    d_ohbu[j * d_imt + i] *
                    d_wka[t * d_jmt * d_imt + j * d_imt + i] *
                    d_viv[t * d_jmt * d_imt + j * d_imt + i];
        }

        wka = (d_wka[k * d_jmt * d_imt + j * d_imt + i] -
                work + d_ub[j * d_imt + i]) *
                d_viv[k * d_jmt * d_imt + j * d_imt + i];
        d_up[k * d_jmt * d_imt + j * d_imt + i] =
                d_afc2 * d_u[k * d_jmt * d_imt + j * d_imt + i] +
                d_afc1 * (d_up[k * d_jmt * d_imt + j * d_imt + i] + wka);
        d_u[k * d_jmt * d_imt + j * d_imt + i] = wka;

        d_utf[k * d_jmt * d_imt + j * d_imt + i] += wka;
        d_vtf[k * d_jmt * d_imt + j * d_imt + i] +=
                d_v[k * d_jmt * d_imt + j * d_imt + i];
    }

    return;
}

__global__ void bclinc6_cu(int d_isc, double *d_viv,
        double *d_epea, double *d_epeb, double *d_epla, double *d_eplb, double *d_up, double *d_vp,
        double *d_dlu, double *d_dlv, double *d_gg, int *d_kmu, double *d_fcor,
        double *d_wka1, double *d_wka2, double *d_dxur, double *d_dyur,
        int d_imt, int d_jmt, int d_km) {
    int i, j, k;

    i = (blockIdx.x) * blockDim.x + threadIdx.x;
    j = (blockIdx.y) * blockDim.y + threadIdx.y;
    k = (blockIdx.z) * blockDim.z + threadIdx.z;

    if (i >= 1 && i < d_imt && j < d_jmt - 1 && k < d_km) {
        int kmb = d_kmu[j * d_imt + i] - 1;
        double gradx, grady, ggu, wk1, wk2;

        gradx = d_fcor[j * d_imt + i] * d_vp[k * d_jmt * d_imt + j * d_imt + i];
        grady = -d_fcor[j * d_imt + i] * d_up[k * d_jmt * d_imt + j * d_imt + i];

        if (kmb >= k) {
            ggu = 0.25 * (d_gg[k * d_jmt * d_imt + j * d_imt + i] +
                          d_gg[k * d_jmt * d_imt + j * d_imt + i - 1] +
                          d_gg[k * d_jmt * d_imt + (j + 1) * d_imt + i] +
                          d_gg[k * d_jmt * d_imt + (j + 1) * d_imt + i - 1]);

            gradx += d_dxur[j * d_imt + i] * 0.5 *
                     (d_wka1[k * d_jmt * d_imt + (j + 1) * d_imt + i] -
                      d_wka1[k * d_jmt * d_imt + j * d_imt + i - 1] -
                      d_wka1[k * d_jmt * d_imt + (j + 1) * d_imt + i - 1] +
                      d_wka1[k * d_jmt * d_imt + j * d_imt + i]);

            grady += d_dyur[j * d_imt + i] * 0.5 *
                     (d_wka1[k * d_jmt * d_imt + (j + 1) * d_imt + i] -
                      d_wka1[k * d_jmt * d_imt + j * d_imt + i - 1] +
                      d_wka1[k * d_jmt * d_imt + (j + 1) * d_imt + i - 1] -
                      d_wka1[k * d_jmt * d_imt + j * d_imt + i]);

            gradx -= ggu * d_dxur[j * d_imt + i] * 0.5 *
                     (d_wka2[k * d_jmt * d_imt + (j + 1) * d_imt + i] -
                      d_wka2[k * d_jmt * d_imt + j * d_imt + i - 1] -
                      d_wka2[k * d_jmt * d_imt + (j + 1) * d_imt + i - 1] +
                      d_wka2[k * d_jmt * d_imt + j * d_imt + i]);

            grady -= ggu * d_dyur[j * d_imt + i] * 0.5 *
                     (d_wka2[k * d_jmt * d_imt + (j + 1) * d_imt + i] -
                      d_wka2[k * d_jmt * d_imt + j * d_imt + i - 1] +
                      d_wka2[k * d_jmt * d_imt + (j + 1) * d_imt + i - 1] -
                      d_wka2[k * d_jmt * d_imt + j * d_imt + i]);
        }

        grady -= d_dlv[k * d_jmt * d_imt + j * d_imt + i];
        gradx -= d_dlu[k * d_jmt * d_imt + j * d_imt + i];

        if (d_isc == 0) {
            wk1 = -(d_epea[j * d_imt + i] * grady +
                    d_epeb[j * d_imt + i] * gradx);
            wk2 = -(d_epea[j * d_imt + i] * gradx -
                    d_epeb[j * d_imt + i] * grady);
        } else {
            wk1 = -(d_epla[j * d_imt + i] * grady +
                    d_eplb[j * d_imt + i] * gradx);
            wk2 = -(d_epla[j * d_imt + i] * gradx -
                    d_eplb[j * d_imt + i] * grady);
        }
        d_dlv[k * d_jmt * d_imt + j * d_imt + i] = wk1 * d_viv[k * d_jmt * d_imt + j * d_imt + i];
        d_dlu[k * d_jmt * d_imt + j * d_imt + i] = wk2 * d_viv[k * d_jmt * d_imt + j * d_imt + i];
    }

    return;
}

__global__ void bclinc71_cu(double *d_viv, double *d_odzt, double *d_odzp, double *d_akmu, double d_dtc,
        double *d_u, double *d_up, double *d_dlu, double *d_sbcx, double *d_bbcx, int *d_kmu,
        double *d_a8, double *d_b8, double *d_c8, double *d_d8, double *d_e8, double *d_f8,
        int d_imt, int d_jmt, int d_km) {

    int i, j, k;

    i = (blockIdx.x) * blockDim.x + threadIdx.x;
    j = (blockIdx.y) * blockDim.y + threadIdx.y;

    if (i >= 2 && i < d_imt - 2 && j >= 2 && j < d_jmt - 2) {
        double g0;
        int kz = d_kmu[j * imt + i];

        for (k = 0; k < d_km; k++) {
            d_u[k * d_jmt * d_imt + j * d_imt + i] =
                    d_up[k * d_jmt * d_imt + j * d_imt + i] +
                    d_dlu[k * d_jmt * d_imt + j * d_imt + i] * d_dtc;
        }

        if (kz > 0) {
            for (k = 1; k < kz; k++) {
                d_a8[k * jmt * imt + j * imt + i] =
                        d_akmu[(k - 1) * jmt * imt + j * imt + i] *
                        d_odzt[k] * d_odzp[k] * d_dtc * 0.5;
                d_d8[k * jmt * imt + j * imt + i] =
                        d_u[k * jmt * imt + j * imt + i];
            }
            for (k = 1; k < kz - 1; k++) {
                d_c8[k * jmt * imt + j * imt + i] =
                        d_akmu[k * jmt * imt + j * imt + i] *
                        d_odzt[k + 1] * d_odzp[k] * d_dtc * 0.5;
                d_b8[k * jmt * imt + j * imt + i] =
                        1.0 + d_a8[k * jmt * imt + j * imt + i] +
                        d_c8[k * jmt * imt + j * imt + i];
                d_e8[k * jmt * imt + j * imt + i] = 0.0;
                d_f8[k * jmt * imt + j * imt + i] = 0.0;
            }

            //  B. C. AT TOP
            k = 0;
            d_a8[k * jmt * imt + j * imt + i] = d_odzp[k] * d_dtc * 0.5;
            d_c8[k * jmt * imt + j * imt + i] =
                    d_akmu[k * jmt * imt + j * imt + i] *
                    d_odzt[k + 1] * d_odzp[k] * d_dtc * 0.5;
            d_b8[k * jmt * imt + j * imt + i] = 1.0 + d_c8[k * jmt * imt + j * imt + i];
            d_d8[k * jmt * imt + j * imt + i] = d_u[k * jmt * imt + j * imt + i];
            d_e8[k * jmt * imt + j * imt + i] = 0.0;
            d_f8[k * jmt * imt + j * imt + i] = 0.0;

            //B. C. AT BOTTOM
            d_b8[(kz - 1) * jmt * imt + j * imt + i] = 1.0 + d_a8[(kz - 1) * jmt * imt + j * imt + i];
            d_c8[(kz - 1) * jmt * imt + j * imt + i] = d_odzp[kz - 1] * d_dtc * 0.5;
            d_e8[kz * jmt * imt + j * imt + i] = 0.0;
            d_f8[kz * jmt * imt + j * imt + i] = 0.0;
            d_d8[(kz - 1) * jmt * imt + j * imt + i] =
                    d_u[(kz - 1) * jmt * imt + j * imt + i] -
                    d_bbcx[j * imt + i] * d_odzp[kz - 1] * d_dtc * 0.5;

            //NOW INVERT
            for (k = kz - 1; k >= 0; k--) {
                g0 = 1.0 / (d_b8[k * jmt * imt + j * imt + i] -
                            d_c8[k * jmt * imt + j * imt + i] *
                            d_e8[(k + 1) * jmt * imt + j * imt + i]);
                d_e8[k * jmt * imt + j * imt + i] =
                        d_a8[k * jmt * imt + j * imt + i] * g0;
                d_f8[k * jmt * imt + j * imt + i] =
                        (d_d8[k * jmt * imt + j * imt + i] +
                         d_c8[k * jmt * imt + j * imt + i] *
                         d_f8[(k + 1) * jmt * imt + j * imt + i]) * g0;
            }

            //B.C. AT SURFACE
            d_u[j * imt + i] =
                    (d_e8[j * imt + i] *
                     d_sbcx[j * imt + i] +
                     d_f8[j * imt + i]) *
                    d_viv[j * imt + i];
            for (k = 1; k < kz; k++) {
                d_u[k * jmt * imt + j * imt + i] =
                        (d_e8[k * jmt * imt + j * imt + i] *
                         d_u[(k - 1) * jmt * imt + j * imt + i] +
                         d_f8[k * jmt * imt + j * imt + i]) *
                        d_viv[k * jmt * imt + j * imt + i];
            }
        }
    }

    return;
}

__global__ void bclinc72_cu(double *d_viv, double *d_odzt, double *d_odzp, double *d_akmu, double d_dtc,
        double *d_v, double *d_vp, double *d_bbcy, double *d_dlv, double *d_sbcy,
        int *d_kmu, double *d_a8, double *d_b8, double *d_c8, double *d_d8, double *d_e8, double *d_f8,
        int d_imt, int d_jmt, int d_km) {

    int i, j, k;

    i = (blockIdx.x) * blockDim.x + threadIdx.x;
    j = (blockIdx.y) * blockDim.y + threadIdx.y;

    if (i >= 2 && i < d_imt - 2 && j >= 2 && j < d_jmt - 2) {
        double g0;
        int kz = d_kmu[j * imt + i];

        for (k = 0; k < d_km; k++) {
            d_v[k * d_jmt * d_imt + j * d_imt + i] =
                    d_vp[k * d_jmt * d_imt + j * d_imt + i] +
                    d_dlv[k * d_jmt * d_imt + j * d_imt + i] * d_dtc;
        }

        if (kz > 0) {
            for (k = 1; k < kz; k++) {
                d_a8[k * jmt * imt + j * imt + i] =
                        d_akmu[(k - 1) * jmt * imt + j * imt + i] *
                        d_odzt[k] * d_odzp[k] * d_dtc * 0.5;
                d_d8[k * jmt * imt + j * imt + i] =
                        d_v[k * jmt * imt + j * imt + i];
            }
            for (k = 1; k < kz - 1; k++) {
                d_c8[k * jmt * imt + j * imt + i] =
                        d_akmu[k * jmt * imt + j * imt + i] *
                        d_odzt[k + 1] * d_odzp[k] * d_dtc * 0.5;
                d_b8[k * jmt * imt + j * imt + i] =
                        1.0 + d_a8[k * jmt * imt + j * imt + i] +
                        d_c8[k * jmt * imt + j * imt + i];
                d_e8[k * jmt * imt + j * imt + i] = 0.0;
                d_f8[k * jmt * imt + j * imt + i] = 0.0;
            }

            //  B. C. AT TOP
            k = 0;
            d_a8[k * jmt * imt + j * imt + i] = d_odzp[k] * d_dtc * 0.5;
            d_c8[k * jmt * imt + j * imt + i] =
                    d_akmu[k * jmt * imt + j * imt + i] *
                    d_odzt[k + 1] * d_odzp[k] * d_dtc * 0.5;
            d_b8[k * jmt * imt + j * imt + i] = 1.0 + d_c8[k * jmt * imt + j * imt + i];
            d_d8[k * jmt * imt + j * imt + i] = d_v[k * jmt * imt + j * imt + i];
            d_e8[k * jmt * imt + j * imt + i] = 0.0;
            d_f8[k * jmt * imt + j * imt + i] = 0.0;

            //B. C. AT BOTTOM
            d_b8[(kz - 1) * jmt * imt + j * imt + i] = 1.0 + d_a8[(kz - 1) * jmt * imt + j * imt + i];
            d_c8[(kz - 1) * jmt * imt + j * imt + i] = d_odzp[kz - 1] * d_dtc * 0.5;
            d_e8[kz * jmt * imt + j * imt + i] = 0.0;
            d_f8[kz * jmt * imt + j * imt + i] = 0.0;
            d_d8[(kz - 1) * jmt * imt + j * imt + i] =
                    d_v[(kz - 1) * jmt * imt + j * imt + i] -
                    d_bbcy[j * imt + i] * d_odzp[kz - 1] * d_dtc * 0.5;

            //NOW INVERT
            for (k = kz - 1; k >= 0; k--) {
                g0 = 1.0 / (d_b8[k * jmt * imt + j * imt + i] -
                            d_c8[k * jmt * imt + j * imt + i] *
                            d_e8[(k + 1) * jmt * imt + j * imt + i]);
                d_e8[k * jmt * imt + j * imt + i] =
                        d_a8[k * jmt * imt + j * imt + i] * g0;
                d_f8[k * jmt * imt + j * imt + i] =
                        (d_d8[k * jmt * imt + j * imt + i] +
                         d_c8[k * jmt * imt + j * imt + i] *
                         d_f8[(k + 1) * jmt * imt + j * imt + i]) * g0;
            }

            //B.C. AT SURFACE
            d_v[j * imt + i] =
                    (d_e8[j * imt + i] *
                     d_sbcy[j * imt + i] +
                     d_f8[j * imt + i]) *
                    d_viv[j * imt + i];
            for (k = 1; k < kz; k++) {
                d_v[k * jmt * imt + j * imt + i] =
                        (d_e8[k * jmt * imt + j * imt + i] *
                         d_v[(k - 1) * jmt * imt + j * imt + i] +
                         d_f8[k * jmt * imt + j * imt + i]) *
                        d_viv[k * jmt * imt + j * imt + i];
            }
        }
    }

    return;
}

__global__ void bclinc81_cu(double *d_viv, double *d_dzp, double *d_u1, double *d_u, double *d_ub, double *d_utf, double *d_ohbu,
        int d_imt, int d_jmt, int d_km) {
    int i, j, k;

    i = (blockIdx.x) * blockDim.x + threadIdx.x;
    j = (blockIdx.y) * blockDim.y + threadIdx.y;
    k = (blockIdx.z) * blockDim.z + threadIdx.z;

    if (i < d_imt && j < d_jmt && k < d_km) {
        double work = 0;
	int t;

        for (t = 0; t < d_km; t++) {
            work += d_dzp[t] * d_ohbu[j * d_imt + i] *
                    d_viv[t * d_jmt * d_imt + j * d_imt + i] *
                    d_u1[t * d_jmt * d_imt + j * d_imt + i];
        }

        d_u[k * d_jmt * d_imt + j * d_imt + i] =
                (d_u1[k * d_jmt * d_imt + j * d_imt + i] -
                 work + d_ub[j * d_imt + i]) *
                d_viv[k * d_jmt * d_imt + j * d_imt + i];

        d_utf[k * d_jmt * d_imt + j * d_imt + i] +=
                d_u[k * d_jmt * d_imt + j * d_imt + i];
    }

    return;
}

__global__ void bclinc82_cu(double *d_viv, double *d_dzp, double *d_v1, double *d_v,double *d_vb, double *d_vtf, double *d_ohbu,
        int d_imt, int d_jmt, int d_km) {
    int i, j, k;

    i = (blockIdx.x) * blockDim.x + threadIdx.x;
    j = (blockIdx.y) * blockDim.y + threadIdx.y;
    k = (blockIdx.z) * blockDim.z + threadIdx.z;

    if (i < d_imt && j < d_jmt && k < d_km) {
        double work = 0;
	int t;
        for (t = 0; t < d_km; t++) {
            work += d_dzp[t] * d_ohbu[j * d_imt + i] *
                    d_viv[t * d_jmt * d_imt + j * d_imt + i] *
                    d_v1[t * d_jmt * d_imt + j * d_imt + i];
        }

        d_v[k * d_jmt * d_imt + j * d_imt + i] =
                (d_v1[k * d_jmt * d_imt + j * d_imt + i] -
                 work + d_vb[j * d_imt + i]) *
                d_viv[k * d_jmt * d_imt + j * d_imt + i];
        d_vtf[k * d_jmt * d_imt + j * d_imt + i] +=
                d_v[k * d_jmt * d_imt + j * d_imt + i];
    }

    return;
}

extern "C" void bclinc_() {
    int errorcode;

    dim3 blockSize(16, 8, 1);
    dim3 gridSize((imt + blockSize.x - 1) / blockSize.x, (jmt + blockSize.y - 1) / blockSize.y, 1);

    dim3 blockSize1(16, 8, 4);
    dim3 gridSize1((imt + blockSize1.x - 1) / blockSize1.x,
                   (jmt + blockSize1.y - 1) / blockSize1.y,
                   (km + blockSize1.z - 1) / blockSize1.z);


        iCheck();
    hipLaunchKernelGGL(bclinc12_cu, dim3(gridSize), dim3(blockSize), 0, 0, pconst_mod_mp_od0_, pconst_mod_mp_c0f_, pconst_mod_mp_cag_,
            d_snlat, pconst_mod_mp_sag_, d_up, d_vp, d_bbcy,
            d_sbcy, d_sbcx, d_bbcx, d_su, d_sv,
            d_kmu, imt, jmt, km);//d_sbcx/d_sbcy/d_bbcx/d_bbcy changed
        iCheck();

    hipLaunchKernelGGL(bclinc11_cu, dim3(gridSize), dim3(blockSize), 0, 0, pconst_mod_mp_od0_, pconst_mod_mp_isc_, pconst_mod_mp_onbb_,
            d_ohbt, d_vit, d_dzp, d_h0bf, d_h0bl, d_gg, d_psa, d_wka, d_wka1, d_wka2, d_zkt,
            imt, jmt, km);//d_h0bf/d_wka/d_wka1/d_wka2 changed
        iCheck();
                //CHECK(hipDeviceSynchronize());

    if (pconst_mod_mp_isc_ < 1) {
        hipLaunchKernelGGL(bclinc6_cu, dim3(gridSize1), dim3(blockSize1 ), 0, 0, pconst_mod_mp_isc_, d_viv,
                d_epea, d_epeb, d_epla, d_eplb, d_up, d_vp,
                d_dlu, d_dlv, d_gg, d_kmu, d_fcor,
                d_wka1, d_wka2, d_dxur, d_dyur, imt, jmt, km);//d_dlu/d_dlv changed
        iCheck();
        //CHECK(hipDeviceSynchronize());

        hipLaunchKernelGGL(bclinc71_cu, dim3(gridSize), dim3(blockSize), 0, 0, d_viv, d_odzt, d_odzp,
                d_akmu, pconst_mod_mp_dtc_, d_u, d_up, d_dlu, d_sbcx, d_bbcx, d_kmu,
                d_a8, d_b8, d_c8, d_d8, d_e8, d_f8,
                imt, jmt, km);//d_u changed
        iCheck();
        CHECK(hipMemcpy(dyn_mod_mp_u_, d_u, dataSize3d, hipMemcpyDeviceToHost));
        hipLaunchKernelGGL(bclinc72_cu, dim3(gridSize), dim3(blockSize), 0, 0, d_viv, d_odzt, d_odzp,
                d_akmu, pconst_mod_mp_dtc_, d_v, d_vp, d_bbcy, d_dlv,
                d_sbcy, d_kmu,
                d_a8_1, d_b8_1, d_c8_1, d_d8_1, d_e8_1, d_f8_1,
                imt, jmt, km);//d_v changed
        iCheck();
        CHECK(hipMemcpy(dyn_mod_mp_v_, d_v, dataSize3d, hipMemcpyDeviceToHost));

        //CHECK(hipDeviceSynchronize());

        pop_haloupdate_bclinc2_(&errorcode);//d_U/d_V

        CHECK(hipMemcpy(d_wka1, dyn_mod_mp_u_, dataSize3d, hipMemcpyHostToDevice));
        hipLaunchKernelGGL(bclinc81_cu, dim3(gridSize1), dim3(blockSize1), 0, 0, d_viv, d_dzp, d_wka1, d_u, d_ub,
                d_utf, d_ohbu, imt, jmt, km);//d_u/d_utf changed
        iCheck();

        CHECK(hipMemcpy(d_wka2, dyn_mod_mp_v_, dataSize3d, hipMemcpyHostToDevice));
        hipLaunchKernelGGL(bclinc82_cu, dim3(gridSize1), dim3(blockSize1), 0, 0, d_viv, d_dzp, d_wka2, d_v, d_vb,
                d_vtf, d_ohbu, imt, jmt, km);//d_v/d_vtf changed
        iCheck();
                //CHECK(hipDeviceSynchronize());
    }
    else {
        hipLaunchKernelGGL(bclinc2_cu, dim3(gridSize1), dim3(blockSize1 ), 0, 0, pconst_mod_mp_isc_,
                d_epea, d_epeb, d_epla, d_eplb, d_up,
                d_dlu, d_gg, d_fcor,
                d_wka1, d_wka2, d_dxur, d_dyur,d_viv, d_akmu, pconst_mod_mp_dtc2_, d_vp, d_bbcy,
                d_dlv, d_sbcy, d_kmu, d_wka, d_odzt, d_odzp, d_a8, d_b8, d_c8, d_d8, d_e8, d_f8,
                imt, jmt, km);//d_wka changed
        iCheck();

        //CHECK(hipDeviceSynchronize());
        hipLaunchKernelGGL(bclinc3_cu, dim3(gridSize), dim3(blockSize ), 0, 0, d_viv, d_akmu, pconst_mod_mp_dtc2_, d_vp, d_bbcy,
                d_dlv, d_sbcy, d_kmu, d_wka, d_odzt, d_odzp, d_a8, d_b8, d_c8, d_d8, d_e8, d_f8,
                imt, jmt, km);//d_wka changed
        iCheck();

        hipMemcpy(work_mod_mp_wka_, d_wka, dataSize3d, hipMemcpyDeviceToHost);
        pop_haloupdate_bclinc3_(&errorcode);//d_wka
        hipMemcpy(d_wka, work_mod_mp_wka_, dataSize3d, hipMemcpyHostToDevice);

        hipLaunchKernelGGL(bclinc4_cu, dim3(gridSize), dim3(blockSize ), 0, 0, d_ohbu, d_viv, d_dzp, d_odzt, d_odzp,
                d_akmu, pconst_mod_mp_afc1_, pconst_mod_mp_afc2_, pconst_mod_mp_dtc2_,
                d_v, d_up, d_vp, d_vb, d_dlu, d_sbcx, d_bbcx,
                d_kmu, d_wka, d_a8, d_b8, d_c8, d_d8, d_e8, d_f8,
                imt, jmt, km);//d_wka/d_vp/d_v changed
        iCheck();

        hipMemcpy(work_mod_mp_wka_, d_wka, dataSize3d, hipMemcpyDeviceToHost);
        pop_haloupdate_bclinc3_(&errorcode);//d_wka
        hipMemcpy(d_wka, work_mod_mp_wka_, dataSize3d, hipMemcpyHostToDevice);

        hipLaunchKernelGGL(bclinc5_cu, dim3(gridSize1), dim3(blockSize1 ), 0, 0, d_viv, d_dzp,
                pconst_mod_mp_afc1_, pconst_mod_mp_afc2_, d_u, d_v,
                d_up, d_ub, d_wka, d_ohbu, d_vtf, d_utf,
                imt, jmt, km);//d_up/d_u/d_v changed
        iCheck();
    } 

    pconst_mod_mp_isc_++;

    return;
}
