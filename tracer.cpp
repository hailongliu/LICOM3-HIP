/*by zly 2018.11*/
// Fix some bugs by Wpf, 2019.1.17
// same result now, 2019.1.25
#include "hip/hip_runtime.h"

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cuda_data.h"
#include "def-undef.h"
#include "precision_mod.h"
#include "param_mod.h"
#include "pconst_mod.h"
#include "constant_mod.h"
#include "dyn_mod.h"
#include "work_mod.h"
#include "domain.h"
#include "blocks.h"
#include "distribution.h"
#include "LICOM_Error_mod.h"
#include "forc_mod.h"
#include "gather_scatter.h"
#include "isopyc_mod.h"
#include "pmix_mod.h"
#include "buf_mod.h"
#include "grid.h"
#include "tracer_mod.h"
#include "hmix_del4.h"
#include "hmix_del2.h"
#include "common.h"
//#include "mpif.h"

extern "C" void isopyc_();
extern "C" void isoflux_(int *);
extern "C" void pop_haloupdate_tracer1_(int *);
extern "C" void pop_haloupdate_tracer2_(int *);
extern "C" void deallocate_tracer1_();
extern "C" void deallocate_tracer2_();
extern "C" void allocate_tracer_();
extern "C" void gather_scatter_tracer_(double(*)[imt_global], double(*)[jmt][imt]);
extern "C" void mpi3_(double *);
//extern "C" void tracer_tdel_(int *);
extern int msg_mod_mp_mpi_comm_ocn_;
//extern "C" void mpi_allreduce_(double *, double *, int *, int *, int *, int *, int *);
//double isopyc_mod_mp_ahisop_[max_blocks_clinic][jmt][imt];
//double (*isopyc_mod_mp_k3_)[3][jmt][km + 1][imt];
__device__ double min1(double a, double b, double c, double d, double e, double f) {
    double i = a;
    if (i > b) i = b;
    if (i > c) i = c;
    if (i > d) i = d;
    if (i > e) i = e;
    if (i > f) i = f;

    return i;
}

__device__ double min2(double a, double b) {
    double i = a;
    if (i > b) i = b;

    return i;
}

__device__ double max1(double a, double b, double c, double d, double e, double f) {
    double i = a;
    if (i < b) i = b;
    if (i < c) i = c;
    if (i < d) i = d;
    if (i < e) i = e;
    if (i < f) i = f;

    return i;
}

__device__ double max2(double a, double b) {
    double i = a;
    if (i < b) i = b;

    return i;
}
__global__ void advection_tracer_1_0(double *d_wkb, double *d_wkd, double *d_u_wface, double *d_v_sface,
                                   double *d_hts, double *d_htw, double *d_dxu, double *d_dyu,
                                   int d_nx, int d_km, int d_jmt, int d_imt, int Str) {

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;
    if (Str == 0 || Str == 2) {
        if (k < d_km && j < d_jmt && 1 <= i && i < (d_imt - 1)) {
            d_v_sface[k * d_jmt * d_imt + j * d_imt + i] =
                    (d_wkb[k * d_jmt * d_imt + j * d_imt + i] +
                     d_wkb[k * d_jmt * d_imt + j * d_imt + i + 1]) *
                    d_hts[j * d_nx + i] * 0.25;
        }

        if (k < d_km && 1 <= j && j < d_jmt - 1 && 0 <= i && i < (d_imt)) {
            d_u_wface[k * d_jmt * d_imt + j * d_imt + i] =
                    (d_wkd[k * d_jmt * d_imt + (j - 1) * d_imt + i] +
                     d_wkd[k * d_jmt * d_imt + j * d_imt + i]) *
                    d_htw[j * d_nx + i] * 0.25;
        }
    }

    if (Str == 1) {
        if (k < d_km && j < d_jmt && 1 <= i && i < (d_imt - 1)) {
            d_v_sface[k * d_jmt * d_imt + j * d_imt + i] =
                    (d_wkb[k * d_jmt * d_imt + j * d_imt + i] * d_dxu[j * d_nx + i] +
                     d_wkb[k * d_jmt * d_imt + j * d_imt + i + 1] * d_dxu[j * d_nx + (i + 1)]) * 0.25;
        }

        if (k < d_km && 1 <= j && j < d_jmt - 1 && 0 <= i && i < (d_imt)) {
            d_u_wface[k * d_jmt * d_imt + j * d_imt + i] =
                    (d_wkd[k * d_jmt * d_imt + (j - 1) * d_imt + i] * d_dyu[(j - 1) * d_nx + i] +
                     d_wkd[k * d_jmt * d_imt + j * d_imt + i] * d_dyu[j * d_nx + i]) * 0.25;
        }
    }
}
__global__ void advection_tracer_1(double *d_u_wface, double *d_v_sface,
                                   double *d_ws, double *d_at, double *d_adv_tt, double *d_tarea_r,
                                   double *d_odzp, double *d_ax,
                                   double *d_ay, double *d_az, int d_nx, int d_km, int d_jmt,
                                   int d_imt, int d_mtracer, int d_nss, int Str) {

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;
    double adv_z1, adv_z2;

    if (Str == 0) {
        if (2 <= j && j < (d_jmt - 2) && 2 <= i && i < (d_imt - 2)) {
            d_adv_tt[j * d_imt + i] =
                    (-d_u_wface[j * d_imt + i] * (d_at[j * d_imt + i] - d_at[j * d_imt + i - 1]) -
                     d_u_wface[j * d_imt + i + 1] * (d_at[j * d_imt + i + 1] - d_at[j * d_imt + i]) -
                     d_v_sface[j * d_imt + i] * (d_at[(j + 1) * d_imt + i] - d_at[j * d_imt + i]) -
                     d_v_sface[(j - 1) * d_imt + i] * (d_at[j * d_imt + i] - d_at[(j - 1) * d_imt + i])) *
                    d_tarea_r[j * d_nx + i];

            d_ax[d_mtracer * d_km * d_jmt * d_imt + j * d_imt + i] =
                    d_ax[d_mtracer * d_km * d_jmt * d_imt + j * d_imt + i] +
                    (-d_u_wface[j * d_imt + i] * (d_at[j * d_imt + i] - d_at[j * d_imt + i - 1]) -
                     d_u_wface[j * d_imt + i + 1] * (d_at[j * d_imt + i + 1] - d_at[j * d_imt + i])) *
                    d_tarea_r[j * d_nx + i] / (double) d_nss;

            d_ay[d_mtracer * d_km * d_jmt * d_imt + j * d_imt + i] =
                    d_ay[d_mtracer * d_km * d_jmt * d_imt + j * d_imt + i] +
                    (-d_v_sface[j * d_imt + i] * (d_at[(j + 1) * d_imt + i] - d_at[j * d_imt + i]) -
                     d_v_sface[(j - 1) * d_imt + i] * (d_at[j * d_imt + i] - d_at[(j - 1) * d_imt + i])) *
                    d_tarea_r[j * d_nx + i] / (double) d_nss;

            adv_z2 = d_ws[1 * d_jmt * d_imt + j * d_imt + i] *
                     (d_at[j * d_imt + i] - d_at[1 * d_jmt * d_imt + j * d_imt + i]);

            d_adv_tt[j * d_imt + i] = d_adv_tt[j * d_imt + i] - p5 * d_odzp[0] * adv_z2;

            d_az[d_mtracer * d_km * d_jmt * d_imt + j * d_imt + i] =
                    d_az[d_mtracer * d_km * d_jmt * d_imt + j * d_imt + i] -
                    p5 * d_odzp[0] * adv_z2 / (double) d_nss;
        }

        if (1 <= k && k < d_km - 1 && 2 <= j && j < (d_jmt - 2) && 2 <= i && i < (d_imt - 2)) {
            d_adv_tt[k * d_jmt * d_imt + j * d_imt + i] =
                    (-d_u_wface[k * d_jmt * d_imt + j * d_imt + i] *
                     (d_at[k * d_jmt * d_imt + j * d_imt + i] - d_at[k * d_jmt * d_imt + j * d_imt + i - 1]) -
                     d_u_wface[k * d_jmt * d_imt + j * d_imt + i + 1] *
                     (d_at[k * d_jmt * d_imt + j * d_imt + i + 1] - d_at[k * d_jmt * d_imt + j * d_imt + i]) -
                     d_v_sface[k * d_jmt * d_imt + j * d_imt + i] *
                     (d_at[k * d_jmt * d_imt + (j + 1) * d_imt + i] - d_at[k * d_jmt * d_imt + j * d_imt + i]) -
                     d_v_sface[k * d_jmt * d_imt + (j - 1) * d_imt + i] *
                     (d_at[k * d_jmt * d_imt + j * d_imt + i] - d_at[k * d_jmt * d_imt + (j - 1) * d_imt + i])) *
                    d_tarea_r[j * d_nx + i];

            d_ax[d_mtracer * d_km * d_jmt * d_imt + k * d_jmt * d_imt + j * d_imt + i] =
                    d_ax[d_mtracer * d_km * d_jmt * d_imt + k * d_jmt * d_imt + j * d_imt + i] +
                    (-d_u_wface[k * d_jmt * d_imt + j * d_imt + i] *
                     (d_at[k * d_jmt * d_imt + j * d_imt + i] - d_at[k * d_jmt * d_imt + j * d_imt + i - 1]) -
                     d_u_wface[k * d_jmt * d_imt + j * d_imt + i + 1] *
                     (d_at[k * d_jmt * d_imt + j * d_imt + i + 1] - d_at[k * d_jmt * d_imt + j * d_imt + i])) *
                    d_tarea_r[j * d_nx + i] / (double) d_nss;

            d_ay[d_mtracer * d_km * d_jmt * d_imt + k * d_jmt * d_imt + j * d_imt + i] =
                    d_ay[d_mtracer * d_km * d_jmt * d_imt + k * d_jmt * d_imt + j * d_imt + i] +
                    (-d_v_sface[k * d_jmt * d_imt + j * d_imt + i] *
                     (d_at[k * d_jmt * d_imt + (j + 1) * d_imt + i] - d_at[k * d_jmt * d_imt + j * d_imt + i]) -
                     d_v_sface[k * d_jmt * d_imt + (j - 1) * d_imt + i] *
                     (d_at[k * d_jmt * d_imt + j * d_imt + i] - d_at[k * d_jmt * d_imt + (j - 1) * d_imt + i])) *
                    d_tarea_r[j * d_nx + i] / (double) d_nss;


            adv_z1 = d_ws[k * d_jmt * d_imt + j * d_imt + i] *
                     (d_at[(k - 1) * d_jmt * d_imt + j * d_imt + i] -
                      d_at[k * d_jmt * d_imt + j * d_imt + i]);

            adv_z2 = d_ws[(k + 1) * d_jmt * d_imt + j * d_imt + i] *
                     (d_at[k * d_jmt * d_imt + j * d_imt + i] -
                      d_at[(k + 1) * d_jmt * d_imt + j * d_imt + i]);

            d_adv_tt[k * d_jmt * d_imt + j * d_imt + i] =
                    d_adv_tt[k * d_jmt * d_imt + j * d_imt + i] -
                    p5 * d_odzp[k] * (adv_z1 + adv_z2);

            d_az[d_mtracer * d_km * d_jmt * d_imt + k * d_jmt * d_imt + j * d_imt + i] =
                    d_az[d_mtracer * d_km * d_jmt * d_imt + k * d_jmt * d_imt + j * d_imt + i] -
                    p5 * d_odzp[k] * (adv_z1 + adv_z2) / (double) d_nss;
        }

        if (2 <= j && j < (d_jmt - 2) && 2 <= i && i < (d_imt - 2)) {
            d_adv_tt[(d_km - 1) * d_jmt * d_imt + j * d_imt + i] =
                    (-d_u_wface[(d_km - 1) * d_jmt * d_imt + j * d_imt + i] *
                     (d_at[(d_km - 1) * d_jmt * d_imt + j * d_imt + i] -
                      d_at[(d_km - 1) * d_jmt * d_imt + j * d_imt + i - 1]) -
                     d_u_wface[(d_km - 1) * d_jmt * d_imt + j * d_imt + i + 1] *
                     (d_at[(d_km - 1) * d_jmt * d_imt + j * d_imt + i + 1] -
                      d_at[(d_km - 1) * d_jmt * d_imt + j * d_imt + i]) -
                     d_v_sface[(d_km - 1) * d_jmt * d_imt + j * d_imt + i] *
                     (d_at[(d_km - 1) * d_jmt * d_imt + (j + 1) * d_imt + i] -
                      d_at[(d_km - 1) * d_jmt * d_imt + j * d_imt + i]) -
                     d_v_sface[(d_km - 1) * d_jmt * d_imt + (j - 1) * d_imt + i] *
                     (d_at[(d_km - 1) * d_jmt * d_imt + j * d_imt + i] -
                      d_at[(d_km - 1) * d_jmt * d_imt + (j - 1) * d_imt + i])) *
                    d_tarea_r[j * d_nx + i];

            d_ax[d_mtracer * d_km * d_jmt * d_imt + (d_km - 1) * d_jmt * d_imt + j * d_imt + i] =
                    d_ax[d_mtracer * d_km * d_jmt * d_imt + (d_km - 1) * d_jmt * d_imt + j * d_imt + i] +
                    (-d_u_wface[(d_km - 1) * d_jmt * d_imt + j * d_imt + i] *
                     (d_at[(d_km - 1) * d_jmt * d_imt + j * d_imt + i] -
                      d_at[(d_km - 1) * d_jmt * d_imt + j * d_imt + i - 1]) -
                     d_u_wface[(d_km - 1) * d_jmt * d_imt + j * d_imt + i + 1] *
                     (d_at[(d_km - 1) * d_jmt * d_imt + j * d_imt + i + 1] -
                      d_at[(d_km - 1) * d_jmt * d_imt + j * d_imt + i])) *
                    d_tarea_r[j * d_nx + i] / (double) d_nss;

            d_ay[d_mtracer * d_km * d_jmt * d_imt + (d_km - 1) * d_jmt * d_imt + j * d_imt + i] =
                    d_ay[d_mtracer * d_km * d_jmt * d_imt + (d_km - 1) * d_jmt * d_imt + j * d_imt + i] +
                    (-d_v_sface[(d_km - 1) * d_jmt * d_imt + j * d_imt + i] *
                     (d_at[(d_km - 1) * d_jmt * d_imt + (j + 1) * d_imt + i] -
                      d_at[(d_km - 1) * d_jmt * d_imt + j * d_imt + i]) -
                     d_v_sface[(d_km - 1) * d_jmt * d_imt + (j - 1) * d_imt + i] *
                     (d_at[(d_km - 1) * d_jmt * d_imt + j * d_imt + i] -
                      d_at[(d_km - 1) * d_jmt * d_imt + (j - 1) * d_imt + i])) *
                    d_tarea_r[j * d_nx + i] / (double) d_nss;


            adv_z1 = d_ws[(d_km - 1) * d_jmt * d_imt + j * d_imt + i] *
                     (d_at[(d_km - 2) * d_jmt * d_imt + j * d_imt + i] -
                      d_at[(d_km - 1) * d_jmt * d_imt + j * d_imt + i]);

            d_adv_tt[(d_km - 1) * d_jmt * d_imt + j * d_imt + i] =
                    d_adv_tt[(d_km - 1) * d_jmt * d_imt + j * d_imt + i] -
                    p5 * d_odzp[d_km - 1] * adv_z1;

            d_az[d_mtracer * d_km * d_jmt * d_imt + (d_km - 1) * d_jmt * d_imt + j * d_imt + i] =
                    d_az[d_mtracer * d_km * d_jmt * d_imt + (d_km - 1) * d_jmt * d_imt + j * d_imt + i] -
                    p5 * d_odzp[d_km - 1] * adv_z1 / (double) d_nss;
        }
    } else if (Str == 1) {
        if (k < d_km && 2 <= j && j < (d_jmt - 2) && 2 <= i && i < (d_imt - 2)) {
            d_adv_tt[k * d_jmt * d_imt + j * d_imt + i] =
                    (-d_u_wface[k * d_jmt * d_imt + j * d_imt + i] *
                     (d_at[k * d_jmt * d_imt + j * d_imt + i] + d_at[k * d_jmt * d_imt + j * d_imt + i - 1]) +
                     d_u_wface[k * d_jmt * d_imt + j * d_imt + i + 1] *
                     (d_at[k * d_jmt * d_imt + j * d_imt + i + 1] + d_at[k * d_jmt * d_imt + j * d_imt + i]) -
                     d_v_sface[k * d_jmt * d_imt + (j - 1) * d_imt + i] *
                     (d_at[k * d_jmt * d_imt + j * d_imt + i] + d_at[k * d_jmt * d_imt + (j - 1) * d_imt + i]) +
                     d_v_sface[k * d_jmt * d_imt + j * d_imt + i] *
                     (d_at[k * d_jmt * d_imt + (j + 1) * d_imt + i] + d_at[k * d_jmt * d_imt + j * d_imt + i])) *
                    d_tarea_r[j * d_nx + i];

            if (k == 0)
                adv_z1 = 0.0e0;
            else
                adv_z1 = d_ws[k * d_jmt * d_imt + j * d_imt + i] *
                         (d_at[(k - 1) * d_jmt * d_imt + j * d_imt + i] +
                          d_at[k * d_jmt * d_imt + j * d_imt + i]) * p5;

            if (k == d_km - 1)
                adv_z2 = 0.0e0;
            else
                adv_z2 = d_ws[(k + 1) * d_jmt * d_imt + j * d_imt + i] *
                         (d_at[k * d_jmt * d_imt + j * d_imt + i] +
                          d_at[(k + 1) * d_jmt * d_imt + j * d_imt + i]) * p5;

            d_adv_tt[k * d_jmt * d_imt + j * d_imt + i] =
                    d_adv_tt[k * d_jmt * d_imt + j * d_imt + i] - d_odzp[k] * (adv_z2 - adv_z1);
        }
    }

}

__global__ void advection_tracer_2(double *d_wkb, double *d_wkd, double *d_u_wface, double *d_v_sface,
                                   double *d_hts, double *d_htw, double *d_dxu, double *d_dyu, double *d_ws,
                                   double *d_at, double *d_adv_tt,
                                   double *d_tarea_r, double *d_odzp, double *d_ax, double *d_ay, double *d_az,
                                   double *d_hun, double *d_hue,
                                   double *d_at00, double *d_atmax, double *d_atmin, double *d_odzt, double *d_vit,
                                   double d_dts,
                                   int d_ny, int d_nx, int d_km, int d_jmt, int d_imt, int d_ntra, int d_mtracer,
                                   int d_nss, int Str) {

    double adv_xy1, adv_xy2, adv_xy3, adv_xy4;
    double wt1, wt2;
    double adv_zz, adv_za, adv_zb1, adv_zb2;
    double adv_x0, adv_y0, adv_xx, adv_yy;
    double adv_c1, adv_c2, adv_zc;

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    if (k < d_km) {
        if (1 <= i && i < (d_imt - 1) && 1 <= j && j < (d_jmt - 1)) {
            adv_x0 = (d_at[k * d_jmt * d_imt + j * d_imt + i + 1] + d_at[k * d_jmt * d_imt + j * d_imt + i]) *
                     d_u_wface[k * d_jmt * d_imt + j * d_imt + i + 1] * d_tarea_r[j * d_nx + i] -
                     (d_at[k * d_jmt * d_imt + j * d_imt + i] + d_at[k * d_jmt * d_imt + j * d_imt + i - 1]) *
                     d_u_wface[k * d_jmt * d_imt + j * d_imt + i] * d_tarea_r[j * d_nx + i];

            adv_y0 = ((d_at[k * d_jmt * d_imt + (j + 1) * d_imt + i] + d_at[k * d_jmt * d_imt + j * d_imt + i]) *
                      d_v_sface[k * d_jmt * d_imt + j * d_imt + i] -
                      (d_at[k * d_jmt * d_imt + j * d_imt + i] + d_at[k * d_jmt * d_imt + (j - 1) * d_imt + i]) *
                      d_v_sface[k * d_jmt * d_imt + (j - 1) * d_imt + i]) * d_tarea_r[j * d_nx + i];

            adv_xy1 = -d_dts *
                      (d_at[k * d_jmt * d_imt + j * d_imt + i + 1] - d_at[k * d_jmt * d_imt + j * d_imt + i]) *
                      2.0 * d_tarea_r[j * d_nx + i] * d_u_wface[k * d_jmt * d_imt + j * d_imt + i + 1] *
                      d_u_wface[k * d_jmt * d_imt + j * d_imt + i + 1] /
                      (d_htw[j * d_nx + i + 1] * d_hun[j * d_nx + i + 1]);

            adv_xy2 = d_dts *
                      (d_at[k * d_jmt * d_imt + j * d_imt + i] - d_at[k * d_jmt * d_imt + j * d_imt + i - 1]) *
                      2.0 * d_tarea_r[j * d_nx + i] * d_u_wface[k * d_jmt * d_imt + j * d_imt + i] *
                      d_u_wface[k * d_jmt * d_imt + j * d_imt + i] / (d_htw[j * d_nx + i] * d_hun[j * d_nx + i]);

            adv_xy3 = -d_dts *
                      (d_at[k * d_jmt * d_imt + (j + 1) * d_imt + i] - d_at[k * d_jmt * d_imt + j * d_imt + i]) *
                      2.0 * d_tarea_r[j * d_nx + i] *
                      d_v_sface[k * d_jmt * d_imt + j * d_imt + i] * d_v_sface[k * d_jmt * d_imt + j * d_imt + i] /
                      (d_hts[j * d_nx + i] * d_hue[j * d_nx + i]);

            adv_xy4 = d_dts *
                      (d_at[k * d_jmt * d_imt + j * d_imt + i] - d_at[k * d_jmt * d_imt + (j - 1) * d_imt + i]) *
                      2.0 * d_tarea_r[j * d_nx + i] * d_v_sface[k * d_jmt * d_imt + (j - 1) * d_imt + i] *
                      d_v_sface[k * d_jmt * d_imt + (j - 1) * d_imt + i] /
                      (d_hts[(j - 1) * d_nx + i] * d_hue[(j - 1) * d_nx + i]);

            adv_c1 = -d_at[k * d_jmt * d_imt + j * d_imt + i] *
                     (d_u_wface[k * d_jmt * d_imt + j * d_imt + i + 1] -
                      d_u_wface[k * d_jmt * d_imt + j * d_imt + i]) * d_tarea_r[j * d_nx + i] * 2.0;

            adv_c2 = -d_at[k * d_jmt * d_imt + j * d_imt + i] *
                     (d_v_sface[k * d_jmt * d_imt + j * d_imt + i] -
                      d_v_sface[k * d_jmt * d_imt + (j - 1) * d_imt + i]) * d_tarea_r[j * d_nx + i] * 2.0;

            if (k == 0) {
                adv_za = -0.5 * d_odzp[0] * d_ws[d_jmt * d_imt + j * d_imt + i] *
                         (d_at[d_jmt * d_imt + j * d_imt + i] + d_at[j * d_imt + i]);

                adv_zb1 = 0.0;
                adv_zb2 = 0.5 * d_odzp[0] *
                          d_ws[d_jmt * d_imt + j * d_imt + i] *
                          d_ws[d_jmt * d_imt + j * d_imt + i] * d_odzt[1] *
                          (d_at[j * d_imt + i] - d_at[d_jmt * d_imt + j * d_imt + i]) * d_dts;

                adv_zc = d_odzp[0] * d_at[j * d_imt + i] * d_ws[d_jmt * d_imt + j * d_imt + i];

            } else if (k == (d_km - 1)) {
                adv_za = 0.5 * d_odzp[d_km - 1] * d_ws[(d_km - 1) * d_jmt * d_imt + j * d_imt + i] *
                         (d_at[(d_km - 1) * d_jmt * d_imt + j * d_imt + i] +
                          d_at[(d_km - 2) * d_jmt * d_imt + j * d_imt + i]);

                adv_zb1 = -0.5 * d_odzp[d_km - 1] *
                          d_ws[(d_km - 1) * d_jmt * d_imt + j * d_imt + i] *
                          d_ws[(d_km - 1) * d_jmt * d_imt + j * d_imt + i] * d_odzt[d_km - 1] *
                          (d_at[(d_km - 2) * d_jmt * d_imt + j * d_imt + i] -
                           d_at[(d_km - 1) * d_jmt * d_imt + j * d_imt + i]) * d_dts;

                adv_zb2 = 0.0;
                adv_zc = -d_odzp[d_km - 1] * d_at[(d_km - 1) * d_jmt * d_imt + j * d_imt + i] *
                         d_ws[(d_km - 1) * d_jmt * d_imt + j * d_imt + i];
            } else {
                adv_za = 0.5 * d_odzp[k] * d_ws[k * d_jmt * d_imt + j * d_imt + i] *
                         (d_at[k * d_jmt * d_imt + j * d_imt + i] + d_at[(k - 1) * d_jmt * d_imt + j * d_imt + i]) -
                         0.5 * d_odzp[k] * d_ws[(k + 1) * d_jmt * d_imt + j * d_imt + i] *
                         (d_at[k * d_jmt * d_imt + j * d_imt + i] +
                          d_at[(k + 1) * d_jmt * d_imt + j * d_imt + i]);

                adv_zb1 = -0.5 * d_odzp[k] *
                          d_ws[k * d_jmt * d_imt + j * d_imt + i] * d_ws[k * d_jmt * d_imt + j * d_imt + i] *
                          d_odzt[k] *
                          (d_at[(k - 1) * d_jmt * d_imt + j * d_imt + i] -
                           d_at[k * d_jmt * d_imt + j * d_imt + i]) * d_dts;

                adv_zb2 = 0.5 * d_odzp[k] *
                          d_ws[(k + 1) * d_jmt * d_imt + j * d_imt + i] *
                          d_ws[(k + 1) * d_jmt * d_imt + j * d_imt + i] * d_odzt[k + 1] *
                          (d_at[k * d_jmt * d_imt + j * d_imt + i] -
                           d_at[(k + 1) * d_jmt * d_imt + j * d_imt + i]) * d_dts;

                adv_zc = -d_odzp[k] * d_at[k * d_jmt * d_imt + j * d_imt + i] *
                         (d_ws[k * d_jmt * d_imt + j * d_imt + i] -
                          d_ws[(k + 1) * d_jmt * d_imt + j * d_imt + i]);
            }

            adv_xx = -(adv_x0 + adv_xy1 + adv_xy2 + adv_c1);
            adv_yy = -(adv_y0 + adv_xy3 + adv_xy4 + adv_c2);
            adv_zz = -(adv_zb1 + adv_zb2 + adv_za + adv_zc);
            d_at00[k * d_jmt * d_imt + j * d_imt + i] =
                    d_at[k * d_jmt * d_imt + j * d_imt + i] +
                    (adv_xx + adv_yy + adv_zz) * d_dts;
        }
    }
    wt1 = -1.0e10;
    wt2 = +1.0e10;

    if (k < d_km) {
        if (1 <= j && j < (d_jmt - 1) && 1 <= i && i < (d_imt - 1)) {
            if (k == 0) {
                d_atmax[k * d_jmt * d_imt + j * d_imt + i] = max1(
                        d_at[k * d_jmt * d_imt + j * d_imt + i] *
                        d_vit[k * d_jmt * d_imt + j * d_imt + i] +
                        (1.0e0 - d_vit[k * d_jmt * d_imt + j * d_imt + i]) * wt1,
                        d_at[k * d_jmt * d_imt + (j - 1) * d_imt + i] *
                        d_vit[k * d_jmt * d_imt + (j - 1) * d_imt + i] +
                        (1.0e0 - d_vit[k * d_jmt * d_imt + (j - 1) * d_imt + i]) * wt1,
                        d_at[k * d_jmt * d_imt + (j + 1) * d_imt + i] *
                        d_vit[k * d_jmt * d_imt + (j + 1) * d_imt + i] +
                        (1.0e0 - d_vit[k * d_jmt * d_imt + (j + 1) * d_imt + i]) * wt1,
                        d_at[k * d_jmt * d_imt + j * d_imt + (i - 1)] *
                        d_vit[k * d_jmt * d_imt + j * d_imt + (i - 1)] +
                        (1.0e0 - d_vit[k * d_jmt * d_imt + j * d_imt + (i - 1)]) * wt1,
                        d_at[k * d_jmt * d_imt + j * d_imt + (i + 1)] *
                        d_vit[k * d_jmt * d_imt + j * d_imt + (i + 1)] +
                        (1.0e0 - d_vit[k * d_jmt * d_imt + j * d_imt + (i + 1)]) * wt1,
                        d_at[(k + 1) * d_jmt * d_imt + j * d_imt + i] *
                        d_vit[(k + 1) * d_jmt * d_imt + j * d_imt + i] +
                        (1.0e0 - d_vit[(k + 1) * d_jmt * d_imt + j * d_imt + i]) * wt1);

                d_atmin[k * d_jmt * d_imt + j * d_imt + i] = min1(
                        d_at[k * d_jmt * d_imt + j * d_imt + i]
                        * d_vit[k * d_jmt * d_imt + j * d_imt + i] +
                        (1.0e0 - d_vit[k * d_jmt * d_imt + j * d_imt + i]) * wt2,
                        d_at[k * d_jmt * d_imt + (j - 1) * d_imt + i] * d_vit[(j - 1) * d_imt + i] +
                        (1.0e0 - d_vit[k * d_jmt * d_imt + (j - 1) * d_imt + i]) * wt2,
                        d_at[k * d_jmt * d_imt + (j + 1) * d_imt + i] *
                        d_vit[k * d_jmt * d_imt + (j + 1) * d_imt + i] +
                        (1.0e0 - d_vit[k * jmt * imt + (j + 1) * imt + i]) * wt2,
                        d_at[k * d_jmt * d_imt + j * d_imt + (i - 1)] *
                        d_vit[k * d_jmt * d_imt + j * d_imt + (i - 1)] +
                        (1.0e0 - d_vit[k * d_jmt * d_imt + j * d_imt + (i - 1)]) * wt2,
                        d_at[k * d_jmt * d_imt + j * d_imt + (i + 1)] *
                        d_vit[k * d_jmt * d_imt + j * d_imt + (i + 1)] +
                        (1.0e0 - d_vit[k * d_jmt * d_imt + j * d_imt + (i + 1)]) * wt2,
                        d_at[(k + 1) * d_jmt * d_imt + j * d_imt + i] *
                        d_vit[(k + 1) * d_jmt * d_imt + j * d_imt + i] +
                        (1.0e0 - d_vit[(k + 1) * d_jmt * d_imt + j * d_imt + i]) * wt2);
            } else if (k == (d_km - 1)) {
                d_atmax[k * d_jmt * d_imt + j * d_imt + i] = max1(
                        d_at[k * d_jmt * d_imt + j * d_imt + i] * d_vit[k * d_jmt * d_imt + j * d_imt + i] +
                        (1.0e0 - d_vit[k * d_jmt * d_imt + j * d_imt + i]) * wt1,
                        d_at[k * d_jmt * d_imt + (j - 1) * d_imt + i] *
                        d_vit[k * d_jmt * d_imt + (j - 1) * d_imt + i] +
                        (1.0e0 - d_vit[k * d_jmt * d_imt + (j - 1) * d_imt + i]) * wt1,
                        d_at[k * d_jmt * d_imt + (j + 1) * d_imt + i] *
                        d_vit[k * d_jmt * d_imt + (j + 1) * d_imt + i] +
                        (1.0e0 - d_vit[k * d_jmt * d_imt + (j + 1) * d_imt + i]) * wt1,
                        d_at[k * d_jmt * d_imt + j * d_imt + (i - 1)] *
                        d_vit[k * d_jmt * d_imt + j * d_imt + (i - 1)] +
                        (1.0e0 - d_vit[k * d_jmt * d_imt + j * d_imt + (i - 1)]) * wt1,
                        d_at[k * d_jmt * d_imt + j * d_imt + (i + 1)] *
                        d_vit[k * d_jmt * d_imt + j * d_imt + (i + 1)] +
                        (1.0e0 - d_vit[k * d_jmt * d_imt + j * d_imt + (i + 1)]) * wt1,
                        d_at[(k - 1) * d_jmt * d_imt + j * d_imt + i] *
                        d_vit[(k - 1) * d_jmt * d_imt + j * d_imt + i] +
                        (1.0e0 - d_vit[(k - 1) * d_jmt * d_imt + j * d_imt + i]) * wt1);

                d_atmin[k * d_jmt * d_imt + j * d_imt + i] = min1(
                        d_at[k * d_jmt * d_imt + j * d_imt + i] * d_vit[k * d_jmt * d_imt + j * d_imt + i] +
                        (1.0e0 - d_vit[k * d_jmt * d_imt + j * d_imt + i]) * wt2,
                        d_at[k * d_jmt * d_imt + (j - 1) * d_imt + i] *
                        d_vit[k * d_jmt * d_imt + (j - 1) * d_imt + i] +
                        (1.0e0 - d_vit[k * d_jmt * d_imt + (j - 1) * d_imt + i]) * wt2,
                        d_at[k * d_jmt * d_imt + (j + 1) * d_imt + i] *
                        d_vit[k * d_jmt * d_imt + (j + 1) * d_imt + i] +
                        (1.0e0 - d_vit[k * d_jmt * d_imt + (j + 1) * d_imt + i]) * wt2,
                        d_at[k * d_jmt * d_imt + j * d_imt + (i - 1)] *
                        d_vit[k * d_jmt * d_imt + j * d_imt + (i - 1)] +
                        (1.0e0 - d_vit[k * d_jmt * d_imt + j * d_imt + (i - 1)]) * wt2,
                        d_at[k * d_jmt * d_imt + j * d_imt + (i + 1)] *
                        d_vit[k * d_jmt * d_imt + j * d_imt + (i + 1)] +
                        (1.0e0 - d_vit[k * d_jmt * d_imt + j * d_imt + (i + 1)]) * wt2,
                        d_at[(k - 1) * d_jmt * d_imt + j * d_imt + i] *
                        d_vit[(k - 1) * d_jmt * d_imt + j * d_imt + i] +
                        (1.0e0 - d_vit[(k - 1) * d_jmt * d_imt + j * d_imt + i]) * wt2);
            } else {
                d_atmax[k * d_jmt * d_imt + j * d_imt + i] = max2(max1(
                        d_at[k * d_jmt * d_imt + j * d_imt + i] *
                        d_vit[k * d_jmt * d_imt + j * d_imt + i] +
                        (1.0e0 - d_vit[k * d_jmt * d_imt + j * d_imt + i]) * wt1,
                        d_at[k * d_jmt * d_imt + (j - 1) * d_imt + i] *
                        d_vit[k * d_jmt * d_imt + (j - 1) * d_imt + i] +
                        (1.0e0 - d_vit[k * d_jmt * d_imt + (j - 1) * d_imt + i]) * wt1,
                        d_at[k * d_jmt * d_imt + (j + 1) * d_imt + i] *
                        d_vit[k * d_jmt * d_imt + (j + 1) * d_imt + i] +
                        (1.0e0 - d_vit[k * d_jmt * d_imt + (j + 1) * d_imt + i]) * wt1,
                        d_at[k * d_jmt * d_imt + j * d_imt + (i - 1)] *
                        d_vit[k * d_jmt * d_imt + j * d_imt + (i - 1)] +
                        (1.0e0 - d_vit[k * d_jmt * d_imt + j * d_imt + (i - 1)]) * wt1,
                        d_at[k * d_jmt * d_imt + j * d_imt + (i + 1)] *
                        d_vit[k * d_jmt * d_imt + j * d_imt + (i + 1)] +
                        (1.0e0 - d_vit[k * d_jmt * d_imt + j * d_imt + (i + 1)]) * wt1,
                        d_at[(k - 1) * d_jmt * d_imt + j * d_imt + i] *
                        d_vit[(k - 1) * d_jmt * d_imt + j * d_imt + i] +
                        (1.0e0 - d_vit[(k - 1) * d_jmt * d_imt + j * d_imt + i]) * wt1),
                                                                  d_at[(k + 1) * d_jmt * d_imt + j * d_imt + i] *
                                                                  d_vit[(k + 1) * d_jmt * d_imt + j * d_imt + i] +
                                                                  (1.0e0 -
                                                                   d_vit[(k + 1) * d_jmt * d_imt + j * d_imt + i]) *
                                                                  wt1);

                d_atmin[k * d_jmt * d_imt + j * d_imt + i] = min2(min1(
                        d_at[k * d_jmt * d_imt + j * d_imt + i] *
                        d_vit[k * d_jmt * d_imt + j * d_imt + i] +
                        (1.0e0 - d_vit[k * d_jmt * d_imt + j * d_imt + i]) * wt2,
                        d_at[k * d_jmt * d_imt + (j - 1) * d_imt + i] *
                        d_vit[k * d_jmt * d_imt + (j - 1) * d_imt + i] +
                        (1.0e0 - d_vit[k * d_jmt * d_imt + (j - 1) * d_imt + i]) * wt2,
                        d_at[k * d_jmt * d_imt + (j + 1) * d_imt + i] *
                        d_vit[k * d_jmt * d_imt + (j + 1) * d_imt + i] +
                        (1.0e0 - d_vit[k * d_jmt * d_imt + (j + 1) * d_imt + i]) * wt2,
                        d_at[k * d_jmt * d_imt + j * d_imt + (i - 1)] *
                        d_vit[k * d_jmt * d_imt + j * d_imt + (i - 1)] +
                        (1.0e0 - d_vit[k * d_jmt * d_imt + j * d_imt + (i - 1)]) * wt2,
                        d_at[k * d_jmt * d_imt + j * d_imt + (i + 1)] *
                        d_vit[k * d_jmt * d_imt + j * d_imt + (i + 1)] +
                        (1.0e0 - d_vit[k * d_jmt * d_imt + j * d_imt + (i + 1)]) * wt2,
                        d_at[(k - 1) * d_jmt * d_imt + j * d_imt + i] *
                        d_vit[(k - 1) * d_jmt * d_imt + j * d_imt + i] +
                        (1.0e0 - d_vit[(k - 1) * d_jmt * d_imt + j * d_imt + i]) * wt2),
                                                                  d_at[(k + 1) * d_jmt * d_imt + j * d_imt + i] *
                                                                  d_vit[(k + 1) * d_jmt * d_imt + j * d_imt + i] +
                                                                  (1.0e0 -
                                                                   d_vit[(k + 1) * d_jmt * d_imt + j * d_imt + i]) *
                                                                  wt2);
            }
        }
    }
}

__global__ void advection_tracer_3(double *d_wkb, double *d_wkd, double *d_u_wface, double *d_v_sface,
                                   double *d_hts, double *d_htw,
                                   double *d_dxu, double *d_dyu, double *d_ws, double *d_at, double *d_adv_tt,
                                   double *d_tarea_r,
                                   double *d_odzp, double *d_ax, double *d_ay, double *d_az, double *d_hun,
                                   double *d_hue, double *d_at00,
                                   double *d_atmax, double *d_atmin, double *d_odzt, double *d_vit, double d_dts,
                                   int d_ny, int d_nx, int d_km,
                                   int d_jmt, int d_imt, int d_ntra, int d_mtracer, int d_nss, int Str) {

    double adv_xy1, adv_xy2, adv_xy3, adv_xy4;
    double adv_zz, adv_za, adv_zb1, adv_zb2;
    double adv_x0, adv_y0, adv_xx, adv_yy;
    double adv_c1, adv_c2, adv_zc;

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    if (k == 0 && 1 <= j && j < (d_jmt - 1) && 1 <= i && i < (d_imt - 1)) {
        if (d_at00[k * d_jmt * d_imt + j * d_imt + i] > d_atmax[k * d_jmt * d_imt + j * d_imt + i] ||
            d_at00[k * d_jmt * d_imt + j * d_imt + i] < d_atmin[k * d_jmt * d_imt + j * d_imt + i] ||
            ((d_at00[k * d_jmt * d_imt + j * d_imt + i + 1] > d_atmax[k * d_jmt * d_imt + j * d_imt + i + 1] ||
              d_at00[k * d_jmt * d_imt + j * d_imt + i + 1] < d_atmin[k * d_jmt * d_imt + j * d_imt + i + 1]) &&
             i <= d_imt - 3)) {
            adv_xy1 = -(d_at[k * d_jmt * d_imt + j * d_imt + i + 1] - d_at[k * d_jmt * d_imt + j * d_imt + i]) *
                      fabs(d_u_wface[k * d_jmt * d_imt + j * d_imt + i + 1]) * d_tarea_r[j * d_nx + i];
        } else {
            adv_xy1 = -d_dts *
                      (d_at[k * d_jmt * d_imt + j * d_imt + i + 1] - d_at[k * d_jmt * d_imt + j * d_imt + i])
                      * 2.0 * d_tarea_r[j * d_nx + i] *
                      d_u_wface[k * d_jmt * d_imt + j * d_imt + i + 1] *
                      d_u_wface[k * d_jmt * d_imt + j * d_imt + i + 1] /
                      (d_htw[j * d_imt + i + 1] * d_hun[j * d_imt + i + 1]);
        }
        if (d_at00[k * d_jmt * d_imt + j * d_imt + i] > d_atmax[k * d_jmt * d_imt + j * d_imt + i] ||
            d_at00[k * d_jmt * d_imt + j * d_imt + i] < d_atmin[k * d_jmt * d_imt + j * d_imt + i] ||
            ((d_at00[k * d_jmt * d_imt + j * d_imt + i - 1] > d_atmax[k * d_jmt * d_imt + j * d_imt + i - 1] ||
              d_at00[k * d_jmt * d_imt + j * d_imt + i - 1] < d_atmin[k * d_jmt * d_imt + j * d_imt + i - 1]) &&
             i >= 2)) {
            adv_xy2 = (d_at[k * d_jmt * d_imt + j * d_imt + i] - d_at[k * d_jmt * d_imt + j * d_imt + i - 1]) *
                      fabs(d_u_wface[k * d_jmt * d_imt + j * d_imt + i]) * d_tarea_r[j * d_nx + i];
        } else {
            adv_xy2 = d_dts *
                      (d_at[k * d_jmt * d_imt + j * d_imt + i] - d_at[k * d_jmt * d_imt + j * d_imt + i - 1])
                      * 2.0 * d_tarea_r[j * d_nx + i] *
                      d_u_wface[k * d_jmt * d_imt + j * d_imt + i] * d_u_wface[k * d_jmt * d_imt + j * d_imt + i] /
                      (d_htw[j * d_imt + i] * d_hun[j * d_imt + i]);
        }
        if (d_at00[k * d_jmt * d_imt + j * d_imt + i] > d_atmax[k * d_jmt * d_imt + j * d_imt + i] ||
            d_at00[k * d_jmt * d_imt + j * d_imt + i] < d_atmin[k * d_jmt * d_imt + j * d_imt + i] ||
            ((d_at00[k * d_jmt * d_imt + (j + 1) * d_imt + i] > d_atmax[k * d_jmt * d_imt + (j + 1) * d_imt + i] ||
              d_at00[k * d_jmt * d_imt + (j + 1) * d_imt + i] < d_atmin[k * d_jmt * d_imt + (j + 1) * d_imt + i]) &&
             j <= d_jmt - 3)) {
            adv_xy3 = -(d_at[k * d_jmt * d_imt + (j + 1) * d_imt + i] - d_at[k * d_jmt * d_imt + j * d_imt + i]) *
                      fabs(d_v_sface[k * d_jmt * d_imt + j * d_imt + i]) * d_tarea_r[j * d_nx + i];
        } else {
            adv_xy3 = -d_dts *
                      (d_at[k * d_jmt * d_imt + (j + 1) * d_imt + i] - d_at[k * d_jmt * d_imt + j * d_imt + i])
                      * 2.0 * d_tarea_r[j * d_nx + i] *
                      d_v_sface[k * d_jmt * d_imt + j * d_imt + i] * d_v_sface[k * d_jmt * d_imt + j * d_imt + i] /
                      (d_hts[j * d_imt + i] * d_hue[j * d_imt + i]);
        }
        if (d_at00[k * d_jmt * d_imt + j * d_imt + i] > d_atmax[k * d_jmt * d_imt + j * d_imt + i] ||
            d_at00[k * d_jmt * d_imt + j * d_imt + i] < d_atmin[k * d_jmt * d_imt + j * d_imt + i] ||
            ((d_at00[k * d_jmt * d_imt + (j - 1) * d_imt + i] > d_atmax[k * d_jmt * d_imt + (j - 1) * d_imt + i] ||
              d_at00[k * d_jmt * d_imt + (j - 1) * d_imt + i] < d_atmin[k * d_jmt * d_imt + (j - 1) * d_imt + i]) &&
             j >= 2)) {
            adv_xy4 = (d_at[k * d_jmt * d_imt + j * d_imt + i] - d_at[k * d_jmt * d_imt + (j - 1) * d_imt + i]) *
                      fabs(d_v_sface[k * d_jmt * d_imt + (j - 1) * d_imt + i]) * d_tarea_r[j * d_nx + i];
        } else {
            adv_xy4 = d_dts *
                      (d_at[k * d_jmt * d_imt + j * d_imt + i] - d_at[k * d_jmt * d_imt + (j - 1) * d_imt + i])
                      * 2.0 * d_tarea_r[j * d_nx + i] *
                      d_v_sface[k * d_jmt * d_imt + (j - 1) * d_imt + i] *
                      d_v_sface[k * d_jmt * d_imt + (j - 1) * d_imt + i] /
                      (d_hts[(j - 1) * d_imt + i] * d_hue[(j - 1) * d_imt + i]);
        }
        adv_zb1 = 0.0;
        if (d_at00[k * d_jmt * d_imt + j * d_imt + i] > d_atmax[k * d_jmt * d_imt + j * d_imt + i] ||
            d_at00[k * d_jmt * d_imt + j * d_imt + i] < d_atmin[k * d_jmt * d_imt + j * d_imt + i] ||
            ((d_at00[(k + 1) * d_jmt * d_imt + j * d_imt + i] > d_atmax[(k + 1) * d_jmt * d_imt + j * d_imt + i] ||
              d_at00[(k + 1) * d_jmt * d_imt + j * d_imt + i] < d_atmin[(k + 1) * d_jmt * d_imt + j * d_imt + i]) &&
             k <= d_km - 2)) {
            adv_zb2 = 0.5 * fabs(d_ws[(k + 1) * d_jmt * d_imt + j * d_imt + i]) * d_odzp[k] *
                      (d_at[k * d_jmt * d_imt + j * d_imt + i] - d_at[(k + 1) * d_jmt * d_imt + j * d_imt + i]);
        } else {
            adv_zb2 = 0.5 * d_odzp[0] *
                      d_ws[1 * d_jmt * d_imt + j * d_imt + i] * d_ws[1 * d_jmt * d_imt + j * d_imt + i] *
                      d_odzt[1] *
                      (d_at[j * d_imt + i] - d_at[1 * d_jmt * d_imt + j * d_imt + i]) * d_dts;
        }
        adv_za = -0.5 * d_odzp[0] * d_ws[1 * d_jmt * d_imt + j * d_imt + i] *
                 (d_at[1 * d_jmt * d_imt + j * d_imt + i] + d_at[j * d_imt + i]);

        adv_zc = d_odzp[0] * d_at[j * d_imt + i] * d_ws[1 * d_jmt * d_imt + j * d_imt + i];
        adv_c1 = -d_at[k * d_jmt * d_imt + j * d_imt + i] *
                 (d_u_wface[k * d_jmt * d_imt + j * d_imt + i + 1] - d_u_wface[k * d_jmt * d_imt + j * d_imt + i])
                 * d_tarea_r[j * d_nx + i] * 2.0;
        adv_c2 = -d_at[k * d_jmt * d_imt + j * d_imt + i] *
                 (d_v_sface[k * d_jmt * d_imt + j * d_imt + i] - d_v_sface[k * d_jmt * d_imt + (j - 1) * d_imt + i])
                 * d_tarea_r[j * d_nx + i] * 2.0;
        adv_x0 = (d_at[k * d_jmt * d_imt + j * d_imt + i + 1] + d_at[k * d_jmt * d_imt + j * d_imt + i]) *
                 d_u_wface[k * d_jmt * d_imt + j * d_imt + i + 1] * d_tarea_r[j * d_nx + i] -
                 (d_at[k * d_jmt * d_imt + j * d_imt + i] + d_at[k * d_jmt * d_imt + j * d_imt + i - 1]) *
                 d_u_wface[k * d_jmt * d_imt + j * d_imt + i] * d_tarea_r[j * d_nx + i];
        adv_y0 = ((d_at[k * d_jmt * d_imt + (j + 1) * d_imt + i] + d_at[k * d_jmt * d_imt + j * d_imt + i]) *
                  d_v_sface[k * d_jmt * d_imt + j * d_imt + i] -
                  (d_at[k * d_jmt * d_imt + j * d_imt + i] + d_at[k * d_jmt * d_imt + (j - 1) * d_imt + i]) *
                  d_v_sface[k * d_jmt * d_imt + (j - 1) * d_imt + i]) * d_tarea_r[j * d_nx + i];
        adv_xx = -(adv_x0 + adv_xy1 + adv_xy2 + adv_c1);
        adv_yy = -(adv_y0 + adv_xy3 + adv_xy4 + adv_c2);
        adv_zz = -(adv_za + adv_zb1 + adv_zb2 + adv_zc);

        d_adv_tt[k * d_jmt * d_imt + j * d_imt + i] = adv_xx + adv_yy + adv_zz;
        d_ax[d_mtracer * d_km * d_jmt * d_imt + k * d_jmt * d_imt + j * d_imt + i] = adv_xx;
        d_ay[d_mtracer * d_km * d_jmt * d_imt + k * d_jmt * d_imt + j * d_imt + i] = adv_yy;
        d_az[d_mtracer * d_km * d_jmt * d_imt + k * d_jmt * d_imt + j * d_imt + i] = adv_zz;
    }

    if (k >= 1 && k < d_km - 1) {
        if (1 <= j && j < (d_jmt - 1) && 1 <= i && i < (d_imt - 1)) {
            if ((d_at00[k * d_jmt * d_imt + j * d_imt + i] > d_atmax[k * d_jmt * d_imt + j * d_imt + i] ||
                 d_at00[k * d_jmt * d_imt + j * d_imt + i] < d_atmin[k * d_jmt * d_imt + j * d_imt + i] ||
                 d_at00[k * d_jmt * d_imt + j * d_imt + i + 1] > d_atmax[k * d_jmt * d_imt + j * d_imt + i + 1] ||
                 d_at00[k * d_jmt * d_imt + j * d_imt + i + 1] < d_atmin[k * d_jmt * d_imt + j * d_imt + i + 1]) &&
                i <= d_imt - 3) {
                adv_xy1 =
                        -(d_at[k * d_jmt * d_imt + j * d_imt + i + 1] - d_at[k * d_jmt * d_imt + j * d_imt + i]) *
                        fabs(d_u_wface[k * d_jmt * d_imt + j * d_imt + i + 1]) * d_tarea_r[j * d_nx + i];
            } else {
                adv_xy1 = -d_dts *
                          (d_at[k * d_jmt * d_imt + j * d_imt + i + 1] - d_at[k * d_jmt * d_imt + j * d_imt + i])
                          * 2.0 * d_tarea_r[j * d_nx + i] *
                          d_u_wface[k * d_jmt * d_imt + j * d_imt + i + 1] *
                          d_u_wface[k * d_jmt * d_imt + j * d_imt + i + 1] /
                          (d_htw[j * d_imt + i + 1] * d_hun[j * d_imt + i + 1]);
            }
            if ((d_at00[k * d_jmt * d_imt + j * d_imt + i] > d_atmax[k * d_jmt * d_imt + j * d_imt + i] ||
                 d_at00[k * d_jmt * d_imt + j * d_imt + i] < d_atmin[k * d_jmt * d_imt + j * d_imt + i] ||
                 d_at00[k * d_jmt * d_imt + j * d_imt + i - 1] > d_atmax[k * d_jmt * d_imt + j * d_imt + i - 1] ||
                 d_at00[k * d_jmt * d_imt + j * d_imt + i - 1] < d_atmin[k * d_jmt * d_imt + j * d_imt + i - 1]) &&
                i >= 2) {
                adv_xy2 =
                        (d_at[k * d_jmt * d_imt + j * d_imt + i] - d_at[k * d_jmt * d_imt + j * d_imt + i - 1]) *
                        fabs(d_u_wface[k * d_jmt * d_imt + j * d_imt + i]) * d_tarea_r[j * d_nx + i];
            } else {
                adv_xy2 = d_dts *
                          (d_at[k * d_jmt * d_imt + j * d_imt + i] - d_at[k * d_jmt * d_imt + j * d_imt + i - 1])
                          * 2.0 * d_tarea_r[j * d_nx + i] *
                          d_u_wface[k * d_jmt * d_imt + j * d_imt + i] *
                          d_u_wface[k * d_jmt * d_imt + j * d_imt + i] /
                          (d_htw[j * d_imt + i] * d_hun[j * d_imt + i]);
            }
            if ((d_at00[k * d_jmt * d_imt + j * d_imt + i] > d_atmax[k * d_jmt * d_imt + j * d_imt + i] ||
                 d_at00[k * d_jmt * d_imt + j * d_imt + i] < d_atmin[k * d_jmt * d_imt + j * d_imt + i] ||
                 d_at00[k * d_jmt * d_imt + (j + 1) * d_imt + i] >
                 d_atmax[k * d_jmt * d_imt + (j + 1) * d_imt + i] ||
                 d_at00[k * d_jmt * d_imt + (j + 1) * d_imt + i] <
                 d_atmin[k * d_jmt * d_imt + (j + 1) * d_imt + i]) && j <= d_jmt - 3) {
                adv_xy3 =
                        -(d_at[k * d_jmt * d_imt + (j + 1) * d_imt + i] -
                          d_at[k * d_jmt * d_imt + j * d_imt + i]) *
                        fabs(d_v_sface[k * d_jmt * d_imt + j * d_imt + i]) * d_tarea_r[j * d_nx + i];
            } else {
                adv_xy3 = -d_dts *
                          (d_at[k * d_jmt * d_imt + (j + 1) * d_imt + i] -
                           d_at[k * d_jmt * d_imt + j * d_imt + i])
                          * 2.0 * d_tarea_r[j * d_nx + i] *
                          d_v_sface[k * d_jmt * d_imt + j * d_imt + i] *
                          d_v_sface[k * d_jmt * d_imt + j * d_imt + i] /
                          (d_hts[j * d_imt + i] * d_hue[j * d_imt + i]);
            }
            if ((d_at00[k * d_jmt * d_imt + j * d_imt + i] > d_atmax[k * d_jmt * d_imt + j * d_imt + i] ||
                 d_at00[k * d_jmt * d_imt + j * d_imt + i] < d_atmin[k * d_jmt * d_imt + j * d_imt + i] ||
                 d_at00[k * d_jmt * d_imt + (j - 1) * d_imt + i] >
                 d_atmax[k * d_jmt * d_imt + (j - 1) * d_imt + i] ||
                 d_at00[k * d_jmt * d_imt + (j - 1) * d_imt + i] <
                 d_atmin[k * d_jmt * d_imt + (j - 1) * d_imt + i]) && j >= 2) {
                adv_xy4 =
                        (d_at[k * d_jmt * d_imt + j * d_imt + i] -
                         d_at[k * d_jmt * d_imt + (j - 1) * d_imt + i]) *
                        fabs(d_v_sface[k * d_jmt * d_imt + (j - 1) * d_imt + i]) * d_tarea_r[j * d_nx + i];
            } else {
                adv_xy4 = d_dts *
                          (d_at[k * d_jmt * d_imt + j * d_imt + i] -
                           d_at[k * d_jmt * d_imt + (j - 1) * d_imt + i])
                          * 2.0 * d_tarea_r[j * d_nx + i] *
                          d_v_sface[k * d_jmt * d_imt + (j - 1) * d_imt + i] *
                          d_v_sface[k * d_jmt * d_imt + (j - 1) * d_imt + i] /
                          (d_hts[(j - 1) * d_imt + i] * d_hue[(j - 1) * d_imt + i]);
            }
            if ((d_at00[k * d_jmt * d_imt + j * d_imt + i] > d_atmax[k * d_jmt * d_imt + j * d_imt + i] ||
                 d_at00[k * d_jmt * d_imt + j * d_imt + i] < d_atmin[k * d_jmt * d_imt + j * d_imt + i] ||
                 d_at00[(k - 1) * d_jmt * d_imt + j * d_imt + i] >
                 d_atmax[(k - 1) * d_jmt * d_imt + j * d_imt + i] ||
                 d_at00[(k - 1) * d_jmt * d_imt + j * d_imt + i] <
                 d_atmin[(k - 1) * d_jmt * d_imt + j * d_imt + i]) && k >= 1) {
                adv_zb1 = -0.5 * fabs(d_ws[k * d_jmt * d_imt + j * d_imt + i]) * d_odzp[k] *
                          (d_at[(k - 1) * d_jmt * d_imt + j * d_imt + i] -
                           d_at[k * d_jmt * d_imt + j * d_imt + i]);
            } else {
                adv_zb1 = -0.5 * d_odzp[k] *
                          d_ws[k * d_jmt * d_imt + j * d_imt + i] * d_ws[k * d_jmt * d_imt + j * d_imt + i] *
                          d_odzt[k] *
                          (d_at[(k - 1) * d_jmt * d_imt + j * d_imt + i] -
                           d_at[k * d_jmt * d_imt + j * d_imt + i])
                          * d_dts;
            }
            if ((d_at00[k * d_jmt * d_imt + j * d_imt + i] > d_atmax[k * d_jmt * d_imt + j * d_imt + i] ||
                 d_at00[k * d_jmt * d_imt + j * d_imt + i] < d_atmin[k * d_jmt * d_imt + j * d_imt + i] ||
                 d_at00[(k + 1) * d_jmt * d_imt + j * d_imt + i] >
                 d_atmax[(k + 1) * d_jmt * d_imt + j * d_imt + i] ||
                 d_at00[(k + 1) * d_jmt * d_imt + j * d_imt + i] <
                 d_atmin[(k + 1) * d_jmt * d_imt + j * d_imt + i]) && k <= d_km - 2) {
                adv_zb2 = 0.5 * fabs(d_ws[(k + 1) * d_jmt * d_imt + j * d_imt + i]) * d_odzp[k] *
                          (d_at[k * d_jmt * d_imt + j * d_imt + i] -
                           d_at[(k + 1) * d_jmt * d_imt + j * d_imt + i]);
            } else {
                adv_zb2 = 0.5 * d_odzp[k] *
                          d_ws[(k + 1) * d_jmt * d_imt + j * d_imt + i] *
                          d_ws[(k + 1) * d_jmt * d_imt + j * d_imt + i] *
                          d_odzt[k + 1] *
                          (d_at[k * d_jmt * d_imt + j * d_imt + i] -
                           d_at[(k + 1) * d_jmt * d_imt + j * d_imt + i])
                          * d_dts;
            }

            adv_c1 = -d_at[k * d_jmt * d_imt + j * d_imt + i] *
                     (d_u_wface[k * d_jmt * d_imt + j * d_imt + i + 1] -
                      d_u_wface[k * d_jmt * d_imt + j * d_imt + i])
                     * d_tarea_r[j * d_nx + i] * 2.0;

            adv_c2 = -d_at[k * d_jmt * d_imt + j * d_imt + i] *
                     (d_v_sface[k * d_jmt * d_imt + j * d_imt + i] -
                      d_v_sface[k * d_jmt * d_imt + (j - 1) * d_imt + i])
                     * d_tarea_r[j * d_nx + i] * 2.0;

            adv_za = 0.5 * d_odzp[k] * d_ws[k * d_jmt * d_imt + j * d_imt + i] *
                     (d_at[k * d_jmt * d_imt + j * d_imt + i] + d_at[(k - 1) * d_jmt * d_imt + j * d_imt + i])
                     - 0.5 * d_odzp[k] * d_ws[(k + 1) * d_jmt * d_imt + j * d_imt + i] *
                       (d_at[k * d_jmt * d_imt + j * d_imt + i] + d_at[(k + 1) * d_jmt * d_imt + j * d_imt + i]);

            adv_zc = -d_odzp[k] * d_at[k * d_jmt * d_imt + j * d_imt + i] *
                     (d_ws[k * d_jmt * d_imt + j * d_imt + i] - d_ws[(k + 1) * d_jmt * d_imt + j * d_imt + i]);

            adv_x0 = (d_at[k * d_jmt * d_imt + j * d_imt + i + 1] + d_at[k * d_jmt * d_imt + j * d_imt + i]) *
                     d_u_wface[k * d_jmt * d_imt + j * d_imt + i + 1] * d_tarea_r[j * d_nx + i] -
                     (d_at[k * d_jmt * d_imt + j * d_imt + i] + d_at[k * d_jmt * d_imt + j * d_imt + i - 1]) *
                     d_u_wface[k * d_jmt * d_imt + j * d_imt + i] * d_tarea_r[j * d_nx + i];

            adv_y0 = ((d_at[k * d_jmt * d_imt + (j + 1) * d_imt + i] + d_at[k * d_jmt * d_imt + j * d_imt + i]) *
                      d_v_sface[k * d_jmt * d_imt + j * d_imt + i] -
                      (d_at[k * d_jmt * d_imt + j * d_imt + i] + d_at[k * d_jmt * d_imt + (j - 1) * d_imt + i]) *
                      d_v_sface[k * d_jmt * d_imt + (j - 1) * d_imt + i]) * d_tarea_r[j * d_nx + i];

            adv_xx = -(adv_x0 + adv_xy1 + adv_xy2 + adv_c1);
            adv_yy = -(adv_y0 + adv_xy3 + adv_xy4 + adv_c2);
            adv_zz = -(adv_za + adv_zb1 + adv_zb2 + adv_zc);
            d_adv_tt[k * d_jmt * d_imt + j * d_imt + i] = adv_xx + adv_yy + adv_zz;
            d_ax[d_mtracer * d_km * d_jmt * d_imt + k * d_jmt * d_imt + j * d_imt + i] = adv_xx;
            d_ay[d_mtracer * d_km * d_jmt * d_imt + k * d_jmt * d_imt + j * d_imt + i] = adv_yy;
            d_az[d_mtracer * d_km * d_jmt * d_imt + k * d_jmt * d_imt + j * d_imt + i] = adv_zz;
        }
    }
    if (k == d_km - 1 && 1 <= j && j < (d_jmt - 1) && 1 <= i && i < (d_imt - 1)) {
        if (d_at00[k * d_jmt * d_imt + j * d_imt + i] > d_atmax[k * d_jmt * d_imt + j * d_imt + i] ||
            d_at00[k * d_jmt * d_imt + j * d_imt + i] < d_atmin[k * d_jmt * d_imt + j * d_imt + i] ||
            ((d_at00[k * d_jmt * d_imt + j * d_imt + i + 1] > d_atmax[k * d_jmt * d_imt + j * d_imt + i + 1] ||
              d_at00[k * d_jmt * d_imt + j * d_imt + i + 1] < d_atmin[k * d_jmt * d_imt + j * d_imt + i + 1]) &&
             i <= d_imt - 3)) {
            adv_xy1 = -
                              (d_at[k * d_jmt * d_imt + j * d_imt + i + 1] -
                               d_at[k * d_jmt * d_imt + j * d_imt + i]) *
                      fabs(d_u_wface[k * d_jmt * d_imt + j * d_imt + i + 1]) * d_tarea_r[j * d_nx + i];
        } else {
            adv_xy1 = -d_dts *
                      (d_at[k * d_jmt * d_imt + j * d_imt + i + 1] - d_at[k * d_jmt * d_imt + j * d_imt + i])
                      * 2.0 * d_tarea_r[j * d_nx + i] *
                      d_u_wface[k * d_jmt * d_imt + j * d_imt + i + 1] *
                      d_u_wface[k * d_jmt * d_imt + j * d_imt + i + 1] /
                      (d_htw[j * d_imt + i + 1] * d_hun[j * d_imt + i + 1]);
        }
        if (d_at00[k * d_jmt * d_imt + j * d_imt + i] > d_atmax[k * d_jmt * d_imt + j * d_imt + i] ||
            d_at00[k * d_jmt * d_imt + j * d_imt + i] < d_atmin[k * d_jmt * d_imt + j * d_imt + i] ||
            ((d_at00[k * d_jmt * d_imt + j * d_imt + i - 1] > d_atmax[k * d_jmt * d_imt + j * d_imt + i - 1] ||
              d_at00[k * d_jmt * d_imt + j * d_imt + i - 1] < d_atmin[k * d_jmt * d_imt + j * d_imt + i - 1]) &&
             i >= 2)) {
            adv_xy2 =
                    (d_at[k * d_jmt * d_imt + j * d_imt + i] - d_at[k * d_jmt * d_imt + j * d_imt + i - 1]) *
                    fabs(d_u_wface[k * d_jmt * d_imt + j * d_imt + i]) * d_tarea_r[j * d_nx + i];
        } else {
            adv_xy2 = d_dts *
                      (d_at[k * d_jmt * d_imt + j * d_imt + i] - d_at[k * d_jmt * d_imt + j * d_imt + i - 1])
                      * 2.0 * d_tarea_r[j * d_nx + i] *
                      d_u_wface[k * d_jmt * d_imt + j * d_imt + i] * d_u_wface[k * d_jmt * d_imt + j * d_imt + i] /
                      (d_htw[j * d_imt + i] * d_hun[j * d_imt + i]);
        }
        if (d_at00[k * d_jmt * d_imt + j * d_imt + i] > d_atmax[k * d_jmt * d_imt + j * d_imt + i] ||
            d_at00[k * d_jmt * d_imt + j * d_imt + i] < d_atmin[k * d_jmt * d_imt + j * d_imt + i] ||
            ((d_at00[k * d_jmt * d_imt + (j + 1) * d_imt + i] > d_atmax[k * d_jmt * d_imt + (j + 1) * d_imt + i] ||
              d_at00[k * d_jmt * d_imt + (j + 1) * d_imt + i] < d_atmin[k * d_jmt * d_imt + (j + 1) * d_imt + i]) &&
             j <= d_jmt - 3)) {
            adv_xy3 = -
                              (d_at[k * d_jmt * d_imt + (j + 1) * d_imt + i] -
                               d_at[k * d_jmt * d_imt + j * d_imt + i]) *
                      fabs(d_v_sface[k * d_jmt * d_imt + j * d_imt + i]) * d_tarea_r[j * d_nx + i];
        } else {
            adv_xy3 = -d_dts *
                      (d_at[k * d_jmt * d_imt + (j + 1) * d_imt + i] - d_at[k * d_jmt * d_imt + j * d_imt + i])
                      * 2.0 * d_tarea_r[j * d_nx + i] *
                      d_v_sface[k * d_jmt * d_imt + j * d_imt + i] * d_v_sface[k * d_jmt * d_imt + j * d_imt + i] /
                      (d_hts[j * d_imt + i] * d_hue[j * d_imt + i]);
        }
        if (d_at00[k * d_jmt * d_imt + j * d_imt + i] > d_atmax[k * d_jmt * d_imt + j * d_imt + i] ||
            d_at00[k * d_jmt * d_imt + j * d_imt + i] < d_atmin[k * d_jmt * d_imt + j * d_imt + i] ||
            ((d_at00[k * d_jmt * d_imt + (j - 1) * d_imt + i] > d_atmax[k * d_jmt * d_imt + (j - 1) * d_imt + i] ||
              d_at00[k * d_jmt * d_imt + (j - 1) * d_imt + i] < d_atmin[k * d_jmt * d_imt + (j - 1) * d_imt + i]) &&
             j >= 2)) {
            adv_xy4 =
                    (d_at[k * d_jmt * d_imt + j * d_imt + i] - d_at[k * d_jmt * d_imt + (j - 1) * d_imt + i]) *
                    fabs(d_v_sface[k * d_jmt * d_imt + (j - 1) * d_imt + i]) * d_tarea_r[j * d_nx + i];
        } else {
            adv_xy4 = d_dts *
                      (d_at[k * d_jmt * d_imt + j * d_imt + i] - d_at[k * d_jmt * d_imt + (j - 1) * d_imt + i])
                      * 2.0 * d_tarea_r[j * d_nx + i] *
                      d_v_sface[k * d_jmt * d_imt + (j - 1) * d_imt + i] *
                      d_v_sface[k * d_jmt * d_imt + (j - 1) * d_imt + i] /
                      (d_hts[(j - 1) * d_imt + i] * d_hue[(j - 1) * d_imt + i]);
        }
        if (d_at00[k * d_jmt * d_imt + j * d_imt + i] > d_atmax[k * d_jmt * d_imt + j * d_imt + i] ||
            d_at00[k * d_jmt * d_imt + j * d_imt + i] < d_atmin[k * d_jmt * d_imt + j * d_imt + i] ||
            d_at00[(k - 1) * d_jmt * d_imt + j * d_imt + i] > d_atmax[(k - 1) * d_jmt * d_imt + j * d_imt + i] ||
            d_at00[(k - 1) * d_jmt * d_imt + j * d_imt + i] < d_atmin[(k - 1) * d_jmt * d_imt + j * d_imt + i]) {
            adv_zb1 = -0.5 * fabs(d_ws[k * d_jmt * d_imt + j * d_imt + i]) * d_odzp[k] *
                      (d_at[(k - 1) * d_jmt * d_imt + j * d_imt + i] - d_at[k * d_jmt * d_imt + j * d_imt + i]);
        } else {
            adv_zb1 = -0.5 * d_odzp[k] *
                      d_ws[k * d_jmt * d_imt + j * d_imt + i] * d_ws[k * d_jmt * d_imt + j * d_imt + i] *
                      d_odzt[k] *
                      (d_at[(k - 1) * d_jmt * d_imt + j * d_imt + i] - d_at[k * d_jmt * d_imt + j * d_imt + i]) *
                      d_dts;
        }
        adv_zb2 = 0.0;
        adv_c1 = -d_at[k * d_jmt * d_imt + j * d_imt + i] *
                 (d_u_wface[k * d_jmt * d_imt + j * d_imt + i + 1] - d_u_wface[k * d_jmt * d_imt + j * d_imt + i])
                 * d_tarea_r[j * d_nx + i] * 2.0;
        adv_c2 = -d_at[k * d_jmt * d_imt + j * d_imt + i] *
                 (d_v_sface[k * d_jmt * d_imt + j * d_imt + i] - d_v_sface[k * d_jmt * d_imt + (j - 1) * d_imt + i])
                 * d_tarea_r[j * d_nx + i] * 2.0;
        adv_za = 0.5 * d_odzp[km - 1] * d_ws[(km - 1) * d_jmt * d_imt + j * d_imt + i] *
                 (d_at[(km - 1) * d_jmt * d_imt + j * d_imt + i] +
                  d_at[(km - 2) * d_jmt * d_imt + j * d_imt + i]);
        adv_zc = -d_odzp[km - 1] * d_at[(km - 1) * d_jmt * d_imt + j * d_imt + i] *
                 d_ws[(km - 1) * d_jmt * d_imt + j * d_imt + i];
        adv_x0 = (d_at[k * d_jmt * d_imt + j * d_imt + i + 1] + d_at[k * d_jmt * d_imt + j * d_imt + i]) *
                 d_u_wface[k * d_jmt * d_imt + j * d_imt + i + 1] * d_tarea_r[j * d_nx + i] -
                 (d_at[k * d_jmt * d_imt + j * d_imt + i] + d_at[k * d_jmt * d_imt + j * d_imt + i - 1]) *
                 d_u_wface[k * d_jmt * d_imt + j * d_imt + i] * d_tarea_r[j * d_nx + i];
        adv_y0 = ((d_at[k * d_jmt * d_imt + (j + 1) * d_imt + i] + d_at[k * d_jmt * d_imt + j * d_imt + i]) *
                  d_v_sface[k * d_jmt * d_imt + j * d_imt + i] -
                  (d_at[k * d_jmt * d_imt + j * d_imt + i] + d_at[k * d_jmt * d_imt + (j - 1) * d_imt + i]) *
                  d_v_sface[k * d_jmt * d_imt + (j - 1) * d_imt + i]) * d_tarea_r[j * d_nx + i];
        adv_xx = -(adv_x0 + adv_xy1 + adv_xy2 + adv_c1);
        adv_yy = -(adv_y0 + adv_xy3 + adv_xy4 + adv_c2);
        adv_zz = -(adv_za + adv_zb1 + adv_zb2 + adv_zc);

        d_adv_tt[k * d_jmt * d_imt + j * d_imt + i] = adv_xx + adv_yy + adv_zz;
        d_ax[d_mtracer * d_km * d_jmt * d_imt + k * d_jmt * d_imt + j * d_imt + i] = adv_xx;
        d_ay[d_mtracer * d_km * d_jmt * d_imt + k * d_jmt * d_imt + j * d_imt + i] = adv_yy;
        d_az[d_mtracer * d_km * d_jmt * d_imt + k * d_jmt * d_imt + j * d_imt + i] = adv_zz;
    }
}

__global__ void invtrit(int d_nx, int d_ny, int d_imt, int d_jmt, int d_km,
                        int *d_kmt, double *d_wkc, double *d_vtl, double *d_stf, double *d_vit,
                        double *d_odzt, double *d_odzp, double d_dts) {
    int i, j, k;
    i = (blockIdx.x) * blockDim.x + threadIdx.x;
    j = (blockIdx.y) * blockDim.y + threadIdx.y;

    double a_8[km], b_8[km], c_8[km], d_8[km];
    double e_8[km + 1], f_8[km + 1];
    double g0;
    int kz;

    if (j >= 2 && j < d_jmt - 2 && i >= 2 && i < d_imt - 2) {
        if (d_kmt[j * d_nx + i] > 0) {
            kz = d_kmt[j * d_nx + i];
            for (k = 1; k < kz; k++) {
                a_8[k] = d_wkc[(k - 1) * d_jmt * d_imt + j * d_imt + i] *
                         d_odzt[k] * d_odzp[k] * d_dts * 0.5e0;

                d_8[k] = d_vtl[k * d_jmt * d_imt + j * d_imt + i];
            }
            for (k = 1; k < kz - 1; k++) {
                c_8[k] = d_wkc[k * d_jmt * d_imt + j * d_imt + i] *
                         d_odzt[k + 1] * d_odzp[k] * d_dts * 0.5e0;

                b_8[k] = 1.0 + a_8[k] + c_8[k];
                e_8[k] = 0.0;
                f_8[k] = 0.0;
            }
            k = 0;
            a_8[k] = d_odzp[k] * d_dts * 0.5e0;
            c_8[k] = d_wkc[k * d_jmt * d_imt + j * d_imt + i] * d_odzt[k + 1] * d_odzp[k] * d_dts * 0.5e0;
            b_8[k] = 1.0 + c_8[k];
            d_8[k] = d_vtl[k * d_jmt * d_imt + j * d_imt + i];
            e_8[k] = 0.0;
            f_8[k] = 0.0;
            b_8[kz - 1] = 1.0 + a_8[kz - 1];
            c_8[kz - 1] = d_odzp[kz - 1] * d_dts * 0.5e0;
            e_8[kz] = 0.0;
            f_8[kz] = 0.0;
            for (k = kz - 1; k >= 0; k--) {
                g0 = 1.0 / (b_8[k] - c_8[k] * e_8[k + 1]);
                e_8[k] = a_8[k] * g0;
                f_8[k] = (d_8[k] + c_8[k] * f_8[k + 1]) * g0;
            }
            d_vtl[j * d_imt + i] = (e_8[0] * d_stf[j * d_imt + i] + f_8[0])
                                   * d_vit[j * d_imt + i];
            for (k = 1; k < kz; k++) {
                d_vtl[k * d_jmt * d_imt + j * d_imt + i] =
                        (e_8[k] * d_vtl[(k - 1) * d_jmt * d_imt + j * d_imt + i] + f_8[k]) *
                        d_vit[k * d_jmt * d_imt + j * d_imt + i];
            }
        }
    }

}

__global__ void operators_div(int d_imt, int d_jmt, int d_km, int d_nx, int d_ny,
                              double *d_wka, int *d_kmt, double *d_uk, double *d_vk, double *d_htw,
                              double *d_hts, double *d_tarea_r) {
    int i, j, k;
    i = (blockIdx.x) * blockDim.x + threadIdx.x;
    j = (blockIdx.y) * blockDim.y + threadIdx.y;
    k = (blockIdx.z) * blockDim.z + threadIdx.z;

    if (k < d_km) {
        if (j < d_ny && i < d_nx) {
            d_wka[k * d_jmt * d_imt + j * d_imt + i] = 0.0e0;
        }
        if (j >= 1 && j < d_ny && i >= 0 && i < d_nx - 1 && k <= (d_kmt[j * d_nx + i] - 1)) {
            d_wka[k * d_jmt * d_imt + j * d_imt + i] =
                    0.5e0 *
                    ((d_uk[k * d_jmt * d_imt + j * d_imt + i + 1] +
                      d_uk[k * d_jmt * d_imt + (j - 1) * d_imt + i + 1]) * d_htw[j * d_nx + i + 1] -
                     (d_uk[k * d_jmt * d_imt + j * d_imt + i] +
                      d_uk[k * d_jmt * d_imt + (j - 1) * d_imt + i]) * d_htw[j * d_nx + i] +
                     (d_vk[k * d_jmt * d_imt + j * d_imt + i + 1] +
                      d_vk[k * d_jmt * d_imt + j * d_imt + i]) * d_hts[j * d_nx + i] -
                     (d_vk[k * d_jmt * d_imt + (j - 1) * d_imt + i + 1] +
                      d_vk[k * d_jmt * d_imt + (j - 1) * d_imt + i]) * d_hts[(j - 1) * d_nx + i]) *
                    d_tarea_r[j * d_nx + i];
        }
    }
}

__global__ void upwell_2(int d_imt, int d_jmt, int d_km,
                         double *d_work, double *d_uk, double *d_vk, double *d_wkd, double *d_wkb,
                         double *d_ohbu) {
    int i, j, k;
    i = (blockIdx.x) * blockDim.x + threadIdx.x;
    j = (blockIdx.y) * blockDim.y + threadIdx.y;
    k = (blockIdx.z) * blockDim.z + threadIdx.z;

    if (k >= 0 && k < d_km) {
        if (j >= 0 && j < d_jmt - 1 && i >= 1 && i < d_imt) {
            d_uk[k * d_jmt * d_imt + j * d_imt + i] =
                    (1.0e0 + d_work[j * d_imt + i] * d_ohbu[j * d_imt + i]) *
                    d_wkd[k * d_jmt * d_imt + j * d_imt + i];

            d_vk[k * d_jmt * d_imt + j * d_imt + i] =
                    (1.0e0 + d_work[j * d_imt + i] * d_ohbu[j * d_imt + i]) *
                    d_wkb[k * d_jmt * d_imt + j * d_imt + i];
        }
    }
}

__global__ void upwell_3(int d_imt, int d_jmt, int d_km,
                         double *d_work, double *d_wka, double *d_vit, double *d_ws, double *d_ohbt,
                         double *d_stf, double *d_dzp) {
    int i, j, k;
    i = (blockIdx.x) * blockDim.x + threadIdx.x;
    j = (blockIdx.y) * blockDim.y + threadIdx.y;

    if (j < d_jmt && i < d_imt) {
        d_work[j * d_imt + i] = 0.0e0;
    }

    for (k = 0; k < d_km; k++) {
        if (j >= 1 && j < d_jmt - 1 && i >= 1 && i < d_imt - 1) {
            d_work[j * d_imt + i] = d_work[j * d_imt + i] -
                                    d_dzp[k] *
                                    d_wka[k * d_jmt * d_imt + j * d_imt + i] *
                                    d_vit[k * d_jmt * d_imt + j * d_imt + i];
        }
    }

    for (k = 1; k < d_km; k++) {
        if (j >= 1 && j < d_jmt - 1 && i >= 1 && i < d_imt - 1) {
            d_ws[k * d_jmt * imt + j * d_imt + i] =
                    d_vit[k * d_jmt * d_imt + j * d_imt + i] *
                    (d_ws[(k - 1) * d_jmt * d_imt + j * d_imt + i] +
                     d_dzp[k - 1] *
                     (d_work[j * d_imt + i] *
                      d_ohbt[j * d_imt + i] +
                      d_wka[(k - 1) * d_jmt * d_imt + j * d_imt + i]));
        }
    }

    if (j >= 1 && j < d_jmt - 1 && i >= 1 && i < d_imt - 1) {
        d_work[j * d_imt + i] = 1.0e0 / (1.0e0 + d_stf[j * d_imt + i] * d_ohbt[j * d_imt + i]);
    }
}

__global__ void upwell_4(int d_imt, int d_jmt, int d_km,
                         double *d_work, double *d_wka, double *d_ws, double *d_ohbt, double *d_stf,
                         double *d_dzp) {
    int i, j, k;
    i = (blockIdx.x) * blockDim.x + threadIdx.x;
    j = (blockIdx.y) * blockDim.y + threadIdx.y;
    k = (blockIdx.z) * blockDim.z + threadIdx.z;

    if (k >= 1 && k < d_km) {
        if (1 <= j && j < d_jmt - 1 && 1 <= i && i < d_imt - 1) {
            d_ws[k * d_jmt * d_imt + j * d_imt + i] =
                    d_ws[k * d_jmt * d_imt + j * d_imt + i] * d_work[j * d_imt + i];
        }
    }
}
__global__ void tracer_1_0(double *d_h0f, double *d_h0l, double *d_stf, double d_onbc,int d_imt, int d_jmt){
    int i, j, k;
    i = (blockIdx.x) * blockDim.x + threadIdx.x;
    j = (blockIdx.y) * blockDim.y + threadIdx.y;
//    k = (blockIdx.z) * blockDim.z + threadIdx.z;

        if (j >= 0 && j < d_jmt && i >= 0 && i < d_imt) {
            d_h0f[j * d_imt + i] = d_h0f[j * d_imt + i] * d_onbc;

            d_stf[j * d_imt + i] = 0.5e0 * d_h0f[j * d_imt + i] +
                                   0.5e0 * d_h0l[j * d_imt + i];
        }
    }


__global__ void tracer_1(int d_nx, int d_ny, int d_imt, int d_jmt, int d_km,
                         double *d_stf, double *d_utf, double *d_vtf, double *d_utl, double *d_vtl,
                         double *d_work, double *d_wka, double *d_wkb, double *d_wkd,
                         double *d_au0, double *d_aus, double *d_auw, double *d_ausw, double d_oncc) {
    int i, j, k;
    i = (blockIdx.x) * blockDim.x + threadIdx.x;
    j = (blockIdx.y) * blockDim.y + threadIdx.y;
    k = (blockIdx.z) * blockDim.z + threadIdx.z;

    if (k >= 0 && k < d_km) {
        if (j >= 0 && j < d_jmt && i >= 0 && i < d_imt) {
            d_utf[k * d_jmt * d_imt + j * d_imt + i] *= d_oncc;
            d_vtf[k * d_jmt * d_imt + j * d_imt + i] *= d_oncc;
        }
    }

    if (k >= 0 && k < d_km) {
        if (j >= 0 && j < d_jmt && i >= 0 && i < d_imt) {
            d_wkd[k * d_jmt * d_imt + j * d_imt + i] = 0.5e0 * d_utf[k * d_jmt * d_imt + j * d_imt + i] +
                                                       0.5e0 * d_utl[k * d_jmt * d_imt + j * d_imt + i];

            d_wkb[k * d_jmt * d_imt + j * d_imt + i] = 0.5e0 * d_vtf[k * d_jmt * d_imt + j * d_imt + i] +
                                                       0.5e0 * d_vtl[k * d_jmt * d_imt + j * d_imt + i];
        }
    }

    //upwell_(work_mod_mp_wkd_, work_mod_mp_wkb_, work_mod_mp_stf_);
    if (k == 0 && j >= 0 && j < d_jmt && i >= 0 && i < d_imt) {
        d_work[j * d_imt + i] = 0.0e0;
    }
    if (k >= 0 && k < d_km) {
        if (j >= 0 && j < d_jmt && i >= 0 && i < d_imt) {
            d_wka[k * d_jmt * d_imt + j * d_imt + i] = 0.0e0;
        }
    }
    //tgrid_to_ugrid
    if (k == 0) {
        if (j >= 0 && j < d_ny - 1 && i >= 1 && i < d_nx) {
            d_work[j * d_imt + i] =
                    d_au0[j * d_nx + i] * d_stf[j * d_imt + i] +
                    d_aus[j * d_nx + i] * d_stf[(j + 1) * d_imt + i] +
                    d_auw[j * d_nx + i] * d_stf[j * d_imt + i - 1] +
                    d_ausw[j * d_nx + i] * d_stf[(j + 1) * d_imt + i - 1];
        }

        if (i >= 0 && i < d_nx && j == 0) {
            d_work[(d_ny - 1) * d_imt + i] = 0.0e0;

        }
        if (j >= 0 && j < d_ny && i == 0) {
            d_work[j * d_imt] = 0.0e0;
        }
    }
    //sync
}

__global__ void tracer_3(int d_imt, int d_jmt, int d_km, int n, double *d_adv_tt, double *d_tf) {
    int i, j, k;
    i = (blockIdx.x) * blockDim.x + threadIdx.x;
    j = (blockIdx.y) * blockDim.y + threadIdx.y;
    k = (blockIdx.z) * blockDim.z + threadIdx.z;

    if (k >= 0 && k < d_km) {
        if (j >= 0 && j < d_jmt && i >= 0 && i < d_imt) {
            d_tf[k * d_jmt * d_imt + j * d_imt + i] = 0.0e0;
            d_adv_tt[k * d_jmt * d_imt + j * d_imt + i] = 0.0;
        }
    }

}

__global__ void tracer_4(int d_imt, int d_jmt, int d_km, int n, double *d_tf, double *d_adv_tt, double *d_vit,
                         double *d_akt, double *d_wkc, double d_dwndmix,double *d_atb,
                         double d_ah,int *d_kmtn, int *d_kmts, int *d_kmte, int *d_kmtw,
                         double *d_dtn, double *d_dts, double *d_dte, double *d_dtw, int *d_kmt,
                         double *d_d2tk, double *d_ahf) {
    int i, j, k;
    i = (blockIdx.x) * blockDim.x + threadIdx.x;
    j = (blockIdx.y) * blockDim.y + threadIdx.y;
    k = (blockIdx.z) * blockDim.z + threadIdx.z;

    if (k >= 0 && k < d_km) {
        if (j >= 2 && j < d_jmt - 2 && i >= 2 && i < d_imt - 2) {
            d_tf[k * d_jmt * d_imt + j * d_imt + i] =
                    d_adv_tt[k * d_jmt * d_imt + j * d_imt + i] *
                    d_vit[k * d_jmt * d_imt + j * d_imt + i];
        }
    }

    if (k >= 0 && k < d_km) {
        if (j >= 1 && j < d_jmt - 1 && i >= 1 && i < d_imt - 1) {
            if (d_akt[n * d_km * d_jmt * d_imt + j * d_imt + i] < d_dwndmix) {
                d_akt[n * d_km * d_jmt * d_imt + j * d_imt + i] = d_dwndmix;
            }

            d_wkc[k * d_jmt * d_imt + j * d_imt + i] =
                    d_akt[n * d_km * d_jmt * d_imt + k * d_jmt * d_imt + j * d_imt + i] ;
                   /* +d_ahisop[j * d_imt + i] *
                    d_k3[2 * d_jmt * (d_km + 1) * d_imt + j * (d_km + 1) * d_imt + (k + 1) * d_imt + i];*/ //ISO
// hdifft_del2
                 double hdtk = 0.0e0;
                 double  cc = 0.0e0;

                if (k <= d_kmt[j * d_imt + i] - 1) {
                     double c;

                    if (k <= d_kmtn[j * d_imt + i] - 1) {
                        c = d_dtn[j * d_imt + i];
                        hdtk += c * d_atb[n * (d_km+1) * d_jmt * d_imt + (k+1) * d_jmt * d_imt+(j - 1) * d_imt + i];
                        cc -= c;
                    }

                    if (k <= d_kmts[j * d_imt + i] - 1) {
                         c = d_dts[j * d_imt + i];
                        hdtk += c * d_atb[n * (d_km+1) * d_jmt * d_imt + (k+1) * d_jmt * d_imt+(j + 1) * d_imt + i];
                         cc -= c;
                     }
 
                     if (k <= d_kmte[j * d_imt + i] - 1) {
                         c = d_dte[j * d_imt + i];
                         hdtk += c * d_atb[n * (d_km+1) * d_jmt * d_imt + (k+1) * d_jmt * d_imt+j * d_imt + i + 1];
                         cc -= c;
                     }
 
                     if (k <= d_kmtw[j * d_imt + i] - 1) {
                       c = d_dtw[j * d_imt + i];
                         hdtk += c * d_atb[n * (d_km+1) * d_jmt * d_imt + (k+1) * d_jmt * d_imt+j * d_imt + i - 1];
                         cc -= c;
                     }
                 }
 
 
                 if (j >= 2 && j < d_jmt-2 && i >=  2 && i < d_imt-2 ) {
                     hdtk = d_ah * (cc * d_atb[n * (d_km+1) * d_jmt * d_imt + (k+1) * d_jmt * d_imt+j * d_imt + i] + hdtk);
                 } else {
                     hdtk = 0;
                   }
           d_tf[k * d_jmt * d_imt + j * d_imt + i] += hdtk;
        }
    }
}

__global__ void tracer_del4_1(int d_imt, int d_jmt, int d_km, int n, double *d_tf, double *d_adv_tt, double *d_vit,
                         double *d_akt, double *d_wkc, double d_dwndmix,double *d_atb,
                         double d_ah,int *d_kmtn, int *d_kmts, int *d_kmte, int *d_kmtw,
                         double *d_dtn, double *d_dts, double *d_dte, double *d_dtw, int *d_kmt,
                         double *d_d2tk, double *d_ahf) {
    int i, j, k;
    i = (blockIdx.x) * blockDim.x + threadIdx.x;
    j = (blockIdx.y) * blockDim.y + threadIdx.y;
    k = (blockIdx.z) * blockDim.z + threadIdx.z;

    if (k >= 0 && k < d_km) {
        if (j >= 1 && j < d_jmt - 1 && i >= 1 && i < d_imt - 1) {
#ifdef  BIHAR
                 double hdtk = 0.0e0;
                 double  cc = 0.0e0;

                if (k <= d_kmt[j * d_imt + i] - 1) {
                     double c;

                    if (k <= d_kmtn[j * d_imt + i] - 1) {
                        c = d_dtn[j * d_imt + i];
                        hdtk += c * d_atb[n * (d_km+1) * d_jmt * d_imt + (k+1) * d_jmt * d_imt+(j - 1) * d_imt + i];
                        cc -= c;
                    }

                    if (k <= d_kmts[j * d_imt + i] - 1) {
                         c = d_dts[j * d_imt + i];
                        hdtk += c * d_atb[n * (d_km+1) * d_jmt * d_imt + (k+1) * d_jmt * d_imt+(j + 1) * d_imt + i];
                         cc -= c;
                     }
 
                     if (k <= d_kmte[j * d_imt + i] - 1) {
                         c = d_dte[j * d_imt + i];
                         hdtk += c * d_atb[n * (d_km+1) * d_jmt * d_imt + (k+1) * d_jmt * d_imt+j * d_imt + i + 1];
                         cc -= c;
                     }
 
                     if (k <= d_kmtw[j * d_imt + i] - 1) {
                       c = d_dtw[j * d_imt + i];
                         hdtk += c * d_atb[n * (d_km+1) * d_jmt * d_imt + (k+1) * d_jmt * d_imt+j * d_imt + i - 1];
                         cc -= c;
                     }
                 }
 
 
           //      if (j >= 2 && j < d_jmt-2 && i >=  2 && i < d_imt-2 ) {
                     hdtk = d_ahf[j * d_imt + i] * (cc * d_atb[n * (d_km+1) * d_jmt * d_imt + (k+1) * d_jmt * d_imt+j * d_imt + i] + hdtk);
             //    } else {
              //       hdtk = 0;
                //   }
           d_d2tk[k * d_jmt * d_imt + j * d_imt + i] = hdtk;
#else
#endif //BIHAR
        }
    }
}

__global__ void tracer_del4_2(int d_imt, int d_jmt, int d_km, int n, double *d_tf, double *d_adv_tt, double *d_vit,
                         double *d_akt, double *d_wkc, double d_dwndmix,double *d_atb,
                         double d_ah,int *d_kmtn, int *d_kmts, int *d_kmte, int *d_kmtw,
                         double *d_dtn, double *d_dts, double *d_dte, double *d_dtw, int *d_kmt,
                         double *d_d2tk, double *d_ahf) {
    int i, j, k;
    i = (blockIdx.x) * blockDim.x + threadIdx.x;
    j = (blockIdx.y) * blockDim.y + threadIdx.y;
    k = (blockIdx.z) * blockDim.z + threadIdx.z;

    if (k >= 0 && k < d_km) {
        if (j >= 2 && j < d_jmt - 2 && i >= 2 && i < d_imt - 2) {
            d_tf[k * d_jmt * d_imt + j * d_imt + i] =
                    d_adv_tt[k * d_jmt * d_imt + j * d_imt + i] *
                    d_vit[k * d_jmt * d_imt + j * d_imt + i];
        }
    }

    if (k >= 0 && k < d_km) {
        if (j >= 1 && j < d_jmt - 1 && i >= 1 && i < d_imt - 1) {
            if (n==0 && d_akt[n * d_km * d_jmt * d_imt + j * d_imt + i] < d_dwndmix) {
                d_akt[n * d_km * d_jmt * d_imt + j * d_imt + i] = d_dwndmix;
            }

            d_wkc[k * d_jmt * d_imt + j * d_imt + i] =
                    d_akt[n * d_km * d_jmt * d_imt + k * d_jmt * d_imt + j * d_imt + i] ;
//second phase of tdel4
                 double hdtk = 0.0e0;
                 double  cc = 0.0e0;

                if (k <= d_kmt[j * d_imt + i] - 1) {
                     double c;

                    if (k <= d_kmtn[j * d_imt + i] - 1) {
                        c = d_dtn[j * d_imt + i];
                        hdtk += c * d_d2tk[ k * d_jmt * d_imt+(j - 1) * d_imt + i];
                        cc -= c;
                    }

                    if (k <= d_kmts[j * d_imt + i] - 1) {
                         c = d_dts[j * d_imt + i];
                        hdtk += c * d_d2tk[ k * d_jmt * d_imt+(j + 1) * d_imt + i];
                         cc -= c;
                     }
 
                     if (k <= d_kmte[j * d_imt + i] - 1) {
                         c = d_dte[j * d_imt + i];
                         hdtk += c * d_d2tk[ k * d_jmt * d_imt+j * d_imt + i + 1];
                         cc -= c;
                     }
 
                     if (k <= d_kmtw[j * d_imt + i] - 1) {
                       c = d_dtw[j * d_imt + i];
                         hdtk += c * d_d2tk[ k * d_jmt * d_imt+j * d_imt + i - 1];
                         cc -= c;
                     }
                 }
 
 
                 if (j >= 2 && j < d_jmt-2 && i >=  2 && i < d_imt-2 ) {
                     hdtk = d_ah * (cc * d_d2tk[ k * d_jmt * d_imt+j * d_imt + i] + hdtk);
                 } else {
                     hdtk = 0;
                   }
           d_tf[k * d_jmt * d_imt + j * d_imt + i] += hdtk;
        }
    }
}

__global__ void tracer_del4_0(int d_imt, int d_jmt, int d_km, int n, double *d_tf, double *d_adv_tt, double *d_vit,
                         double *d_akt, double *d_wkc, double d_dwndmix,double *d_atb,
                         double d_ah,int *d_kmtn, int *d_kmts, int *d_kmte, int *d_kmtw,
                         double *d_dtn, double *d_dts, double *d_dte, double *d_dtw, int *d_kmt,
                         double *d_d2tk, double *d_ahf) {
    int i, j, k;
    i = (blockIdx.x) * blockDim.x + threadIdx.x;
    j = (blockIdx.y) * blockDim.y + threadIdx.y;
    k = (blockIdx.z) * blockDim.z + threadIdx.z;

    if (k >= 0 && k < d_km) {
        if (j >= 2 && j < d_jmt - 2 && i >= 2 && i < d_imt - 2) {
            d_tf[k * d_jmt * d_imt + j * d_imt + i] =
                    d_adv_tt[k * d_jmt * d_imt + j * d_imt + i] *
                    d_vit[k * d_jmt * d_imt + j * d_imt + i];
        }
    }

    if (k >= 0 && k < d_km) {
        if (j >= 1 && j < d_jmt - 1 && i >= 1 && i < d_imt - 1) {
            if (d_akt[n * d_km * d_jmt * d_imt + j * d_imt + i] < d_dwndmix) {
                d_akt[n * d_km * d_jmt * d_imt + j * d_imt + i] = d_dwndmix;
            }

            d_wkc[k * d_jmt * d_imt + j * d_imt + i] =
                    d_akt[n * d_km * d_jmt * d_imt + k * d_jmt * d_imt + j * d_imt + i] ;
//pre phase of tdel4
        }
    }
}

__global__ void tracer_5(int d_imt, int d_jmt, int d_km, int n, int d_nx,
                         double *d_swv, double *d_pen, double *d_vit, double *d_tf, double *d_penetrate, double *d_odzp,
                         double *d_odzt,
                         double *d_wkc, double *d_atb, double *d_dz, double *d_stf, double *d_ssf, double *d_sss,
                         double *d_ifrac, double *d_net, double *d_temp11, double *d_tarea,
                         int *d_kmt, double d_boundary_restore, int d_jst, double d_gamma, double *d_dt_restore,
                         double *d_fresh, double *d_seaice, double d_od0) {
    int i, j, k;
    double wt1, wt2;
    double aidif = 0.5e0;
    i = (blockIdx.x) * blockDim.x + threadIdx.x;
    j = (blockIdx.y) * blockDim.y + threadIdx.y;
    k = (blockIdx.z) * blockDim.z + threadIdx.z;

    if (n == 0) {
        if (k >= 1 && k < d_km - 1) {
            if (j >= 2 && j < d_jmt - 2 && i >= 2 && i < d_imt - 2) {
                wt1 = d_swv[j * d_imt + i] * d_pen[k - 1] * d_vit[k * d_jmt * d_imt + j * d_imt + i];

                wt2 = d_swv[j * d_imt + i] * d_pen[k] * d_vit[(k + 1) * d_jmt * d_imt + j * d_imt + i];
                d_tf[k * d_jmt * d_imt + j * d_imt + i] = d_tf[k * d_jmt * d_imt + j * d_imt + i] +
                                                          (wt1 - wt2) * d_odzp[k];
                d_penetrate[k * d_jmt * d_imt + j * d_imt + i] = (wt1 - wt2) * d_odzp[k];
            }
        }

        if (k == 0 && j >= 2 && j < d_jmt - 2 && i >= 2 && i < d_imt - 2) {
            wt1 = d_swv[j * d_imt + i] * d_pen[0] * d_vit[d_jmt * d_imt + j * d_imt + i];

            wt2 = d_swv[j * d_imt + i] * d_pen[d_km - 2] * d_vit[(d_km - 1) * d_jmt * d_imt + j * d_imt + i];
            d_tf[j * d_imt + i] = d_tf[j * d_imt + i] - d_odzp[0] * wt1;
            d_tf[(d_km - 1) * d_jmt * d_imt + j * d_imt + i] =
                    d_tf[(d_km - 1) * d_jmt * d_imt + j * d_imt + i] +
                    d_odzp[d_km - 1] * wt2;

            d_penetrate[j * d_imt + i] = -wt1 * d_odzp[0];

            d_penetrate[(d_km - 1) * d_jmt * d_imt + j * d_imt + i] = wt2 * d_odzp[d_km - 1];
        }
    }

    if (k >= 1 && k < d_km - 1) {
        if (j >= 2 && j < d_jmt - 2 && i >= 2 && i < d_imt - 2) {
            wt1 = d_wkc[(k - 1) * d_jmt * d_imt + j * d_imt + i] *
                  (d_atb[n * (d_km + 1) * d_jmt * d_imt + k * d_jmt * d_imt + j * d_imt + i] -
                   d_atb[n * (d_km + 1) * d_jmt * d_imt + (k + 1) * d_jmt * d_imt + j * d_imt + i]) *
                  d_odzt[k] * d_vit[k * d_jmt * d_imt + j * d_imt + i];

            wt2 = d_wkc[k * d_jmt * d_imt + j * d_imt + i] *
                  (d_atb[n * (d_km + 1) * d_jmt * d_imt + (k + 1) * d_jmt * d_imt + j * d_imt + i] -
                   d_atb[n * (d_km + 1) * d_jmt * d_imt + (k + 2) * d_jmt * d_imt + j * d_imt + i]) *
                  d_odzt[k + 1] * d_vit[(k + 1) * d_jmt * d_imt + j * d_imt + i];

            d_tf[k * d_jmt * d_imt + j * d_imt + i] = d_tf[k * d_jmt * d_imt + j * d_imt + i] +
                                                      d_odzp[k] * (wt1 - wt2) * 0.5e0;
            d_dz[n * d_km * d_jmt * d_imt + k * d_jmt * d_imt + j * d_imt + i] =
                    d_odzp[k] * (wt1 - wt2) * 0.5e0;
        }
    }

    if (k == 0 && j >= 2 && j < d_jmt - 2 && i >= 2 && i < d_imt - 2) {
        wt1 = d_wkc[j * d_imt + i] *
              (d_atb[n * (d_km + 1) * d_jmt * d_imt + d_jmt * d_imt + j * d_imt + i] -
               d_atb[n * (d_km + 1) * d_jmt * d_imt + 2 * d_jmt * d_imt + j * d_imt + i]) *
              d_odzt[1] * d_vit[d_jmt * d_imt + j * d_imt + i];

        wt2 = d_wkc[(d_km - 2) * d_jmt * d_imt + j * d_imt + i] *
              (d_atb[n * (d_km + 1) * d_jmt * d_imt + (d_km - 1) * d_jmt * d_imt + j * d_imt + i] -
               d_atb[n * (d_km + 1) * d_jmt * d_imt + d_km * d_jmt * d_imt + j * d_imt + i]) *
              d_odzt[d_km - 1] * d_vit[(d_km - 1) * d_jmt * d_imt + j * d_imt + i];

        d_tf[j * d_imt + i] = d_tf[j * d_imt + i] - d_odzp[0] * wt1 * 0.5e0;

        d_tf[(d_km - 1) * d_jmt * d_imt + j * d_imt + i] = d_tf[(d_km - 1) * d_jmt * d_imt + j * d_imt + i] +
                                                           d_odzp[d_km - 1] * wt2 * (1.0e0 - aidif);

        d_dz[n * d_km * d_jmt * d_imt + j * d_imt + i] = -d_odzp[0] * wt1 * 0.5e0;

        d_dz[n * d_km * d_jmt * d_imt + (d_km - 1) * d_jmt * d_imt + j * d_imt + i] =
                d_odzp[d_km - 1] * wt2 * 0.5e0;
    }

    if (k == 0 && n == 1) {
        if (j >= 0 && j < d_jmt && i >= 0 && i < d_imt) {
            if (d_kmt[j * d_nx + i] > 0) {
#ifdef COUP	
wpf,0409, the new version need change here

                if (d_boundary_restore == 2) {
                    d_stf[j * d_imt + i] =
                            d_ssf[j * d_imt + i] / d_odzp[0] +
                             d_gamma * 
                            (d_dt_restore[n * d_jmt * d_imt + j * d_imt + i] -
                             d_atb[(d_km + 1) * d_jmt * d_imt + d_jmt * d_imt + j * d_imt + i]) / d_odzp[0];
                } else {
                    d_stf[j * d_imt + i] = d_ssf[j * d_imt + i] / d_odzp[0];
                }
#else
	#ifdef FRC_CORE
//		work_mod_mp_stf_[iblock][j][i] = (forc_mod_mp_fresh_[j][i] * 34.7 *OD0 / 1000.0
//			+ pconst_mod_mp_gamma_ * (forc_mod_mp_sss_[iblock][j][i] - tracer_mod_mp_atb_[iblock][1][1][j][i])*
//			forc_mod_mp_seaice_[iblock][j][i] / pconst_mod_mp_odzp_[0]
//			+ pconst_mod_mp_gamma_ * 30.0 /365.0/4.0 * 50.0 *
//			(forc_mod_mp_sss_[iblock][j][i] - tracer_mod_mp_atb_[iblock][1][1][j][i]) /
//			 pconst_mod_mp_odzp_[0] * (1.0 - forc_mod_mp_seaice_[iblock][j][i]));
                    d_stf[j * d_imt + i] = d_fresh[j * d_imt + i] * 34.7 *d_od0 / 1000.0 +d_gamma *
                            (d_sss[j * d_imt + i] - d_atb[(d_km + 1) * d_jmt * d_imt + d_jmt * d_imt + j * d_imt + i]) *
                            d_seaice[j * d_imt + i]/ d_odzp[0] +
                            d_gamma * 30.0 /365.0/4.0 * 50.0 *
                            (d_sss[j * d_imt + i] - d_atb[(d_km + 1) * d_jmt * d_imt + d_jmt * d_imt + j * d_imt + i]) /
                            d_odzp[0] * (1.0 - d_seaice[j * d_imt + i]);
	#else
                    d_stf[j * d_imt + i] = d_gamma *
                            (d_sss[j * d_imt + i] -
                             d_atb[(d_km + 1) * d_jmt * d_imt + d_jmt * d_imt + j * d_imt + i]) / d_odzp[0];
	#endif
#endif
                d_tf[j * d_imt + i] = d_tf[j * d_imt + i] + d_stf[j * d_imt + i] * (1.0 - aidif) * d_odzp[0];
                d_net[d_jmt * d_imt + j * d_imt + i] = d_stf[j * d_imt + i] * d_odzp[0];
            }
        }

#ifdef SSSNORM
        if (k == 0 && j >= 0 && j < d_jmt && i >= 0 && i < d_imt) {
            d_temp11[j * d_imt + i] = d_tarea[j * d_imt + i] *
                                      d_net[d_jmt * d_imt + j * d_imt + i] *
                                      d_vit[j * d_imt + i];
        }
#endif //SSSNORM
    }
}

__global__ void tracer_6(int d_imt, int d_jmt, int d_km, int n, int d_nx,
                         double *d_tf, double *d_vit, double *d_stf, double *d_tsf, double *d_net,
                         double *d_odzp, int *d_kmt, int d_jst, double d_err_norm2,
                         double *d_swv, double *d_nswv,double *d_sst,double *d_atb,double *d_dqdt,double d_od0cp, 
                         double *d_seaice, double d_gamma) {
    int i, j;
    i = (blockIdx.x) * blockDim.x + threadIdx.x;
    j = (blockIdx.y) * blockDim.y + threadIdx.y;

    if (n == 1) {
#ifdef SSSNORM
        if (j >= 0 && j < d_jmt && i >= 0 && i < d_imt) {
            d_tf[j * d_imt + i] = d_tf[j * d_imt + i] + d_err_norm2 * d_vit[j * d_imt + i];
        }
#endif //SSSNORM
    } else {
        if (j >= 0 && j < d_jmt && i >= 0 && i < d_imt) {
            if (d_kmt[j * d_nx + i] > 0) {
#ifdef COUP
                d_stf[j * d_nx + i] = d_tsf[j * d_nx + i];
#else
#ifdef FRC_CORE
                d_stf[j * d_nx + i] = ((d_swv[j * d_imt + i] + d_nswv[j * d_imt + i])*d_od0cp +
            				d_seaice[j * d_imt + i] * d_gamma *
            				(d_sst[j * d_imt + i] - d_atb[d_jmt * d_imt+ j * d_imt + i])/
            				d_odzp[0]);
//        STF (I,J,IBLOCK) = ((SWV(I,J,IBLOCK)+NSWV(I,J,IBLOCK))*OD0CP+ & 
//        SEAICE(I,J,IBLOCK)*GAMMA*(SST(I,J,IBLOCK)-ATB(I,J,1,1,IBLOCK))/ODZP(1))
#else
                d_stf[j * d_nx + i] = (d_swv[j * d_imt + i] + d_nswv[j * d_imt + i] - 
            				d_dqdt[j * d_imt + i] * 
            				(d_sst[j * d_imt + i] - d_atb[d_jmt * d_imt+ j * d_imt + i])
            				)*d_od0cp;
//        STF (I,J,IBLOCK) = (SWV(I,J,IBLOCK)+NSWV(I,J,IBLOCK)-DQDT(I,J,IBLOCK)*
//                (SST (I,J,IBLOCK) - ATB (I,J,1,1,IBLOCK)))*OD0CP
#endif //FRC_CORE
#endif //COUP
                d_tf[j * d_nx + i] = d_tf[j * d_nx + i] + d_stf[j * d_nx + i] * d_odzp[0] * 0.5e0;
                d_net[j * d_nx + i] = d_stf[j * d_nx + i] * d_odzp[0];
            }
        }
    }
}

__global__ void tracer_7(int d_imt, int d_jmt, int d_km, int n,
                         double *d_tf, double *d_vit, double *d_atb,
                         double *d_restore, double *d_restore_at, double *d_dt_restore,
                         int d_boundary_restore, int d_simple_assm, int d_ktv) {
    int i, j, k;
    double lamda[km];
    double lamda1[km];
    i = (blockIdx.x) * blockDim.x + threadIdx.x;
    j = (blockIdx.y) * blockDim.y + threadIdx.y;
    k = (blockIdx.z) * blockDim.z + threadIdx.z;

//    if (n == 0) {
        if (k >= 0 && k < d_km) {
            lamda[k] = 1.0e0 / (15.e0 * 86400.e0);
        }

        if (k >= 0 && k < 45 && k < d_km) {
            lamda1[k] = 1.0e0 / (5.e0 * 86400.e0);
        }
        if (k >= 45 && k < d_km) {
            lamda1[k] = 1.0e0 / (5.e0 * 86400.e0) / ((k + 1) * (k + 1) / 5.0 / 5.0);
        }
//    }

    if (d_simple_assm) {
        if (k >= 0 && k < d_ktv && k < d_km) {
            if (j >= 2 && j < d_jmt - 2 && i >= 2 && i < d_imt - 2) {
                d_tf[k * d_jmt * d_imt + j * d_imt + i] =
                        d_tf[k * d_jmt * d_imt + j * d_imt + i] +
                        d_vit[k * d_jmt * d_imt + j * d_imt + i] *
                        (d_restore_at[n * d_km * d_jmt * d_imt + k * d_jmt * d_imt + j * d_imt + i] -
                         d_atb[n * (d_km + 1) * d_jmt * d_imt + (k + 1) * d_jmt * d_imt + j * d_imt + i]) *
                        lamda1[k];

                d_dt_restore[n * d_km * d_jmt * d_imt + k * d_jmt * d_imt + j * d_imt + i] =
                        d_vit[k * d_jmt * d_imt + j * d_imt + i] *
                        (d_restore_at[n * d_km * d_jmt * d_imt + k * d_jmt * d_imt + j * d_imt + i] -
                         d_atb[n * (d_km + 1) * d_jmt * d_imt + (k + 1) * d_jmt * d_imt + j * d_imt + i]) *
                        lamda1[k];
            }
        }
    }

//    lamda[k] = 1.0e0 / (15.e0 * 86400.e0);
    if (d_boundary_restore == 1) {
        if (k >= 1 && k < d_km) {
            if (j >= 2 && j < d_jmt - 2 && i >= 2 && i < d_imt - 2) {
                d_tf[k * d_jmt * d_imt + j * d_imt + i] =
                        d_tf[k * d_jmt * d_imt + j * d_imt + i] +
                        d_vit[k * d_jmt * d_imt + j * d_imt + i] *
                        (d_restore[n * d_km * d_jmt * d_imt + k * d_jmt * d_imt + j * d_imt + i] -
                         d_atb[n * (d_km + 1) * d_jmt * d_imt + (k + 1) * d_jmt * d_imt + j * d_imt + i]) *
                        lamda[k];
            }
        }
    }
}

__global__ void tracer_8(int d_imt, int d_jmt, int d_km, int n,
                         double *d_ori, double *d_vtl, double *d_at, double *d_tf, double d_dts) {
    int i, j, k;
    i = (blockIdx.x) * blockDim.x + threadIdx.x;
    j = (blockIdx.y) * blockDim.y + threadIdx.y;
    k = (blockIdx.z) * blockDim.z + threadIdx.z;

    if (k >= 0 && k < d_km) {
        if (j >= 2 && j < d_jmt - 2 && i >= 2 && i < d_imt - 2) {
            d_vtl[k * d_jmt * d_imt + j * d_imt + i] =
                    d_at[n * d_km * d_jmt * d_imt + k * d_jmt * d_imt + j * d_imt + i] +
                    d_dts * d_tf[k * d_jmt * d_imt + j * d_imt + i];

        }
    }

//wpf, 0321, sd version unuse this code

    if (k >= 0 && k < d_km) {
        if (j >= 0 && j < d_jmt && i >= 0 && i < d_imt) {
            d_ori[k * d_jmt * d_imt + j * d_imt + i] = d_vtl[k * d_jmt * d_imt + j * d_imt + i];
        }
    }

}

__global__ void tracer_9(int d_imt, int d_jmt, int d_km, int n,
                         double *d_at, double *d_vtl, double *d_atb, double *d_tend, double *d_vit,
                         double *d_dt_diff, double *d_ori, double *d_stf, double *d_odzp,
                         double d_dts) {
    int i, j, k;
    i = (blockIdx.x) * blockDim.x + threadIdx.x;
    j = (blockIdx.y) * blockDim.y + threadIdx.y;
    k = (blockIdx.z) * blockDim.z + threadIdx.z;

//wpf, 0321, sd version unuse this code

    if (j >= 0 && j < d_jmt && i >= 0 && i < d_imt) {
        d_dt_diff[n * d_km * d_jmt * d_imt + j * d_imt + i] =
                (d_vtl[j * d_imt + i] - d_ori[j * d_imt + i]) /
                d_dts * d_vit[j * d_imt + i] - d_stf[j * d_imt + i] * 0.5e0 * d_odzp[0];
    }

    if (k >= 1 && k < d_km) {
        if (j >= 0 && j < d_jmt && i >= 0 && i < d_imt) {
            d_dt_diff[n * d_km * d_jmt * d_imt + k * d_jmt * d_imt + j * d_imt + i] =
                    (d_vtl[k * d_jmt * d_imt + j * d_imt + i] -
                     d_ori[k * d_jmt * d_imt + j * d_imt + i]) / d_dts *
                    d_vit[k * d_jmt * d_imt + j * d_imt + i];
        }
    }

    if (k >= 0 && k < d_km) {
        if (j >= 0 && j < d_jmt && i >= 0 && i < d_imt) {
            d_tend[n * d_km * d_jmt * d_imt + k * d_jmt * d_imt + j * d_imt + i] =
                    (d_vtl[k * d_jmt * d_imt + j * d_imt + i]-
                     d_at[n * d_km * d_jmt * d_imt + k * d_jmt * d_imt + j * d_imt + i]) 
                     / d_dts * d_vit[k * d_jmt * d_imt + j * d_imt + i];
        }
    }

    if (k >= 0 && k < d_km) {
        if (j >= 0 && j < d_jmt && i >= 0 && i < d_imt) {
            d_at[n * d_km * d_jmt * d_imt + k * d_jmt * d_imt + j * d_imt + i] =
                    d_vtl[k * d_jmt * d_imt + j * d_imt + i];

            d_atb[n * (d_km + 1) * d_jmt * d_imt + (k + 1) * d_jmt * d_imt + j * d_imt + i] =
                    d_at[n * d_km * d_jmt * d_imt + k * d_jmt * d_imt + j * d_imt + i];
        }
    }
}

extern "C" void tracer_(int ii) {
    int i, j, n, errorcode;
    int n_tmpt,tmpt3;
    int tmpt1 = 1;
    int tmpt2 = 0;

    double err_norm1, err_norm2;
    double temp11[jmt][imt];
    double temp12[jmt_global][imt_global];
    int Str;

    if (strncmp(&pconst_mod_mp_adv_tracer_[0], "centered", 8) == 0) Str = 0;
    else if (strncmp(&pconst_mod_mp_adv_tracer_[0], "flux", 4) == 0) Str = 1;
    else if (strncmp(&pconst_mod_mp_adv_tracer_[0], "tspas", 5) == 0) Str = 2;

    //deallocate_tracer1_();
    //allocate_tracer_();//stf,tf,wkb,wkc,wkd,restore_at

//sd    if(pconst_mod_mp_simple_assm_ == false){
//sd    	hipMemset(d_restore_at, 0.0, dataSize4d);
//sd    }

    CHECK(hipMemset(d_at00, 0.0, dataSize3d));
    CHECK(hipMemset(d_atmax, 0.0, dataSize3d));
    CHECK(hipMemset(d_atmin, 0.0, dataSize3d));
    CHECK(hipMemset(d_uk, c0, dataSize3d));
    CHECK(hipMemset(d_vk, c0, dataSize3d));

    dim3 threadsPerBlock_2d(8, 8, 1);
    dim3 numBlocks_2d((imt + threadsPerBlock_2d.x - 1) / threadsPerBlock_2d.x,
                      (jmt + threadsPerBlock_2d.y - 1) / threadsPerBlock_2d.y);

    dim3 threadsPerBlock_3d(8, 8, 4);
    dim3 numBlocks_3d((imt + threadsPerBlock_3d.x - 1) / threadsPerBlock_3d.x,
                      (jmt + threadsPerBlock_3d.y - 1) / threadsPerBlock_3d.y,
                      (km + threadsPerBlock_3d.z - 1) / threadsPerBlock_3d.z);
    hipLaunchKernelGGL(tracer_1_0, dim3(numBlocks_2d), dim3(threadsPerBlock_2d ), 0, 0, d_h0f, d_h0l, d_stf, pconst_mod_mp_onbc_,imt, jmt);
    hipLaunchKernelGGL(tracer_1, dim3(numBlocks_3d), dim3(threadsPerBlock_3d ), 0, 0, nx_block, ny_block, imt, jmt, km,
            d_stf, d_utf, d_vtf, d_utl, d_vtl,d_work, d_wka, d_wkb, d_wkd, d_au0, d_aus, d_auw, d_ausw, pconst_mod_mp_oncc_);
                   iCheck();

    hipLaunchKernelGGL(upwell_2, dim3(numBlocks_3d), dim3(threadsPerBlock_3d ), 0, 0, imt, jmt, km, d_work
            , d_uk, d_vk, d_wkd, d_wkb, d_ohbu);
                    iCheck();

    hipLaunchKernelGGL(operators_div, dim3(numBlocks_3d), dim3(threadsPerBlock_3d ), 0, 0, imt, jmt, km, nx_block, ny_block,
            d_wka, d_kmt, d_uk, d_vk, d_htw, d_hts, d_tarea_r);
                    iCheck();

    hipLaunchKernelGGL(upwell_3, dim3(numBlocks_2d), dim3(threadsPerBlock_2d ), 0, 0, imt, jmt, km, d_work, d_wka,
            d_vit, d_ws, d_ohbt, d_stf, d_dzp);
                    iCheck();

    hipLaunchKernelGGL(upwell_4, dim3(numBlocks_3d), dim3(threadsPerBlock_3d ), 0, 0, imt, jmt, km, d_work, d_wka,
            d_ws, d_ohbt, d_stf, d_dzp);
                    iCheck();

#if (defined ISO)
    CHECK(hipMemcpy(pmix_mod_mp_rict_replace_, d_rict_replace, dataSize3333d, hipMemcpyDeviceToHost));
    if(ii!=0){
    	CHECK(hipMemcpy(tracer_mod_mp_atb_, d_atb, dataSize44d, hipMemcpyDeviceToHost));
    }
    isopyc_();
#endif

    pconst_mod_mp_ist_ += 1;

    //CHECK(hipMemcpy(d_ahisop, isopyc_mod_mp_ahisop_, dataSize2d, hipMemcpyHostToDevice));
    //CHECK(hipMemcpy(d_k3, isopyc_mod_mp_k3_[0], dataSize444d, hipMemcpyHostToDevice));

                    iCheck();
    for (n = 0; n < ntra; n++) { //ntra = 2
        hipLaunchKernelGGL(tracer_3, dim3(numBlocks_3d), dim3(threadsPerBlock_3d ), 0, 0, imt, jmt, km, n, d_adv_tt, d_tf);
                    iCheck();

        n_tmpt = n + 1;
        hipLaunchKernelGGL(advection_tracer_1_0, dim3(numBlocks_3d), dim3(threadsPerBlock_3d ), 0, 0,
                                   d_wkb, d_wkd, d_u_wface, d_v_sface, d_hts, d_htw, d_dxu, d_dyu,
                                   nx_block, km, jmt, imt, Str) ;
        hipLaunchKernelGGL(advection_tracer_1, dim3(numBlocks_3d), dim3(threadsPerBlock_3d ), 0, 0, d_u_wface, d_v_sface, d_ws, d_at + n * km * imt * jmt, d_adv_tt, d_tarea_r, d_odzp, d_ax, d_ay, d_az, nx_block, km, jmt, imt, n, pconst_mod_mp_nss_, Str);
 
                   iCheck();

        if (Str == 2) {
            hipLaunchKernelGGL(advection_tracer_2, dim3(numBlocks_3d), dim3(threadsPerBlock_3d ), 0, 0, d_wkb, d_wkd, d_u_wface, d_v_sface, d_hts,
                    d_htw, d_dxu, d_dyu, d_ws, d_at + n * km * imt * jmt, d_adv_tt, d_tarea_r, d_odzp,
                    d_ax, d_ay, d_az, d_hun, d_hue, d_at00, d_atmax, d_atmin, d_odzt, d_vit, pconst_mod_mp_dts_,
                    ny_block, nx_block, km, jmt, imt, ntra, n, pconst_mod_mp_nss_, Str);
                    iCheck();

            hipLaunchKernelGGL(advection_tracer_3, dim3(numBlocks_3d), dim3(threadsPerBlock_3d ), 0, 0, d_wkb, d_wkd, d_u_wface, d_v_sface, d_hts,
                    d_htw, d_dxu, d_dyu, d_ws, d_at + n * km * imt * jmt, d_adv_tt, d_tarea_r, d_odzp,
                    d_ax, d_ay, d_az, d_hun, d_hue, d_at00, d_atmax, d_atmin, d_odzt, d_vit, pconst_mod_mp_dts_,
                    ny_block, nx_block, km, jmt, imt, ntra, n, pconst_mod_mp_nss_, Str);
                    iCheck();
        }

#ifdef BIHAR

        CHECK(hipMemset(d_d2tk, c0, dataSize3d));
        hipLaunchKernelGGL(tracer_del4_1, dim3(numBlocks_3d), dim3(threadsPerBlock_3d ), 0, 0, imt, jmt, km, n, d_tf, d_adv_tt, d_vit,
                d_akt, d_wkc, pconst_mod_mp_dwndmix_,d_atb,hmix_del4_mp_ah_, 
                d_kmtn, d_kmts, d_kmte, d_kmtw, d_dtn, d_dts, d_dte, d_dtw, d_kmt, d_d2tk, d_ahf);
        hipLaunchKernelGGL(tracer_del4_2, dim3(numBlocks_3d), dim3(threadsPerBlock_3d ), 0, 0, imt, jmt, km, n, d_tf, d_adv_tt, d_vit,
                d_akt, d_wkc, pconst_mod_mp_dwndmix_,d_atb,hmix_del4_mp_ah_, 
                d_kmtn, d_kmts, d_kmte, d_kmtw, d_dtn, d_dts, d_dte, d_dtw, d_kmt, d_d2tk, d_ahf);
/*
//use F90 to do difft_del4, when HIP code ready,change here
        hipLaunchKernelGGL(tracer_del4_0, dim3(numBlocks_3d), dim3(threadsPerBlock_3d ), 0, 0, imt, jmt, km, n, d_tf, d_adv_tt, d_vit,
                d_akt, d_wkc, pconst_mod_mp_dwndmix_,d_atb,hmix_del4_mp_ah_, 
                d_kmtn, d_kmts, d_kmte, d_kmtw, d_dtn, d_dts, d_dte, d_dtw, d_kmt, d_d2tk, d_ahf);
//        CHECK(hipMemcpy(tracer_mod_mp_atb_, d_atb, dataSize44d, hipMemcpyDeviceToHost));
        CHECK(hipMemcpy(tracer_mod_mp_atb_, &d_atb[n * (km + 1) * jmt * imt] , dataSize44d/2, hipMemcpyDeviceToHost));
        CHECK(hipMemcpy(work_mod_mp_tf_, d_tf, dataSize3d, hipMemcpyDeviceToHost));
        tracer_tdel_(&n_tmpt);
        CHECK(hipMemcpy(d_tf, work_mod_mp_tf_, dataSize3d, hipMemcpyHostToDevice));
//get tf back
*/
#else
        hipLaunchKernelGGL(tracer_4, dim3(numBlocks_3d), dim3(threadsPerBlock_3d ), 0, 0, imt, jmt, km, n, d_tf, d_adv_tt, d_vit,
                d_akt, d_wkc, pconst_mod_mp_dwndmix_,d_atb,hmix_del2_mp_ah_, 
                d_kmtn, d_kmts, d_kmte, d_kmtw, d_dtn, d_dts, d_dte, d_dtw, d_kmt, d_d2tk, d_ahf);
#endif //BIHAR

#if (defined ISO)
        CHECK(hipMemcpy(work_mod_mp_tf_, d_tf, dataSize3d, hipMemcpyDeviceToHost));
        //work_mod_mp_tf_ will be used in function isoflus; if cancel this word, it will cause big error.
       if(n!=0){
	        hipMemcpy(tracer_mod_mp_atb_, d_atb, dataSize44d, hipMemcpyDeviceToHost);
	}

        isoflux_(&n_tmpt);

        CHECK(hipMemcpy(d_tf, work_mod_mp_tf_, dataSize3d, hipMemcpyHostToDevice));
        //work_mod_mp_tf_ will be used in function isoflus; if cancel this word, it will cause big error.
#endif

        CHECK(hipMemcpy(d_net, tracer_mod_mp_net_, dataSize33d, hipMemcpyHostToDevice));
        //function pop_haloupdate_tracer1 need net; cancel this word will cause error.

//        CHECK(hipMemcpy(d_temp11, temp11, dataSize2d, hipMemcpyHostToDevice));

        hipLaunchKernelGGL(tracer_5, dim3(numBlocks_3d), dim3(threadsPerBlock_3d ), 0, 0, imt, jmt, km, n, nx_block,
                d_swv, d_pen, d_vit, d_tf, d_penetrate, d_odzp, d_odzt,
                d_wkc, d_atb, d_dz, d_stf, d_ssf, d_sss, d_ifrac, d_net, d_temp11, d_tarea,
                d_kmt, pconst_mod_mp_boundary_restore_, jst, pconst_mod_mp_gamma_, d_dt_restore,
                d_fresh, d_seaice, pconst_mod_mp_od0_);
                    iCheck();
#ifdef SSSNORM
        CHECK(hipMemcpy(temp11, d_temp11, dataSize2d, hipMemcpyDeviceToHost));

             err_norm2=0.0e0;
        if (n == 1) {
                err_norm1 = 0.0e0;
                for (j = 2; j < jmt-2; j++) {
                    for (i = 2; i < imt-2; i++) {
                        err_norm1 = err_norm1 + temp11[j][i];
                    }
                }
            mpi3_(&err_norm1);
//            if (param_mod_mp_mytid_ == 0) printf("err_norm1= %le\n",err_norm1);

            err_norm2 = -err_norm1 / grid_mp_area_t_;

            tracer_mod_mp_fw_norm2_ = err_norm2;

        }
#endif //SSSNORM
        hipLaunchKernelGGL(tracer_6, dim3(numBlocks_2d), dim3(threadsPerBlock_2d ), 0, 0, imt, jmt, km, n, nx_block, d_tf,
                d_vit, d_stf, d_tsf, d_net, d_odzp, d_kmt, jst, err_norm2,
                d_swv, d_nswv, d_sst, d_atb, d_dqdt, pconst_mod_mp_od0cp_, d_seaice, pconst_mod_mp_gamma_);
                    iCheck();

        CHECK(hipMemcpy(tracer_mod_mp_net_, d_net, dataSize33d, hipMemcpyDeviceToHost));
        pop_haloupdate_tracer1_(&errorcode);//d_net

        hipLaunchKernelGGL(tracer_7, dim3(numBlocks_3d), dim3(threadsPerBlock_3d ), 0, 0, imt, jmt, km, n, d_tf, d_vit, d_atb,
                d_restore, d_restore_at, d_dt_restore,
                pconst_mod_mp_boundary_restore_, tmpt2, tmpt2);
                    iCheck();


        hipLaunchKernelGGL(tracer_8, dim3(numBlocks_3d), dim3(threadsPerBlock_3d ), 0, 0, imt, jmt, km, n, d_ori, d_vtl,
                d_at, d_tf, pconst_mod_mp_dts_);
                    iCheck();

        hipLaunchKernelGGL(invtrit, dim3(numBlocks_2d), dim3(threadsPerBlock_2d ), 0, 0, nx_block, ny_block, imt, jmt, km,
                d_kmt, d_wkc, d_vtl, d_stf, d_vit, d_odzt, d_odzp, pconst_mod_mp_dts_);
                    iCheck();

        CHECK(hipMemcpy(dyn_mod_mp_vtl_, d_vtl, dataSize3d, hipMemcpyDeviceToHost));
        pop_haloupdate_tracer2_(&errorcode);//d_vtl
        CHECK(hipMemcpy(d_vtl, dyn_mod_mp_vtl_, dataSize3d, hipMemcpyHostToDevice));

        hipLaunchKernelGGL(tracer_9, dim3(numBlocks_3d), dim3(threadsPerBlock_3d ), 0, 0, imt, jmt, km, n, d_at, d_vtl,
                d_atb, d_tend, d_vit, d_dt_diff, d_ori, d_stf, d_odzp, pconst_mod_mp_dts_);
                    iCheck();
    }

    //deallocate_tracer2_();

    return;
}
