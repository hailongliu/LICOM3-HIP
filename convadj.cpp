#include "hip/hip_runtime.h"
#include "cuda_data.h"

#include "precision_mod.h"
#include "param_mod.h"
#include "pconst_mod.h"
#include "constant_mod.h"
#include "tracer_mod.h"
#include "domain.h"
#include "grid.h"

#include "string.h"
#include "stdlib.h"
#include "stdio.h"

#if (defined LOWRES)
#include "output_mod.h"
#endif

//#include"strcmp.cu"
//#include"dens_cu.cu"

//extern __device__ int strcmp(char, const char);
//
//extern __device__ double dens_c_cu(double *, double *, int *, double *);
//
//extern __device__ double dens_f_cu(double *, double *, int *, double *);

__global__ void convadj_cu(double *d_to, double *d_so, double d_pconst_dts,int *d_kmt, double *d_vit, double *d_dzp,
                           double *d_at, double *d_tend, double *d_dt_conv, double *d_atb,
                           double *d_icmon, double *d_c) {

    double rhoup[km] = {0.0e0};
    double rholo[km] = {0.0e0};

    double trasum[2];
    double tup, sup, tlo, slo, dztsum, tramix;
    int kcon, lctot, lcven, l1, l, lcon, lcona, lconb, lmix, n;

    double c2dtts;
    double tmpt1, tmpt2;

    int i, j, k;
    i = (blockIdx.x) * blockDim.x + threadIdx.x;
    j = (blockIdx.y) * blockDim.y + threadIdx.y;

    if (i < imt && j < jmt) {
        kcon = d_kmt[j * imt + i];
        lctot = 0;
        lcven = 0;

        //if (kcon == 0) goto III;
        if (kcon != 0) {
            lcven = 1;
            lcon = 0;

            for (l = 0; l < km - 1; l++) {
                l1 = l + 1;
                tup = d_at[l1 * jmt * imt + j * imt + i] - d_to[l1];
                sup = d_at[1 * km * jmt * imt + l1 * jmt * imt + j * imt + i] - d_so[l1];
                tlo = d_at[l * jmt * imt + j * imt + i] - d_to[l1];
                slo = d_at[1 * km * jmt * imt + l * jmt * imt + j * imt + i] - d_so[l1];
//                rhoup[l1] = dens_c_cu(&tup, &sup, &l1, d_c);//????
                rhoup[l1] = (d_c[l1] +
                             (d_c[3 * km + l1] + d_c[6 * km + l1] * sup) * sup +
                             (d_c[2 * km + l1] + d_c[7 * km + l1] * sup + d_c[5 * km + l1] * tup) *
                             tup) * tup
                            + (d_c[1 * km + l1] + (d_c[4 * km + l1] + d_c[8 * km + l1] * sup) * sup) *
                              sup;
//                rholo[l] = dens_c_cu(&tlo, &slo, &l1, d_c);//????
                rholo[l] = (d_c[l1] +
                            (d_c[3 * km + l1] + d_c[6 * km + l1] * slo) * slo +
                            (d_c[2 * km + l1] + d_c[7 * km + l1] * slo + d_c[5 * km + l1] * tlo) *
                            tlo) * tlo
                           + (d_c[1 * km + l1] + (d_c[4 * km + l1] + d_c[8 * km + l1] * slo) * slo) *
                             slo;
            }

            for (k = kcon - 1; k >= 1; k--) {
                if (rholo[k - 1] > rhoup[k]) lcon = k;
            }

            if (lcon != 0) {

                //  CONV_1 : DO
                for (;;) {
                    lcona = lcon;
                    lconb = lcon + 1;

                    dztsum = d_dzp[lcona - 1] + d_dzp[lconb - 1];

                    for (n = 0; n < 2; n++) {
                        trasum[n] = d_at[n * km * jmt * imt + (lcona - 1) * jmt * imt + j * imt + i] *
                                    d_dzp[lcona - 1] +
                                    d_at[n * km * jmt * imt + (lconb - 1) * jmt * imt + j * imt + i] *
                                    d_dzp[lconb - 1];
                        tramix = trasum[n] / dztsum;
                        d_at[n * km * jmt * imt + (lcona - 1) * jmt * imt + j * imt + i] = tramix;
                        d_at[n * km * jmt * imt + (lconb - 1) * jmt * imt + j * imt + i] = tramix;
                    }

                    //  CONV_2 : DO
                    for (;;) {

                        if (lconb != kcon) {
                            l1 = lconb + 1;
                            tmpt1 = d_at[(lconb - 1) * jmt * imt + j * imt + i] - d_to[l1 - 1];
                            tmpt2 = d_at[1 * km * jmt * imt + (lconb - 1) * jmt * imt + j * imt + i] -
                                    d_so[l1 - 1];
//                            rholo[lconb - 1] = dens_f_cu(&tmpt1, &tmpt2, &l1, d_c);//????
                            rholo[lconb - 1] = (d_c[l1 - 1] +
                                                (d_c[3 * km + l1 - 1] + d_c[6 * km + l1 - 1] * tmpt2) * tmpt2 +
                                                (d_c[2 * km + l1 - 1] + d_c[7 * km + l1 - 1] * tmpt2 +
                                                 d_c[5 * km + l1 - 1] * tmpt1) * tmpt1) * tmpt1
                                               + (d_c[1 * km + l1 - 1] +
                                                  (d_c[4 * km + l1 - 1] + d_c[8 * km + l1 - 1] * tmpt2) * tmpt2) *
                                                 tmpt2;

                            if (rholo[lconb - 1] > rhoup[l1 - 1]) {
                                lconb = lconb + 1;
                                dztsum = dztsum + d_dzp[lconb - 1];

                                for (n = 0; n < 2; n++) {
                                    trasum[n] = trasum[n] +
                                                d_at[n * km * jmt * imt + (lconb - 1) * jmt * imt +
                                                     j * imt + i] * d_dzp[lconb - 1];
                                    tramix = trasum[n] / dztsum;

                                    for (lmix = lcona - 1; lmix < lconb; lmix++) {
                                        d_at[n * km * jmt * imt + lmix * jmt * imt + j * imt + i] = tramix;
                                    }
                                }

                                //  CYCLE CONV_2
                                goto CONV_2;
                                //continue;
                            }
                        }

                        if (lcona > 1) {
                            l1 = lcona - 1;
                            tmpt1 = d_at[(l1 - 1) * jmt * imt + j * imt + i] - d_to[lcona - 1];
                            tmpt2 = d_at[1 * km * jmt * imt + (l1 - 1) * jmt * imt + j * imt + i] -
                                    d_so[lcona - 1];
//                            rholo[l1 - 1] = dens_f_cu(&tmpt1, &tmpt2, &lcona, d_c);//????
                            rholo[l1 - 1] = (d_c[lcona - 1] +
                                             (d_c[3 * km + lcona - 1] + d_c[6 * km + lcona - 1] * tmpt2) * tmpt2 +
                                             (d_c[2 * km + lcona - 1] + d_c[7 * km + lcona - 1] * tmpt2 +
                                              d_c[5 * km + lcona - 1] * tmpt1) * tmpt1) * tmpt1
                                            + (d_c[1 * km + lcona - 1] +
                                               (d_c[4 * km + lcona - 1] + d_c[8 * km + lcona - 1] * tmpt2) * tmpt2) *
                                              tmpt2;

                            tmpt1 = d_at[(lcona - 1) * jmt * imt + j * imt + i] - d_to[lcona - 1];
                            tmpt2 = d_at[1 * km * jmt * imt + (lcona - 1) * jmt * imt + j * imt + i] -
                                    d_so[lcona - 1];
//                            rhoup[lcona - 1] = dens_f_cu(&tmpt1, &tmpt2, &lcona, d_c);//????
                            rhoup[lcona - 1] = (d_c[lcona - 1] +
                                             (d_c[3 * km + lcona - 1] + d_c[6 * km + lcona - 1] * tmpt2) * tmpt2 +
                                             (d_c[2 * km + lcona - 1] + d_c[7 * km + lcona - 1] * tmpt2 +
                                              d_c[5 * km + lcona - 1] * tmpt1) * tmpt1) * tmpt1
                                            + (d_c[1 * km + lcona - 1] +
                                               (d_c[4 * km + lcona - 1] + d_c[8 * km + lcona - 1] * tmpt2) * tmpt2) *
                                              tmpt2;

                            if (rholo[lcona - 2] > rhoup[lcona - 1]) {
                                lcona = lcona - 1;
                                dztsum = dztsum + d_dzp[lcona - 1];
                                for (n = 0; n < 2; n++) {
                                    trasum[n] = trasum[n] +
                                                d_at[n * km * jmt * imt + (lcona - 1) * jmt * imt +
                                                     j * imt + i] * d_dzp[lcona - 1];
                                    tramix = trasum[n] / dztsum;
                                    for (lmix = lcona; lmix <= lconb; lmix++) {
                                        d_at[n * km * jmt * imt + (lmix - 1) * jmt * imt + j * imt +
                                             i] = tramix;
                                    }
                                }
                                // CYCLE CONV_2
                                goto CONV_2;
                                //continue;
                            }
                        }

                        // EXIT CONV_2
                        goto EXIT_CONV_2;
                        //break;

                        CONV_2:;
                    }

                    EXIT_CONV_2:;

                    lctot = lctot + lconb - lcona + 1;

                    if (lcona == 1) lcven = lconb - lcona + 1;

                    if (lconb == kcon) {
#if (defined LOWRES)
                        d_icmon[j * imt + i] += lctot;
                        d_icmon[1 * jmt * imt + j * imt + i] += lcven;
#endif
                        // CYCLE III
                        goto III;
                        //break;
                    }

                    lcon = lconb;

                    //  CONV_3 : DO
                    for (;;) {
                        lcon = lcon + 1;
                        if (lcon == kcon) {
#if (defined LOWRES)
                            d_icmon[j * imt + i] += lctot;
                            d_icmon[j * imt + i] += lctot;
                            d_icmon[1 * jmt * imt + j * imt + i] += lcven;
#endif
                            // CYCLE III
                            goto III; //?????
                            //break;
                        }

                        if (rholo[lcon - 1] <= rhoup[lcon])
                            //CYCLE CONV_3
                            goto CONV_3;
                            //continue;

                        // EXIT CONV_3
                        goto EXIT_CONV_3;
                        //break;

                        CONV_3:;
                    }// END DO CONV_3

                    //if (lcon == kcon) break;

                    EXIT_CONV_3:;
                }// END DO CONV_1
III:;
            }
        }

//wpf, 0409, the follow code is reuse in new 5k version

//        if (strcmp(d_adv_tracer, "centered") == 0) {//no
//            if (*d_ist >= 1) {
//                c2dtts = *d_pconst_dts * 2.0e0;
//            } else {
//                c2dtts = *d_pconst_dts;
//            }
//        } else if (strcmp(d_adv_tracer, "tspas") == 0) {//yes
        c2dtts = d_pconst_dts;
//        }
        //else {
        // if(mytid==0) write(16,*)'error in convadj'
        //output_convadj_();
        // call exit_licom(sigAbort,'The false advection option for tracer in convadj')
        //}

        for (n = 0; n < ntra; n++) {
            for (k = 0; k < km; k++) {
                d_dt_conv[n * km * jmt * imt + k * jmt * imt + j * imt + i] =
                        (d_at[n * km * jmt * imt + k * jmt * imt + j * imt + i] -
                         d_atb[n * (km + 1) * jmt * imt + (k + 1) * jmt * imt + j * imt + i]) / c2dtts *
                        d_vit[k * jmt * imt + j * imt + i]; //for output dt diffusion //LPF20160823
                d_tend[n * km * jmt * imt + k * jmt * imt + j * imt + i] += d_dt_conv[
                        n * km * jmt * imt + k * jmt * imt + j * imt + i];//total tendency

//                if (strcmp(d_adv_tracer, "tspas") == 0) {//yes
                d_atb[n * (km + 1) * jmt * imt + (k + 1) * jmt * imt + j * imt + i] = 
            	    d_at[n * km * jmt * imt + k * jmt * imt + j * imt + i];
//                }
            }
        }
    }
/*
//  if (strcmp(d_adv_tracer, "tspas") == 0) {//yes
        for (n = 0; n < ntra; n++) {
            for (k = 0; k < km; k++) {
                d_atb[n * (km + 1) * jmt * imt + (k + 1) * jmt * imt + j * imt + i] = 
            	    d_at[n * km * jmt * imt + k * jmt * imt + j * imt + i];
            }
        }
   }
*/
    return;
}

extern "C" void convadj_() {

    dim3 blockSize(8, 16, 1);
    dim3 gridSize((imt + blockSize.x - 1) / blockSize.x, (jmt + blockSize.y - 1) / blockSize.y, 1);

    hipLaunchKernelGGL(convadj_cu, dim3(gridSize), dim3(blockSize ), 0, 0, d_to, d_so, pconst_mod_mp_dts_, d_kmt, d_vit, d_dzp, d_at, d_tend, d_dt_conv, d_atb, d_icmon, d_c);

    return;
}
