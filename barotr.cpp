/*by zly 2018.11*/
// Fix some bugs by Wpf, 2019.1.22
#include "hip/hip_runtime.h"
#include "cuda_data.h"

#include "param_mod.h"
#include "work_mod.h"
#include "pconst_mod.h"
#include "dyn_mod.h"
#include "grid.h"
#include "precision_mod.h"
#include "hmix_del2.h"
#include "operators.h"
#include "blocks.h"
#include "constant_mod.h"
#include "domain.h"
#include "POP_HaloMod.h"
#include "POP_GridHorzMod.h"
#include "hmix_del4.h"
#include<stdio.h>
#include "common.h"

extern "C" void pop_haloupdate_barotr1_(int *);

extern "C" void pop_haloupdate_barotr2_(int *, int *);
//extern "C" void pop_haloupdate2d8rmy_(int *, double *);

extern "C" void deallocate_barotr_();
//extern "C" void barotr_f_();
//extern "C" void barotr_tdel_();

__global__ void barotr1_cu(double *d_work, double *d_wka, double *d_dzph,double *d_h0,
                           double *d_ub, double *d_vb, double *d_au0, double *d_ausw, double *d_aus,
                           double *d_auw, int d_imt, int d_jmt) {
    int i, j;

    i = (blockIdx.x) * blockDim.x + threadIdx.x;
    j = (blockIdx.y) * blockDim.y + threadIdx.y;

    if (j < d_jmt && i < d_imt) {
        if (j == d_jmt - 1 || i == 0) {
            d_work[j * d_imt + i] = 0.0e0;
        } else {
            d_work[j * d_imt + i] =
                    d_au0[j * d_imt + i] * d_h0[j * d_imt + i] +
                    d_aus[j * d_imt + i] * d_h0[(j + 1) * d_imt + i] +
                    d_auw[j * d_imt + i] * d_h0[j * d_imt + i - 1] +
                    d_ausw[j * d_imt + i] * d_h0[(j + 1) * d_imt + i - 1];
        }

        if (i < d_imt - 1 && j > 0) {
            d_wka[j * d_imt + i] = d_ub[j * d_imt + i] * (d_dzph[j * d_imt + i] + d_work[j * d_imt + i]);
            d_wka[1 * d_jmt * d_imt + j * d_imt + i] =
                    d_vb[j * d_imt + i] * (d_dzph[j * d_imt + i] + d_work[j * d_imt + i]);
        }
    }
     return;
}

__global__ void barotr1_cu1(double *d_work, double *d_wka, double *d_vit, double *d_htw, double *d_hts, int *d_kmt,
                           double *d_tarea_r, int d_imt, int d_jmt) {
    int i, j;

    i = (blockIdx.x) * blockDim.x + threadIdx.x;
    j = (blockIdx.y) * blockDim.y + threadIdx.y;


        if (i < d_imt - 2 && i > 1 && j > 1 && j < d_jmt - 2) {
            if (0 <= d_kmt[j * d_imt + i] - 1) {
                d_work[j * d_imt + i] =
                        0.5 * ((d_wka[j * d_imt + i + 1] + d_wka[(j - 1) * d_imt + i + 1]) * d_htw[j * d_imt + i + 1] -
                               (d_wka[j * d_imt + i] + d_wka[(j - 1) * d_imt + i]) * d_htw[j * d_imt + i] +
                               (d_wka[1 * d_imt * d_jmt + j * d_imt + i + 1] +
                                d_wka[1 * d_imt * d_jmt + j * d_imt + i]) *
                               d_hts[j * d_imt + i] -
                               (d_wka[1 * d_imt * d_jmt + (j - 1) * d_imt + i + 1] +
                                d_wka[1 * d_imt * d_jmt + (j - 1) * d_imt + i]) * d_hts[(j - 1) * d_imt + i]) *
                        d_tarea_r[j * d_imt + i] * d_vit[j * d_imt + i] * -0.250;
            } else {
                d_work[j * d_imt + i] = 0;
            }

        }
         return;
    }
__global__ void barotr2_cu1(double d_dtb, double *d_h0, double *d_h0p, double *d_work, int d_imt, int d_jmt) {
    int i, j;

    i = (blockIdx.x) * blockDim.x + threadIdx.x;
    j = (blockIdx.y) * blockDim.y + threadIdx.y;

    if (i < d_imt && j < d_jmt) {
       d_h0[j * d_imt + i] = d_h0p[j * d_imt + i]+ d_work[j * d_imt + i] * d_dtb;
   }
    return;
}


//        d_h0p[j * d_imt + i] = d_h0[j * d_imt + i];
//
//

#ifndef BIHAR
__global__ void barotr2_cu(double d_dtb, double *d_viv, double *d_ebea, double *d_ebeb, double *d_fcor,
                           double *d_h0, double *d_h0p, double *d_ubp, double *d_vbp,
                           double *d_dlub, double *d_dlvb, double *d_wka, double *d_work, double *d_wgp,
                           double *d_pax, double *d_pay, double *d_pxb, double *d_pyb, double *d_whx, double *d_why,
                           double d_am, double *d_duc, double *d_dum, double *d_dun,
                           double *d_dus, double *d_due, double *d_duw, double *d_dmc,
                           double *d_dmn, double *d_dms, double *d_dme, double *d_dmw,
                           double *d_dxur, double *d_dyur, int *d_kmu,
                           double *d_au0, double *d_ausw, double *d_aus, double *d_auw,
                           int d_jb, int d_je, int d_ib, int d_ie, int d_imt, int d_jmt, int d_km) {
    int i, j;

    i = (blockIdx.x) * blockDim.x + threadIdx.x;
    j = (blockIdx.y) * blockDim.y + threadIdx.y;

//---------------------------------------------------------------------
//     COMPUTE THE "ARTIFICIAL" HORIZONTAL VISCOSITY
//---------------------------------------------------------------------
/*
		tmpt1 = 1;
		for (iblock = 0; iblock < domain_mp_nblocks_clinic_; iblock++) {
			// get_thisblock_f(&this_block, &iblock);
//			iblock = iblock+1;  //if function used below is written in c ,decomment this line
//			hmix_del4_mp_hdiffu_del4_(&tmpt1, hduk, hdvk, dyn_mod_mp_ubp_[iblock], dyn_mod_mp_vbp_[iblock], &jb, &je, &ib, &ie, &iblock);
			this_block = blocks_mp_all_blocks_[domain_mp_blocks_clinic_[iblock]-1];
			tmpt1 = 1;
			this_block.local_id=iblock+1;
			hmix_del4_mp_hdiffu_del4_(&tmpt1, hduk, hdvk, dyn_mod_mp_ubp_[iblock], dyn_mod_mp_vbp_[iblock], &this_block);
			for (j = 2; j < jmt - 2; j++) {
				for (i = 2; i < imt - 2; i++) {
					work_mod_mp_wka_[iblock][4][j][i] = hduk[j][i];
					work_mod_mp_wka_[iblock][5][j][i] = hdvk[j][i];
				}
			}
		}
*/
    if (i < d_imt && j < d_jmt) {
        double hduk = 0, hdvk = 0;
//            d_wka[j * d_imt + i] = 0; //jjrtest


//        d_h0p[j * d_imt + i] = d_h0[j * d_imt + i];

//        if (i > 1 && i < d_imt - 2 && j > 1 && j < d_jmt - 2) {
            if (d_kmu[j * d_imt + i] - 1 >= 0) {
                if (j >= d_jb - 1 && j < d_je && i >= d_ib - 1 && i < d_ie) {
                    hduk = d_am *
                           (((d_duc[j * d_imt + i] + d_dum[j * d_imt + i]) * d_ubp[j * d_imt + i] +
                             d_dun[j * d_imt + i] * d_ubp[(j - 1) * d_imt + i] +
                             d_dus[j * d_imt + i] * d_ubp[(j + 1) * d_imt + i] +
                             d_due[j * d_imt + i] * d_ubp[j * d_imt + i + 1] +
                             d_duw[j * d_imt + i] * d_ubp[j * d_imt + i - 1]) +
                            (d_dmc[j * d_imt + i] * d_vbp[j * d_imt + i] +
                             d_dmn[j * d_imt + i] * d_vbp[(j - 1) * d_imt + i] +
                             d_dms[j * d_imt + i] * d_vbp[(j + 1) * d_imt + i] +
                             d_dme[j * d_imt + i] * d_vbp[j * d_imt + i + 1] +
                             d_dmw[j * d_imt + i] * d_vbp[j * d_imt + i - 1])) * d_viv[j * d_imt + i];

                    hdvk = d_am *
                           (((d_duc[j * d_imt + i] + d_dum[j * d_imt + i]) * d_vbp[j * d_imt + i] +
                             d_dun[j * d_imt + i] * d_vbp[(j - 1) * d_imt + i] +
                             d_dus[j * d_imt + i] * d_vbp[(j + 1) * d_imt + i] +
                             d_due[j * d_imt + i] * d_vbp[j * d_imt + i + 1] +
                             d_duw[j * d_imt + i] * d_vbp[j * d_imt + i - 1]) -
                            (d_dmc[j * d_imt + i] * d_ubp[j * d_imt + i] +
                             d_dmn[j * d_imt + i] * d_ubp[(j - 1) * d_imt + i] +
                             d_dms[j * d_imt + i] * d_ubp[(j + 1) * d_imt + i] +
                             d_dme[j * d_imt + i] * d_ubp[j * d_imt + i + 1] +
                             d_dmw[j * d_imt + i] * d_ubp[j * d_imt + i - 1])) * d_viv[j * d_imt + i]; //zf change d_ubp, d_vbp
                }

            }

  //      }
           

        double gradx = 0.0e0, grady = 0.0e0;
        if (i > 0 && i < d_imt - 1 && j > 1 && j < d_jmt ) { //jjr bug
            if (d_kmu[j * d_imt + i] - 1 >= 0) {
                gradx = d_dxur[j * d_imt + i] * 0.5 *
                        (d_h0[(j + 1) * d_imt + i] -
                         d_h0[j * d_imt + i - 1] -
                         d_h0[(j + 1) * d_imt + i - 1] +
                         d_h0[j * d_imt + i]);

                grady = d_dyur[j * d_imt + i] * 0.5 *
                        (d_h0[(j + 1) * d_imt + i] -
                         d_h0[j * d_imt + i - 1] +
                         d_h0[(j + 1) * d_imt + i - 1] -
                          d_h0[j * d_imt + i]);
                        }
            double gstar = (d_wgp[j * d_imt + i] - 1.0) * g;
            d_wka[j * d_imt + i] = hduk + gstar * gradx;
            d_wka[1 * d_jmt * d_imt + j * d_imt + i] = hdvk + gstar * grady;

            d_wka[4 * d_jmt * d_imt + j * d_imt + i] = hduk;
            d_wka[5 * d_jmt * d_imt + j * d_imt + i] = hdvk;

                }

//        d_h0[j * d_imt + i] = d_h0p[j * d_imt + i]+ d_work[j * d_imt + i] * d_dtb;
        if (j == d_jmt - 1 || i == 0) {
            d_work[j * d_imt + i] = 0.0e0;
        } else {
            d_work[j * d_imt + i] =
                    d_au0[j * d_imt + i] * d_h0[j * d_imt + i] +
                    d_aus[j * d_imt + i] * d_h0[(j + 1) * d_imt + i] +
                    d_auw[j * d_imt + i] * d_h0[j * d_imt + i - 1] +
                    d_ausw[j * d_imt + i] * d_h0[(j + 1) * d_imt + i - 1];
        }

        if (i > 1 && i < d_imt - 2 && j > 1 && j < d_jmt - 2) {
            d_wka[j * d_imt + i] = d_viv[j * d_imt + i] *
                                   (d_wka[j * d_imt + i] + d_dlub[j * d_imt + i] -
                                    d_fcor[j * d_imt + i] * d_vbp[j * d_imt + i] +
                                    d_pax[j * d_imt + i] + d_pxb[j * d_imt + i] -
                                    d_work[j * d_imt + i] * d_whx[j * d_imt + i]);

//zf change d_vb to d_vbp
            d_wka[1 * d_jmt * d_imt + j * d_imt + i] =
                    d_viv[j * d_imt + i] *
                    (d_wka[1 * d_jmt * d_imt + j * d_imt + i] + d_dlvb[j * d_imt + i] +
                     d_fcor[j * d_imt + i] * d_ubp[j * d_imt + i] + 
                     d_pay[j * d_imt + i] + d_pyb[j * d_imt + i] -
                     d_work[j * d_imt + i] * d_why[j * d_imt + i]);
//zf change d_ub to d_ubp

            d_wka[2 * d_jmt * d_imt + j * d_imt + i] =
                    d_ebea[j * d_imt + i] * d_wka[j * d_imt + i] -
                    d_ebeb[j * d_imt + i] * d_wka[1 * d_jmt * d_imt + j * d_imt + i];

            d_wka[3 * d_jmt * d_imt + j * d_imt + i] =
                    d_ebea[j * d_imt + i] * d_wka[1 * d_jmt * d_imt + j * d_imt + i] +
                    d_ebeb[j * d_imt + i] * d_wka[j * d_imt + i];
        }
//
    }

    return;
}
#endif

#ifdef BIHAR

__global__ void bhdel4_zcurl(double *d_dxu, double *d_dyu,
                double *d_ubp,  double *d_vbp, double *d_wka,int *d_kmt) {
        int i, j, k;
        i = blockIdx.x * blockDim.x + threadIdx.x;
        j = blockIdx.y * blockDim.y + threadIdx.y;

        unsigned int tid = j * imt + i;

        if (i < imt && j < jmt) {
             d_wka[ tid] =0.0e0;
        if ( i >= 1 &&  j >= 1) {
                     if(0<=d_kmt[tid] - 1) {
                        d_wka[ tid] =
                                0.5 *
                                ((d_vbp[ tid] * d_dyu[tid] +
                                  d_vbp[ tid - imt] * d_dyu[tid - imt] -
                                  d_vbp[ tid - 1] * d_dyu [tid - 1] -
                                  d_vbp[ tid - imt - 1] * d_dyu [tid - imt - 1] -
                                  d_ubp[ tid] * d_dxu[tid] -
                                  d_ubp[ tid - 1] * d_dxu[tid - 1 ] +
                                  d_ubp[ tid -imt ] * d_dxu[tid - imt] +
                                  d_ubp[ tid -imt -1 ] * d_dxu[tid - imt - 1]));
                }
      }
  }
 return;
}

__global__ void bhdel4_div (double *d_hts, double *d_htw,
                double *d_tarea_r, double *d_ubp,  double *d_vbp, double *d_wka,int *d_kmt) {
        int i, j, k;
        i = blockIdx.x * blockDim.x + threadIdx.x;
        j = blockIdx.y * blockDim.y + threadIdx.y;

        unsigned int tid = j * imt + i;

        if (i < imt && j < jmt ) {
             d_wka[ tid] =0.0e0;
        if (i < imt - 1 && j >= 1) {
                     if(0<=d_kmt[tid] - 1) {
                        d_wka[tid] =
                                0.5 *
                                (( d_ubp[tid + 1] +
                                  d_ubp[tid - imt + 1]) * d_htw[tid + 1] -
                                 (d_ubp[ tid] +
                                  d_ubp[ tid - imt]) * d_htw[tid] +
                                 (d_vbp[tid + 1] +
                                  d_vbp[ tid]) * d_hts[tid] -
                                 (d_vbp[tid - imt + 1] +
                                  d_vbp[tid - imt]) * d_hts[tid - imt]) *
                                d_tarea_r[tid];
                }
      }
  }
 return;
}

__global__ void bhdel4_grad(double *d_wka,double *d_dxur, double *d_dyur, int *d_kmu,double *d_am_factor){ //d_pp=d_wka

    int i, j, k;
    double gradx,grady;

    i = blockIdx.x * blockDim.x + threadIdx.x;
    j = blockIdx.y * blockDim.y + threadIdx.y;

    if (j < jmt && i < imt) {


        gradx = 0.0;//c0
        grady= 0.0;//c0

        if ((1 <= i) && (j < (jmt - 1))) {
            if (0 <= d_kmu[j * imt + i] - 1) {
                gradx =
                        d_dxur[j * imt + i] * p5 *
                        (d_wka[ (j + 1) * imt + i] -
                         d_wka[j * imt + i - 1] -
                         d_wka[ (j + 1) * imt + i - 1] +
                         d_wka[j * imt + i]);

                grady=
                        d_dyur[j * imt + i] * p5 *
                        (d_wka[(j + 1) * imt + i] -
                         d_wka[ j * imt + i - 1] +
                         d_wka[ (j + 1) * imt + i - 1] -
                         d_wka[ j * imt + i]);
            }
        }
         d_am_factor[ j * imt + i]+=gradx*gradx+grady*grady;
    }

    return;
}


__global__ void barotr2_udel4_0( double *d_ubp, double *d_vbp,
                           double *d_dlub, double *d_dlvb, double *d_uarea_r, double *d_work, double *d_wgp,
                           double *d_pax, double *d_pay, double *d_pxb, double *d_pyb, double *d_whx, double *d_why,
                           double d_am, double *d_duc, double *d_dum, double *d_dun,
                           double *d_dus, double *d_due, double *d_duw, double *d_dmc,
                           double *d_dmn, double *d_dms, double *d_dme, double *d_dmw,
                           double *d_dxur, double *d_dyur, int *d_kmu,
                           double *d_au0, double *d_ausw, double *d_aus, double *d_auw,
                           int d_jb, int d_je, int d_ib, int d_ie, int d_imt, int d_jmt, int d_km,
                           double *d_amf, double *d_am_factor, double *d_d2uk, double *d_d2vk) {

    int i, j;

    i = (blockIdx.x) * blockDim.x + threadIdx.x;
    j = (blockIdx.y) * blockDim.y + threadIdx.y;

    if (i < d_imt && j < d_jmt) {
        double hduk = 0, hdvk = 0;
             double am_f,dxdy,t_area;
              am_f=1.0e0;

//        if (i > 1 && i < d_imt - 2 && j > 1 && j < d_jmt - 2) {
                if (j >= d_jb - 2 && j <= d_je && i >= d_ib - 2 && i <= d_ie) {
            if (d_kmu[j * d_imt + i] - 1 >= 0) {
                          t_area=sqrt(1.0e0/d_uarea_r[j * d_imt + i]);
                         dxdy=t_area*t_area*t_area*t_area*t_area*45.0e0;
                          am_f= sqrt (d_am_factor[j * d_imt + i]) * dxdy /abs(d_am*d_amf[j * d_imt + i]);
           //                am_f=  dxdy /abs(d_am*d_amf[j * d_imt + i]);
                          am_f=fmin(40.0e0,am_f);
                          am_f=fmax(1.0e0,am_f);

            
                    hduk =  (((d_duc[j * d_imt + i] + d_dum[j * d_imt + i]) * d_ubp[j * d_imt + i] +
                             d_dun[j * d_imt + i] * d_ubp[(j - 1) * d_imt + i] +
                             d_dus[j * d_imt + i] * d_ubp[(j + 1) * d_imt + i] +
                             d_due[j * d_imt + i] * d_ubp[j * d_imt + i + 1] +
                             d_duw[j * d_imt + i] * d_ubp[j * d_imt + i - 1]) +
                            (d_dmc[j * d_imt + i] * d_vbp[j * d_imt + i] +
                             d_dmn[j * d_imt + i] * d_vbp[(j - 1) * d_imt + i] +
                             d_dms[j * d_imt + i] * d_vbp[(j + 1) * d_imt + i] +
                             d_dme[j * d_imt + i] * d_vbp[j * d_imt + i + 1] +
                             d_dmw[j * d_imt + i] * d_vbp[j * d_imt + i - 1]));

                    hdvk = (((d_duc[j * d_imt + i] + d_dum[j * d_imt + i]) * d_vbp[j * d_imt + i] +
                             d_dun[j * d_imt + i] * d_vbp[(j - 1) * d_imt + i] +
                             d_dus[j * d_imt + i] * d_vbp[(j + 1) * d_imt + i] +
                             d_due[j * d_imt + i] * d_vbp[j * d_imt + i + 1] +
                             d_duw[j * d_imt + i] * d_vbp[j * d_imt + i - 1]) -
                            (d_dmc[j * d_imt + i] * d_ubp[j * d_imt + i] +
                             d_dmn[j * d_imt + i] * d_ubp[(j - 1) * d_imt + i] +
                             d_dms[j * d_imt + i] * d_ubp[(j + 1) * d_imt + i] +
                             d_dme[j * d_imt + i] * d_ubp[j * d_imt + i + 1] +
                             d_dmw[j * d_imt + i] * d_ubp[j * d_imt + i - 1]));
                }
            }
      //  }
        d_d2uk[j * d_imt + i] =am_f * d_amf[j * d_imt + i] * hduk;
        d_d2vk[j * d_imt + i] = am_f * d_amf[j * d_imt + i] * hdvk;
    }

    return;
}


__global__ void barotr2_udel4_1(double d_dtb, double *d_viv, double *d_ebea, double *d_ebeb, double *d_fcor,
                           double *d_h0, double *d_h0p, double *d_ubp, double *d_vbp,
                           double *d_dlub, double *d_dlvb, double *d_wka, double *d_work, double *d_wgp,
                           double *d_pax, double *d_pay, double *d_pxb, double *d_pyb, double *d_whx, double *d_why,
                           double d_am, double *d_duc, double *d_dum, double *d_dun,
                           double *d_dus, double *d_due, double *d_duw, double *d_dmc,
                           double *d_dmn, double *d_dms, double *d_dme, double *d_dmw,
                           double *d_dxur, double *d_dyur, int *d_kmu,
                           double *d_au0, double *d_ausw, double *d_aus, double *d_auw,
                           int d_jb, int d_je, int d_ib, int d_ie, int d_imt, int d_jmt, int d_km,
                           double *d_amf, double *d_am_factor, double *d_d2uk, double *d_d2vk) {
    int i, j;

    i = (blockIdx.x) * blockDim.x + threadIdx.x;
    j = (blockIdx.y) * blockDim.y + threadIdx.y;

    if (i < d_imt && j < d_jmt) {
        double hduk = 0, hdvk = 0;

//        if (i > 1 && i < d_imt - 2 && j > 1 && j < d_jmt - 2) {
                if (j >= d_jb - 1 && j < d_je && i >= d_ib - 1 && i < d_ie) {
            if (d_kmu[j * d_imt + i] - 1 >= 0) {
                    hduk = d_am *
                           (((d_duc[j * d_imt + i] + d_dum[j * d_imt + i]) * d_d2uk[j * d_imt + i] +
                             d_dun[j * d_imt + i] * d_d2uk[(j - 1) * d_imt + i] +
                             d_dus[j * d_imt + i] * d_d2uk[(j + 1) * d_imt + i] +
                             d_due[j * d_imt + i] * d_d2uk[j * d_imt + i + 1] +
                             d_duw[j * d_imt + i] * d_d2uk[j * d_imt + i - 1]) +
                            (d_dmc[j * d_imt + i] * d_d2vk[j * d_imt + i] +
                             d_dmn[j * d_imt + i] * d_d2vk[(j - 1) * d_imt + i] +
                             d_dms[j * d_imt + i] * d_d2vk[(j + 1) * d_imt + i] +
                             d_dme[j * d_imt + i] * d_d2vk[j * d_imt + i + 1] +
                             d_dmw[j * d_imt + i] * d_d2vk[j * d_imt + i - 1]));

                    hdvk = d_am *
                           (((d_duc[j * d_imt + i] + d_dum[j * d_imt + i]) * d_d2vk[j * d_imt + i] +
                             d_dun[j * d_imt + i] * d_d2vk[(j - 1) * d_imt + i] +
                             d_dus[j * d_imt + i] * d_d2vk[(j + 1) * d_imt + i] +
                             d_due[j * d_imt + i] * d_d2vk[j * d_imt + i + 1] +
                             d_duw[j * d_imt + i] * d_d2vk[j * d_imt + i - 1]) -
                            (d_dmc[j * d_imt + i] * d_d2uk[j * d_imt + i] +
                             d_dmn[j * d_imt + i] * d_d2uk[(j - 1) * d_imt + i] +
                             d_dms[j * d_imt + i] * d_d2uk[(j + 1) * d_imt + i] +
                             d_dme[j * d_imt + i] * d_d2uk[j * d_imt + i + 1] +
                             d_dmw[j * d_imt + i] * d_d2uk[j * d_imt + i - 1]));
                }
            }
      //  }

        double gradx = 0.0e0, grady = 0.0e0;
        if (i > 1 && i < d_imt - 2 && j > 1 && j < d_jmt - 2) {
            if (d_kmu[j * d_imt + i] - 1 >= 0) {
                            gradx = d_dxur[j * d_imt + i] * 0.5 *
                        (d_h0[(j + 1) * d_imt + i] -
                         d_h0[j * d_imt + i - 1] -
                         d_h0[(j + 1) * d_imt + i - 1] +
                         d_h0[j * d_imt + i]);

                grady = d_dyur[j * d_imt + i] * 0.5 *
                        (d_h0[(j + 1) * d_imt + i] -
                         d_h0[j * d_imt + i - 1] +
                         d_h0[(j + 1) * d_imt + i - 1] -
                          d_h0[j * d_imt + i]);
                        }
            double gstar = (d_wgp[j * d_imt + i] - 1.0) * g;
            d_wka[j * d_imt + i] = hduk + gstar * gradx;
            d_wka[1 * d_jmt * d_imt + j * d_imt + i] = hdvk + gstar * grady;

            d_wka[4 * d_jmt * d_imt + j * d_imt + i] = hduk;
            d_wka[5 * d_jmt * d_imt + j * d_imt + i] = hdvk;

                }

//        d_h0[j * d_imt + i] = d_h0p[j * d_imt + i]+ d_work[j * d_imt + i] * d_dtb;
        if (j == d_jmt - 1 || i == 0) {
            d_work[j * d_imt + i] = 0.0e0;
        } else {
            d_work[j * d_imt + i] =
                    d_au0[j * d_imt + i] * d_h0[j * d_imt + i] +
                    d_aus[j * d_imt + i] * d_h0[(j + 1) * d_imt + i] +
                    d_auw[j * d_imt + i] * d_h0[j * d_imt + i - 1] +
                    d_ausw[j * d_imt + i] * d_h0[(j + 1) * d_imt + i - 1];
        }

        if (i > 1 && i < d_imt - 2 && j > 1 && j < d_jmt - 2) {
            d_wka[j * d_imt + i] = d_viv[j * d_imt + i] *
                                   (d_wka[j * d_imt + i] + d_dlub[j * d_imt + i] -
                                    d_fcor[j * d_imt + i] * d_vbp[j * d_imt + i] +
                                    d_pax[j * d_imt + i] + d_pxb[j * d_imt + i] -
                                    d_work[j * d_imt + i] * d_whx[j * d_imt + i]);

//zf change d_vb to d_vbp
            d_wka[1 * d_jmt * d_imt + j * d_imt + i] =
                    d_viv[j * d_imt + i] *
                    (d_wka[1 * d_jmt * d_imt + j * d_imt + i] + d_dlvb[j * d_imt + i] +
                     d_fcor[j * d_imt + i] * d_ubp[j * d_imt + i] + 
                     d_pay[j * d_imt + i] + d_pyb[j * d_imt + i] -
                     d_work[j * d_imt + i] * d_why[j * d_imt + i]);
//zf change d_ub to d_ubp

            d_wka[2 * d_jmt * d_imt + j * d_imt + i] =
                    d_ebea[j * d_imt + i] * d_wka[j * d_imt + i] -
                    d_ebeb[j * d_imt + i] * d_wka[1 * d_jmt * d_imt + j * d_imt + i];

            d_wka[3 * d_jmt * d_imt + j * d_imt + i] =
                    d_ebea[j * d_imt + i] * d_wka[1 * d_jmt * d_imt + j * d_imt + i] +
                    d_ebeb[j * d_imt + i] * d_wka[j * d_imt + i];
        }

    }

    return;
}

__global__ void barotr2_udel4_2(double d_dtb, double *d_viv, double *d_ebea, double *d_ebeb, double *d_fcor,
                           double *d_h0, double *d_h0p, double *d_ubp, double *d_vbp,
                           double *d_dlub, double *d_dlvb, double *d_wka, double *d_work, double *d_wgp,
                           double *d_pax, double *d_pay, double *d_pxb, double *d_pyb, double *d_whx, double *d_why,
                           double d_am, double *d_duc, double *d_dum, double *d_dun,
                           double *d_dus, double *d_due, double *d_duw, double *d_dmc,
                           double *d_dmn, double *d_dms, double *d_dme, double *d_dmw,
                           double *d_dxur, double *d_dyur, int *d_kmu,
                           double *d_au0, double *d_ausw, double *d_aus, double *d_auw,
                           int d_jb, int d_je, int d_ib, int d_ie, int d_imt, int d_jmt, int d_km) {
    int i, j;

    i = (blockIdx.x) * blockDim.x + threadIdx.x;
    j = (blockIdx.y) * blockDim.y + threadIdx.y;

//---------------------------------------------------------------------
//     COMPUTE THE "ARTIFICIAL" HORIZONTAL VISCOSITY
//---------------------------------------------------------------------
    if (i < d_imt && j < d_jmt) {
        double hduk = 0, hdvk = 0;

        if (i > 1 && i < d_imt - 2 && j > 1 && j < d_jmt - 2) {
            if (d_kmu[j * d_imt + i] - 1 >= 0) {
                if (j >= d_jb - 1 && j < d_je && i >= d_ib - 1 && i < d_ie) {
//#ifdef BIHAR
		hduk=d_wka[4 * d_jmt * d_imt + j * d_imt + i];
		hdvk=d_wka[5 * d_jmt * d_imt + j * d_imt + i];
//#endif
                }
            }
        }

        double gradx = 0.0e0, grady = 0.0e0;
        if (i > 1 && i < d_imt - 2 && j > 1 && j < d_jmt - 2) {
            if (d_kmu[j * d_imt + i] - 1 >= 0) {
                            gradx = d_dxur[j * d_imt + i] * 0.5 *
                        (d_h0[(j + 1) * d_imt + i] -
                         d_h0[j * d_imt + i - 1] -
                         d_h0[(j + 1) * d_imt + i - 1] +
                         d_h0[j * d_imt + i]);

                grady = d_dyur[j * d_imt + i] * 0.5 *
                        (d_h0[(j + 1) * d_imt + i] -
                         d_h0[j * d_imt + i - 1] +
                         d_h0[(j + 1) * d_imt + i - 1] -
                          d_h0[j * d_imt + i]);
                        }
            double gstar = (d_wgp[j * d_imt + i] - 1.0) * g;
            d_wka[j * d_imt + i] = hduk + gstar * gradx;
            d_wka[1 * d_jmt * d_imt + j * d_imt + i] = hdvk + gstar * grady;

                }

//        d_h0[j * d_imt + i] = d_h0p[j * d_imt + i]+ d_work[j * d_imt + i] * d_dtb;
        if (j == d_jmt - 1 || i == 0) {
            d_work[j * d_imt + i] = 0.0e0;
        } else {
            d_work[j * d_imt + i] =
                    d_au0[j * d_imt + i] * d_h0[j * d_imt + i] +
                    d_aus[j * d_imt + i] * d_h0[(j + 1) * d_imt + i] +
                    d_auw[j * d_imt + i] * d_h0[j * d_imt + i - 1] +
                    d_ausw[j * d_imt + i] * d_h0[(j + 1) * d_imt + i - 1];
        }

        if (i > 1 && i < d_imt - 2 && j > 1 && j < d_jmt - 2) {
            d_wka[j * d_imt + i] = d_viv[j * d_imt + i] *
                                   (d_wka[j * d_imt + i] + d_dlub[j * d_imt + i] -
                                    d_fcor[j * d_imt + i] * d_vbp[j * d_imt + i] +
                                    d_pax[j * d_imt + i] + d_pxb[j * d_imt + i] -
                                    d_work[j * d_imt + i] * d_whx[j * d_imt + i]);

//zf change d_vb to d_vbp
            d_wka[1 * d_jmt * d_imt + j * d_imt + i] =
                    d_viv[j * d_imt + i] *
                    (d_wka[1 * d_jmt * d_imt + j * d_imt + i] + d_dlvb[j * d_imt + i] +
                     d_fcor[j * d_imt + i] * d_ubp[j * d_imt + i] + 
                     d_pay[j * d_imt + i] + d_pyb[j * d_imt + i] -
                     d_work[j * d_imt + i] * d_why[j * d_imt + i]);
//zf change d_ub to d_ubp

            d_wka[2 * d_jmt * d_imt + j * d_imt + i] =
                    d_ebea[j * d_imt + i] * d_wka[j * d_imt + i] -
                    d_ebeb[j * d_imt + i] * d_wka[1 * d_jmt * d_imt + j * d_imt + i];

            d_wka[3 * d_jmt * d_imt + j * d_imt + i] =
                    d_ebea[j * d_imt + i] * d_wka[1 * d_jmt * d_imt + j * d_imt + i] +
                    d_ebeb[j * d_imt + i] * d_wka[j * d_imt + i];
        }
//
    }

    return;
}
#endif //BIHAR

__global__ void barotr3_cu(double d_dtb, double *d_dzph, double *d_ub, double *d_vb, double *d_ubp, double *d_vbp, double *d_wka, double *d_work,
                            int d_imt, int d_jmt) {

    int i, j;

    i = (blockIdx.x) * blockDim.x + threadIdx.x;
    j = (blockIdx.y) * blockDim.y + threadIdx.y;

    if (i < d_imt && j < d_jmt) {

        d_ub[j * d_imt + i] =d_ubp[j*d_imt+i]+ d_wka[2 * d_jmt * d_imt + j * d_imt + i] * d_dtb;
        d_vb[j * d_imt + i] =d_vbp[j*d_imt+i]+ d_wka[3 * d_jmt * d_imt + j * d_imt + i] * d_dtb;
        //d_ub[j * d_imt + i] += d_wka[2 * d_jmt * d_imt + j * d_imt + i] * d_dtb;
        //d_vb[j * d_imt + i] += d_wka[3 * d_jmt * d_imt + j * d_imt + i] * d_dtb;
// zf changed
        if (j >= 1 && i < d_imt - 1) {
            d_wka[j * d_imt + i] = d_ub[j * d_imt + i] * (d_dzph[j * d_imt + i] + d_work[j * d_imt + i]);
            d_wka[1 * d_jmt * d_imt + j * d_imt + i] =
                    d_vb[j * d_imt + i] * (d_dzph[j * d_imt + i] + d_work[j * d_imt + i]);
        }
}
  return;
}

__global__ void barotr3_cu1(double *d_vit, double *d_h0p,
                           double *d_ub, double *d_vb, double *d_ubp, double *d_vbp, double *d_wka, double *d_work,
                           double *d_htw, double *d_hts, int *d_kmt, double *d_tarea_r,
                           double d_ah, int *d_kmtn, int *d_kmts, int *d_kmte, int *d_kmtw,
                           double *d_dtn, double *d_dts, double *d_dte, double *d_dtw,
                           int d_jb, int d_je, int d_ib, int d_ie, int d_nc, int d_imt, int d_jmt, int d_km,
                           double *d_d2tk, double *d_ahf) {

    int i, j;

    i = (blockIdx.x) * blockDim.x + threadIdx.x;
    j = (blockIdx.y) * blockDim.y + threadIdx.y;

        if (i < d_imt - 1 && j >= 1 && j < d_jmt) {
            double div_out = 0, hdtk = 0;
            if (0 <= d_kmt[j * d_imt + i] - 1) {
                div_out =
                        0.5 * ((d_wka[j * d_imt + i + 1] + d_wka[(j - 1) * d_imt + i + 1]) * d_htw[j * d_imt + i + 1] -
                               (d_wka[j * d_imt + i] + d_wka[(j - 1) * d_imt + i]) * d_htw[j * d_imt + i] +
                               (d_wka[1 * d_imt * d_jmt + j * d_imt + i + 1] +
                                d_wka[1 * d_imt * d_jmt + j * d_imt + i]) *
                               d_hts[j * d_imt + i] -
                               (d_wka[1 * d_imt * d_jmt + (j - 1) * d_imt + i + 1] +
                                d_wka[1 * d_imt * d_jmt + (j - 1) * d_imt + i]) * d_hts[(j - 1) * d_imt + i]) *
                        d_tarea_r[j * d_imt + i];
            }

//wpf 0309, sd mod(nc,4); old version mod(nc,2)
            if ((d_nc + 1) % 4 == 0) { //4 or 2
#ifdef  BIHAR
//				hmix_del4_hdifft_del4_(&tmpt1, hdtk, dyn_mod_mp_h0p_[iblock], &jb, &je, &ib, &ie, &iblock);
                double cc = 0;
                if (j >= d_jb - 1 && j < d_je && i >= d_ib - 1 && i < d_ie) {

                if (0 <= d_kmt[j * d_imt + i] - 1) {
                    double c;

                    if (0 <= d_kmtn[j * d_imt + i] - 1) {
                        c = d_dtn[j * d_imt + i];
                        hdtk += c * d_d2tk[(j - 1) * d_imt + i];
                        cc -= c;
                    }


                    if (0 <= d_kmts[j * d_imt + i] - 1) {
                        c = d_dts[j * d_imt + i];
                        hdtk += c * d_d2tk[(j + 1) * d_imt + i];
                        cc -= c;
                    }

                    if (0 <= d_kmte[j * d_imt + i] - 1) {
                        c = d_dte[j * d_imt + i];
                        hdtk += c * d_d2tk[j * d_imt + i + 1];
                        cc -= c;
                    }

                    if (0 <= d_kmtw[j * d_imt + i] - 1) {
                        c = d_dtw[j * d_imt + i];
                        hdtk += c * d_d2tk[j * d_imt + i - 1];
                        cc -= c;
                    }
                }


                    hdtk = d_ah * (cc * d_d2tk[j * d_imt + i] + hdtk);
                } else {
                    hdtk = 0;
                }
#else
                double cc = 0;

                if (j >= d_jb - 1 && j < d_je && i >= d_ib - 1 && i < d_ie) {
                if (0 <= d_kmt[j * d_imt + i] - 1) {
                    double c;

                    if (0 <= d_kmtn[j * d_imt + i] - 1) {
                        c = d_dtn[j * d_imt + i];
                        hdtk += c * d_h0p[(j - 1) * d_imt + i];
                        cc -= c;
                    }


                    if (0 <= d_kmts[j * d_imt + i] - 1) {
                        c = d_dts[j * d_imt + i];
                        hdtk += c * d_h0p[(j + 1) * d_imt + i];
                        cc -= c;
                    }

                    if (0 <= d_kmte[j * d_imt + i] - 1) {
                        c = d_dte[j * d_imt + i];
                        hdtk += c * d_h0p[j * d_imt + i + 1];
                        cc -= c;
                    }

                    if (0 <= d_kmtw[j * d_imt + i] - 1) {
                        c = d_dtw[j * d_imt + i];
                        hdtk += c * d_h0p[j * d_imt + i - 1];
                        cc -= c;
                    }
                }


                    hdtk = d_ah * (cc * d_h0p[j * d_imt + i] + hdtk);
                } else {
                    hdtk = 0;
                }
#endif
            }

            d_work[j * d_imt + i] = (hdtk*1.0e0 - div_out) * d_vit[j * d_imt + i];
             //test
        }
    return;
}

__global__ void barotr3_tdel4_1(double *d_vit, double *d_h0p,
                           double *d_ub, double *d_vb, double *d_ubp, double *d_vbp, double *d_wka, double *d_work,
                           double *d_htw, double *d_hts, int *d_kmt, double *d_tarea_r,
                           double d_ah, int *d_kmtn, int *d_kmts, int *d_kmte, int *d_kmtw,
                           double *d_dtn, double *d_dts, double *d_dte, double *d_dtw,
                           int d_jb, int d_je, int d_ib, int d_ie, int d_nc, int d_imt, int d_jmt, int d_km,
                           double *d_d2tk, double *d_ahf) {

    int i, j;

    i = (blockIdx.x) * blockDim.x + threadIdx.x;
    j = (blockIdx.y) * blockDim.y + threadIdx.y;

        if (i < d_imt && j < d_jmt) {
            d_d2tk[j * d_imt + i] = 0.0e0;
        }
        if (i < d_imt - 1 && j >= 1 && j < d_jmt ) {
            double div_out = 0, hdtk = 0;
            if ((d_nc + 1) % 4 == 0) {
#ifdef  BIHAR
                double cc = 0;
                if (j >= d_jb - 2 && j <= d_je && i >= d_ib - 2 && i <= d_ie) {

                if (0 <= d_kmt[j * d_imt + i] - 1) {
                    double c;

                    if (0 <= d_kmtn[j * d_imt + i] - 1) {
                        c = d_dtn[j * d_imt + i];
                        hdtk += c * d_h0p[(j - 1) * d_imt + i];
                        cc -= c;
                    }


                    if (0 <= d_kmts[j * d_imt + i] - 1) {
                        c = d_dts[j * d_imt + i];
                        hdtk += c * d_h0p[(j + 1) * d_imt + i];
                        cc -= c;
                    }

                    if (0 <= d_kmte[j * d_imt + i] - 1) {
                        c = d_dte[j * d_imt + i];
                        hdtk += c * d_h0p[j * d_imt + i + 1];
                        cc -= c;
                    }

                    if (0 <= d_kmtw[j * d_imt + i] - 1) {
                        c = d_dtw[j * d_imt + i];
                        hdtk += c * d_h0p[j * d_imt + i - 1];
                        cc -= c;
                    }
                }


                    hdtk = d_ahf[j * d_imt + i] * (cc * d_h0p[j * d_imt + i] + hdtk);
                } else {
                    hdtk = 0;
                }
#endif
            }
            d_d2tk[j * d_imt + i] = hdtk;

        }
    return;
}

__global__ void barotr3_tdel4_2(double *d_vit, double *d_h0p,
                           double *d_ub, double *d_vb, double *d_ubp, double *d_vbp, double *d_wka, double *d_work,
                           double *d_htw, double *d_hts, int *d_kmt, double *d_tarea_r,
                           double d_ah, int *d_kmtn, int *d_kmts, int *d_kmte, int *d_kmtw,
                           double *d_dtn, double *d_dts, double *d_dte, double *d_dtw,
                           int d_jb, int d_je, int d_ib, int d_ie, int d_nc, int d_imt, int d_jmt, int d_km,
                           double *d_hdtk, double *d_ahf) {

    int i, j;

    i = (blockIdx.x) * blockDim.x + threadIdx.x;
    j = (blockIdx.y) * blockDim.y + threadIdx.y;

        if (i >= 2 && i < d_imt - 2 && j >= 1 && j < d_jmt - 2) {
            double div_out = 0, hdtk = 0;
            if (0 <= d_kmt[j * d_imt + i] - 1) {
                div_out =
                        0.5 * ((d_wka[j * d_imt + i + 1] + d_wka[(j - 1) * d_imt + i + 1]) * d_htw[j * d_imt + i + 1] -
                               (d_wka[j * d_imt + i] + d_wka[(j - 1) * d_imt + i]) * d_htw[j * d_imt + i] +
                               (d_wka[1 * d_imt * d_jmt + j * d_imt + i + 1] +
                                d_wka[1 * d_imt * d_jmt + j * d_imt + i]) *
                               d_hts[j * d_imt + i] -
                               (d_wka[1 * d_imt * d_jmt + (j - 1) * d_imt + i + 1] +
                                d_wka[1 * d_imt * d_jmt + (j - 1) * d_imt + i]) * d_hts[(j - 1) * d_imt + i]) *
                        d_tarea_r[j * d_imt + i];
            }

//wpf 0309, sd mod(nc,4); old version mod(nc,2)
            if ((d_nc + 1) % 4 == 0) {
        	hdtk=d_hdtk[j * d_imt + i];
            }

            d_work[j * d_imt + i] = (hdtk*1.0e0 - div_out) * d_vit[j * d_imt + i];
        }
    return;
}


__global__ void barotr4_cu(double d_dtb, double *d_h0, double *d_h0p, double *d_ub,
                           double *d_vb, double *d_ubp, double *d_vbp, double *d_h0f,
                           double *d_h0bf, double *d_work, int d_imt, int d_jmt) {

    int i, j;

    i = (blockIdx.x) * blockDim.x + threadIdx.x;
    j = (blockIdx.y) * blockDim.y + threadIdx.y;

    if (i < d_imt && j < d_jmt) {
        double h0 = d_h0p[j * d_imt + i] + (d_work[j * d_imt + i] * d_dtb);
//          d_h0[j * d_imt + i]= d_h0p[j * d_imt + i];
        d_h0p[j * d_imt + i] = h0;

        d_ubp[j * d_imt + i] = d_ub[j * d_imt + i];
        d_vbp[j * d_imt + i] = d_vb[j * d_imt + i];

        d_h0f[j * d_imt + i] += h0;
        d_h0bf[j * d_imt + i] += h0;

        d_h0[j * d_imt + i] = h0;
  //      d_h0[j * d_imt + i] = d_h0p[j * d_imt + i];
    }

    return;
}

extern "C" void barotr_() {
    struct block this_block;
    this_block = blocks_mp_all_blocks_[domain_mp_blocks_clinic_[0] - 1];
    int jb = this_block.jb;
    int je = this_block.je;
    int ib = this_block.ib;
    int ie = this_block.ie;

    int nc, i, j;
    int errorcode1, errorcode2;

    dim3 blockSize(16, 8, 1);
    dim3 gridSize((imt + blockSize.x - 1) / blockSize.x, (jmt + blockSize.y - 1) / blockSize.y, 1);
        iCheck();

   CHECK(hipMemset(d_wka, c0, dataSize3d));
   for (nc = 0; nc < pconst_mod_mp_nbb_; nc++) { //jjrtest
   //for (nc = 0; nc < 1; nc++) {
   CHECK(hipMemset(d_wka, c0, dataSize3d));
   CHECK(hipMemset(d_work, c0, dataSize2d));
//    for (nc = 0; nc < 1; nc++) {
        hipLaunchKernelGGL(barotr1_cu, dim3(gridSize), dim3(blockSize ), 0, 0, d_work, d_wka, d_dzph, d_h0, d_ub, d_vb,
                d_au0, d_ausw, d_aus, d_auw,imt, jmt);//d_work/d_wka changed
        iCheck();

        hipLaunchKernelGGL(barotr1_cu1, dim3(gridSize), dim3(blockSize ), 0, 0, d_work, d_wka, d_vit, d_htw, d_hts, d_kmt, d_tarea_r,
                imt, jmt);//d_work/d_wka changed
        iCheck();
       
        hipMemcpy(work_mod_mp_work_, d_work, dataSize2d, hipMemcpyDeviceToHost);
        pop_haloupdate_barotr1_(&errorcode1);//work
        hipMemcpy(d_work, work_mod_mp_work_, dataSize2d, hipMemcpyHostToDevice);

      hipLaunchKernelGGL(barotr2_cu1, dim3(gridSize), dim3(blockSize ), 0, 0, pconst_mod_mp_dtb_,  d_h0,
                d_h0p,  d_work, imt,jmt);//d_work/d_wka/d_h0 changed

#ifdef BIHAR

       hipLaunchKernelGGL(bhdel4_div, dim3(gridSize), dim3(blockSize), 0, 0, d_hts, d_htw, d_tarea_r,
               d_ubp,d_vbp,d_wka,d_kmt) ;
        hipMemset(d_am_factor, c0, dataSize3d);
       hipLaunchKernelGGL(bhdel4_grad, dim3(gridSize), dim3(blockSize), 0, 0, d_wka, d_dxur, d_dyur,
              d_kmu,d_am_factor) ;
       hipLaunchKernelGGL(bhdel4_zcurl, dim3(gridSize), dim3(blockSize), 0, 0, d_dxu, d_dyu,
               d_ubp,d_vbp,d_wka,d_kmt) ;

       hipLaunchKernelGGL(bhdel4_grad, dim3(gridSize), dim3(blockSize ), 0, 0, d_wka, d_dxur, d_dyur,
               d_kmu,d_am_factor) ;

        hipLaunchKernelGGL(barotr2_udel4_0, dim3(gridSize), dim3(blockSize ), 0, 0, 
                d_ubp, d_vbp, d_dlub, d_dlvb, d_uarea_r, d_work, d_wgp, d_pax, d_pay, d_pxb, 
                d_pyb, d_whx, d_why, hmix_del4_mp_am_, d_duc, d_dum, d_dun, d_dus, d_due,
                d_duw, d_dmc, d_dmn, d_dms, d_dme, d_dmw, d_dxur, d_dyur,
                d_kmu, d_au0, d_ausw, d_aus, d_auw, jb, je, ib, ie, imt, jmt, km, 
                d_amf, d_am_factor, d_d2uk, d_d2vk);
        iCheck();
        hipLaunchKernelGGL(barotr2_udel4_1, dim3(gridSize), dim3(blockSize ), 0, 0, pconst_mod_mp_dtb_, d_viv, d_ebea, d_ebeb, d_fcor, d_h0,
                d_h0p, d_ubp, d_vbp, d_dlub, d_dlvb, d_wka, d_work, d_wgp, d_pax, d_pay, d_pxb, 
                d_pyb, d_whx, d_why, hmix_del4_mp_am_, d_duc, d_dum, d_dun, d_dus, d_due,
                d_duw, d_dmc, d_dmn, d_dms, d_dme, d_dmw, d_dxur, d_dyur,
                d_kmu, d_au0, d_ausw, d_aus, d_auw, jb, je, ib, ie, imt, jmt, km, 
                d_amf, d_am_factor, d_d2uk, d_d2vk);
        iCheck();

/*
//use F90 to do diffu_del4, when HIP code ready,change here
        CHECK(hipMemcpy(dyn_mod_mp_ubp_, d_ubp, dataSize2d, hipMemcpyDeviceToHost));
        CHECK(hipMemcpy(dyn_mod_mp_vbp_, d_vbp, dataSize2d, hipMemcpyDeviceToHost));
        barotr_f_();
        hipMemcpy(&d_wka[4 * imt * jmt], work_mod_mp_wka_[0][4], dataSize2d * 2, hipMemcpyHostToDevice);
//put hduk,hdvk in wka[4],[5]
        hipLaunchKernelGGL(barotr2_udel4_2, dim3(gridSize), dim3(blockSize ), 0, 0, pconst_mod_mp_dtb_, d_viv, d_ebea, d_ebeb, d_fcor, d_h0,
                d_h0p, d_ubp, d_vbp, d_dlub, d_dlvb, d_wka, d_work, d_wgp, d_pax, d_pay, d_pxb, 
                d_pyb, d_whx, d_why, hmix_del4_mp_am_, d_duc, d_dum, d_dun, d_dus, d_due,
                d_duw, d_dmc, d_dmn, d_dms, d_dme, d_dmw, d_dxur, d_dyur,
                d_kmu, d_au0, d_ausw, d_aus, d_auw, jb, je, ib, ie, imt, jmt, km);//d_work/d_wka/d_h0 changed
        iCheck();
*/
#else
        hipLaunchKernelGGL(barotr2_cu, dim3(gridSize), dim3(blockSize ), 0, 0, pconst_mod_mp_dtb_, d_viv, d_ebea, d_ebeb, d_fcor, d_h0,
                d_h0p, d_ubp, d_vbp, d_dlub, d_dlvb, d_wka, d_work, d_wgp, d_pax, d_pay, d_pxb, //zf change d_ub, d_vb to d_ubp, d_vbp
                d_pyb, d_whx, d_why, hmix_del2_mp_am_, d_duc, d_dum, d_dun, d_dus, d_due,
                d_duw, d_dmc, d_dmn, d_dms, d_dme, d_dmw, d_dxur, d_dyur,
                d_kmu, d_au0, d_ausw, d_aus, d_auw, jb, je, ib, ie, imt, jmt, km);//d_work/d_wka/d_h0 changed
        iCheck();
#endif
        hipMemcpy(work_mod_mp_wka_[0][2], &d_wka[2 * imt * jmt], dataSize2d * 2, hipMemcpyDeviceToHost);
        pop_haloupdate_barotr2_(&errorcode1, &errorcode2);//wka[:,:,2,:], wka[:,:,3,:]
        hipMemcpy(&d_wka[2 * imt * jmt], work_mod_mp_wka_[0][2], dataSize2d * 2, hipMemcpyHostToDevice);

        hipLaunchKernelGGL(barotr3_cu, dim3(gridSize), dim3(blockSize ), 0, 0, pconst_mod_mp_dtb_, d_dzph, d_ub, d_vb,
                d_ubp, d_vbp, d_wka, d_work,imt, jmt);//d_work/d_wka/d_ub/d_vb changed

        iCheck();

#ifdef BIHAR

        CHECK(hipMemset(d_d2tk, c0, dataSize2d));
        hipLaunchKernelGGL(barotr3_tdel4_1, dim3(gridSize), dim3(blockSize ), 0, 0, d_vit, d_h0p, d_ub, d_vb,
                d_ubp, d_vbp, d_wka, d_work, d_htw, d_hts, d_kmt, d_tarea_r,
                hmix_del4_mp_ah_, d_kmtn, d_kmts, d_kmte, d_kmtw, d_dtn, d_dts, d_dte, d_dtw,
                jb, je, ib, ie, nc, imt, jmt, km, d_d2tk, d_ahf);//d_d2tk changed
        iCheck();
        hipLaunchKernelGGL(barotr3_cu1, dim3(gridSize), dim3(blockSize ), 0, 0, d_vit, d_h0p, d_ub, d_vb,
                d_ubp, d_vbp, d_wka, d_work, d_htw, d_hts, d_kmt, d_tarea_r,
                hmix_del4_mp_ah_, d_kmtn, d_kmts, d_kmte, d_kmtw, d_dtn, d_dts, d_dte, d_dtw,
                jb, je, ib, ie, nc, imt, jmt, km, d_d2tk, d_ahf);//d_work/d_wka/d_ub/d_vb changed
        iCheck();
/*
//use F90 to do difft_del4, when HIP code ready,change here
        CHECK(hipMemcpy(dyn_mod_mp_h0p_, d_h0p, dataSize2d, hipMemcpyDeviceToHost));
        barotr_tdel_();
        hipMemcpy(d_hdtk, work_mod_mp_work_, dataSize2d, hipMemcpyHostToDevice);
//get hdtk back
        hipLaunchKernelGGL(barotr3_tdel4_2, dim3(gridSize), dim3(blockSize ), 0, 0, d_vit, d_h0p, d_ub, d_vb,
                d_ubp, d_vbp, d_wka, d_work, d_htw, d_hts, d_kmt, d_tarea_r,
                hmix_del4_mp_ah_, d_kmtn, d_kmts, d_kmte, d_kmtw, d_dtn, d_dts, d_dte, d_dtw,
                jb, je, ib, ie, nc, imt, jmt, km, d_hdtk, d_ahf);//d_d2tk changed
        iCheck();
*/
#else
        hipLaunchKernelGGL(barotr3_cu1, dim3(gridSize), dim3(blockSize ), 0, 0, d_vit, d_h0p, d_ub, d_vb,
                d_ubp, d_vbp, d_wka, d_work, d_htw, d_hts, d_kmt, d_tarea_r,
                hmix_del2_mp_ah_, d_kmtn, d_kmts, d_kmte, d_kmtw, d_dtn, d_dts, d_dte, d_dtw,
                jb, je, ib, ie, nc, imt, jmt, km, d_d2tk, d_ahf);//d_work/d_wka/d_ub/d_vb changed
        iCheck();
#endif //BIHAR


        hipMemcpy(work_mod_mp_work_, d_work, dataSize2d, hipMemcpyDeviceToHost);
        pop_haloupdate_barotr1_(&errorcode1);//d_work
        hipMemcpy(d_work, work_mod_mp_work_, dataSize2d, hipMemcpyHostToDevice);

       hipLaunchKernelGGL(barotr4_cu, dim3(gridSize), dim3(blockSize ), 0, 0, pconst_mod_mp_dtb_, d_h0, d_h0p, d_ub, d_vb, d_ubp,
                d_vbp, d_h0f, d_h0bf, d_work, imt, jmt);//d_h0/d_ubp/d_vbp/d_h0p/d_h0f/d_h0bf changed
//        hipMemcpy(d_h0, &d_wka[4 * imt * jmt], dataSize2d , hipMemcpyDeviceToDevice);

        iCheck();

        pconst_mod_mp_isb_++;
    }

    //deallocate_barotr_();//dlub,dlvb
    return;
}
