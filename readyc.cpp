/*by tsb 2018.11*/
// Fix some bugs by Wpf, 2019.1.16
// dmax1 -> fmax, dmin1 -> fmin, get_block, ...
//#include <math_functions.h>
#include "hip/hip_runtime.h"
#include <math.h>
#include "cuda_data.h"
#include <stdio.h>
#include "blocks.h"
#include "constant_mod.h"
#include "domain.h"
#include "dyn_mod.h"
#include "forc_mod.h"
#include "grid.h"
#include "hmix_del2.h"
#include "hmix_del4.h"
#include "LICOM_Error_mod.h"
#include "param_mod.h"
#include "pconst_mod.h"
#include "pmix_mod.h"
#include "precision_mod.h"
#include "work_mod.h"
#include "tracer_mod.h"
#include <hip/hip_runtime.h>
#include <stdio.h>
#include <float.h>
#include <math_functions.h>
#include "constant_mod.h"
#include "precision_mod.h"
#include "canuto_mod.h"
#include "tracer_mod.h"
#include "param_mod.h"
#include "pconst_mod.h"
#include "cuda_data.h"
#include "turb.h"
#include "common.h"
//#define sign_int(x,y) ((y)>=0?abs(x):-abs(x))

//extern "C" void allocate_readyc_();
//extern "C" void readycf_();
extern "C" void pop_haloupdate_readyc_(int *);
//extern "C" void readyc_udel4_();

__global__ void readyc1(double *d_h0bl, double *d_h0bf, double *d_h0, double *d_up,
		double *d_vp, double *d_viv, double *d_vit, double *d_wp12, double *d_wp13, double *d_at0,
		double *d_atn, double *d_ate, double *d_atne) {
	int i, j, k;
	i = blockIdx.x * blockDim.x + threadIdx.x;
	j = blockIdx.y * blockDim.y + threadIdx.y;
	k = blockIdx.z * blockDim.z + threadIdx.z;

	if (i < imt && j < jmt && k < km) {
		if (jst - 1 <= j && j < jet && k == 0) {
			d_h0bl[j * imt + i] = d_h0bf[j * imt + i];
			d_h0bf[j * imt + i] = d_h0[j * imt + i];
		}

		if (1 <= j && i < (imt - 1)) {
			d_wp12[k * jmt * imt + j * imt + i] =
				d_vit[k * jmt * imt + j * imt + i] *
				(d_at0[j * imt + i] *
				 d_up[k * jmt * imt + j * imt + i] *
				 d_viv[k * jmt * imt + j * imt + i] +
				 d_atn[j * imt + i] *
				 d_up[k * jmt * imt + (j - 1) * imt + i] *
				 d_viv[k * jmt * imt + (j - 1) * imt + i] +
				 d_ate[j * imt + i] *
				 d_up[k * jmt * imt + j * imt + i + 1] *
				 d_viv[k * jmt * imt + j * imt + i + 1] +
				 d_atne[j * imt + i] *
				 d_up[k * jmt * imt + (j - 1) * imt + i + 1] *
				 d_viv[k * jmt * imt + (j - 1) * imt + i + 1]) /
				(d_viv[k * jmt * imt + j * imt + i] *
				 d_at0[j * imt + i] +
				 d_viv[k * jmt * imt + (j - 1) * imt + i] *
				 d_atn[j * imt + i] +
				 d_viv[k * jmt * imt + j * imt + i + 1] *
				 d_ate[j * imt + i] +
				 d_viv[k * jmt * imt + (j - 1) * imt + i + 1] *
				 d_atne[j * imt + i] + 1.0e-25);

			d_wp13[k * jmt * imt + j * imt + i] =
				d_vit[k * jmt * imt + j * imt + i] *
				(d_at0[j * imt + i] *
				 d_vp[k * jmt * imt + j * imt + i] *
				 d_viv[k * jmt * imt + j * imt + i] +
				 d_atn[j * imt + i] *
				 d_vp[k * jmt * imt + (j - 1) * imt + i] *
				 d_viv[k * jmt * imt + (j - 1) * imt + i] +
				 d_ate[j * imt + i] * d_vp[k * jmt * imt + j * imt + i + 1] *
				 d_viv[k * jmt * imt + j * imt + i + 1] +
				 d_atne[j * imt + i] *
				 d_vp[k * jmt * imt + (j - 1) * imt + i + 1] *
				 d_viv[k * jmt * imt + (j - 1) * imt + i + 1]) /
				(d_viv[k * jmt * imt + j * imt + i] *
				 d_at0[j * imt + i] +
				 d_viv[k * jmt * imt + (j - 1) * imt + i] *
				 d_atn[j * imt + i] +
				 d_viv[k * jmt * imt + j * imt + i + 1] *
				 d_ate[j * imt + i] +
				 d_viv[k * jmt * imt + (j - 1) * imt + i + 1] *
				 d_atne[j * imt + i] + 1.0e-25);
		}
	}

	return;
}

__global__ void readyc1_2(double *d_wp12, double *d_wp13, double *d_vit, double *d_s2t,
		double *d_riv1, double *d_riv2,
		double *d_odzt, double *d_rit,
		double *d_rict, double *d_up, double *d_vp,
		double *d_s2u, double *d_viv, double *d_riu, double *d_ric) {

	int i, j, k;
	i = blockIdx.x * blockDim.x + threadIdx.x;
	j = blockIdx.y * blockDim.y + threadIdx.y;
	k = blockIdx.z * blockDim.z + threadIdx.z;

	if (k < kmm1 && j < jmt && i < imt) {
		double epsln = 1.0e-25;

		if (1 <= j && i < imt - 1) {
			d_riv1[k * jmt * imt + j * imt + i] =
				d_wp12[k * jmt * imt + j * imt + i] *
				d_vit[k * jmt * imt + j * imt + i] -
				d_wp12[(k + 1) * jmt * imt + j * imt + i] *
				d_vit[(k + 1) * jmt * imt + j * imt + i];

			d_riv2[k * jmt * imt + j * imt + i] =
				d_wp13[k * jmt * imt + j * imt + i] *
				d_vit[k * jmt * imt + j * imt + i] -
				d_wp13[(k + 1) * jmt * imt + j * imt + i] *
				d_vit[(k + 1) * jmt * imt + j * imt + i];

			d_s2t[k * jmt * imt + j * imt + i] =
				d_vit[(k + 1) * jmt * imt + j * imt + i] *
				(d_riv1[k * jmt * imt + j * imt + i] *
				 d_riv1[k * jmt * imt + j * imt + i] +
				 d_riv2[k * jmt * imt + j * imt + i] *
				 d_riv2[k * jmt * imt + j * imt + i]) *
				d_odzt[k + 1] * d_odzt[k + 1];

			d_rit[k * jmt * imt + j * imt + i] =
				d_vit[(k + 1) * jmt * imt + j * imt + i] *
				d_rict[k * jmt * imt + j * imt + i] /
				(d_s2t[k * jmt * imt + j * imt + i] + epsln);
		}

		d_riv1[k * jmt * imt + j * imt + i] =
			d_up[k * jmt * imt + j * imt + i] -
			d_up[(k + 1) * jmt * imt + j * imt + i];

		d_riv2[k * jmt * imt + j * imt + i] =
			d_vp[k * jmt * imt + j * imt + i] -
			d_vp[(k + 1) * jmt * imt + j * imt + i];

		        d_s2u[k * jmt * imt + j * imt + i] =
		                d_viv[(k + 1) * jmt * imt + j * imt + i] *
		                (d_riv1[k * jmt * imt + j * imt + i] *
		                 d_riv1[k * jmt * imt + j * imt + i] +
		                 d_riv2[k * jmt * imt + j * imt + i] *
		                 d_riv2[k * jmt * imt + j * imt + i]) *
		                d_odzt[k + 1] * d_odzt[k + 1];

		        d_riu[(k + 1) * jmt * imt + j * imt + i] =
		                d_viv[(k + 1) * jmt * imt + j * imt + i] *
		                d_ric[k * jmt * imt + j * imt + i] /
		                (d_s2u[k * jmt * imt + j * imt + i] + epsln);
	}

	return;
}

__global__ void readyc1_3(double *d_diff_back, double *d_diff_back_sh, double *d_diff_back_nh,
		double d_diff_back_coef_max, double d_diff_back_eq,
		double d_diff_back_coef, double *d_tlat) {

	int i, j;
	i = blockIdx.x * blockDim.x + threadIdx.x;
	j = blockIdx.y * blockDim.y + threadIdx.y;

	if (j < jmt && i < imt) {
		d_diff_back[j * imt + i] = 0.0e0;
		d_diff_back_sh[j * imt + i] = 0.0e0;
		d_diff_back_nh[j * imt + i] = 0.0e0;

		if (d_tlat[j * imt + i] < 0.0) {
			d_diff_back_sh[j * imt + i] =
				d_diff_back_coef_max *
				exp(-(0.4e0 * (d_tlat[j * imt + i] / DEGtoRAD + 28.9))
						* (0.4e0 * (d_tlat[j * imt + i] / DEGtoRAD + 28.9)));
		} else {
			d_diff_back_nh[j * imt + i] =
				d_diff_back_coef_max *
				exp(-(0.4e0 * (d_tlat[j * imt + i] / DEGtoRAD - 28.9))
						* (0.4e0 * (d_tlat[j * imt + i] / DEGtoRAD - 28.9)));
		}

		d_diff_back[j * imt + i] = d_diff_back_eq +
			d_diff_back_sh[j * imt + i] +
			d_diff_back_nh[j * imt + i];

		if (d_tlat[j * imt + i] < -10.0 * DEGtoRAD) {
			d_diff_back[j * imt + i] = d_diff_back[j * imt + i] +
				d_diff_back_coef;
		} else if (fabs(d_tlat[j * imt + i]) <= 10.0 * DEGtoRAD) {
			d_diff_back[j * imt + i] =
				d_diff_back[j * imt + i] +
				d_diff_back_coef *
				((fabs(d_tlat[j * imt + i]) / DEGtoRAD / 10.0) *
				 (fabs(d_tlat[j * imt + i]) / DEGtoRAD / 10.0));
		} else {
			d_diff_back[j * imt + i] =
				d_diff_back[j * imt + i] +
				d_diff_back_coef;
		}
	}

	return;
}

__global__ void readyc2(double *d_akmu, double *d_au0, double *d_akmt, double *d_aus, double *d_auw,
		double *d_ausw, double *d_viv) {

	int i, j, k;

	i = blockIdx.x * blockDim.x + threadIdx.x;
	j = blockIdx.y * blockDim.y + threadIdx.y;
	k = blockIdx.z * blockDim.z + threadIdx.z;

	if (i > 0 && j < jmt - 1 && k < kmm1) {
		d_akmu[k * jmt * imt + j * imt + i] =
			(d_au0[j * imt + i] *
			 d_akmt[k * jmt * imt + j * imt + i] +
			 d_aus[j * imt + i] *
			 d_akmt[k * jmt * imt + (j + 1) * imt + i] +
			 d_auw[j * imt + i] *
			 d_akmt[k * jmt * imt + j * imt + i - 1] +
			 d_ausw[j * imt + i] *
			 d_akmt[k * jmt * imt + (j + 1) * imt + i - 1]) *
			d_viv[(k + 1) * jmt * imt + j * imt + i];
	}

	return;
}

__global__ void readyc3(double *d_wka, double *d_ws, double *d_au0, double *d_ausw, double *d_aus, double *d_auw, double *d_u_wface, double *d_v_sface, double *d_u, double *d_v, double *d_hue, double *d_hun) {
	int i, j, k;

	i = (blockIdx.x) * blockDim.x + threadIdx.x;
	j = (blockIdx.y) * blockDim.y + threadIdx.y;
	k = blockIdx.z * blockDim.z + threadIdx.z;

	unsigned int tid = j * imt + i;
	unsigned int ttid = k * jmt * imt + j * imt + i;
	
	if (j < jmt - 1 && 1 <= i && i < imt && k < km) {
		d_wka[ttid] =
			d_au0[tid] * (d_ws[ttid]) +
			d_aus[tid] * (d_ws[ttid + imt]) +
			d_auw[tid] * (d_ws[ttid - 1]) +
			d_ausw[tid] * (d_ws[ttid + imt - 1]);
	}

	if (i > 0 && i < imt && j < jmt - 1 && k < km) {
		d_u_wface[ttid] = (d_u[ttid - 1] + d_u[ttid]) * p25 * d_hue[tid - 1];
		d_v_sface[ttid] = (d_v[ttid] + d_v[ttid + imt]) * p25 * d_hun[(j + 1) * imt + i];
	}
	return;
}

__global__ void hdiffu_del2_cu(double *d_hduk2, double *d_hdvk2, int d_jb, int d_je,
		int d_ib, int d_ie, double *d_wka, double d_c0f, double *d_up, double *d_vp,
		double *d_duc, double *d_dum, double *d_dun, double *d_dus, double *d_due, double *d_duw,
		double *d_dmc, double *d_dmn, double *d_dms, double *d_dme, double *d_dmw, double d_am,
		double *d_dlu, double *d_dlv, double *d_viv, int *d_kmu) {
	int i, j, k;

	i = (blockIdx.x) * blockDim.x + threadIdx.x;
	j = (blockIdx.y) * blockDim.y + threadIdx.y;
	k = blockIdx.z * blockDim.z + threadIdx.z;

	if (j < jmt && i < imt && k < km) {
		d_hduk2[k * jmt * imt + j * imt + i] = 0.0;
		d_hdvk2[k * jmt * imt + j * imt + i] = 0.0;

		if (j >= d_jb - 1 && j < d_je && i >= d_ib - 1 && i < d_ie) {
			d_hduk2[k * jmt * imt + j * imt + i] =
				d_am * (((d_duc[j * imt + i] +
								d_dum[j * imt + i]) * d_up[k * jmt * imt + j * imt + i] +
							d_dun[j * imt + i] * d_up[k * jmt * imt + (j - 1) * imt + i] +
							d_dus[j * imt + i] * d_up[k * jmt * imt + (j + 1) * imt + i] +
							d_due[j * imt + i] * d_up[k * jmt * imt + j * imt + i + 1] +
							d_duw[j * imt + i] * d_up[k * jmt * imt + j * imt + i - 1]) +
						(d_dmc[j * imt + i] * d_vp[k * jmt * imt + j * imt + i] +
						 d_dmn[j * imt + i] * d_vp[k * jmt * imt + (j - 1) * imt + i] +
						 d_dms[j * imt + i] * d_vp[k * jmt * imt + (j + 1) * imt + i] +
						 d_dme[j * imt + i] * d_vp[k * jmt * imt + j * imt + i + 1] +
						 d_dmw[j * imt + i] * d_vp[k * jmt * imt + j * imt + i - 1])) *
				d_viv[k * jmt * imt + j * imt + i];

			d_hdvk2[k * jmt * imt + j * imt + i] =
				d_am * (((d_duc[j * imt + i] +
								d_dum[j * imt + i]) * d_vp[k * jmt * imt + j * imt + i] +
							d_dun[j * imt + i] * d_vp[k * jmt * imt + (j - 1) * imt + i] +
							d_dus[j * imt + i] * d_vp[k * jmt * imt + (j + 1) * imt + i] +
							d_due[j * imt + i] * d_vp[k * jmt * imt + j * imt + i + 1] +
							d_duw[j * imt + i] * d_vp[k * jmt * imt + j * imt + i - 1]) -
						(d_dmc[j * imt + i] * d_up[k * jmt * imt + j * imt + i] +
						 d_dmn[j * imt + i] * d_up[k * jmt * imt + (j - 1) * imt + i] +
						 d_dms[j * imt + i] * d_up[k * jmt * imt + (j + 1) * imt + i] +
						 d_dme[j * imt + i] * d_up[k * jmt * imt + j * imt + i + 1] +
						 d_dmw[j * imt + i] * d_up[k * jmt * imt + j * imt + i - 1])) *
				d_viv[k * jmt * imt + j * imt + i];
		}

		if ((k > d_kmu[j * imt + i] - 1)) {
			d_hduk2[k * jmt * imt + j * imt + i] = 0.0;
			d_hdvk2[k * jmt * imt + j * imt + i] = 0.0;
		}

		if (1 <= j && j < jmt - 1 && 1 <= i && i < imt - 1)  {//jjrbug
			d_dlv[k * jmt * imt + j * imt + i] =
				d_dlv[k * jmt * imt + j * imt + i] +
				d_hdvk2[k * jmt * imt + j * imt + i];

			d_dlu[k * jmt * imt + j * imt + i] =
				d_dlu[k * jmt * imt + j * imt + i] +
				d_hduk2[k * jmt * imt + j * imt + i];
		}

		d_wka[k * jmt * imt + j * imt + i] =
			d_c0f * sqrt(d_up[k * jmt * imt + j * imt + i] *
					d_up[k * jmt * imt + j * imt + i] +
					d_vp[k * jmt * imt + j * imt + i] *
					d_vp[k * jmt * imt + j * imt + i]);
	}

	return;
}

__global__ void hdel4_zcurl(double *d_dxu, double *d_dyu,
                double *d_up,  double *d_vp, double *d_wka,int *d_kmt) {
        int i, j, k;
        i = blockIdx.x * blockDim.x + threadIdx.x;
        j = blockIdx.y * blockDim.y + threadIdx.y;
        k = blockIdx.z * blockDim.z + threadIdx.z;

        unsigned int tid = j * imt + i;

        if (i < imt && j < jmt && k <km) {
             d_wka[k * jmt * imt + tid] =0.0e0;
        if ( i >= 1 &&  j >= 1) {
                     if(k<=d_kmt[tid] - 1) {
                        d_wka[k * jmt * imt + tid] =
                                0.5 *
                                ((d_vp[k * jmt * imt + tid] * d_dyu[tid] +
                                  d_vp[k * jmt * imt + tid - imt] * d_dyu[tid - imt] - 
                                  d_vp[k * jmt * imt + tid - 1] * d_dyu [tid - 1] -
                                  d_vp[k * jmt * imt + tid - imt - 1] * d_dyu [tid - imt - 1] -
                                  d_up[k * jmt * imt + tid] * d_dxu[tid] -
                                  d_up[k * jmt * imt + tid - 1] * d_dxu[tid - 1 ] +
                                  d_up[k * jmt * imt + tid -imt ] * d_dxu[tid - imt] +
                                  d_up[k * jmt * imt + tid -imt -1 ] * d_dxu[tid - imt - 1]));
                }
      }
  }
 return;
}

__global__ void hdel4_div (double *d_hts, double *d_htw,
                double *d_tarea_r, double *d_up,  double *d_vp, double *d_wka,int *d_kmt) {
        int i, j, k;
        i = blockIdx.x * blockDim.x + threadIdx.x;
        j = blockIdx.y * blockDim.y + threadIdx.y;
        k = blockIdx.z * blockDim.z + threadIdx.z;

        unsigned int tid = j * imt + i;

        if (i < imt && j < jmt && k <km) {
             d_wka[k * jmt * imt + tid] =0.0e0;
        if (i < imt - 1 && j >= 1) {
                     if(k<=d_kmt[tid] - 1) {
                        d_wka[k * jmt * imt + tid] =
                                0.5 *
                                ((d_up[k * jmt * imt + tid + 1] +
                                  d_up[k * jmt * imt + tid - imt + 1]) * d_htw[tid + 1] -
                                 (d_up[k * jmt * imt + tid] +
                                  d_up[k * jmt * imt + tid - imt]) * d_htw[tid] +
                                 (d_vp[k * jmt * imt + tid + 1] +
                                  d_vp[k * jmt * imt + tid]) * d_hts[tid] -
                                 (d_vp[k * jmt * imt + tid - imt + 1] +
                                  d_vp[k * jmt * imt + tid - imt]) * d_hts[tid - imt]) *
                                d_tarea_r[tid];
                }
      }
  }
 return;
}




__global__ void hdel4_grad(double *d_wka,double *d_dxur, double *d_dyur, int *d_kmu,double *d_am_factor){ //d_pp=d_wka

    int i, j, k;
    double gradx,grady;

    i = blockIdx.x * blockDim.x + threadIdx.x;
    j = blockIdx.y * blockDim.y + threadIdx.y;
    k = blockIdx.z * blockDim.z + threadIdx.z;

    if (k < km && j < jmt && i < imt) {


        gradx = 0.0;//c0
        grady= 0.0;//c0

        if ((1 <= i) && (j < (jmt - 1))) {
            if (k <= d_kmu[j * imt + i] - 1) {
                gradx =
                        d_dxur[+j * imt + i] * p5 *
                        (d_wka[k * jmt * imt + (j + 1) * imt + i] -
                         d_wka[k * jmt * imt + j * imt + i - 1] -
                         d_wka[k * jmt * imt + (j + 1) * imt + i - 1] +
                         d_wka[k * jmt * imt + j * imt + i]);

                grady=
                        d_dyur[+j * imt + i] * p5 *
                        (d_wka[k * jmt * imt + (j + 1) * imt + i] -
                         d_wka[k * jmt * imt + j * imt + i - 1] +
                         d_wka[k * jmt * imt + (j + 1) * imt + i - 1] -
                         d_wka[k * jmt * imt + j * imt + i]);
            }
        }
         d_am_factor[k * jmt * imt + j * imt + i]+=gradx*gradx+grady*grady;
    }

    return;
}


__global__ void hdiffu_del4_0(double *d_hduk2, double *d_hdvk2, int d_jb, int d_je,
		int d_ib, int d_ie, double *d_uarea_r, double d_c0f, double *d_up, double *d_vp,
		double *d_duc, double *d_dum, double *d_dun, double *d_dus, double *d_due, double *d_duw,
		double *d_dmc, double *d_dmn, double *d_dms, double *d_dme, double *d_dmw, double d_am,
		double *d_dlu, double *d_dlv, double *d_viv, int *d_kmu,
		double *d_amf, double *d_am_factor, double *d_d2uk, double *d_d2vk) {
	int i, j, k;

	i = (blockIdx.x) * blockDim.x + threadIdx.x;
	j = (blockIdx.y) * blockDim.y + threadIdx.y;
	k = blockIdx.z * blockDim.z + threadIdx.z;

	if (j < jmt && i < imt && k < km) {
             double am_f,dxdy,t_area;
              am_f=1.0e0;
              
		d_hduk2[k * jmt * imt + j * imt + i] = 0.0;
		d_hdvk2[k * jmt * imt + j * imt + i] = 0.0;

		if (j >= d_jb - 2 && j <= d_je && i >= d_ib - 2 && i <= d_ie) {
                          if ((k <= d_kmu[j * imt + i] - 1)) {
                          t_area=sqrt(1.0e0/d_uarea_r[j * imt + i]);
                          dxdy=t_area*t_area*t_area*t_area*t_area*45.0e0;
                           am_f= sqrt (d_am_factor[k * jmt * imt + j * imt + i] )* dxdy /abs(d_am*d_amf[j * imt + i]); 
                          am_f=fmin(40.0e0,am_f);
                          am_f=fmax(1.0e0,am_f);
                     
			d_hduk2[k * jmt * imt + j * imt + i] = (((d_duc[j * imt + i] +
								d_dum[j * imt + i]) * d_up[k * jmt * imt + j * imt + i] +
							d_dun[j * imt + i] * d_up[k * jmt * imt + (j - 1) * imt + i] +
							d_dus[j * imt + i] * d_up[k * jmt * imt + (j + 1) * imt + i] +
							d_due[j * imt + i] * d_up[k * jmt * imt + j * imt + i + 1] +
							d_duw[j * imt + i] * d_up[k * jmt * imt + j * imt + i - 1]) +
						(d_dmc[j * imt + i] * d_vp[k * jmt * imt + j * imt + i] +
						 d_dmn[j * imt + i] * d_vp[k * jmt * imt + (j - 1) * imt + i] +
						 d_dms[j * imt + i] * d_vp[k * jmt * imt + (j + 1) * imt + i] +
						 d_dme[j * imt + i] * d_vp[k * jmt * imt + j * imt + i + 1] +
						 d_dmw[j * imt + i] * d_vp[k * jmt * imt + j * imt + i - 1]));

			d_hdvk2[k * jmt * imt + j * imt + i] = (((d_duc[j * imt + i] +
								d_dum[j * imt + i]) * d_vp[k * jmt * imt + j * imt + i] +
							d_dun[j * imt + i] * d_vp[k * jmt * imt + (j - 1) * imt + i] +
							d_dus[j * imt + i] * d_vp[k * jmt * imt + (j + 1) * imt + i] +
							d_due[j * imt + i] * d_vp[k * jmt * imt + j * imt + i + 1] +
							d_duw[j * imt + i] * d_vp[k * jmt * imt + j * imt + i - 1]) -
						(d_dmc[j * imt + i] * d_up[k * jmt * imt + j * imt + i] +
						 d_dmn[j * imt + i] * d_up[k * jmt * imt + (j - 1) * imt + i] +
						 d_dms[j * imt + i] * d_up[k * jmt * imt + (j + 1) * imt + i] +
						 d_dme[j * imt + i] * d_up[k * jmt * imt + j * imt + i + 1] +
						 d_dmw[j * imt + i] * d_up[k * jmt * imt + j * imt + i - 1]));
		}
        }
/*
		if ((k > d_kmu[j * imt + i] - 1)) {
			d_hduk2[k * jmt * imt + j * imt + i] = 0.0;
			d_hdvk2[k * jmt * imt + j * imt + i] = 0.0;
		}
  */
		d_d2uk[k * jmt * imt + j * imt + i] = am_f * d_amf[j * imt + i] * d_hduk2[k * jmt * imt + j * imt + i];
		d_d2vk[k * jmt * imt + j * imt + i] = am_f * d_amf[j * imt + i] * d_hdvk2[k * jmt * imt + j * imt + i];

  }
	return;
}

__global__ void hdiffu_del4_1(double *d_hduk2, double *d_hdvk2, int d_jb, int d_je,
		int d_ib, int d_ie, double *d_wka, double d_c0f, double *d_up, double *d_vp,
		double *d_duc, double *d_dum, double *d_dun, double *d_dus, double *d_due, double *d_duw,
		double *d_dmc, double *d_dmn, double *d_dms, double *d_dme, double *d_dmw, double d_am,
		double *d_dlu, double *d_dlv, double *d_viv, int *d_kmu,
		double *d_amf, double *d_am_factor, double *d_d2uk, double *d_d2vk) {
	int i, j, k;

	i = (blockIdx.x) * blockDim.x + threadIdx.x;
	j = (blockIdx.y) * blockDim.y + threadIdx.y;
	k = blockIdx.z * blockDim.z + threadIdx.z;

	if (j < jmt && i < imt && k < km) {
		d_hduk2[k * jmt * imt + j * imt + i] = 0.0;
		d_hdvk2[k * jmt * imt + j * imt + i] = 0.0;

		if (j >= d_jb - 1 && j < d_je && i >= d_ib - 1 && i < d_ie) {
			d_hduk2[k * jmt * imt + j * imt + i] =
				d_am * (((d_duc[j * imt + i] +
								d_dum[j * imt + i]) * d_d2uk[k * jmt * imt + j * imt + i] +
							d_dun[j * imt + i] * d_d2uk[k * jmt * imt + (j - 1) * imt + i] +
							d_dus[j * imt + i] * d_d2uk[k * jmt * imt + (j + 1) * imt + i] +
							d_due[j * imt + i] * d_d2uk[k * jmt * imt + j * imt + i + 1] +
							d_duw[j * imt + i] * d_d2uk[k * jmt * imt + j * imt + i - 1]) +
						(d_dmc[j * imt + i] * d_d2vk[k * jmt * imt + j * imt + i] +
						 d_dmn[j * imt + i] * d_d2vk[k * jmt * imt + (j - 1) * imt + i] +
						 d_dms[j * imt + i] * d_d2vk[k * jmt * imt + (j + 1) * imt + i] +
						 d_dme[j * imt + i] * d_d2vk[k * jmt * imt + j * imt + i + 1] +
						 d_dmw[j * imt + i] * d_d2vk[k * jmt * imt + j * imt + i - 1]));

			d_hdvk2[k * jmt * imt + j * imt + i] =
				d_am * (((d_duc[j * imt + i] +
								d_dum[j * imt + i]) * d_d2vk[k * jmt * imt + j * imt + i] +
							d_dun[j * imt + i] * d_d2vk[k * jmt * imt + (j - 1) * imt + i] +
							d_dus[j * imt + i] * d_d2vk[k * jmt * imt + (j + 1) * imt + i] +
							d_due[j * imt + i] * d_d2vk[k * jmt * imt + j * imt + i + 1] +
							d_duw[j * imt + i] * d_d2vk[k * jmt * imt + j * imt + i - 1]) -
						(d_dmc[j * imt + i] * d_d2uk[k * jmt * imt + j * imt + i] +
						 d_dmn[j * imt + i] * d_d2uk[k * jmt * imt + (j - 1) * imt + i] +
						 d_dms[j * imt + i] * d_d2uk[k * jmt * imt + (j + 1) * imt + i] +
						 d_dme[j * imt + i] * d_d2uk[k * jmt * imt + j * imt + i + 1] +
						 d_dmw[j * imt + i] * d_d2uk[k * jmt * imt + j * imt + i - 1]));
		}

		if ((k > d_kmu[j * imt + i] - 1)) {
			d_hduk2[k * jmt * imt + j * imt + i] = 0.0;
			d_hdvk2[k * jmt * imt + j * imt + i] = 0.0;
		}

		if (1 <= j && j < jmt - 1 && 1 <= i && i < imt - 1) { //jjr bug
			d_dlv[k * jmt * imt + j * imt + i] =
				d_dlv[k * jmt * imt + j * imt + i] +
				d_hdvk2[k * jmt * imt + j * imt + i];

			d_dlu[k * jmt * imt + j * imt + i] =
				d_dlu[k * jmt * imt + j * imt + i] +
				d_hduk2[k * jmt * imt + j * imt + i];
		}

		d_wka[k * jmt * imt + j * imt + i] =
			d_c0f * sqrt(d_up[k * jmt * imt + j * imt + i] *
					d_up[k * jmt * imt + j * imt + i] +
					d_vp[k * jmt * imt + j * imt + i] *
					d_vp[k * jmt * imt + j * imt + i]);
	}

	return;
}

__global__ void hdiffu_del4_2(double *d_hduk2, double *d_hdvk2, int d_jb, int d_je,
		int d_ib, int d_ie, double *d_wka, double d_c0f, double *d_up, double *d_vp,
		double *d_duc, double *d_dum, double *d_dun, double *d_dus, double *d_due, double *d_duw,
		double *d_dmc, double *d_dmn, double *d_dms, double *d_dme, double *d_dmw, double d_am,
		double *d_dlu, double *d_dlv, double *d_viv, int *d_kmu) {
	int i, j, k;

	i = (blockIdx.x) * blockDim.x + threadIdx.x;
	j = (blockIdx.y) * blockDim.y + threadIdx.y;
	k = blockIdx.z * blockDim.z + threadIdx.z;

	if (j < jmt && i < imt && k < km) {
		d_wka[k * jmt * imt + j * imt + i] =
			d_c0f * sqrt(d_up[k * jmt * imt + j * imt + i] *
					d_up[k * jmt * imt + j * imt + i] +
					d_vp[k * jmt * imt + j * imt + i] *
					d_vp[k * jmt * imt + j * imt + i]);
	}

	return;
}
__global__ void upwell_cu_0(double *d_au0, double *d_aus, double *d_auw, double *d_ausw,
		double *d_h0, double *d_ohbu,  double *d_u, double *d_uk, double *d_v,
		double *d_vk, double *d_work) {
	int i, j, k;
	i = blockIdx.x * blockDim.x + threadIdx.x;
	j = blockIdx.y * blockDim.y + threadIdx.y;

	unsigned int tid = j * imt + i;
	if (i < imt && j < jmt) {

		if (j < jmt - 1 && i >= 1) {
			d_work[tid] =
				d_au0[tid] * d_h0[tid] +
				d_aus[tid] * d_h0[(j + 1) * imt + i] +
				d_auw[tid] * d_h0[tid - 1] +
				d_ausw[tid] * d_h0[(j + 1) * imt + i - 1];

			for (k = 0; k < km; k++) {
				d_uk[k * jmt * imt + tid] = (1.0e0 + d_work[tid] * d_ohbu[tid]) *
					d_u[k * jmt * imt + tid];
				d_vk[k * jmt * imt + tid] = (1.0e0 + d_work[tid] * d_ohbu[tid]) *
					d_v[k * jmt * imt + tid];
			}
		}
	}
          return; 
}

__global__ void upwell_cu(double *d_dzp, double *d_h0, double *d_hts, double *d_htw, double *d_ohbt,
		double *d_tarea_r, double *d_uk, double *d_vit, double *d_vk, double *d_wka, double *d_work, double *d_ws, int *d_kmt) {
	int i, j, k;
	i = blockIdx.x * blockDim.x + threadIdx.x;
	j = blockIdx.y * blockDim.y + threadIdx.y;

	unsigned int tid = j * imt + i;

	d_work[tid] = 0.0e0;
	if (i < imt && j < jmt) {
	if (i < imt - 1 && j >= 1) {
		for (k = 0; k < km; k++) {
                     d_wka[k * jmt * imt + tid] =0.0e0;
                     if(k<=d_kmt[tid] - 1) {
			d_wka[k * jmt * imt + tid] =
				0.5 *
				((d_uk[k * jmt * imt + tid + 1] +
				  d_uk[k * jmt * imt + tid - imt + 1]) * d_htw[tid + 1] -
				 (d_uk[k * jmt * imt + tid] +
				  d_uk[k * jmt * imt + tid - imt]) * d_htw[tid] +
				 (d_vk[k * jmt * imt + tid + 1] +
				  d_vk[k * jmt * imt + tid]) * d_hts[tid] -
				 (d_vk[k * jmt * imt + tid - imt + 1] +
				  d_vk[k * jmt * imt + tid - imt]) * d_hts[tid - imt]) *
				d_tarea_r[tid];
		}
	}
      }
      d_work[tid] = 0.0e0;
     }
//	if (j == 0 || j >= jmt - 1 || i == 0 || i >= imt - 1) return;
      if (j>0 && j < jmt - 1 && 1 <= i && i<imt - 1) {
		for (k = 0; k < km; k++) {
			d_work[tid] -= d_dzp[k] *
				d_wka[k * jmt * imt + tid] *
				d_vit[k * jmt * imt + tid];
		}

	for (k = 1; k < km; k++) {
		d_ws[k * jmt * imt + tid] =
			d_vit[k * jmt * imt + tid] *
			(d_ws[(k - 1) * jmt * imt + tid] +
			 d_dzp[k - 1] *
			 (d_work[tid] * d_ohbt[tid] +
			  d_wka[(k - 1) * jmt * imt + tid]));
	}

	d_work[tid] = 1.0e0 / (1.0e0 + d_h0[tid] * d_ohbt[tid]);

	for (k = 1; k < km; k++) {
		d_ws[k * jmt * imt + tid] *= d_work[tid];
	}


}
           return;
}
__global__ void advection_momentum_cu(double *d_dlu, double *d_dlv, 
		double *d_odzp, double *d_u, double *d_uarea_r, double *d_u_wface, double *d_v,
		double *d_v_sface, double *d_wka) {
	int i, j, k;
	i = blockIdx.x * blockDim.x + threadIdx.x;
	j = blockIdx.y * blockDim.y + threadIdx.y;
	k = blockIdx.z * blockDim.z + threadIdx.z;

	if (i > 0 && i < imt - 1 && j > 0 && j < jmt - 1 && k < km) { //jjrbug
		double adv_z1 = 0.0e0, adv_z2 = 0.0e0, adv_z3 = 0.0e0, adv_z4 = 0.0e0;
		unsigned int tid = j * imt + i;
		unsigned int ttid = k * jmt * imt + j * imt + i;

		d_dlu[ttid] = (-d_u_wface[ttid] * (d_u[ttid] - d_u[ttid - 1]) -
				d_u_wface[ttid + 1] * (d_u[ttid + 1] - d_u[ttid]) -
				d_v_sface[ttid] * (d_u[k * jmt * imt + (j + 1) * imt + i] - d_u[ttid]) -
				d_v_sface[ttid - imt] * (d_u[ttid] - d_u[ttid - imt])) * d_uarea_r[tid];
		d_dlv[ttid] = (-d_u_wface[ttid] * (d_v[ttid] - d_v[ttid - 1]) -
				d_u_wface[ttid + 1] * (d_v[ttid + 1] - d_v[ttid]) -
				d_v_sface[ttid] * (d_v[k * jmt * imt + (j + 1) * imt + i] - d_v[ttid]) -
				d_v_sface[ttid - imt] * (d_v[ttid] - d_v[ttid - imt])) * d_uarea_r[tid];
		if (k > 0) {
			adv_z1 = d_wka[ttid] * (d_u[(k - 1) * jmt * imt + tid] - d_u[ttid]);
			adv_z3 = d_wka[ttid] * (d_v[(k - 1) * jmt * imt + tid] - d_v[ttid]);
		}

		if (k < km - 1) {
			adv_z2 = d_wka[jmt * imt + ttid] * (d_u[ttid] - d_u[jmt * imt + ttid]);
			adv_z4 = d_wka[jmt * imt + ttid] * (d_v[ttid] - d_v[jmt * imt + ttid]);
		}
		d_dlu[ttid] += -p5 * d_odzp[k] * (adv_z1 + adv_z2);
		d_dlv[ttid] += -p5 * d_odzp[k] * (adv_z3 + adv_z4);
	}

	return;
}

__global__ void readyc4(double *d_sbcx, double *d_sbcy, double *d_su, double *d_sv, double *d_bbcx,
		double d_od0, double d_c0f, double *d_up, double *d_vp, double *d_snlat,
		double d_cag, double d_sag, double *d_dlub, double *d_dlvb, double *d_dzp,
		double *d_ohbu, double *d_dlu, double *d_viv, double *d_dlv, int *d_kmu, double *d_bbcy) {
	int i, j, k, kmb;
	i = (blockIdx.x) * blockDim.x + threadIdx.x;
	j = (blockIdx.y) * blockDim.y + threadIdx.y;

	if (i < imt && j < jmt) {
		d_dlub[j * imt + i] = 0.0;
		d_dlvb[j * imt + i] = 0.0;

		for (k = 0; k < km; k++) {
			d_dlub[j * imt + i] += d_dzp[k] *
				d_ohbu[j * imt + i] *
				d_dlu[k * jmt * imt + j * imt + i] *
				d_viv[k * jmt * imt + j * imt + i];

			d_dlvb[j * imt + i] += d_dzp[k] *
				d_ohbu[j * imt + i] *
				d_dlv[k * jmt * imt + j * imt + i] *
				d_viv[k * jmt * imt + j * imt + i];
		}

		if (1 <= j && j < jmt - 1 && 1 <= i && i < imt - 1) {
			kmb = d_kmu[j * imt + i] - 1;

			if (kmb >= 0) {
				d_sbcx[j * imt + i] = d_su[j * imt + i] * d_od0;

				d_sbcy[j * imt + i] = d_sv[j * imt + i] * d_od0;

				d_bbcx[j * imt + i] =
					d_c0f * sqrt(d_up[kmb * jmt * imt + j * imt + i] *
							d_up[kmb * jmt * imt + j * imt + i] +
							d_vp[kmb * jmt * imt + j * imt + i] *
							d_vp[kmb * jmt * imt + j * imt + i]) *
					(d_up[kmb * jmt * imt + j * imt + i] * d_cag +
					 d_snlat[j * imt + i] *
					 d_vp[kmb * jmt * imt + j * imt + i] * d_sag);

				d_bbcy[j * imt + i] =
					d_c0f * sqrt(d_up[kmb * jmt * imt + j * imt + i] *
							d_up[kmb * jmt * imt + j * imt + i] +
							d_vp[kmb * jmt * imt + j * imt + i] *
							d_vp[kmb * jmt * imt + j * imt + i]) *
					(-d_snlat[j * imt + i] *
					 d_up[kmb * jmt * imt + j * imt + i] * d_sag +
					 d_vp[kmb * jmt * imt + j * imt + i] * d_cag);
			} else {
				d_sbcx[j * imt + i] = 0.0e0;
				d_sbcy[j * imt + i] = 0.0e0;
				d_bbcx[j * imt + i] = 0.0e0;
				d_bbcy[j * imt + i] = 0.0e0;
			}

			d_dlub[j * imt + i] = d_dlub[j * imt + i] +
				(d_sbcx[j * imt + i] -
				 d_bbcx[j * imt + i]) *
				d_ohbu[j * imt + i];

			d_dlvb[j * imt + i] = d_dlvb[j * imt + i] +
				(d_sbcy[j * imt + i] -
				 d_bbcy[j * imt + i]) *
				d_ohbu[j * imt + i];
		}
	}

	return;
}


__global__ void readyc4_2(double *d_su, double *d_sv, double d_od0, double *d_up, double *d_vp,
		double d_sag, double d_cag, double *d_akmu, double *d_odzt, double *d_viv,
		double *d_wka, double *d_snlat, double *d_dlv, double *d_odzp, double *d_dlu) {
	int i, j, k;
	i = (blockIdx.x) * blockDim.x + threadIdx.x;
	j = (blockIdx.y) * blockDim.y + threadIdx.y;
	k = (blockIdx.z) * blockDim.z + threadIdx.z;

	if (k < km && 1 <= j && j < jmt - 1 && 1 <= i && i < imt - 1) {
		double diff_u1, diff_u2, diff_v1, diff_v2;
		double aidif = 0.5;

		if (k == 0) {
			diff_v1 = d_sv[j * imt + i] * d_od0 * (1 - aidif);
			diff_u1 = d_su[j * imt + i] * d_od0 * (1 - aidif);
		} else {
			diff_v1 = d_akmu[(k - 1) * jmt * imt + j * imt + i] *
				(1 - aidif) *
				(d_vp[(k - 1) * jmt * imt + j * imt + i] -
				 d_vp[k * jmt * imt + j * imt + i]) *
				d_odzt[k] *
				d_viv[k * jmt * imt + j * imt + i] +
				(1.0e0 - d_viv[k * jmt * imt + j * imt + i]) *
				d_wka[(k - 1) * jmt * imt + j * imt + i] *
				(1 - aidif) *
				(-d_snlat[j * imt + i] *
				 d_up[(k - 1) * jmt * imt + j * imt + i] * d_sag +
				 d_vp[(k - 1) * jmt * imt + j * imt + i] * d_cag);

			diff_u1 = d_akmu[(k - 1) * jmt * imt + j * imt + i] *
				(1 - aidif) *
				(d_up[(k - 1) * jmt * imt + j * imt + i] -
				 d_up[k * jmt * imt + j * imt + i]) * d_odzt[k] *
				d_viv[k * jmt * imt + j * imt + i] +
				(1.0e0 - d_viv[k * jmt * imt + j * imt + i]) *
				d_wka[(k - 1) * jmt * imt + j * imt + i] * (1 - aidif) *
				(d_up[(k - 1) * jmt * imt + j * imt + i] * d_cag +
				 d_snlat[j * imt + i] *
				 d_vp[(k - 1) * jmt * imt + j * imt + i] * d_sag);
		}

		if (k == km - 1) {
			diff_v2 = d_wka[(km - 1) * jmt * imt + j * imt + i] *
				(-d_snlat[j * imt + i] *
				 d_up[(km - 1) * jmt * imt + j * imt + i] * d_sag +
				 d_vp[(km - 1) * jmt * imt + j * imt + i] * d_cag) *
				(1 - aidif);

			diff_u2 = d_wka[(km - 1) * jmt * imt + j * imt + i] *
				(d_up[(km - 1) * jmt * imt + j * imt + i] * d_cag +
				 d_snlat[j * imt + i] *
				 d_vp[(km - 1) * jmt * imt + j * imt + i] * d_sag) *
				(1 - aidif);
		} else {
			diff_v2 = d_akmu[k * jmt * imt + j * imt + i] * (1 - aidif) *
				(d_vp[k * jmt * imt + j * imt + i] -
				 d_vp[(k + 1) * jmt * imt + j * imt + i]) *
				d_odzt[k + 1] * d_viv[(k + 1) * jmt * imt + j * imt + i] +
				(1.0e0 - d_viv[(k + 1) * jmt * imt + j * imt + i]) *
				d_wka[k * jmt * imt + j * imt + i] *
				(1 - aidif) *
				(-d_snlat[j * imt + i] *
				 d_up[k * jmt * imt + j * imt + i] * d_sag +
				 d_vp[k * jmt * imt + j * imt + i] * d_cag);

			diff_u2 = d_akmu[k * jmt * imt + j * imt + i] * (1 - aidif) *
				(d_up[k * jmt * imt + j * imt + i] -
				 d_up[(k + 1) * jmt * imt + j * imt + i]) *
				d_odzt[k + 1] * d_viv[(k + 1) * jmt * imt + j * imt + i] +
				(1.0e0 - d_viv[(k + 1) * jmt * imt + j * imt + i]) *
				d_wka[k * jmt * imt + j * imt + i] *
				(1 - aidif) *
				(d_up[k * jmt * imt + j * imt + i] * d_cag +
				 d_snlat[j * imt + i] *
				 d_vp[k * jmt * imt + j * imt + i] * d_sag);
		}

		d_dlv[k * jmt * imt + j * imt + i] = d_dlv[k * jmt * imt + j * imt + i] +
			d_odzp[k] * (diff_v1 - diff_v2);

		d_dlu[k * jmt * imt + j * imt + i] = d_dlu[k * jmt * imt + j * imt + i] +
			d_odzp[k] * (diff_u1 - diff_u2);
	}

	return;
}

__inline__ __device__ double wavelat(double xf, double yn) { return xf * acosh(yn / xf); }

__device__ double eplatidepend_fun(double f, double an) {

	double an0, f_30;
	an0 = 5.24e-3;
	double anum, den;
	double omega_tmp;

	omega_tmp = 4.0e0*atan(1.0e0) / 43082.0e0;

	f_30 = omega_tmp;

	den = wavelat(f_30, an0);
	anum = wavelat(f, an);
	return anum / den;
}

__inline__ __device__ int sign_int(int x, int y) { return (y) >= 0 ? abs(x) : -abs(x); }

__inline__ __device__ double sign(double x, double y) { return (y) >= 0.0e0 ? fabs(x) : -fabs(x); }

__device__ void formld(int n, double *z, double *t, double *amld) {
	int k;
	double tm;
	for (k = 0; k < n; k++) {
		if (fabs(t[k] - t[0]) > 0.1) {
			tm = t[0] - sign(0.1e0, t[0] - t[k]);
			*amld = z[k] + (z[k - 1] - z[k]) * (tm - t[k]) / (t[k - 1] - t[k] + 1.0e-20); //jjr bug
			return;
		}
	}
	*amld = z[n - 1];
}

__device__ void interp1d_expabs(double x, double *x_1, double *slq2_1, double *sm_1, double *sh_1, double *ss_1,
		double *slq2, double *sm, double *sh, double *ss, double delta, double rat) {

	double deltax, deltaxta, dsh_x, dslq2_x, dsm_x, dss_x, tabindx;
	int lx0, lx1;

	if (x > x_1[2 * mt]) {
		x = x_1[2 * mt];
	} else if (x < x_1[0]) {
		x = x_1[0];
	}

	if (fabs(x) >= x_1[mt + mt0] && (fabs(x) >= x_1[2 * mt])) {
		lx0 = (int) sign((double) mt, x);
		lx1 = lx0;
	} else {
		if (fabs(x) < x_1[mt + mt0]) {
			lx1 = (int) (x / delta) + (int) sign((double) 1, x);
		} else {
			tabindx = sign((double) (mt0) + ((log(fabs(x)) - log(x_1[mt + mt0])) / log(rat)), x);
			lx1 = (int) (tabindx) + (int) sign((double) 1, x);
		}

		if (fabs(x_1[mt + lx1]) < (fabs(x))) {
			lx1 += sign_int(1, lx1);
		} else if (fabs(x_1[mt + lx1 - sign_int(1, lx1)]) > fabs(x)) {
			lx1 += -sign_int(1, lx1);
		}
		lx0 = lx1 - (int) sign((double) 1, x);
		if (fabs(x) <= 0.0e0) { lx1 = 1; }
	}

	if ((x > 0.0e0 && (x < x_1[mt + lx0] || x > x_1[mt + lx1])) ||
			(x < 0.0e0 && (x > x_1[mt + lx0] || x < x_1[mt + lx1]))) {
		//printf("exit_122");
		return;
	}

	deltaxta = x_1[mt + lx1] - x_1[mt + lx0];
	deltax = x - x_1[mt + lx0];

	if (lx1 == lx0) {
		dslq2_x = 0.0e0;
		dsm_x = 0.0e0;
		dsh_x = 0.0e0;
		dss_x = 0.0e0;
	} else {
		dslq2_x = (slq2_1[mt + lx1] - slq2_1[mt + lx0]) / deltaxta;
		dsm_x = (sm_1[mt + lx1] - sm_1[mt + lx0]) / deltaxta;
		dsh_x = (sh_1[mt + lx1] - sh_1[mt + lx0]) / deltaxta;
		dss_x = (ss_1[mt + lx1] - ss_1[mt + lx0]) / deltaxta;
	}
	*slq2 = slq2_1[mt + lx0] + dslq2_x * deltax;
	*sm = sm_1[mt + lx0] + dsm_x * deltax;
	*sh = sh_1[mt + lx0] + dsh_x * deltax;
	*ss = ss_1[mt + lx0] + dss_x * deltax;

	return;
}

__device__ void interp2d_expabs(double *ri, double *rid, double *slq2, double *sm, double *sh, double *ss,
		double *d_rib, double *d_shb, double *d_slq2b, double *d_smb, double *d_ssb,
		int *d_irimax, double d_dri, double d_rri) {
	double deltari, deltarid, deltaridta, deltarita, dsh_ri, dsh_rid, dslq2_rid, dsm_ri, dsm_rid, dss_ri, dss_rid, tabindrid, dslq2_ri, tabindri;
	int lri0, lri1, lrid0, lrid1;

	if (*ri > d_rib[2 * mt]) {
		if (fabs(*rid) <= *ri) {
			*rid = d_rib[2 * mt] * (*rid / *ri);
			*ri = d_rib[2 * mt];
		} else if (*rid > *ri) {
			*ri = d_rib[2 * mt] * (*ri / *rid);
			*rid = d_rib[2 * mt];
		} else if (*rid < -*ri) {
			*ri = d_rib[0] * (*ri / *rid);
			*rid = d_rib[0];
		}
	} else if (*ri < d_rib[0]) {
		if (fabs(*rid) <= -*ri) {
			*rid = d_rib[0] * (*rid / *ri);
			*ri = d_rib[0];
		} else if (*rid > -*ri) {
			*ri = d_rib[2 * mt] * (*ri / *rid);
			*rid = d_rib[2 * mt];
		} else if (*rid < *ri) {
			*ri = d_rib[0] * (*ri / *rid);
			*rid = d_rib[0];
		}
	} else if (*rid > d_rib[2 * mt]) {
		*ri = d_rib[2 * mt] * (*ri / *rid);
		*rid = d_rib[2 * mt];
	} else if (*rid < d_rib[0]) {
		*ri = d_rib[0] * (*ri / *rid);
		*rid = d_rib[0];
	}//y
	double ri_tmp=*ri;
	double rid_tmp=*rid;
	double tmp=d_rib[mt + mt0];
	double tmp1=d_rib[2 * mt];
	if (fabs(rid_tmp) >= tmp && (fabs(rid_tmp)) >= tmp1) {
		lrid0 = (int) sign((double) (mt), rid_tmp);
		lrid1 = lrid0;
	} else {
		if (fabs(rid_tmp) < tmp) {
			lrid1 = (int) (rid_tmp / d_dri) + (int) sign((double) 1, rid_tmp);
		} else {
			tabindrid = sign((double) (mt0) + ((log(fabs(rid_tmp)) - log(tmp)) / log(d_rri)), rid_tmp);
			lrid1 = (int) (tabindrid) + (int) sign((double) 1, rid_tmp);
		}
		//
		if ((fabs(d_rib[mt + lrid1])) < (fabs(rid_tmp))) {
			lrid1 += sign_int(1,lrid1);
		} else if (fabs(d_rib[mt+lrid1-sign_int(1, lrid1)]) > fabs(rid_tmp)) {
			lrid1 += -sign_int(1,lrid1);
		}
		lrid0 = lrid1 - (int) sign((double) 1, rid_tmp);
		if (rid_tmp == 0.0e0) lrid1 = 1;
	}//y
	//tmp=d_rib[mt+lrid0];
	//tmp1=d_rib[mt+lrid1];
	if ((rid_tmp > 0.0e0 && (rid_tmp < d_rib[mt+lrid0] || rid_tmp > d_rib[mt+lrid1] ))||
			(rid_tmp < 0.0e0 && (rid_tmp > d_rib[mt+lrid0] || rid_tmp < d_rib[mt+lrid1] ))) { //jjr bug lt
		//printf("exit_204");
		return;
	}//n
	if (ri_tmp > fmin(d_rib[mt + d_irimax[mt + lrid0]], d_rib[mt + d_irimax[mt + lrid1]])) {
		*slq2 = 0.0e0;
		*sm = 0.0e0;
		*sh = 0.0e0;
		*ss = 0.0e0;
		return;
	}//y

	if (fabs(ri_tmp) >= tmp && fabs(ri_tmp) >= tmp1) {
		lri0 = (int) sign((double) (mt), ri_tmp);
		lri1 = lri0;
	} else {//y
		if (fabs(ri_tmp) < tmp) {
			lri1 = (int) (ri_tmp / d_dri) + (int) sign((double) 1, ri_tmp);
		} else {
			tabindri = sign((double) (mt0) + ((log(fabs(ri_tmp)) - log(tmp)) / log(d_rri)), ri_tmp);
			lri1 = (int) (tabindri) + (int) sign((double) 1, ri_tmp);
		}//y
		if ((fabs(d_rib[mt + lri1])) < (fabs(ri_tmp))) {
			lri1 +=  sign_int(1, lri1);
		}//y
		else if ( fabs(d_rib[mt+lri1- sign_int(1,lri1)]) > fabs(ri_tmp)) {//y
			lri1 -=  sign_int(1, lri1);
		}//n
		lri0 = lri1 - (int) sign((double) 1, ri_tmp);
		if (ri_tmp == 0.0e0) lri1 = 1;
	}

	if ((ri_tmp > 0.0e0 && (ri_tmp < d_rib[mt + lri0] || ri_tmp > d_rib[mt + lri1])) ||
			(ri_tmp < 0.0e0 && (ri_tmp > d_rib[mt + lri0] || ri_tmp < d_rib[mt + lri1]))) {
		//printf("exit_236");
		return;
	}

	deltaridta = d_rib[mt + lrid1] - d_rib[mt + lrid0];
	deltarita = d_rib[mt + lri1] - d_rib[mt + lri0];
	deltarid = rid_tmp - d_rib[mt + lrid0];
	deltari = ri_tmp - d_rib[mt + lri0];

	if (lrid1 == lrid0) {
		dslq2_rid = 0.0e0;
		dslq2_ri = 0.0e0;
		dsm_rid = 0.0e0;
		dsm_ri = 0.0e0;
		dsh_rid = 0.e0;
		dsh_ri = 0.e0;
		dss_rid = 0.e0;
		dss_ri = 0.e0;
	} else {
		dslq2_rid =
			(d_slq2b[(2 * mt + 1) * (lrid1 + mt) + lri0 + mt] - d_slq2b[(2 * mt + 1) * (lrid0 + mt) + lri0 + mt]) /
			deltaridta;
		dslq2_ri = (d_slq2b[(2 * mt + 1) * (lrid0 + mt) + lri1 + mt] - d_slq2b[(2 * mt + 1) * (lrid0 + mt) + lri0 + mt]) /
			deltarita;
		dsm_rid = (d_smb[(2 * mt + 1) * (lrid1 + mt) + lri0 + mt] - d_smb[(2 * mt + 1) * (lrid0 + mt) + lri0 + mt]) /
			deltaridta;
		dsm_ri = (d_smb[(2 * mt + 1) * (lrid0 + mt) + lri1 + mt] - d_smb[(2 * mt + 1) * (lrid0 + mt) + lri0 + mt]) /
			deltarita;
		dsh_rid = (d_shb[(2 * mt + 1) * (lrid1 + mt) + lri0 + mt] - d_shb[(2 * mt + 1) * (lrid0 + mt) + lri0 + mt]) /
			deltaridta;
		dsh_ri = (d_shb[(2 * mt + 1) * (lrid0 + mt) + lri1 + mt] - d_shb[(2 * mt + 1) * (lrid0 + mt) + lri0 + mt]) /
			deltarita;
		dss_rid = (d_ssb[(2 * mt + 1) * (lrid1 + mt) + lri0 + mt] - d_ssb[(2 * mt + 1) * (lrid0 + mt) + lri0 + mt]) /
			deltaridta;
		dss_ri = (d_ssb[(2 * mt + 1) * (lrid0 + mt) + lri1 + mt] - d_ssb[(2 * mt + 1) * (lrid0 + mt) + lri0 + mt]) /
			deltarita;
	}

	*slq2 = d_slq2b[(2 * mt + 1) * (lrid0 + mt) + lri0 + mt] + dslq2_ri * deltari + dslq2_rid * deltarid;
	*sm = d_smb[(2 * mt + 1) * (lrid0 + mt) + lri0 + mt] + dsm_ri * deltari + dsm_rid * deltarid;
	*sh = d_shb[(2 * mt + 1) * (lrid0 + mt) + lri0 + mt] + dsh_ri * deltari + dsh_rid * deltarid;
	*ss = d_ssb[(2 * mt + 1) * (lrid0 + mt) + lri0 + mt] + dss_ri * deltari + dss_rid * deltarid;

	return;
}

__device__ void turb_device(double *z, double *t, double *rh, double *ri, double *rid, double *s2, double fricmx,
		double wndmix, double *an2, double buoytur, double buoysol, double coriol, double *amld,
		double *akm, double *akh, double *aks, int n, int na, int nmax, double *d_amtaun2a1,
		double *d_and2on2a1, double *d_back_ra_r, double *d_rib, double *d_sha1, double *d_shb,
		double *d_sh_r1, double *d_slq2b, double *d_slq2_r1, double *d_sma1, double *d_smb,
		double *d_sm_r1, double *d_ssa1, double *d_ssb, double *d_ss_r1, double *d_turb_param,
		int *d_irimax) {
	//int i = blockIdx.x * blockDim.x + threadIdx.x;
	//int j = blockIdx.y * blockDim.y + threadIdx.y;
	int k;
	double v_back, t_back, s_back;
	double visc_cbu_limit1, diff_cbt_limit1;
	bool lifupper, lifepsy;

	double al2, tmp, slq2, sm, ss, sh;
	double aldeep;
	int ilmldpoint = 0, ilmldpointneg = 0, ipoint = 0, icall = 0;

	double amtaun2;
	double deltheta_r_1, d_rri = d_turb_param[0], d_rnd2on2 = d_turb_param[1], d_dri = d_turb_param[2], d_deltheta_r = d_turb_param[3], d_b1 = d_turb_param[4], d_theta_rcrp = d_turb_param[5], d_theta_rcrn = d_turb_param[6], d_dand2on2 = d_turb_param[10];
	double anlq2, back_ra_r1, back_ri1 = 0.0e0, back_ric1, back_rid1, back_rit1;
	double dback_ra_r_o_dtheta, delback_ra_r, delsh_back, delslq2_back, delsm_back, delss_back, dsh_back_o_dtheta, dslq2_back_o_dtheta, dsm_back_o_dtheta, dss_back_o_dtheta, eplatidepend, epson2, epson2_;
	double epsy;
	int ifpureshear, ifrafglt, inegproblem = 0, iproblem = 0, itheta_r0, itheta_r1, jtheta_r0;
	double ra_r, ric, rit, sh_back, slq2_back, sm_back, ss_back, theta_r, theta_r0, theta_r1, tmp_back;

	int nb = n > 0 ? n : 0;

	visc_cbu_limit1 = fricmx;
	diff_cbt_limit1 = fricmx;

	double buoytot = buoytur + buoysol;
	if (n > 0) {
		formld(n, z, t, amld);
	} else if (n == 0) {
		*amld = z[0];
	} else {
		*amld = 0.0e0;
	}

	double al0 = 0.17 * *amld;
	double ifbelow = 0;
	double ifnofsmall, an;
	double ri1, rid1;
	for (k = 0; k < n; k++) {
		ifnofsmall = 0;
		if (an2[k] >= 0.0e0) an = sqrt(an2[k]);
		if (an < fabs(coriol)) {
			ifnofsmall = 1;
		}
		ri1 = ri[k];
		double and2, and2on2;
		rid1 = rid[k];

		and2 = rid[k] / (ri[k] + 1.0e-25) * an2[k];
		and2on2 = and2 / (an2[k] + 1.0e-25);

		if (an2[k] < d_rib[0] * s2[k]) {//(ifzeroshear)&&
			ifpureshear = 1;
		} else {
			ifpureshear = 0;
		}
		//double sm1,sh1,ss1,slq21;
		if (ifpureshear == 1) {
			interp1d_expabs(and2on2, d_and2on2a1, d_amtaun2a1, d_sma1, d_sha1, d_ssa1, &amtaun2, &sm, &sh, &ss,
					d_dand2on2, d_rnd2on2);//dand2on2=dri
			slq2 = (-amtaun2) / ((d_b1 * d_b1) * ri[k] + 1.0e-25);
		}else {
			//TODO the next line is need, but exit addressing conflict 
			interp2d_expabs(&ri1, &rid1, &slq2, &sm, &sh, &ss, d_rib, d_shb, d_slq2b, d_smb, d_ssb, d_irimax, d_dri,d_rri);
			if (slq2 < 0.0e0) {
				//printf("exit_245");
				return;
			}

			if (slq2 == 0.0e0) ifbelow = 1;
		}


		lifupper = ((ifbelow == 0 || ifnofsmall == 1 || ((ri1 < 0.0e0) && (k <= 1))) &&
				(slq2 > 0.0e0));//ifepson2<2||ifdeeplat=1&&
		if (lifupper) {
			++ilmldpoint;
			if (ri1 < 0.0e0) ++ilmldpointneg;
		}
		++ipoint;
		if (k == 0) ++icall;
		//n
		epsy = 0.0e0;
		lifepsy = false;
		if (n > 0) {//(isurfuse==1)
			if (slq2 == 0.0e0) {
				epsy = 0.0e0;
			} else {
				epsy = -buoytot / ((1.0e0 / ((d_b1 * d_b1) * (slq2))) - 0.5e0 * sm);
			}

			lifepsy = ((epsy >= 0.0e0) && lifupper);

			if ((epsy < 0.0e0) && lifupper) {
				++iproblem;
				if (ri1 < 0.0e0) ++inegproblem;
			}
		}
		//y
		double akz = 0.4 * z[k];
		double al = akz * al0 / (al0 + akz);
		al2 = al * al;

		if (ifbelow != 1 && !lifepsy) {//(ifepson2==2)&&
			if (ri1 > 0.0e0) {
				anlq2 = slq2 * ri1;
				if (anlq2 > 0.281e0) { //icondear=0
					al2 *= 0.281e0 / anlq2;
					slq2 = 0.281e0 / (ri1 + 1.0e-20);
				}
			}
		}
		if (an2[k] < 0.0e0) {
			epson2_ = epson2__;
		} else {
			if (ifnofsmall == 1) {
				eplatidepend = 0.0e0;
			} else {
				eplatidepend = eplatidepend_fun(fabs(coriol), an);
			}
			eplatidepend = fmax(eplatidepend, eplatidependmin);
			epson2_ = epson2__ * eplatidepend;
		}
		epson2 = epson2_;

		if (ri[k] > 0.0e0) {
			rit = (ri[k] + rid[k]) / 2.0e0;
			ric = (ri[k] - rid[k]) / 2.0e0;
			ra_r = sqrt((rit * rit) + (ric * ric));
			if (fabs(rit) <= 0.0e0) {
				if (fabs(ric) <= 0.0e0) {
					theta_r = atan(1.0e0) ;
				} else {
					theta_r = 2.0e0*atan(1.0e0) ;
				}
			} else {
				theta_r = atan(ric / rit);
			}
			if (fabs(theta_r)  > 2.0e0*atan(1.0e0) ) {
				//printf("exit_430");
				return;
			}
			if (theta_r  < -atan(1.0e0)) theta_r += 4.0e0*atan(1.0e0);

			jtheta_r0 = (int) ((theta_r + atan(1.0e0)) / d_deltheta_r);
			//jtheta_r1 = jtheta_r0+1;

			itheta_r0 = jtheta_r0 - n_theta_r_oct;
			itheta_r1 = itheta_r0 + 1;

			theta_r0 = itheta_r0 * d_deltheta_r;
			theta_r1 = itheta_r1 * d_deltheta_r;
			if ((theta_r0 <= d_theta_rcrp) && (theta_r > d_theta_rcrp)) {
				theta_r = theta_r1;
				theta_r0 = theta_r1;
				itheta_r0 = itheta_r1;
				++itheta_r1;
				theta_r1 += d_deltheta_r;
			} else if ((theta_r1 >= d_theta_rcrn) && (theta_r < d_theta_rcrn)) {
				theta_r = theta_r0;
				theta_r1 = theta_r0;
				itheta_r1 = itheta_r0;
				--itheta_r0;
				theta_r0 += -d_deltheta_r;
			}

			//theta_r_deg = theta_r * 180.0e0/pi_const;

			if ((itheta_r1 > 3 * n_theta_r_oct) || (itheta_r0 < -n_theta_r_oct)) {
				//printf("exit_459");
				return;
			}
			deltheta_r_1 = theta_r - theta_r0;

			delback_ra_r = d_back_ra_r[n_theta_r_oct + itheta_r1] - d_back_ra_r[n_theta_r_oct + itheta_r0];

			dback_ra_r_o_dtheta = delback_ra_r / d_deltheta_r;
			back_ra_r1 = d_back_ra_r[n_theta_r_oct + itheta_r0] + deltheta_r_1 * dback_ra_r_o_dtheta;

			ifrafglt = 0;
			if (ifrafgmax == 1) {
				if ((theta_r <= d_theta_rcrp) || (theta_r >= d_theta_rcrn)) {
					if (back_ra_r1 > ra_r) {
						ifrafglt = 1;
						back_ra_r1 = ra_r;
					}
				}
			}
			if (back_ra_r1 < 0.0e0) {
				//printf("exit_480");
				return;
			}

			back_rit1 = cos(theta_r) * back_ra_r1;
			back_ric1 = sin(theta_r) * back_ra_r1;
			back_ri1 = back_rit1 + back_ric1;
			back_rid1 = back_rit1 - back_ric1;


			if (ifrafglt == 1) {//ifbg_theta_interp=1
				interp2d_expabs(&back_ri1, &back_rid1, &slq2_back, &sm_back, &sh_back, &ss_back, d_rib, d_shb, d_slq2b,d_smb, d_ssb, d_irimax, d_dri, d_rri);
			} else {//ifbg_theta_interp
				deltheta_r_1 = theta_r - itheta_r0 * d_deltheta_r;
				delsm_back = d_sm_r1[n_theta_r_oct + itheta_r1] - d_sm_r1[n_theta_r_oct + itheta_r0];
				dsm_back_o_dtheta = delsm_back / d_deltheta_r;
				sm_back = d_sm_r1[n_theta_r_oct + itheta_r0] + deltheta_r_1 * dsm_back_o_dtheta;
				delsh_back = d_sh_r1[n_theta_r_oct + itheta_r1] - d_sh_r1[n_theta_r_oct + itheta_r0];
				dsh_back_o_dtheta = delsh_back / d_deltheta_r;

				sh_back = d_sh_r1[n_theta_r_oct + itheta_r0] + deltheta_r_1 * dsh_back_o_dtheta;
				delss_back = d_ss_r1[n_theta_r_oct + itheta_r1] - d_ss_r1[n_theta_r_oct + itheta_r0];
				dss_back_o_dtheta = delss_back / d_deltheta_r;
				ss_back = d_ss_r1[n_theta_r_oct + itheta_r0] + deltheta_r_1 * dss_back_o_dtheta;
				delslq2_back = d_slq2_r1[n_theta_r_oct + itheta_r1] - d_slq2_r1[n_theta_r_oct + itheta_r0];
				dslq2_back_o_dtheta = delslq2_back / d_deltheta_r;
				slq2_back = d_slq2_r1[n_theta_r_oct + itheta_r0] + deltheta_r_1 * dslq2_back_o_dtheta;
			}

			if (slq2_back < 0.0e0) {
				//printf("exit_513");
				return;
			}
			//s2_back = (ri1/back_ri1) * s2[k];
		}
		if (ri[k] <= 0.0e0) {//ifsali
			back_ra_r1 = 0.0e0;
			back_rit1 = 0.0e0;
			back_ric1 = 0.0e0;
			back_ri1 = 0.0e0;
			back_rid1 = 0.0e0;
		}

		//if(ri1<=0.0e0) s2_back = 0.0e0;
		if (ri1 <= 0.0e0) {
			sm_back = 0.0e0;
			sh_back = 0.0e0;
			ss_back = 0.0e0;
		}

		if ((sm_back < 0.0e0) || (sh_back < 0.0e0) || (ss_back < 0.0e0)) {
			//printf("exit_531");
			return;
		}
		tmp_back = 0.5e0 * d_b1 * d_b1 * back_ri1 * slq2_back * epson2;
		v_back = tmp_back * sm_back;
		t_back = tmp_back * sh_back;
		s_back = tmp_back * ss_back;

		if ((v_back < 0.0e0) || (t_back < 0.0e0) || (s_back < 0.0e0)) {
			//printf("exit_544");
			return;
		}

		if ((ri[k] > 0.0e0) && ((fabs(v_back) <= 0.0e0) || (fabs(t_back) <= 0.0e0) || (fabs(s_back) < 0.0e0))) {
			//printf("exit_549");
			return;
		}
		aldeep = 0.0e0;
		double delz, delrh, del2rh;
		if (ifbelow == 1) {//ifepson2
			if (ri1 >= 0.0e0) {
				tmp = 0.5e0 * d_b1 * d_b1 * ri1 * slq2 * epson2;
			} else if (k >= 2) {
				if (k == n - 1) {
					delz = z[k] - z[k - 1];
					delrh = rh[k] - rh[k - 1];
					del2rh = rh[k] - 2.0e0 * rh[k - 1] + rh[k - 2];
				} else {
					delz = z[k + 1] - z[k - 1];
					delrh = rh[k + 1] - rh[k - 1];
					del2rh = rh[k + 1] - 2.0e0 * rh[k] + rh[k - 1];
				}
				double dzrh = delrh / delz;
				double d2zrh = 4.0e0 * del2rh / (delz * delz);

				double rdzlndzrh = dzrh / d2zrh;
				double al0deep = 0.17e0 * fabs(rdzlndzrh);
				double akz = 0.4e0 * z[k];
				aldeep = akz * al0deep / (al0deep + akz);
				al2 = aldeep * aldeep;

				if (ifpureshear != 1) {
					if (ifshearmin) {
						s2[k] = fmax(s2[k], s2min);
					}
					tmp = 0.5e0 * d_b1 * al2 * sqrt(s2[k] / (slq2 + 1.0e-40));
				}
			} else {
				if (ifpureshear != 1) {
					s2[k] = fmax(s2[k], s2min);
					if (lifepsy) {
						tmp = 0.5e0 * epsy / (s2[k] + 1.0e-40);
					} else {
						tmp = 0.5e0 * d_b1 * al2 * sqrt(s2[k] / (slq2 + 1.0e-40));
					}
				}
			}
		} else {
			if (ifpureshear != 1) {
				s2[k] = fmax(s2[k], s2min);
				if (lifepsy) {
					tmp = 0.5e0 * epsy / (s2[k] + 1.0e-40);
				} else {
					tmp = 0.5e0 * d_b1 * al2 * sqrt(s2[k] / (slq2 + 1.0e-40));
				}
			}
		}
		if (ifpureshear == 1) tmp = 0.5e0 * (d_b1 * d_b1) * al2 * sqrt(-an2[k] / amtaun2);
		akm[k] = fmin(tmp * sm + v_back, visc_cbu_limit1);
		akh[k] = fmin(tmp * sh + t_back, diff_cbt_limit1);
		aks[k] = fmin(tmp * ss + s_back, diff_cbt_limit1);
	}

	for (k = nb; k < nmax; k++) {
		akm[k] = 0.0e0;
		akh[k] = 0.0e0;
		aks[k] = 0.0e0;
	}

	if (n > 0) {
		if (akm[0] < wndmix) akm[0] = wndmix;
		if (akh[0] < wndmix) akh[0] = wndmix;
		if (aks[0] < wndmix) aks[0] = wndmix;
	}

	for (k = 0; k < n; k++) {
		if ((akm[k] < 0.0e0) || (akh[k] < 0.0e0) || (aks[k] < 0.0e0)) {
			//printf("exit_622");
			return;
		}
	}
	return;
}

__global__ void turb_cu(double *d_akmt, double *d_akt, double *d_ak_tide, double *d_amld,
		double *d_amtaun2a1, double *d_and2on2a1, double *d_at, double *d_back_ra_r, double *d_buoysol,
		double *d_buoytur, double *d_dzp, double *d_fcort, double *d_fztidal, double *d_fz_tide,
		double *d_pdensity, double *d_rib, double *d_ricdttms, double *d_richardson, double *d_rict,
		double *d_rit, double *d_s2t, double *d_sha1, double *d_shb, double *d_sh_r1, double *d_slq2b,
		double *d_slq2_r1, double *d_sma1, double *d_smb, double *d_sm_r1, double *d_ssa1,
		double *d_ssb, double *d_ss_r1, double *d_turb_param, double *d_vit,
		double *d_wave_dis, double *d_wp3_tidal, double *d_zkp, int *d_irimax, int *d_kmt,
		double *d_diff_back,double *d_wk1, double *d_wk2, double* d_wk3,double *d_wp1, double *d_wp3, double *d_wp7, double *d_wp8) {
	int i, j, k;
	i = blockIdx.x * blockDim.x + threadIdx.x;
	j = blockIdx.y * blockDim.y + threadIdx.y;


	if (i > 0 && i < imt - 1 && j > 0 && j < jmt - 1 ){
		unsigned int tid = j * imt + i;
		if(d_vit[tid] > 0.5) {
		//	double wk1[km], wk2[km], wk3[km];
		//	double wp1[km], wp3[km];
			double wp4[km], wp5[km], wp6[km];
		//	double wp7[km], wp8[km];
			double mldtmp, wp10, wp11;
			int d_ncc = d_turb_param[7];
			double dfricmx = d_turb_param[8], dwndmix = d_turb_param[9];
			int iwk = d_kmt[tid] - 1;
			int iwk1 = (d_kmt[tid] > 2) ? iwk : 1;

			for (k = 0; k < km; k++) {
                                d_wp1[i * km * jmt+j*km+ k] = 0.0e0;
                                d_wp3[i * km * jmt+j*km+ k] = 0.0e0;
                                d_wp7[i * km * jmt+j*km+ k] = 0.0e0;
                                d_wp8[i * km * jmt+j*km+ k] = 0.0e0;

				wp4[k] = 0.0e0;
				wp5[k] = 0.0e0;
				wp6[k] = 0.0e0;
			}
			for (k = 0; k < iwk; k++) {
			       d_wp8[i * km * jmt+j*km+ k] = -d_vit[(k + 1) * imt * jmt + tid] * d_zkp[k + 1] * 1.0e+2;
			       d_wp1[i * km * jmt+j*km+ k]= d_vit[(k + 1) * imt * jmt + tid] * (d_at[k * imt * jmt + tid] -
						(d_at[k * imt * jmt + tid] -
						 d_at[(k + 1) * imt * jmt + tid]) *
						d_dzp[k] / (d_dzp[k] + d_dzp[k + 1]));
			}
			for (k = 0; k < iwk; k++) {
				wp4[k] = d_vit[(k + 1) * imt * jmt + tid] * d_rit[k * imt * jmt + tid];
				wp5[k] = d_vit[(k + 1) * imt * jmt + tid] * d_ricdttms[k * imt * jmt + tid];
				wp6[k] = d_vit[(k + 1) * imt * jmt + tid] * d_s2t[k * imt * jmt + tid];
			 d_wp7[i * km * jmt+j*km+ k]  = d_vit[(k + 1) * imt * jmt + tid] * d_rict[k * imt * jmt + tid];
				//tq = wp1 - d_ncc;
				//sq = (wp2-35.0e0)*1.0e-3 - d_turb_param[8];
			 d_wp3[i * km * jmt+j*km+ k]  = d_vit[(k + 1) * imt * jmt + tid] * (d_pdensity[k * imt * jmt + tid] +
						(d_pdensity[(k + 1) * imt * jmt + tid] -
						 d_pdensity[k * imt * jmt + tid]) * d_dzp[k] /
						(d_dzp[k] + d_dzp[k + 1])) * 1.0e-3;
			}
			//wp9=d_vit[tid]*d_ustar[tid]*1.0e+2;
			wp10 = d_vit[tid] * d_buoytur[tid] * 1.0e+4;
			wp11 = d_vit[tid] * d_buoysol[tid] * 1.0e+4;
                        turb_device(&d_wp8[i * km * jmt+j*km],&d_wp1[i * km * jmt+j*km],&d_wp3[i * km * jmt+j*km], wp4, wp5, wp6, 
                                          dfricmx * 1.0e+4, dwndmix * 1.0e+4, &d_wp7[i * km * jmt+j*km], wp10, wp11, d_fcort[tid],
                                        &mldtmp, &d_wk1[i * km * jmt+j*km], &d_wk2[i * km * jmt + j*km ], &d_wk3[ i * km * jmt + j*km],
                                        iwk, iwk1, km - 1, d_amtaun2a1, d_and2on2a1, d_back_ra_r, d_rib, d_sha1,
                                        d_shb, d_sh_r1, d_slq2b, d_slq2_r1, d_sma1, d_smb, d_sm_r1, d_ssa1, d_ssb, d_ss_r1,
                                        d_turb_param,d_irimax);



			d_amld[tid] = mldtmp * 1.0e-2;

#if ( defined TIDEMIX )
			//                    for(k=0 ;  k<km; k++)  {
			//                        pconst_mod_mp_ak_tide_[iblock][k][j][i]=0.0;
			//                    }
			if(isnan(d_wave_dis[tid])) printf("%s:%d,%d,%d",__FILE__,__LINE__,i,j);
			for (k = 0; k < iwk; k++) {
				d_ak_tide[k * imt * jmt + tid] = back_tidalmixing +
					mixing_ef * local_mixing_fraction * d_wave_dis[tid] *
					d_fz_tide[k * imt * jmt + tid] /
					(fmax(d_rict[k * imt * jmt + tid], 1.0e-8) *
					 ( d_wp3[i * km * jmt+j*km+ k]* 1000.0));

				if (d_ak_tide[k * imt * jmt + tid] > max_tidalmixing) d_ak_tide[k * imt * jmt + tid] = max_tidalmixing;
				d_richardson[k * imt * jmt + tid] = d_rict[k * imt * jmt + tid];
				d_fztidal[k * imt * jmt + tid] = d_fz_tide[k * imt * jmt + tid];
				d_wp3_tidal[k * imt * jmt + tid] = d_wp3[i * km * jmt+j*km+ k];
			}
#if ( defined CANUTOMIXOUT )
			for(k=0 ; k<grid_mp_kmt_[iblock][j][i]-1;  k++) {
				pconst_mod_mp_wp1_canuto_[iblock][k][j][i]=wp1[k];     //yuzp-2016/12/4
				pconst_mod_mp_wp2_canuto_[iblock][k][j][i]=wp2[k];    //yuzp-2016/12/4
				pconst_mod_mp_wp3_canuto_[iblock][k][j][i]=wp3[k];     //yuzp-2016/12/4
				pconst_mod_mp_wp4_canuto_[iblock][k][j][i]=wp4[k];   //yuzp-2016/12/4
				pconst_mod_mp_wp5_canuto_[iblock][k][j][i]=wp5[k];     //yuzp-2016/12/4
				pconst_mod_mp_wp6_canuto_[iblock][k][j][i]=wp6[k];    //yuzp-2016/12/4
				pconst_mod_mp_wp7_canuto_[iblock][k][j][i]=wp7[k];     //yuzp-2016/12/4
				pconst_mod_mp_wp8_canuto_[iblock][k][j][i]=wp8[k];    //yuzp-2016/12/4

				pconst_mod_mp_wp12_canuto_[iblock][k][j][i]=wp12[iblock][k][j][i];   //yuzp-2016/12/4
				pconst_mod_mp_wp13_canuto_[iblock][k][j][i]=wp13[iblock][k][j][i];    //yuzp-2016/12/4
				pconst_mod_mp_wk4_canuto_[iblock][k][j][i]=wk4[k];     //yuzp-2016/12/4
			}
#endif
			for (k = iwk - 2; k >= 0; k--) {
				d_ak_tide[k * imt * jmt + tid] = fmin(d_ak_tide[k * imt * jmt + tid],
						d_ak_tide[(k + 1) * imt * jmt + tid]);
			}
#endif  //TIDEMIX

			for (k = 0; k < km - 1; k++) {
#if ( defined TIDEMIX )
                       d_akmt[k * imt * jmt + tid] = d_wk1[i * km * jmt + j * km + k] * 1.0e-4 +
                        d_ak_tide[k * imt * jmt + tid] * 5.0 ;
                       d_akt[k * imt * jmt + tid] = (d_wk2[i * km * jmt + j * km + k] * 1.0e-4 +
                        d_ak_tide[k * imt * jmt + tid] ) / (float) d_ncc ;
                       d_akt[km * imt * jmt + k * imt * jmt + tid] = (d_wk3[i * km * jmt +j * km + k] * 1.0e-4 +
                        d_ak_tide[k * imt * jmt + tid] )/ (float) d_ncc ;

//				d_akmt[k * imt * jmt + tid] = wk1[k] * 1.0e-4 + d_ak_tide[k * imt * jmt + tid] * 5.0 ;
//				d_akt[k * imt * jmt + tid] += (wk2[k] * 1.0e-4 + d_ak_tide[k * imt * jmt + tid]) / (float) d_ncc;
//				d_akt[km * imt * jmt + k * imt * jmt + tid] += (wk3[k] * 1.0e-4 + d_ak_tide[k * imt * jmt + tid] ) / (float) d_ncc;
#if ( defined CANUTOMIXOUT )
				pconst_mod_mp_wk1_canuto_[iblock][k][j][i]=wk1[k];     //yuzp-2016/12/4
				pconst_mod_mp_wk2_canuto_[iblock][k][j][i]=wk2[k];     //yuzp-2016/12/4
				pconst_mod_mp_wk3_canuto_[iblock][k][j][i]=wk3[k];     //yuzp-2016/12/4
#endif
#ifdef BCKMEX
				//                            pconst_mod_mp_akmt_[iblock][k][j][i]=pconst_mod_mp_akmt_[iblock][k][j][i]+diff_back[iblock][j][i]*10.0*1.e-4; //10--Pr_number
				//                            pconst_mod_mp_akt_[iblock][0][k][j][i]=pconst_mod_mp_akt_[iblock][0][k][j][i]+diff_back[iblock][j][i]/(float)pconst_mod_mp_ncc_*1.e-4;
				//                            pconst_mod_mp_akt_[iblock][1][k][j][i]=pconst_mod_mp_akt_[iblock][1][k][j][i]+diff_back[iblock][j][i]/(float)pconst_mod_mp_ncc_*1.e-4;
				d_akmt[k * imt * jmt + tid] += d_diff_back[tid] * 10.0 * 1.0e-4;
				d_akt[k * imt * jmt + tid]  += d_diff_back[tid] * 1.0e-4 / (float) d_ncc;
				d_akt[km * imt * jmt + k * imt * jmt + tid] += d_diff_back[tid] * 1.0e-4 / (float) d_ncc;
#endif
#else  //TIDEMIX
				//                            pconst_mod_mp_akmt_[iblock][k][j][i]=wk1[k]*1.0e-4;
				//                            pconst_mod_mp_akt_[iblock][0][k][j][i]=pconst_mod_mp_akt_[iblock][0][k][j][i]+ wk2[k]*1.0e-4/(float)pconst_mod_mp_ncc_;
				//                            pconst_mod_mp_akt_[iblock][1][k][j][i]=pconst_mod_mp_akt_[iblock][1][k][j][i]+ WK3(K)*1.0e-4/(float)pconst_mod_mp_ncc_;
#ifdef BCKMEX
				//                            pconst_mod_mp_akmt_[iblock][k][j][i]=pconst_mod_mp_akmt_[iblock][k][j][i]+diff_back[iblock][j][i]*10.0*1.e-4; //10--Pr_number
				//                            pconst_mod_mp_akt_[iblock][0][k][j][i]=pconst_mod_mp_akt_[iblock][0][k][j][i]+diff_back[iblock][j][i]/(float)pconst_mod_mp_ncc_*1.e-4;
				//                            pconst_mod_mp_akt_[iblock][1][k][j][i]=pconst_mod_mp_akt_[iblock][1][k][j][i]+diff_back[iblock][j][i]/(float)pconst_mod_mp_ncc_*1.e-4;
				d_akmt[k * imt * jmt + tid] += d_diff_back[tid] * 10.0 * 1.0e-4;
				d_akt[k * imt * jmt + tid]  += d_diff_back[tid] * 1.0e-4 / (float) d_ncc;
				d_akt[km * imt * jmt + k * imt * jmt + tid] += d_diff_back[tid] * 1.0e-4 / (float) d_ncc;
#endif
#endif //TIDEMIX
// wpf, for test model without turb2, lhl suggest the following value
//            d_akmt[k * imt * jmt + tid] = 1.0e-4;
//            d_akt[k * imt * jmt + tid]  = 1.0e-5;
//            d_akt[km * imt * jmt + k * imt * jmt + tid] = 1.0e-5;
//        d_amld[tid] = 1.0e-2;
			}
		}}

	return;
}

extern "C" void readyc_() {
	//allocate_readyc_();

	int errorcode;

	struct block this_block;
	this_block = blocks_mp_all_blocks_[domain_mp_blocks_clinic_[0] - 1];
	int jb = this_block.jb;
	int je = this_block.je;
	int ib = this_block.ib;
	int ie = this_block.ie;

	dim3
		blockSize(8, 16, 1);
	dim3
		gridSize((imt + blockSize.x - 1) / blockSize.x,
				(jmt + blockSize.y - 1) / blockSize.y);
	dim3
		blockSize3d(8, 16, 4);
	dim3
		gridSize3d((imt + blockSize3d.x - 1) / blockSize3d.x,
				(jmt + blockSize3d.y - 1) / blockSize3d.y,
				(km + blockSize3d.z - 1) / blockSize3d.z);

	hipMemset(d_s2t, c0, dataSize3d);
	hipMemset(d_riu, c0, dataSize3d); //jjr test bug
	hipMemset(d_wp12, c0, dataSize3d);
	hipMemset(d_wp13, c0, dataSize3d);
	CHECK(hipMemset(d_amld, c0, dataSize2d));
	hipMemset(d_riv1, c0, dataSize3d);
	hipMemset(d_riv2, c0, dataSize3d);
	//hipMemset(d_diff_back_sh, c0, dataSize2d);
	//hipMemset(d_diff_back_nh, c0, dataSize2d);
	CHECK(hipMemset(d_akmt, c0, dataSize3d));

	iCheck();
	hipLaunchKernelGGL(readyc1, dim3(gridSize3d), dim3(blockSize3d ), 0, 0, d_h0bl, d_h0bf, d_h0,
			d_up, d_vp, d_viv, d_vit, d_wp12, d_wp13,
			d_at0, d_atn, d_ate, d_atne);
	iCheck();

	hipLaunchKernelGGL(readyc1_2, dim3(gridSize3d), dim3(blockSize3d ), 0, 0, d_wp12, d_wp13, d_vit, d_s2t,
			d_riv1, d_riv2, d_odzt, d_rit, d_rict, d_up, d_vp,
			d_s2u, d_viv, d_riu, d_ric);
	iCheck();

#ifdef BCKMEX
	hipLaunchKernelGGL(readyc1_3, dim3(gridSize), dim3(blockSize ), 0, 0, d_diff_back, d_diff_back_sh,
			d_diff_back_nh, pconst_mod_mp_diff_back_coef_max_,
			pconst_mod_mp_diff_back_eq_, pconst_mod_mp_diff_back_coef_, d_tlat);
	CHECK(hipDeviceSynchronize()); 
	iCheck();
#endif
	CHECK(hipMemset(d_akmu, c0, dataSize3d));
	CHECK(hipMemset(d_ak_tide, c0, dataSize3d));
        CHECK(hipMemset(d_wk1, c0, dataSize3d));
        CHECK(hipMemset(d_wk2, c0, dataSize3d));
        CHECK(hipMemset(d_wk3, c0, dataSize3d));


	hipLaunchKernelGGL(turb_cu, dim3(gridSize), dim3(blockSize ), 0, 0, d_akmt, d_akt, d_ak_tide, d_amld,
			d_amtaun2a1, d_and2on2a1, d_at, d_back_ra_r, d_buoysol, d_buoytur, d_dzp, d_fcort, d_fztidal,
			d_fz_tide, d_pdensity, d_rib, d_ricdttms, d_richardson, d_rict, d_rit, d_s2t, d_sha1, d_shb,
			d_sh_r1, d_slq2b, d_slq2_r1, d_sma1, d_smb, d_sm_r1, d_ssa1, d_ssb, d_ss_r1, d_turb_param, d_vit,
			d_wave_dis, d_wp3_tidal, d_zkp, d_irimax, d_kmt, d_diff_back, d_wk1,d_wk2,d_wk3,d_wp1,d_wp3,d_wp7,d_wp8);
	iCheck();
//	CHECK(hipMemcpy(tracer_mod_mp_amld_, d_amld, dataSize2d, hipMemcpyDeviceToHost));
//	pop_haloupdate_readyc_(&errorcode);//amld

	hipLaunchKernelGGL(readyc2, dim3(gridSize3d), dim3(blockSize3d ), 0, 0, d_akmu, d_au0, d_akmt, d_aus, d_auw, d_ausw, d_viv);
	iCheck();

/*	hipMemcpy(pmix_mod_mp_rict_, d_rict, dataSize3333d, hipMemcpyDeviceToHost);
	int i, j, k;
	for (j = 0; j < jmt; j++) {
		for (i = 0; i < imt; i++) {
			pmix_mod_mp_rict_ref_[0][j][i] = pmix_mod_mp_rict_[0][13][j][i];  //why 14layer
			for (k = 0; k < kmm1; k++) {
				if (fabs(pconst_mod_mp_zkp_[k + 1]) > tracer_mod_mp_amld_[0][j][i]) {
					pmix_mod_mp_rict_ref_[0][j][i] = pmix_mod_mp_rict_[0][k][j][i];
					break;
				}
			}
		}
	}
*/
	CHECK(hipMemset(d_uk, c0, dataSize3d));
	CHECK(hipMemset(d_vk, c0, dataSize3d));
	CHECK(hipMemset(d_wka, c0, dataSize3d));
	CHECK(hipMemset(d_work, c0, dataSize2d));

	hipLaunchKernelGGL(upwell_cu_0, dim3(gridSize), dim3(blockSize ), 0, 0, d_au0, d_aus, d_auw, d_ausw,
			d_h0, d_ohbu,  d_u, d_uk, d_v, d_vk, d_work);

	hipLaunchKernelGGL(upwell_cu, dim3(gridSize), dim3(blockSize ), 0, 0, d_dzp, d_h0, d_hts, d_htw, d_ohbt,
			d_tarea_r, d_uk, d_vit, d_vk, d_wka, d_work, d_ws, d_kmt);
	iCheck();

	CHECK(hipMemset(d_dlu, c0, dataSize3d));
	CHECK(hipMemset(d_dlv, c0, dataSize3d));
	CHECK(hipMemset(d_u_wface, c0, dataSize3d));
	CHECK(hipMemset(d_v_sface, c0, dataSize3d));
	CHECK(hipMemset(d_wka, c0, dataSize3d));

	hipLaunchKernelGGL(readyc3, dim3(gridSize3d), dim3(blockSize3d ), 0, 0, d_wka, d_ws, d_au0, d_ausw, d_aus, d_auw, d_u_wface, d_v_sface, d_u, d_v, d_hue, d_hun);
	iCheck();


	CHECK(hipDeviceSynchronize());
	hipLaunchKernelGGL(advection_momentum_cu, dim3(gridSize3d), dim3(blockSize3d), 0, 0, d_dlu, d_dlv, d_odzp,
			d_u, d_uarea_r, d_u_wface, d_v, d_v_sface, d_wka);
	iCheck();

/*	
	CHECK(hipMemcpy(dyn_mod_mp_u_, d_u, dataSize3d, hipMemcpyDeviceToHost));
	CHECK(hipMemcpy(dyn_mod_mp_v_, d_v, dataSize3d, hipMemcpyDeviceToHost));
	CHECK(hipMemcpy(work_mod_mp_wka_, d_wka, dataSize3d, hipMemcpyDeviceToHost));
	readycf_();
	CHECK(hipMemcpy(d_dlu, dyn_mod_mp_dlu_, dataSize3d, hipMemcpyHostToDevice));
	CHECK(hipMemcpy(d_dlv, dyn_mod_mp_dlv_, dataSize3d, hipMemcpyHostToDevice));
*/

#ifdef BIHAR
/*
//use F90 to do diffu_del4, when HIP code ready,change here
        CHECK(hipMemcpy(dyn_mod_mp_up_, d_up, dataSize3d, hipMemcpyDeviceToHost));
        CHECK(hipMemcpy(dyn_mod_mp_vp_, d_vp, dataSize3d, hipMemcpyDeviceToHost));
        CHECK(hipMemcpy(dyn_mod_mp_dlu_, d_dlu, dataSize3d, hipMemcpyDeviceToHost));
        CHECK(hipMemcpy(dyn_mod_mp_dlv_, d_dlv, dataSize3d, hipMemcpyDeviceToHost));
        readyc_udel4_();
	CHECK(hipMemcpy(d_dlu, dyn_mod_mp_dlu_, dataSize3d, hipMemcpyHostToDevice));
	CHECK(hipMemcpy(d_dlv, dyn_mod_mp_dlv_, dataSize3d, hipMemcpyHostToDevice));
//get dlu,dlv back
	hipLaunchKernelGGL(hdiffu_del4_2, dim3(gridSize3d), dim3(blockSize3d ), 0, 0, d_hduk2, d_hdvk2, jb, je, ib, ie, d_wka,
			pconst_mod_mp_c0f_, d_up, d_vp, d_duc, d_dum, d_dun, d_dus, d_due,
			d_duw, d_dmc, d_dmn, d_dms, d_dme, d_dmw, hmix_del4_mp_am_, d_dlu, d_dlv, d_viv, d_kmu);
	iCheck();
*/


       hipLaunchKernelGGL(hdel4_div, dim3(gridSize3d), dim3(blockSize3d ), 0, 0, d_hts, d_htw, d_tarea_r,
               d_up,d_vp,d_wka,d_kmt) ;

	hipMemset(d_am_factor, c0, dataSize3d);
       hipLaunchKernelGGL(hdel4_grad, dim3(gridSize3d), dim3(blockSize3d ), 0, 0, d_wka, d_dxur, d_dyur,
               d_kmu,d_am_factor) ;
       hipLaunchKernelGGL(hdel4_zcurl, dim3(gridSize3d), dim3(blockSize3d ), 0, 0, d_dxu, d_dyu, 
               d_up,d_vp,d_wka,d_kmt) ;

       hipLaunchKernelGGL(hdel4_grad, dim3(gridSize3d), dim3(blockSize3d ), 0, 0, d_wka, d_dxur, d_dyur,
               d_kmu,d_am_factor) ;

	hipLaunchKernelGGL(hdiffu_del4_0, dim3(gridSize3d), dim3(blockSize3d ), 0, 0, d_hduk2, d_hdvk2, jb, je, ib, ie, d_uarea_r,
			pconst_mod_mp_c0f_, d_up, d_vp, d_duc, d_dum, d_dun, d_dus, d_due,
			d_duw, d_dmc, d_dmn, d_dms, d_dme, d_dmw, hmix_del4_mp_am_, d_dlu, d_dlv, d_viv, d_kmu,
			d_amf, d_am_factor, d_d2uk, d_d2vk);
	iCheck();
	hipLaunchKernelGGL(hdiffu_del4_1, dim3(gridSize3d), dim3(blockSize3d ), 0, 0, d_hduk2, d_hdvk2, jb, je, ib, ie, d_wka,
			pconst_mod_mp_c0f_, d_up, d_vp, d_duc, d_dum, d_dun, d_dus, d_due,
			d_duw, d_dmc, d_dmn, d_dms, d_dme, d_dmw, hmix_del4_mp_am_, d_dlu, d_dlv, d_viv, d_kmu,
			d_amf, d_am_factor, d_d2uk, d_d2vk);
	iCheck();

#else
	hipLaunchKernelGGL(hdiffu_del2_cu, dim3(gridSize3d), dim3(blockSize3d ), 0, 0, d_hduk2, d_hdvk2, jb, je, ib, ie, d_wka,
			pconst_mod_mp_c0f_, d_up, d_vp, d_duc, d_dum, d_dun, d_dus, d_due,
			d_duw, d_dmc, d_dmn, d_dms, d_dme, d_dmw, hmix_del2_mp_am_, d_dlu, d_dlv, d_viv, d_kmu);
	iCheck();
#endif BIHAR

	hipLaunchKernelGGL(readyc4, dim3(gridSize), dim3(blockSize ), 0, 0, d_sbcx, d_sbcy, d_su, d_sv, d_bbcx, pconst_mod_mp_od0_,
			pconst_mod_mp_c0f_, d_up, d_vp, d_snlat, pconst_mod_mp_cag_, pconst_mod_mp_sag_, d_dlub, d_dlvb,
			d_dzp, d_ohbu, d_dlu, d_viv, d_dlv, d_kmu, d_bbcy);
	iCheck();

	hipLaunchKernelGGL(readyc4_2, dim3(gridSize3d), dim3(blockSize3d ), 0, 0, d_su, d_sv, pconst_mod_mp_od0_, d_up, d_vp,
			pconst_mod_mp_sag_, pconst_mod_mp_cag_, d_akmu, d_odzt, d_viv, d_wka, d_snlat, d_dlv, d_odzp, d_dlu);
	iCheck();
}

