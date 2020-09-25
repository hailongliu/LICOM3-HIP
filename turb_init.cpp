#include <stdio.h>
#include "hip/hip_runtime.h"
#include "canuto_mod.h"
#include "turb.h"
#include "param_mod.h"
#include "grid.h"
#include "forc_mod.h"
#include "cuda_data.h"
#include "pconst_mod.h"
#include "canuto_mod.h"
#include "turb.h"
#include "common.h"

size_t dataSize1d_c1 = (2*mt+1) * sizeof(double);
size_t dataSize1d_c2 = (4*n_theta_r_oct+1) * sizeof(double);
double *d_sma1 = NULL;
double *d_sha1 = NULL;
double *d_ssa1 = NULL;
double *d_back_ra_r = NULL;
double *d_sm_r1 = NULL;
double *d_sh_r1 = NULL;
double *d_ss_r1 = NULL;
double *d_slq2_r1 = NULL;
double *d_and2on2a1 = NULL;
double *d_amtaun2a1 = NULL;
double *d_turb_param = NULL;
double *d_bclinc_param = NULL;

int *d_irimax = NULL;
double *d_rib = NULL;
double *d_shb = NULL;
double *d_slq2b = NULL;
double *d_smb = NULL;
double *d_ssb = NULL;
double *d_fcort = NULL;

extern "C" void turbinit(){
	
	CHECK(hipMalloc((void **)&d_sma1, dataSize1d_c1));
	CHECK(hipMalloc((void **)&d_sha1, dataSize1d_c1));
	CHECK(hipMalloc((void **)&d_ssa1, dataSize1d_c1));
	CHECK(hipMalloc((void **)&d_back_ra_r, dataSize1d_c2));
	CHECK(hipMalloc((void **)&d_sm_r1, dataSize1d_c2));
	CHECK(hipMalloc((void **)&d_sh_r1, dataSize1d_c2));
	CHECK(hipMalloc((void **)&d_ss_r1, dataSize1d_c2));
	CHECK(hipMalloc((void **)&d_slq2_r1, dataSize1d_c2));
	CHECK(hipMalloc((void **)&d_and2on2a1, dataSize1d_c1));
	CHECK(hipMalloc((void **)&d_amtaun2a1, dataSize1d_c1));
	CHECK(hipMalloc((void **)&d_turb_param, 11*sizeof(double)));
	CHECK(hipMalloc((void **)&d_bclinc_param, 10*sizeof(double)));
	
	CHECK(hipMalloc((void **)&d_irimax, (2*mt+1)*sizeof(int)));
	CHECK(hipMalloc((void **)&d_rib, (2*mt+1)*sizeof(double)));
	CHECK(hipMalloc((void **)&d_shb, (2*mt+1)*(2*mt+1)*sizeof(double)));
	CHECK(hipMalloc((void **)&d_slq2b, (2*mt+1)*(2*mt+1)*sizeof(double)));
	CHECK(hipMalloc((void **)&d_smb, (2*mt+1)*(2*mt+1)*sizeof(double)));
	CHECK(hipMalloc((void **)&d_ssb, (2*mt+1)*(2*mt+1)*sizeof(double)));
	CHECK(hipMalloc((void **)&d_fcort, imt*jmt*sizeof(double)));

	CHECK(hipMemcpy(d_irimax, canuto_mod_mp_irimax_, (2*mt+1)*sizeof(int), hipMemcpyHostToDevice));
	CHECK(hipMemcpy(d_rib, canuto_mod_mp_rib_, (2*mt+1)*sizeof(double), hipMemcpyHostToDevice));
	CHECK(hipMemcpy(d_slq2b, canuto_mod_mp_slq2b_, (2*mt+1)*(2*mt+1)*sizeof(double), hipMemcpyHostToDevice));
	CHECK(hipMemcpy(d_smb, canuto_mod_mp_smb_, (2*mt+1)*(2*mt+1)*sizeof(double), hipMemcpyHostToDevice));
	CHECK(hipMemcpy(d_shb, canuto_mod_mp_shb_, (2*mt+1)*(2*mt+1)*sizeof(double), hipMemcpyHostToDevice));
	CHECK(hipMemcpy(d_ssb, canuto_mod_mp_ssb_, (2*mt+1)*(2*mt+1)*sizeof(double), hipMemcpyHostToDevice));
	CHECK(hipMemcpy(d_and2on2a1, canuto_mod_mp_and2on2a1_,(2*mt+1)*sizeof(double), hipMemcpyHostToDevice));
	CHECK(hipMemcpy(d_amtaun2a1, canuto_mod_mp_amtaun2a1_, (2*mt+1)*sizeof(double), hipMemcpyHostToDevice));
	CHECK(hipMemcpy(d_fcort,grid_mp_fcort_, imt*jmt*sizeof(double), hipMemcpyHostToDevice));

	CHECK(hipMemcpy(d_ustar, forc_mod_mp_ustar_, dataSize2d, hipMemcpyHostToDevice));
	CHECK(hipMemcpy(d_wave_dis, forc_mod_mp_wave_dis_, dataSize2d, hipMemcpyHostToDevice));
//	CHECK(hipMemcpy(d_dzp, pconst_mod_mp_dzp_, dataSize1d, hipMemcpyHostToDevice));
	CHECK(hipMemcpy(d_zkp, pconst_mod_mp_zkp_, dataSize1112d, hipMemcpyHostToDevice));
	CHECK(hipMemcpy(d_sma1, canuto_mod_mp_sma1_, dataSize1d_c1, hipMemcpyHostToDevice));
	CHECK(hipMemcpy(d_sha1, canuto_mod_mp_sha1_, dataSize1d_c1, hipMemcpyHostToDevice));
	CHECK(hipMemcpy(d_ssa1, canuto_mod_mp_ssa1_, dataSize1d_c1, hipMemcpyHostToDevice));
	CHECK(hipMemcpy(d_back_ra_r, canuto_mod_mp_back_ra_r_, dataSize1d_c2, hipMemcpyHostToDevice));
	CHECK(hipMemcpy(d_sm_r1, canuto_mod_mp_sm_r1_, dataSize1d_c2, hipMemcpyHostToDevice));
	CHECK(hipMemcpy(d_sh_r1, canuto_mod_mp_sh_r1_, dataSize1d_c2, hipMemcpyHostToDevice));
	CHECK(hipMemcpy(d_ss_r1, canuto_mod_mp_ss_r1_, dataSize1d_c2, hipMemcpyHostToDevice));
	CHECK(hipMemcpy(d_slq2_r1, canuto_mod_mp_slq2_r1_, dataSize1d_c2, hipMemcpyHostToDevice));

        double h_turb_param[11]={canuto_mod_mp_rri_,canuto_mod_mp_rnd2on2_,canuto_mod_mp_dri_, canuto_mod_mp_deltheta_r_,canuto_mod_mp_b1_,canuto_mod_mp_theta_rcrp_, canuto_mod_mp_theta_rcrn_,(double)pconst_mod_mp_ncc_,pconst_mod_mp_dfricmx_,pconst_mod_mp_dwndmix_,canuto_mod_mp_dand2on2_};
        double h_bclinc_param[10]={pconst_mod_mp_od0_,pconst_mod_mp_c0f_,pconst_mod_mp_cag_,pconst_mod_mp_sag_,pconst_mod_mp_onbb_,pconst_mod_mp_afc1_,pconst_mod_mp_afc2_,pconst_mod_mp_dtc_,pconst_mod_mp_dtc2_};
        CHECK(hipMemcpy(d_turb_param, h_turb_param, 11*sizeof(double), hipMemcpyHostToDevice));
        CHECK(hipMemcpy(d_bclinc_param, h_bclinc_param, 9*sizeof(double), hipMemcpyHostToDevice));

}
