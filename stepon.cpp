#include <stdio.h>
#include "hip/hip_runtime.h"
#include "param_mod.h"
#include "pconst_mod.h"
#include "blocks.h"
#include "domain.h"
#include "cuda_data.h"
#include "dyn_mod.h"
#include "tracer_mod.h"

extern "C" void readyt_();
extern "C" void energy();
extern "C" void readyc_();
extern "C" void barotr_();
extern "C" void bclinc_();
extern "C" void tracer_(int);
extern "C" void icesnow_();
extern "C" void convadj_();
extern "C" void accumm_();

extern "C" void cudainit();
extern "C" void cudafinal();
extern "C" void steponinit();
extern "C" void turbinit();

extern "C" void energyMemcpy();
int init_flag=0;
extern "C" void stepon_(int *pconst_mod_mp_nss_, int *num_cpl, int *ncc, int *dts_accum) {
    int ii, jj;

    struct block this_block;
    this_block = blocks_mp_all_blocks_[domain_mp_blocks_clinic_[0] - 1];
    int jb = this_block.jb;
    int je = this_block.je;
    int ib = this_block.ib;
    int ie = this_block.ie;

    if ( init_flag==0) {
        int deviceCount;
        hipGetDeviceCount(&deviceCount);
/*
        if (deviceCount == 1) {
	    hipSetDevice(0);
            printf("myid=%d,deviceCount=%d\n", param_mod_mp_mytid_, deviceCount);
        }else{
		 hipSetDevice(param_mod_mp_mytid_ % deviceCount);
            printf("myid=%d, deviceId=%d\n", param_mod_mp_mytid_, param_mod_mp_mytid_ % deviceCount);
	}
*/
       cudainit();
       turbinit();
       init_flag=1;
    }

    steponinit();

    for (ii = 0; ii < *pconst_mod_mp_nss_ / (*num_cpl); ii++) {
//        if(param_mod_mp_mytid_==0)printf("in GPU step: %d , proc %03d\n",ii, param_mod_mp_mytid_);
        readyt_();

        for (jj = 0; jj < *ncc; jj++) {
            readyc_();

            barotr_();

            bclinc_();
        }

        tracer_(ii);

        icesnow_();
//    energyMemcpy();
  //  energy();



        if(ii==0) convadj_();

    }

    
//wpf, copy back to do addps, nextstep
    hipMemcpy(tracer_mod_mp_at_, d_at, dataSize4d, hipMemcpyDeviceToHost);
    hipMemcpy(dyn_mod_mp_h0_, d_h0, dataSize2d, hipMemcpyDeviceToHost);
    CHECK(hipMemcpy(dyn_mod_mp_u_, d_u, dataSize3d, hipMemcpyDeviceToHost));
    CHECK(hipMemcpy(dyn_mod_mp_v_, d_v, dataSize3d, hipMemcpyDeviceToHost));
    CHECK(hipMemcpy(dyn_mod_mp_ws_, d_ws, dataSize3d, hipMemcpyDeviceToHost));

    return;
}
