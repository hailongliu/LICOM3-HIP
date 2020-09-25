#include "hip/hip_runtime.h"
#include "cuda_data.h"

#include "forc_mod.h"
#include "tracer_mod.h"
#include "pconst_mod.h"
#include "buf_mod.h"
#include "dyn_mod.h"

extern "C" void steponinit(){
    CHECK(hipMemcpy(d_ssf, forc_mod_mp_ssf_, dataSize2d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_sss, forc_mod_mp_sss_, dataSize2d, hipMemcpyHostToDevice));
//    CHECK(hipMemcpy(d_sst, forc_mod_mp_sst_, dataSize2d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_seaice, forc_mod_mp_seaice_, dataSize2d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_fresh, forc_mod_mp_fresh_, dataSize2d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_tsf, forc_mod_mp_tsf_, dataSize2d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_su, forc_mod_mp_su_, dataSize2d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_sv, forc_mod_mp_sv_, dataSize2d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_nswv, forc_mod_mp_nswv_, dataSize2d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_swv, forc_mod_mp_swv_, dataSize2d, hipMemcpyHostToDevice));

//    CHECK(hipMemcpy(d_licomqice, tracer_mod_mp_licomqice_, dataSize2d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_atb, tracer_mod_mp_atb_, dataSize44d, hipMemcpyHostToDevice));
//    CHECK(hipMemcpy(d_at, tracer_mod_mp_at_, dataSize4d, hipMemcpyHostToDevice));
//wpf 0309, sd delete this code
//    if(pconst_mod_mp_simple_assm_ == true){
//    	hipMemcpy(d_restore_at, tracer_mod_mp_restore_at_, dataSize4d, hipMemcpyHostToDevice);
//    }
//    CHECK(hipMemcpy(d_ifrac, buf_mod_mp_ifrac_, dataSize2d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_h0,dyn_mod_mp_h0_, dataSize2d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_up, dyn_mod_mp_up_, dataSize3d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_vp, dyn_mod_mp_vp_, dataSize3d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_utf, dyn_mod_mp_utf_, dataSize3d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_vtf, dyn_mod_mp_vtf_, dataSize3d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_h0p, dyn_mod_mp_h0p_, dataSize2d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_h0f, dyn_mod_mp_h0f_, dataSize2d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_h0bf, dyn_mod_mp_h0bf_, dataSize2d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_vbp, dyn_mod_mp_vbp_, dataSize2d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_ubp, dyn_mod_mp_ubp_, dataSize2d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_vb, dyn_mod_mp_vb_, dataSize2d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_ub, dyn_mod_mp_ub_, dataSize2d, hipMemcpyHostToDevice));
//netstep


   

    return;
}
