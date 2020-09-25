#include "hip/hip_runtime.h"
#include "cuda_data.h"

#include "dyn_mod.h"
#include "tracer_mod.h"
#include "common.h"

extern "C" void energyMemcpy() {
    CHECK(hipMemcpy(dyn_mod_mp_u_, d_u, dataSize3d, hipMemcpyDeviceToHost));
    CHECK(hipMemcpy(dyn_mod_mp_v_, d_v, dataSize3d, hipMemcpyDeviceToHost));
   CHECK(hipMemcpy(tracer_mod_mp_at_, d_at, dataSize4d, hipMemcpyDeviceToHost));
//    CHECK(hipMemcpy(d_h0, dyn_mod_mp_h0_, dataSize2d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy( dyn_mod_mp_h0_,d_h0, dataSize2d, hipMemcpyDeviceToHost));
//    hipMemcpy( dyn_mod_mp_vb_,d_vb,dataSize2d, hipMemcpyDeviceToHost);
 //   hipMemcpy(dyn_mod_mp_ub_,d_ub, dataSize2d, hipMemcpyDeviceToHost);
}
