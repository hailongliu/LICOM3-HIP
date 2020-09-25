#include <sys/time.h>
#include <stdio.h>
#include "param_mod.h"
#ifndef _COMMON_H
#define _COMMON_H
#define iCheck()                                                                \
{                                                                              \
    const hipError_t error = hipGetLastError();                               \
    if (error != hipSuccess )                       \
    {                                                                          \
        fprintf(stderr,"Error: %s:%d, ", __FILE__, __LINE__);                 \
        fprintf(stderr,"code: %d, reason: %s tid=%d\n", error,                       \
                hipGetErrorString(error),param_mod_mp_mytid_);                                    \
        exit(1);                                                               \
    }                                                                          \
}

/*
#define CHECK(call)                                                            \
{                                                                              \
    const hipError_t error = call;                                            \
    if (error != hipSuccess )                       \
    {                                                                          \
        fprintf(stderr,"Error: %s:%d, ", __FILE__, __LINE__);                 \
        fprintf(stderr,"code: %d, reason: %s\n", error,                       \
                hipGetErrorString(error));                                    \
        exit(1);                                                               \
    }                                                                          \
}
*/
#define CHECK_CUBLAS(call)                                                     \
{                                                                              \
    hipblasStatus_t err;                                                        \
    if ((err = (call)) != HIPBLAS_STATUS_SUCCESS)                               \
    {                                                                          \
        fprintf(stderr, "Got CUBLAS error %d at %s:%d\n", err, __FILE__,       \
                __LINE__);                                                     \
        exit(1);                                                               \
    }                                                                          \
}

#define CHECK_CURAND(call)                                                     \
{                                                                              \
    hiprandStatus_t err;                                                        \
    if ((err = (call)) != HIPRAND_STATUS_SUCCESS)                               \
    {                                                                          \
        fprintf(stderr, "Got CURAND error %d at %s:%d\n", err, __FILE__,       \
                __LINE__);                                                     \
        exit(1);                                                               \
    }                                                                          \
}

#define CHECK_CUFFT(call)                                                      \
{                                                                              \
    hipfftResult err;                                                           \
    if ( (err = (call)) != HIPFFT_SUCCESS)                                      \
    {                                                                          \
        fprintf(stderr, "Got CUFFT error %d at %s:%d\n", err, __FILE__,        \
                __LINE__);                                                     \
        exit(1);                                                               \
    }                                                                          \
}

#define CHECK_CUSPARSE(call)                                                   \
{                                                                              \
    hipsparseStatus_t err;                                                      \
    if ((err = (call)) != HIPSPARSE_STATUS_SUCCESS)                             \
    {                                                                          \
        fprintf(stderr, "Got error %d at %s:%d\n", err, __FILE__, __LINE__);   \
        hipError_t cuda_err = hipGetLastError();                             \
        if (cuda_err != hipSuccess)                                           \
        {                                                                      \
            fprintf(stderr, "  CUDA error \"%s\" also detected\n",             \
                    hipGetErrorString(cuda_err));                             \
        }                                                                      \
        exit(1);                                                               \
    }                                                                          \
}

inline double seconds()
{
    struct timeval tp;
    struct timezone tzp;
    int i = gettimeofday(&tp, &tzp);
    return ((double)tp.tv_sec + (double)tp.tv_usec * 1.e-6);
}

#endif // _COMMON_H
