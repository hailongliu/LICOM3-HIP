#ifndef INCLUDE_CONSTANT
#define INCLUDE_CONSTANT

#include <stdlib.h>
#include <math.h>

// #include "precision_mod.h"
// #include "param_mod.h"
// #include "pconst_mod.h"
// #include "pmix_mod.h"
// #include "msg_mod.h"
// #include "shr_const_mod.h"
// #include "shr_msg_mod.h"
// #include "shr_kind_mod.h"
// #include "shr_sys_mod.h"
// #include "netcdf.h"

  //----- constants -----
#define latvap  (double)SHR_CONST_LATVAP // latent heat of evap   ~ J/kg
#define Tzro    (double)SHR_CONST_TKFRZ  // 0 degrees C                       ~ kelvin
#define Tfrz    ((double)Tzro   - 1.8)    // temp of saltwater freezing ~ kelvin
// real(SHR_KIND_R8), parameter ::  pi      (double)SHR_CONST_PI     // a famous math constant
#define pi      (double)4.0*atan(1.0)
// real(SHR_KIND_R8), parameter ::  omega   (double)SHR_CONST_OMEGA  // earth's rotation  ~ rad/sec
#define omega   (double)0.7292e-4
// real(SHR_KIND_R8), parameter ::  g       (double)SHR_CONST_G      // gravity ~ m/s^2
#define g       (double)9.806e0
#define DEGtoRAD (double)(pi/180.0)         // pi/180
#define RADIUS  (double)6371000e0
#define very_small  (double)  1.e-15
#define karman (double)0.4e0


#define c0 (double)0.0e0
#define  c1          1.0     
#define  c2          2.0     
#define  c3          3.0     
#define  c4          4.0     
#define  c5          5.0     
#define  c8          8.0     
#define  c10        10.0     
#define  c16        16.0     
#define  c1000   1000.0     
#define  c10000  10000.0     
#define  c1p5       1.5     
#define  p33     (c1/c3)        
#define  p5      0.500      
#define  p25     0.250      
#define  p125    0.125      
#define  p001    0.001      
#define  eps     1.0e-10    
#define  eps2    1.0e-20    
#define  bignum  1.0e+30    
#define  pi2  c2*pi


#define    undefined     -12345.0
#define    undefined_nf  NF90_FILL_DOUBLE
#define    undefined_nf_int  NF90_FILL_INT
#define    undefined_nf_r4   NF90_FILL_FLOAT
#define    undefined_nf_r8   NF90_FILL_DOUBLE

//   !*** location of fields for staggered grids
#define   field_loc_unknown    0
#define   field_loc_noupdate  -1
#define   field_loc_center     1
#define   field_loc_SWcorner   2
#define   field_loc_Sface      3
#define   field_loc_Wface      4

//    !*** field type attribute - necessary for handling
//    !*** changes of direction across tripole boundary

//    integer (int_kind), parameter, public ::   &
#define    field_type_unknown    0
#define    field_type_noupdate  -1
#define    field_type_scalar     1
#define    field_type_vector     2
#define    field_type_angle      3

//   character (5), parameter, public :: &
#define   blank_fmt  "(' ')"

extern void constant_mod_mp_const_();

#endif // //INCLUDE_CONSTANT
