SUBROUTINE pop_haloupdate_barotr1(errorCode)
use work_mod
use tracer_mod
use dyn_mod
use domain
use POP_HaloMod
! use POP_FieldMod
use POP_GridHorzMod
use precision_mod
use param_mod
!LPF20160515
integer::errorCode
         call POP_HaloUpdate(work , POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)
END SUBROUTINE