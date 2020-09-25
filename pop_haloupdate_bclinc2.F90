SUBROUTINE pop_haloupdate_bclinc2(errorCode)
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
         call POP_HaloUpdate(U, POP_haloClinic, POP_gridHorzLocSWcorner , &
                       POP_fieldKindVector, errorCode,  0.0)
!
         call POP_HaloUpdate(V, POP_haloClinic, POP_gridHorzLocSWcorner , &
                       POP_fieldKindVector, errorCode,  0.0)
END SUBROUTINE