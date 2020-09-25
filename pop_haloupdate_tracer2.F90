SUBROUTINE pop_haloupdate_tracer2(errorCode)
use tracer_mod
use dyn_mod
use domain
use POP_HaloMod
! use POP_FieldMod
use POP_GridHorzMod
!LPF20160515
integer::errorCode
call POP_HaloUpdate(VTL , POP_haloClinic, POP_gridHorzLocCenter,&
    POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)
END SUBROUTINE