SUBROUTINE pop_haloupdate_barotr2(errorCode1, errorCode2)
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

integer::errorCode1
integer::errorCode2
!
    call POP_HaloUpdate(wka(:,:,3,:), POP_haloClinic, POP_gridHorzLocSWcorner , &
                        POP_fieldKindVector, errorCode1, fillValue = 0.0_r8)
!
    call POP_HaloUpdate(wka(:,:,4,:), POP_haloClinic, POP_gridHorzLocSWcorner , &
                        POP_fieldKindVector, errorCode2, fillValue = 0.0_r8)

END SUBROUTINE
