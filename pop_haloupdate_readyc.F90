SUBROUTINE pop_haloupdate_readyc(errorCode)
use tracer_mod
use domain
use POP_HaloMod
! use POP_FieldMod
use POP_GridHorzMod
!LPF20160515
integer::errorCode
!    call POP_HaloUpdate2DR4(amld , POP_haloClinic, POP_gridHorzLocCenter,&
!                          POP_fieldKindScalar, errorCode, 0.0)
    call POP_HaloUpdate(amld , POP_haloClinic, POP_gridHorzLocCenter,&
                         POP_fieldKindScalar, errorCode, 0.0)
END SUBROUTINE