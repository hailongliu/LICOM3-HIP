SUBROUTINE pop_haloupdate_smuv_3d(X, errorCode)
use domain
use precision_mod
use param_mod
use POP_HaloMod
use POP_GridHorzMod
integer::errorCode
REAL(r8)    :: X (IMT,JMT,KM,max_blocks_clinic)
         call POP_HaloUpdate(X, POP_haloClinic, POP_gridHorzLocSWcorner , &
                       POP_fieldKindVector, errorCode, fillValue = 0.0_r8)
END SUBROUTINE