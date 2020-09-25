SUBROUTINE pop_haloupdate_smts(X, errorCode)
use domain
use precision_mod
use param_mod
use POP_HaloMod
use POP_GridHorzMod
integer::errorCode
REAL(r8)    :: X (IMT,JMT,KM,max_blocks_clinic)
         call POP_HaloUpdate(X , POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)
END SUBROUTINE