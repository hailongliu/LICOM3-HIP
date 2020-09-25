!  CVS: $Id: yy00.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     ===============
      SUBROUTINE YY00
!     ===============
 
#include <def-undef.h>
use param_mod
use output_mod

      IMPLICIT NONE

#if (defined LOWRES)
            ERR_norm2mon=0.0D0 
            Z0MON = 0.0
            HIMON = 0.0
            HDMON = 0.0
            ICMON = 0.0
            ICMON = 0.0
 
               TSMON = 0.0
               SSMON = 0.0
               USMON = 0.0
               VSMON = 0.0
               WSMON = 0.0
#if (defined SMAG_OUT)
               AM3MON = 0.0
#endif

               tendmon  =0.0
               axmon    = 0.0
               aymon    = 0.0
               azmon    = 0.0
               dxmon    = 0.0
               dymon    = 0.0
               dzmon    = 0.0
               dt_diffmon=0.0
               dt_convmon=0.0
!              ddymon   = 0.0
#ifdef ISO
               axmon_iso  = 0.0
               aymon_iso  = 0.0
               azmon_iso  = 0.0
               dxmon_iso  = 0.0
               dymon_iso  = 0.0
               dzmon_iso  = 0.0
!              aaymon_iso = 0.0
!              ddymon_iso = 0.0
#if ( defined ISOOUT )
               vntisomon = 0.0
               vetisomon = 0.0
               vbtisomon = 0.0
#endif
#endif

 
     netmon = 0.0
     penmon = 0.0
     mldmon = 0.0
     ifracmon = 0.0
     akmmon = 0.0
     aktmon = 0.0
!lhl1204
#if ( defined TIDEMIX )
               aktidemon = 0.0
               wavedismon = 0.0

               richardsonmon=0.0     !yuzp-2016/11/13
               fztidalmon=0.0     !yuzp-2016/11/13
               wp3_tidalmon=0.0     !yuzp-2016/11/13
               ak_tide1mon=0.0     !yuzp-2016/11/19
#endif
#if ( defined CANUTOMIXOUT )
               wp1_canutomon=0.0     !yuzp-2016/12/4
               wp2_canutomon=0.0     !yuzp-2016/12/4
               wp3_canutomon=0.0     !yuzp-2016/12/4
               wp4_canutomon=0.0     !yuzp-2016/12/4
               wp5_canutomon=0.0     !yuzp-2016/12/4
               wp6_canutomon=0.0     !yuzp-2016/12/4
               wp7_canutomon=0.0     !yuzp-2016/12/4
               wp8_canutomon=0.0     !yuzp-2016/12/4
               wk1_canutomon=0.0     !yuzp-2016/12/4
               wk2_canutomon=0.0     !yuzp-2016/12/4
               wk3_canutomon=0.0     !yuzp-2016/12/4
               wk4_canutomon=0.0     !yuzp-2016/12/4
               wp10_canutomon=0.0     !yuzp-2016/12/4
               wp11_canutomon=0.0     !yuzp-2016/12/4
               fcor_canutomon=0.0     !yuzp-2016/12/4
               fcort_canutomon=0.0     !yuzp-2016/12/4
               wp12_canutomon=0.0     !yuzp-2016/12/4
               wp13_canutomon=0.0     !yuzp-2016/12/4
               alpha_canutomon=0.0     !yuzp-2016/12/4
               beta_canutomon=0.0     !yuzp-2016/12/4
#endif

#if (defined ISO_TYPE_BF)
               athkdfmon = 0.0
#endif

#endif
 
      RETURN
      END SUBROUTINE YY00
 
 
