!  CVS: $Id: addps.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =================
      SUBROUTINE ADDPS
!     ================
!     COMPENSATING THE LOSS OF GROSS MASS
 
#include <def-undef.h>
use precision_mod
use param_mod
use pconst_mod
use dyn_mod
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn
use domain
use grid, only: tarea, area_t
      IMPLICIT NONE
      REAL(r8)    :: ERROR,DH00,error0
      integer     :: iblock, int_dh00
     
      ERROR = 0.0D0
 
!$OMP PARALLEL DO PRIVATE(IBLOCK,J,I), REDUCTION(+:error)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO J = 3, jmt-2
         DO I = 3, imt-2
            ERROR = ERROR + H0 (I,J,Iblock)* VIT (I,J,1,iblock)*tarea(i,j,iblock)
         END DO
      END DO
   END DO
       call mpi_reduce(error,error0,1,MPI_PR,mpi_sum,0,mpi_comm_ocn,ierr)
       call mpi_bcast(error0,1,MPI_PR,0,mpi_comm_ocn,ierr)
      dh00 = - error0/ area_t
      int_dh00 = int(dh00*1.0D22)
      if ( mytid==0 ) write(*,*) "Mass Imbalance : " , int_dh00,int_dh00*1.0D-22

!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO J = 1, JMT
         DO I = 1,IMT
            H0 (I,J,IBLOCK)= (H0 (I,J,IBLOCK) + dble(int_dh00)*1.0D-22)* VIT (I,J,1,IBLOCK)
         END DO
      END DO
   END DO
 
      RETURN
      END SUBROUTINE ADDPS
 
