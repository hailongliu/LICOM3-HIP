!  CVS: $Id: energy.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =================
      SUBROUTINE ENERGY
!     =================

#include <def-undef.h>
use precision_mod
use param_mod
use pconst_mod
use tracer_mod
use dyn_mod
use work_mod
use grid
use constant_mod
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn
use shr_sys_mod
use domain
use blocks
!use constant_mod

      IMPLICIT NONE

      real (r8) :: ek0,ea0,eb0,et0,es0,volume, volume0, nnp, nnp0
      REAL (r8) :: EK,EA,EB,ET,ES, dxx
      real(r8)::t0,t1,clock0f,nst
      save t0
      integer :: nnn,iblock !,month1 !LPF20200315

      type (block) :: this_block

!---------------------------------------------------------------------
!     TOTAL K.E.
!---------------------------------------------------------------------

      EK = 0.0
      month=(iyfm-1)*12+mon0

      EB = 0.0
      nst=0

      EK = 0.0
!      month1=(iyfm-1)*12+mon0
!      if (mytid==0) write(*,*)'month1=,',month1,iyfm,mon0

!!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I),reduction(+:EK)
      do iblock= 1, nblocks_clinic
      this_block = get_block(blocks_clinic(iblock),iblock)
      DO K = 1,KM
         DO J = this_block%jb, this_block%je
            DO I = this_block%ib, this_block%ie
               EK = EK + DZP (K)* uarea(i,j,iblock)* VIV (I,J,K,iblock) &
                    * (U (I,J,K,iblock)* U (I,J,K,iblock) + V (I,J,K,iblock)* V (I,J,K,iblock))
            END DO
         END DO
      END DO
      END DO
      call mpi_reduce(ek,ek0,1,MPI_PR,mpi_sum,0,mpi_comm_ocn,ierr)
      ek0= ek0*0.5D0* d0

!---------------------------------------------------------------------
!     AVAILABLE P.E.
!---------------------------------------------------------------------

      EA = 0.0
!!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I),reduction(+:EA)
      do iblock= 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)
         DO J = this_block%jb, this_block%je
            DO I = this_block%ib, this_block%ie
            EA = EA + tarea(i,j,iblock)* VIT (I,J,1,iblock)* H0 (I,J,iblock)* H0 (I,J,iblock)
         END DO
      END DO
      END DO

      call mpi_reduce(ea,ea0,1,MPI_PR,mpi_sum,0,mpi_comm_ocn,ierr)
      ea0= ea0*0.5D0*d0*g

      EB = 0.0
      nst=0
!!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I),reduction(+:EB)
      do iblock= 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)
         DO J = this_block%jb, this_block%je
            DO I = this_block%ib, this_block%ie
            EB = EB + tarea(i,j,iblock)* VIT (I,J,1,iblock)* AT (I,J,1,1,iblock)
         END DO
      END DO
      end do

      call mpi_reduce(eb,eb0,1,MPI_PR,mpi_sum,0,mpi_comm_ocn,ierr)
      eb0= eb0/area_t

!---------------------------------------------------------------------
!     GLOBAL MEAN TEMPERATURE & SALINITY
!---------------------------------------------------------------------

      ET = 0.0D0
      ES = 0.0D0
!!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I),reduction(+:ET,ES)
      DO iblock= 1, nblocks_clinic
      this_block = get_block(blocks_clinic(iblock),iblock)
      DO K = 1,km
         DO J = this_block%jb, this_block%je
            DO I = this_block%ib, this_block%ie
               ET = ET + DZP (K)*tarea(i,j,iblock)* VIT (I,J,K,iblock)* AT(I,J,K,1,iblock)
               ES = ES + DZP (K)*tarea(i,j,iblock)* VIT (I,J,K,iblock)* AT(I,J,K,2,iblock)
            END DO
         END DO
      END DO
      END DO

      call mpi_reduce(et,et0,1,MPI_PR,mpi_sum,0,mpi_comm_ocn,ierr)
      call mpi_reduce(es,es0,1,MPI_PR,mpi_sum,0,mpi_comm_ocn,ierr)
      et0= et0/ volume_t
      es0= es0/ volume_t *1.0D3

!---------------------------------------------------------------------
!     Calculating CPU time for per-day.
!---------------------------------------------------------------------
      if (mytid==0) then
         t1=clock0f()
      endif

      if (mytid==0 )then
      WRITE (16,FMT='(I5,I3,6D25.15)') MONTH,IDAY,EK0,EA0,EB0,ET0,ES0
!      call flush_(6)
      end if

      if (mytid==0) then
         t0=t1
      endif
!
      RETURN
      END SUBROUTINE ENERGY


      SUBROUTINE chk_var3d(var,ek0,a,kk)
!     =======================

#include <def-undef.h>
use precision_mod
use param_mod
use pconst_mod
use tracer_mod
use dyn_mod
use work_mod
use grid
use constant_mod
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn
use shr_sys_mod
use domain
use blocks

      IMPLICIT NONE

      real (r8) :: ek0,ea0,eb0,et0,es0,volume, volume0, nnp, nnp0
      REAL (r8) :: EK,EA,EB,ET,ES
      real(r8)::t0,t1,clock0f,nst
      save t0
      integer :: nnn,iblock, a,kk
      real(r8) :: var(imt,jmt,kk,max_blocks_clinic), var1(imt,jmt,kk,max_blocks_clinic)
      type (block) :: this_block

      ek = 0.0D0
!!$OMP PARALLEL DO PRIVATE (iblock,k,J,I),reduction(+:EK)
      do iblock= 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)
         DO K = 1, kk
         DO J = this_block%jb, this_block%je
            DO I = this_block%ib, this_block%ie
              if ( a ==1 ) then
              if ( this_block%i_glob(i) > imt_global/2  .or. this_block%j_glob(j) /= 1 ) then
                ek = ek+tarea(i,j,iblock)*var(I,J,k,iblock)*dzp(k)*vit(i,j,k,iblock)
              end if
              else
                ek = ek+uarea(i,j,iblock)*var(I,J,k,iblock)*var(I,J,k,iblock)*dzp(k)*viv(i,j,k,iblock)
              end if
            END DO
          END DO
          END DO
      END DO

      call mpi_reduce(ek,ek0,1,MPI_PR,mpi_sum,0,mpi_comm_ocn,ierr)
      ek0= ek0/volume_t
!

      RETURN
      END SUBROUTINE chk_var3d

      SUBROUTINE chk_var2d(var,ek0,a)
!     =======================
#include <def-undef.h>
use precision_mod
use param_mod
use pconst_mod
use tracer_mod
use dyn_mod
use work_mod
use grid
use constant_mod
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn
use shr_sys_mod
use domain
use blocks

      IMPLICIT NONE

      real (r8) :: ek0,ea0,eb0,et0,es0,volume, volume0, nnp, nnp0
      REAL (r8) :: EK,EA,EB,ET,ES
      real(r8)::t0,t1,clock0f,nst
      real(r8):: var(imt,jmt,max_blocks_clinic)
      save t0
      integer :: nnn,iblock, a

      type (block) :: this_block


      ek = 0.0D0
!!$OMP PARALLEL DO PRIVATE (J,I),reduction(+:EK)
      do iblock= 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)
         DO J = this_block%jb, this_block%je
            DO I = this_block%ib, this_block%ie
              if ( a == 1) then
                ek = ek + tarea(i,j,iblock)* VIT (I,J,1,iblock)* var (I,J,iblock)* var (I,J,iblock)
              else
                ek = ek + uarea(i,j,iblock)* VIV (I,J,1,iblock)* var (I,J,iblock)* var (I,J,iblock)
              end if
         END DO
      END DO
      END DO

      call mpi_reduce(ek,ek0,1,MPI_PR,mpi_sum,0,mpi_comm_ocn,ierr)
!
      RETURN
      END SUBROUTINE chk_var2d

      subroutine write_global(var,nfile)
!
      use precision_mod
      use param_mod
      use gather_scatter
      use domain
      use work_mod
!
      real(r8) :: var(imt,jmt,max_blocks_clinic)
      integer :: nfile
!
      allocate( work1_g(imt_global,jmt_global))
         call gather_global(work1_g, var, master_task,distrb_clinic)
         if (mytid ==0) write(nfile) ((work1_g(i,j), i=1,imt_global), j=4,jmt_global)
      deallocate( work1_g)
!
      end subroutine write_global
