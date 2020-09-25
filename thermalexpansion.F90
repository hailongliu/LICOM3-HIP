       subroutine thermal(tt,ss,pp,aa,bb,mask)
!
! calculate the thermal expansion coefficient and
! the saline contraction coefficient using McDougall (1987) equations
!
! input: tt  potential temperature in C
!        ss  salinity in (s-35)/1000
!        p  pressure in Pa
! output:alpha thermal expansion in 1/C
!        beta saline contraction in 1/psu
! 
       use precision_mod
       use param_mod
       use pconst_mod
       use domain
      IMPLICIT NONE
       real(r8),dimension(imt,jmt,km,max_blocks_clinic) :: tt,ss,pp,tmp,aa,bb,mask
       integer :: iblock
!
!  transform the unit
!  ss from (s-35.0)/1000 to s-35.0
!  pp from Pa to db
!
!       ss=ss*1000.0
!       p=p/10000.0
!
      DO IBLOCK  = 1, NBLOCKS_CLINIC
!
      DO K = 1,KM
         DO J = 1, JMT
            DO I = 2,IMT-1
            tmp(i,j,k,iblock)=-pp(i,j,k,iblock)/OD0/10000.D0*mask(i,j,k,iblock)
            END DO
         END DO   
      END DO
!
      DO K = 1,KM
         DO J = 1, JMT
            DO I = 2,IMT-1
       bb(i,j,k,iblock)=(0.785567d-3-0.301985d-5*tt(i,j,k,iblock)+0.555579d-7*tt(i,j,k,iblock)**2 &
           -0.415613d-9*tt(i,j,k,iblock)**3+(ss(i,j,k,iblock)*1000.0d0)*(-0.356603d-6+0.788212d-8*tt(i,j,k,iblock)&
           +0.408195d-10*tmp(i,j,k,iblock)-0.602281d-15*tmp(i,j,k,iblock)**2) &
           +(ss(i,j,k,iblock)*1000.0d0)**2*(0.515032d-8)+tmp(i,j,k,iblock)*(-0.121555d-7 & 
           +0.192867d-9*tt(i,j,k,iblock)-0.2131127d-11*tt(i,j,k,iblock)**2) &
           +tmp(i,j,k,iblock)**2*(0.176621d-12-0.175379d-14*tt(i,j,k,iblock)) &
           +tmp(i,j,k,iblock)**3*(0.12155d-17))*mask(i,j,k,iblock)
!
       aa(i,j,k,iblock)=(0.665157d-1+0.170907d-1*tt(i,j,k,iblock)-0.203814d-3*tt(i,j,k,iblock)**2&
            +0.298357d-5*tt(i,j,k,iblock)**3-0.255019d-7*tt(i,j,k,iblock)**4  & 
            +(ss(i,j,k,iblock)*1000.0d0)*(0.378110d-2 -0.846960d-4*tt(i,j,k,iblock) & 
            -0.164759d-6*tmp(i,j,k,iblock) - 0.251520d-11*tmp(i,j,k,iblock)**2)&
            +(ss(i,j,k,iblock)*1000.0d0)**2*(-0.678662d-5)+tmp(i,j,k,iblock)*(0.380374d-4  & 
            -0.933746d-6*tt(i,j,k,iblock) +0.791325d-8*tt(i,j,k,iblock)**2) &  
            +0.512857d-12*tmp(i,j,k,iblock)**2*tt(i,j,k,iblock)**2&
            -0.302285d-13*tmp(i,j,k,iblock)**3)*bb(i,j,k,iblock)*mask(i,j,k,iblock)
!
            END DO
         END DO   
      END DO
!
      END DO

!
       return
       end subroutine thermal
