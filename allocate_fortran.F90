subroutine allocate_fortran()
use pmix_mod
use dyn_mod
use work_mod
use tracer_mod
use forc_mod

      if (.not. allocated(stf)) then
       allocate(stf(imt,jmt,max_blocks_clinic))
       stf=0.
      end if
      if (.not. allocated(wave_dis)) then
       allocate(wave_dis(imt,jmt,max_blocks_clinic))
       wave_dis=0.
      end if
      if (.not. allocated(tf)) then
       allocate(tf(imt,jmt,km,max_blocks_clinic))
       tf=0.
      end if
      if (.not. allocated(wkb)) then
       allocate(wkb(imt,jmt,km,max_blocks_clinic))
       wkb=0.
      end if
      if (.not. allocated(wkc)) then
       allocate(wkc(imt,jmt,km,max_blocks_clinic))
       wkc=0.
      end if
      if (.not. allocated(wkd)) then
       allocate(wkd(imt,jmt,km,max_blocks_clinic))
       wkd=0.
      end if


    if(.not. allocated(dlu)) then
     allocate(dlu(imt,jmt,km,max_blocks_clinic))
     dlu=0.
    end if
    if(.not. allocated(dlv)) then
     allocate(dlv(imt,jmt,km,max_blocks_clinic))
     dlv=0.
    end if
if(.not. allocated(gg)) allocate(gg(imt,jmt,km,max_blocks_clinic))
if(.not. allocated(u)) allocate(u(imt,jmt,km,max_blocks_clinic))
if(.not. allocated(v)) allocate(v(imt,jmt,km,max_blocks_clinic))
if(.not. allocated(rit)) allocate(rit(imt,jmt,kmm1,max_blocks_clinic))
if(.not. allocated(ric)) allocate(ric(imt,jmt,kmm1,max_blocks_clinic))
if(.not. allocated(rict)) allocate(rict(imt,jmt,kmm1,max_blocks_clinic))
if(.not. allocated(rict_replace)) allocate(rict_replace(imt,jmt,kmm1,max_blocks_clinic))

if(.not. allocated(riu)) allocate(riu(imt,jmt,0:km,max_blocks_clinic),stat=ierr) 
if(.not. allocated(rict_ref)) allocate(rict_ref(imt,jmt,max_blocks_clinic))
if(.not. allocated(dlub)) allocate(dlub(imt,jmt,max_blocks_clinic),stat=ierr)
if(.not. allocated(dlvb)) allocate(dlvb(imt,jmt,max_blocks_clinic),stat=ierr)

end subroutine allocate_fortran
