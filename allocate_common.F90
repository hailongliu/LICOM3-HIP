
subroutine allocate_common
use forc_mod
use tracer_mod
use dyn_mod
use work_mod
use grid
use pconst_mod
use canuto_mod

!      PARAMETER (mt=1462)
! forc_mod
      allocate(tsf(imt,jmt,max_blocks_clinic))
      tsf=0.
      allocate(ssf(imt,jmt,max_blocks_clinic))
      ssf=0.
      allocate(sss(imt,jmt,max_blocks_clinic))
      sss=0.
      allocate(su(imt,jmt,max_blocks_clinic))
      su=0.
      allocate(sv(imt,jmt,max_blocks_clinic))
      sv=0.
      allocate(swv(imt,jmt,max_blocks_clinic))
      swv=0.
      allocate(nswv(imt,jmt,max_blocks_clinic))
      nswv=0.
      allocate(psa(imt,jmt,max_blocks_clinic))
      psa=0.
      allocate(ustar(imt,jmt,max_blocks_clinic))
      ustar=0.
      allocate(restore(imt,jmt,km,ntra,max_blocks_clinic))
      restore=0.
!tracer_mod
      allocate(licomqice(imt,jmt,max_blocks_clinic))
      licomqice=0.
      allocate(net(imt,jmt,NTRA,max_blocks_clinic))
      net=0.
      allocate(atb(imt,jmt,0:km,NTRA,max_blocks_clinic))
      atb=0.
      allocate(amld(imt,jmt,max_blocks_clinic))
      amld=0.
!dyn_mod
      allocate(utl(imt,jmt,km,max_blocks_clinic))
      utl=0.
      allocate(vtl(imt,jmt,km,max_blocks_clinic))
      vtl=0.
      allocate(utf(imt,jmt,km,max_blocks_clinic))
      utf=0.
      allocate(vtf(imt,jmt,km,max_blocks_clinic))
      vtf=0.
      allocate(up(imt,jmt,km,max_blocks_clinic))
      up=0.
      allocate(vp(imt,jmt,km,max_blocks_clinic))
      vp=0.
      allocate(ws(imt,jmt,kmp1,max_blocks_clinic))
      ws=0.

      allocate(ub(imt,jmt,max_blocks_clinic))
      ub=0.
      allocate(vb(imt,jmt,max_blocks_clinic))
      vb=0.
      allocate(ubp(imt,jmt,max_blocks_clinic))
      ubp=0.
      allocate(vbp(imt,jmt,max_blocks_clinic))
      vbp=0.

      allocate(h0p(imt,jmt,max_blocks_clinic))
      h0p=0.
      allocate(h0l(imt,jmt,max_blocks_clinic))
      h0l=0.
      allocate(h0f(imt,jmt,max_blocks_clinic))
      h0f=0.
      allocate(h0bl(imt,jmt,max_blocks_clinic))
      h0bl=0.
      allocate(h0bf(imt,jmt,max_blocks_clinic))
      h0bf=0.
      allocate(SBCX(imt,jmt,max_blocks_clinic))
      SBCX=0.
      allocate(BBCX(imt,jmt,max_blocks_clinic))
      BBCX=0.
      allocate(SBCY(imt,jmt,max_blocks_clinic))
      SBCY=0.
      allocate(BBCY(imt,jmt,max_blocks_clinic))
      BBCY=0.
!work_mod
      allocate(work(imt,jmt,max_blocks_clinic))
      work=0.
      allocate(wka(imt,jmt,km,max_blocks_clinic))
      wka=0.
!grid_mod
      allocate(dxu(nx_block,ny_block,max_blocks_clinic))
      dxu=0.
      allocate(dyu(nx_block,ny_block,max_blocks_clinic))
      dyu=0.
      allocate(dxt(nx_block,ny_block,max_blocks_clinic))
      dxt=0.
      allocate(dyt(nx_block,ny_block,max_blocks_clinic))
      dyt=0.
      allocate(dxur(nx_block,ny_block,max_blocks_clinic))
      dxur=0.
      allocate(dyur(nx_block,ny_block,max_blocks_clinic))
      dyur=0.
      allocate(dxtr(nx_block,ny_block,max_blocks_clinic))
      dxtr=0.
      allocate(dytr(nx_block,ny_block,max_blocks_clinic))
      dytr=0.
      allocate(hts(nx_block,ny_block,max_blocks_clinic))
      hts=0.
      allocate(htw(nx_block,ny_block,max_blocks_clinic))
      htw=0.
      allocate(hun(nx_block,ny_block,max_blocks_clinic))
      hun=0.
      allocate(hue(nx_block,ny_block,max_blocks_clinic))
      hue=0.
      allocate(ulat(nx_block,ny_block,max_blocks_clinic))
      ulat=0.
      allocate(ulon(nx_block,ny_block,max_blocks_clinic))
      ulon=0.
      allocate(tlat(nx_block,ny_block,max_blocks_clinic))
      tlat=0.
      allocate(tlon(nx_block,ny_block,max_blocks_clinic))
      tlon=0.
      allocate(angle(nx_block,ny_block,max_blocks_clinic))
      angle=0.
      allocate(anglet(nx_block,ny_block,max_blocks_clinic))
      anglet=0.
      allocate(fcor(nx_block,ny_block,max_blocks_clinic))
      fcor=0.
      allocate(fcort(nx_block,ny_block,max_blocks_clinic))
      fcort=0.
      allocate(uarea(nx_block,ny_block,max_blocks_clinic))
      uarea=0.
      allocate(tarea(nx_block,ny_block,max_blocks_clinic))
      tarea=0.
      allocate(uarea_r(nx_block,ny_block,max_blocks_clinic))
      uarea_r=0.
      allocate(tarea_r(nx_block,ny_block,max_blocks_clinic))
      tarea_r=0.
      allocate(ht(nx_block,ny_block,max_blocks_clinic))
      ht=0.
      allocate(hu(nx_block,ny_block,max_blocks_clinic))
      hu=0.
      allocate(hur(nx_block,ny_block,max_blocks_clinic))
      hur=0.
      allocate(kmt(nx_block,ny_block,max_blocks_clinic))
      kmt=0
      allocate(kmu(nx_block,ny_block,max_blocks_clinic))
      kmu=0
      allocate(kmtold(nx_block,ny_block,max_blocks_clinic))
      kmtold=0.

      allocate(kmtn(nx_block,ny_block,max_blocks_clinic))
      kmtn=0.
      allocate(kmts(nx_block,ny_block,max_blocks_clinic))
      kmts=0.
      allocate(kmte(nx_block,ny_block,max_blocks_clinic))
      kmte=0.
      allocate(kmtw(nx_block,ny_block,max_blocks_clinic))
      kmtw=0.
      allocate(kmun(nx_block,ny_block,max_blocks_clinic))
      kmun=0.
      allocate(kmus(nx_block,ny_block,max_blocks_clinic))
      kmus=0.
      allocate(kmue(nx_block,ny_block,max_blocks_clinic))
      kmue=0.
      allocate(kmuw(nx_block,ny_block,max_blocks_clinic))
      kmuw=0.

      allocate(at0(nx_block,ny_block,max_blocks_clinic))
      at0=0.
      allocate(atn(nx_block,ny_block,max_blocks_clinic))
      atn=0.
      allocate(ate(nx_block,ny_block,max_blocks_clinic))
      ate=0.
      allocate(atne(nx_block,ny_block,max_blocks_clinic))
      atne=0.
      allocate(au0(nx_block,ny_block,max_blocks_clinic))
      au0=0.
      allocate(aus(nx_block,ny_block,max_blocks_clinic))
      aus=0.
      allocate(auw(nx_block,ny_block,max_blocks_clinic))
      auw=0.
      allocate(ausw(nx_block,ny_block,max_blocks_clinic))
      ausw=0.

!pconst_mod
      allocate(ebea(imt,jmt,max_blocks_clinic))
      ebea=0.
      allocate(ebeb(imt,jmt,max_blocks_clinic))
      ebeb=0.
      allocate(ebla(imt,jmt,max_blocks_clinic))
      ebla=0.
      allocate(eblb(imt,jmt,max_blocks_clinic))
      eblb=0.
      allocate(epea(imt,jmt,max_blocks_clinic))
      epea=0.
      allocate(epeb(imt,jmt,max_blocks_clinic))
      epeb=0.
      allocate(epla(imt,jmt,max_blocks_clinic))
      epla=0.
      allocate(eplb(imt,jmt,max_blocks_clinic))
      eplb=0.
      allocate(ohbt(imt,jmt,max_blocks_clinic))
      ohbt=0.
      allocate(ohbu(imt,jmt,max_blocks_clinic))
      ohbu=0.
      allocate(dzph(imt,jmt,max_blocks_clinic))
      dzph=0.
      allocate(hbx(imt,jmt,max_blocks_clinic))
      hbx=0.
      allocate(hby(imt,jmt,max_blocks_clinic))
      hby=0.
      allocate(snlat(imt,jmt,max_blocks_clinic))
      snlat=0.
      allocate(vit(imt,jmt,km,max_blocks_clinic))
      vit=0.
      allocate(viv(imt,jmt,km,max_blocks_clinic))
      viv=0.
      allocate(fz_tide(imt,jmt,km,max_blocks_clinic))
      fz_tide=0.
      allocate(ak_tide(imt,jmt,km,max_blocks_clinic))
      ak_tide=0.
!canuto_mod
!/* when larger than 2GB error occurs
!      allocate(slq2b_h(2*mt+1,2*mt+1))
!      allocate(smb_h(2*mt+1,2*mt+1))
!      allocate(shb_h(2*mt+1,2*mt+1))
!      allocate(ssb_h(2*mt+1,2*mt+1))
!*/
end subroutine allocate_common
