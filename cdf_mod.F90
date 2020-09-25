!  CVS: $Id: cdf_mod.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =================
module cdf_mod
#include <def-undef.h>
use precision_mod
use param_mod
!     this head file includes variables related to Netcdf output.
!     written by liu hai long 2001, jun 
#if (defined NETCDF) || (defined ALL)
!     error status return
      integer:: iret

!     file id
      integer:: ncid

!     dimension id
      integer:: lon_dim,lat_dim,lev1_dim,lev_dim,time_dim, basin_dim,y_dim, lat_aux_dim,x_dim, tracer_dim!ZWP2013-10-17

!     dimension lenth
      integer, parameter :: lon_len=imt_global,lat_len=jmt_global,lev_len=km,lev1_len=km+1,time_len=1, basin_len=2 ! ZWP2013-10-17

!     variable id
#if (defined SMAG_OUT)
      integer, parameter :: am_rank=4
      integer:: am_id,am_rank,am_dims(am_rank)
#endif
      integer::  lat_id,lon_id,ulon_id,ulat_id,lev1_id,lev_id,time_id,z0_id,hi_id,hd_id, lat_aux_id, &
               ic1_id,ic2_id, net1_id, net2_id, ts_id,  ss_id,us_id,vs_id,ws_id, psi_euler_id, psi_eddy_id, bsf_id,&
               mth_adv_id,mth_adv_iso_id, mth_dif_id, mld_id,akm_id,akt_id,aks_id &
              ,qice_id,su_id,sv_id,sshf_id,lthf_id,lwv_id,swv_id, fresh_id, runoff_id,icefrac_id &
             ,h0_id,u_id,v_id,t_id,s_id,h1_id,h2_id,h3_id,month_id, day_id, basin_id, area_id,masksurf_id
#ifdef ISOOUT

     integer:: isox_id,isoy_id,isoz_id
     integer, parameter ::isox_rank=4 ,isoy_rank=4 ,isoz_rank=4
     integer:: isox_dims(isox_rank) , isoy_dims(isoy_rank) , isoz_dims(isoz_rank)
#endif

!#ifdef BUDGETOUT
     integer:: dt_convid,dtdiff_id,xadv_id,yadv_id,zadv_id,xadviso_id,yadviso_id,zadviso_id,&
               pen_id,tend_id,xdifiso_id,ydifiso_id,zdifiso_id,xydif_id,zdif_id
     integer:: dt_convssid,dtdiff_ss_id,sss_id,xadv_ss_id,yadv_ss_id,zadv_ss_id,xadviso_ss_id,yadviso_ss_id,&
              zadviso_ss_id,tend_ss_id,xdifiso_ss_id,ydifiso_ss_id,zdifiso_ss_id,xydif_ss_id,zdif_ss_id
     integer, parameter ::xadv_rank=4,yadv_rank=4,zadv_rank=4,tend_rank=4,&
                          xadviso_rank=4 ,yadviso_rank=4,zadviso_rank=4,pen_rank=4,&
                         xdifiso_rank=4,ydifiso_rank=4,zdifiso_rank=4,xydif_rank=4,zdif_rank=4
     integer:: xadv_dims(xadv_rank) , yadv_dims(yadv_rank),zadv_dims(zadv_rank),tend_dims(tend_rank),&
               xadviso_dims(xadviso_rank),yadviso_dims(yadviso_rank) , zadviso_dims(zadviso_rank),pen_dims(pen_rank),&
              xdifiso_dims(xdifiso_rank),ydifiso_dims(ydifiso_rank),zdifiso_dims(zdifiso_rank),xydif_dims(xydif_rank),zdif_dims(zdif_rank) 
!#endif

!     variable rank
!YU
!ZWP      integer, parameter ::lat_rank=1,lon_rank=1,lev_rank=1,time_rank=1,z0_rank=3,hi_rank=3, &
      integer, parameter ::lev_rank=1,time_rank=1,z0_rank=3,hi_rank=3,qice_rank=3,  &
                            hd_rank=3,ic1_rank=3,ic2_rank=3, net1_rank=3,net2_rank=3,ts_rank=4,ss_rank=4, &
                            us_rank=4, vs_rank=4, ws_rank=4, psi_euler_rank=4,psi_eddy_rank=4 , lev1_rank=1, bsf_rank=3,&
               mth_adv_rank=4,mth_adv_iso_rank=4, mth_dif_rank=4, mld_rank=3,akm_rank=4,akt_rank=4,aks_rank=4 &
              ,su_rank=3,sv_rank=3,sshf_rank=3,lthf_rank=3,lwv_rank=3,swv_rank=3&
             ,h0_rank=2,u_rank=3,v_rank=3,t_rank=3,s_rank=3,h1_rank=2,h2_rank=2,h3_rank=2,month_rank=1, day_rank=1,&
              basin_rank=1,ulat_rank=2,ulon_rank=2,lat_rank=2,lon_rank=2, lat_aux_rank=1, &
              area_rank=2,masksurf_rank=2, fresh_rank=3, runoff_rank=3,icefrac_rank=3 !zwp2013-10-17
      integer, parameter ::richardson_rank=4,fztidal_rank=4,wp3_tidal_rank=4     !yuzp-2016/11/13

      integer, parameter ::ak_tide1_rank=4     !yuzp-2016/11/19

      integer, parameter ::wp1_canuto_rank=4,wp2_canuto_rank=4,wp3_canuto_rank=4,wp4_canuto_rank=4,wp5_canuto_rank=4,&
                           wp6_canuto_rank=4,wp7_canuto_rank=4,wp8_canuto_rank=4     !yuzp-2016/12/4
      integer, parameter ::wk1_canuto_rank=4,wk2_canuto_rank=4,wk3_canuto_rank=4,wk4_canuto_rank=4     !yuzp-2016/12/4
      integer, parameter ::fcor_canuto_rank=3,fcort_canuto_rank=3     !yuzp-2016/12/4
      integer, parameter ::wp12_canuto_rank=4,wp13_canuto_rank=4     !yuzp-2016/12/4
      integer, parameter ::wp10_canuto_rank=3,wp11_canuto_rank=3    !yuzp-2016/12/8
      integer, parameter ::alpha_canuto_rank=4,beta_canuto_rank=4     !yuzp-2016/12/4

!
!     variable shapes
!ZWP      integer::  lat_dims(lat_rank),lon_dims(lon_rank),lev_dims(lev_rank),time_dims(time_rank), &
      integer::  lat_dims(lat_rank),lon_dims(lon_rank),ulat_dims(ulat_rank),ulon_dims(ulon_rank),&
                 lev_dims(lev_rank),time_dims(time_rank), &
                z0_dims( z0_rank), hi_dims( hi_rank), hd_dims( hd_rank), ic1_dims( ic1_rank), &
               ic2_dims(ic2_rank),net1_dims(net1_rank), net2_dims(net2_rank),ts_dims( ts_rank), ss_dims( ss_rank), &
                us_dims( us_rank), vs_dims( vs_rank), ws_dims( ws_rank),psi_euler_dims(psi_euler_rank), psi_eddy_dims(psi_eddy_rank), &
              lev1_dims(lev1_rank), bsf_dims(bsf_rank),qice_dims(qice_rank),lat_aux_dims(lat_aux_rank),   &
              mth_adv_dims(mth_adv_rank), mth_adv_iso_dims(mth_adv_iso_rank), mth_dif_dims(mth_dif_rank), mld_dims(mld_rank),akm_dims(akm_rank), &
             akt_dims(akt_rank),aks_dims(aks_rank), day_dims(day_rank), month_dims(month_rank) &
             ,su_dims(su_rank),sv_dims(sv_rank), area_dims(area_rank),masksurf_dims(masksurf_rank), fresh_dims(fresh_rank), runoff_dims(runoff_rank),icefrac_dims(icefrac_rank) &
             ,sshf_dims(sshf_rank),lthf_dims(lthf_rank),lwv_dims(lwv_rank),swv_dims(swv_rank)&
             ,h0_dims(h0_rank),u_dims(u_rank),v_dims(v_rank),t_dims(t_rank),s_dims(s_rank)&
             ,h1_dims(h1_rank),h2_dims(h2_rank),h3_dims(h3_rank),basin_dims(basin_rank)!ZWP2013-10-17
      integer::richardson_dims(richardson_rank),fztidal_dims(fztidal_rank),wp3_tidal_dims(wp3_tidal_rank)     !yuzp-2016/11/13
      integer::ak_tide1_dims(ak_tide1_rank)     !yuzp-2016/11/19

      integer::wp1_canuto_dims(wp1_canuto_rank),wp2_canuto_dims(wp2_canuto_rank),&
               wp3_canuto_dims(wp3_canuto_rank),wp4_canuto_dims(wp4_canuto_rank),&
               wp5_canuto_dims(wp5_canuto_rank),wp6_canuto_dims(wp6_canuto_rank),&
               wp7_canuto_dims(wp7_canuto_rank),wp8_canuto_dims(wp8_canuto_rank)     !yuzp-2016/12/4

      integer::wk1_canuto_dims(wk1_canuto_rank),wk2_canuto_dims(wk2_canuto_rank),&
               wk3_canuto_dims(wk3_canuto_rank),wk4_canuto_dims(wk4_canuto_rank)     !yuzp-2016/12/4

      integer::wp10_canuto_dims(wp10_canuto_rank),wp11_canuto_dims(wp11_canuto_rank)     !yuzp-2016/12/8
      integer::wp12_canuto_dims(wp12_canuto_rank),wp13_canuto_dims(wp13_canuto_rank)     !yuzp-2016/12/4

      integer::alpha_canuto_dims(alpha_canuto_rank),beta_canuto_dims(beta_canuto_rank)     !yuzp-2016/12/4

      integer::fcor_canuto_dims(fcor_canuto_rank),fcort_canuto_dims(fcort_canuto_rank)     !yuzp-2016/12/4

      real(r8) t0_cdf 
      real(r4),allocatable,dimension(:,:,:):: buffer_r4_local !linpf 2012Jul27
      real(r4),allocatable,dimension(:,:):: buffer_r4_global !linpf 2012Jul27
!     variables
!    start and count
      integer:: start1(1),count1(1)
      integer:: start2(2),count2(2)
      integer:: start3(3),count3(3)
      integer:: start4(4),count4(4)
#endif
 
    integer:: tideeg_id
    integer:: richardson_id,fztidal_id,wp3_tidal_id     !yuzp-2016/11/13
    integer:: ak_tide1_id     !yuzp-2016/11/19

    integer:: wp1_canuto_id,wp2_canuto_id,wp3_canuto_id,wp4_canuto_id,&
              wp5_canuto_id,wp6_canuto_id,wp7_canuto_id,wp8_canuto_id     !yuzp-2016/12/4

    integer:: wk1_canuto_id,wk2_canuto_id,wk3_canuto_id,wk4_canuto_id     !yuzp-2016/12/4

    integer:: fcor_canuto_id,fcort_canuto_id     !yuzp-2016/12/4

    integer:: wp10_canuto_id,wp11_canuto_id     !yuzp-2016/12/4

    integer:: wp12_canuto_id,wp13_canuto_id     !yuzp-2016/12/4

    integer:: alpha_canuto_id,beta_canuto_id     !yuzp-2016/12/4

    integer, parameter :: tideeg_rank=2
    integer:: tideeg_dims(tideeg_rank)


#if (defined TIDEMIX)
      integer:: aktide_id
      integer, parameter ::aktide_rank=4
      integer:: aktide_dims(aktide_rank)
#endif
#if (defined ISO_TYPE_BF)
      integer:: athkdf_id
      integer, parameter ::athkdf_rank=4
      integer:: athkdf_dims(athkdf_rank)
#endif
end module cdf_mod
