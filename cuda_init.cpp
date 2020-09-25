#include "hip/hip_runtime.h"
#include "cuda_data.h"
#include "param_mod.h"

#include "dyn_mod.h"
#include "grid.h"
#ifndef BIHAR
#include "hmix_del2.h"
#else
#include "hmix_del4.h"
#endif
#include "pmix_mod.h"
#include "forc_mod.h"
#include "tracer_mod.h"
#include "pconst_mod.h"
#include "buf_mod.h"
#include "work_mod.h"
#include "common.h"
#include "isopyc_mod.h"

//double isopyc_mod_mp_ahisop_[max_blocks_clinic][jmt][imt];
//double (*isopyc_mod_mp_k3_)[3][jmt][km + 1][imt];

extern "C" void allocate_fortran_();
size_t dataSizeI = sizeof(int);
size_t dataSizeD = sizeof(double);
size_t dataSizeC = sizeof(char);
size_t dataSizeB = sizeof(bool);
size_t dataSize1d = km * sizeof(double);
size_t dataSize11d = (km + 1) * sizeof(double);
size_t dataSize12d = (km - 1) * sizeof(double);
size_t dataSize111d = kmm1 * sizeof(double);
size_t dataSize1112d = kmp1 * sizeof(double);
size_t dataSize2d = jmt * imt * sizeof(double);
size_t dataSize22d = ny_block * nx_block * sizeof(double);
size_t dataSize222d = jmt_global * imt_global * sizeof(double);
size_t dataSize3d = km * jmt * imt * sizeof(double);
size_t dataSize33d = ntra * jmt * imt * sizeof(double);
size_t dataSize332d = (km + 1) * jmt * imt * sizeof(double);
size_t dataSize333d = kmp1 * jmt * imt * sizeof(double);
size_t dataSize3333d = kmm1 * jmt * imt * sizeof(double);
size_t dataSize4d = ntra * km * jmt * imt * sizeof(double);
size_t dataSize44d = ntra * (km + 1) * jmt * imt * sizeof(double);
size_t dataSize444d = 3 * jmt * (km + 1) * imt * sizeof(double);

size_t dataSize2i = ny_block * nx_block * sizeof(int);
//size_t dataSize1block = sizeof(struct block);

//struct block *d_thisblock = NULL;
//readyc turb
double *d_wk1 = NULL; //jjrtest
double *d_wk2 = NULL;
double *d_wk3 = NULL;
double *d_wp1 = NULL; //jjrtest
double *d_wp3 = NULL; //jjrtest
double *d_wp7 = NULL; //jjrtest
double *d_wp8 = NULL; //jjrtest

//advection_tracer
double *d_at00 = NULL;
double *d_atmax = NULL;
double *d_atmin = NULL;

//tracer
double *d_adv_tt = NULL;
double *d_ori = NULL;
double *d_temp11 = NULL;

//readyc
double *d_wp12 = NULL;
double *d_wp13 = NULL;
double *d_riv1 = NULL;
double *d_riv2 = NULL;
double *d_diff_back = NULL;
double *d_diff_back_sh = NULL;
double *d_diff_back_nh = NULL;

double *d_uk = NULL; //temp
double *d_vk = NULL; //temp
double *d_u_wface = NULL;//temp
double *d_v_sface = NULL;//temp

double *d_hdvk2 = NULL; //tsb
double *d_hduk2 = NULL; //tsb

//readyt
double *d_pp = NULL;
double *d_ppa = NULL;
double *d_ppb = NULL;
double *d_ppc = NULL;
double *d_alpha = NULL;
double *d_beta = NULL;

//barotr2
//operators_div_cu
double *d_div_out = NULL;
//hmix_del2_hdiffu_del2_cu
double *d_hduk = NULL;
double *d_hdvk = NULL;
//operators_grad_cu
double *d_gradx = NULL;
double *d_grady = NULL;
//hmix_del2_hdifft_del2_cu
double *d_hdtk = NULL;

//bclinc
//invtriu_device
double *d_a8 = NULL;
double *d_b8 = NULL;
double *d_c8 = NULL;
double *d_d8 = NULL;
double *d_e8 = NULL;
double *d_f8 = NULL;

double *d_a8_1 = NULL;
double *d_b8_1 = NULL;
double *d_c8_1 = NULL;
double *d_d8_1 = NULL;
double *d_e8_1 = NULL;
double *d_f8_1 = NULL;

//pconst
double *d_ohbu = NULL;
double *d_viv = NULL;
double *d_vit = NULL;
double *d_dzp = NULL;
float *d_f_dzp = NULL;
double *d_dzph = NULL;

double *d_to = NULL;
double *d_so = NULL;
double *d_pconst_dts = NULL;
double *d_c = NULL;

char *d_adv_tracer = NULL;

int *d_ist = NULL;
//int *d_isb = NULL;
int *d_kvt = NULL;

double *d_ebea = NULL;
double *d_ebeb = NULL;

double *d_snlat = NULL;

double *d_ohbt = NULL;
double *d_zkt = NULL;
double *d_epea = NULL;
double *d_epeb = NULL;
double *d_epla = NULL;
double *d_eplb = NULL;
double *d_akmu = NULL;

int *d_ncc = NULL;
double *d_ahv = NULL;
double *d_akt = NULL;
double *d_onbc = NULL;
double *d_oncc = NULL;
double *d_odzp = NULL;
double *d_odzt = NULL;
double *d_gamma = NULL;
double *d_dwndmix = NULL;
int *d_boundary_restore = NULL;

char *d_adv_momentum = NULL;

int *d_nss = NULL;

double *d_po = NULL;
double *d_tbice = NULL;
double *d_cp = NULL;

double *d_hbx = NULL;
double *d_hby = NULL;

double *d_zkp = NULL;
double *d_ak_tide = NULL;
double *d_fz_tide = NULL;
double *d_fztidal = NULL;
double *d_wp3_tidal = NULL;
double *d_akmt = NULL;
double *d_richardson = NULL;

//grid
double *d_au0 = NULL;
double *d_aus = NULL;
double *d_auw = NULL;
double *d_ausw = NULL;
double *d_at0 = NULL;
double *d_atn = NULL;
double *d_ate = NULL;
double *d_atne = NULL;
int *d_kmtn = NULL;
int *d_kmts = NULL;
int *d_kmte = NULL;
int *d_kmtw = NULL;
double *d_ulat = NULL;
double *d_uarea = NULL;
double *d_tarea = NULL;
double *d_htw = NULL;
double *d_hts = NULL;
double *d_tarea_r = NULL;
double *d_fcor = NULL;
double *d_dxur = NULL;
double *d_dyur = NULL;

int *d_kmu = NULL;
int *d_kmt = NULL;

double *d_area_t = NULL;

double *d_dxu = NULL;
double *d_dyu = NULL;
double *d_hun = NULL;
double *d_hue = NULL;
double *d_uarea_r = NULL;
double *d_tlat = NULL;
//hmix_del2
double *d_duc = NULL;
double *d_dum = NULL;
double *d_dun = NULL;
double *d_dus = NULL;
double *d_due = NULL;
double *d_duw = NULL;
double *d_dmc = NULL;
double *d_dmn = NULL;
double *d_dms = NULL;
double *d_dme = NULL;
double *d_dmw = NULL;
double *d_dtn = NULL;
double *d_dts = NULL;
double *d_dte = NULL;
double *d_dtw = NULL;
double *d_ahf = NULL;
double *d_amf = NULL;
double *d_am_factor = NULL;
double *d_d2tk = NULL;
double *d_d2uk = NULL;
double *d_d2vk = NULL;

//dyn
double *d_u = NULL;
double *d_v = NULL;
double *d_h0 = NULL;
double *d_h0p = NULL;
double *d_h0f = NULL;
double *d_h0bf = NULL;
double *d_ub = NULL;
double *d_dlub = NULL;
double *d_ubp = NULL;
double *d_vb = NULL;
double *d_dlvb = NULL;
double *d_vbp = NULL;

double *d_up = NULL;
double *d_vp = NULL;
double *d_bbcy = NULL;
double *d_dlu = NULL;
double *d_dlv = NULL;
double *d_h0bl = NULL;
double *d_gg = NULL;
double *d_sbcy = NULL;
double *d_utf = NULL;
double *d_vtf = NULL;
double *d_sbcx = NULL;
double *d_bbcx = NULL;

double *d_h0l = NULL;
double *d_utl = NULL;
double *d_vtl = NULL;

double *d_ws = NULL;


//tracer
double *d_at = NULL;
double *d_atb = NULL;
double *d_tend = NULL;
double *d_dt_conv = NULL;

double *d_dz = NULL;
double *d_net = NULL;
double *d_dt_diff = NULL;
//double *d_fw_norm2 = NULL;
double *d_penetrate = NULL;
double *d_restore_at = NULL;
double *d_dt_restore = NULL;
double *d_amld = NULL;

double *d_az = NULL;
double *d_ay = NULL;
double *d_ax = NULL;

double *d_pdensity = NULL;

double *d_licomqice = NULL;


//output
double *d_icmon = NULL;

//work
double *d_wka = NULL;
double *h_wka=NULL;
double *d_wka1 = NULL;
double *d_wka2 = NULL;
double *d_work = NULL;
double *d_work1 = NULL;
double *d_work2 = NULL;
double *h_work=NULL;
double *d_wgp = NULL;
double *d_pax = NULL;
double *d_pay = NULL;
double *d_pxb = NULL;
double *d_pyb = NULL;
double *d_whx = NULL;
double *d_why = NULL;
//double *d_wkk = NULL;

double *d_stf = NULL;
double *d_wkd = NULL;
double *d_wkb = NULL;
double *d_wkc = NULL;
double *d_tf = NULL;

//forc
double *d_su = NULL;
double *d_sv = NULL;
double *d_psa = NULL;

double *d_swv = NULL;
double *d_tsf = NULL;
double *d_ssf = NULL;
double *d_sss = NULL;
double *d_sst = NULL;
double *d_dqdt = NULL;
double *d_seaice = NULL;
double *d_fresh = NULL;
double *d_restore = NULL;

double *d_buoysol = NULL;
double *d_buoytur = NULL;
double *d_nswv = NULL;

double *d_wave_dis = NULL;
double *d_ustar = NULL;


//isopyc
double *d_ahisop = NULL;
double *d_k3 = NULL;

//pmix
double *d_pen = NULL;
double *d_rit = NULL; //MR.WANG TAO CANCEL THIS CODE 20190724 10.15 A.M.
double *d_ric = NULL;
double *d_rict = NULL;
double *d_ricdttms = NULL;
double *d_ricdt = NULL;
double *d_rict_replace = NULL;

double *d_s2u = NULL;
double *d_s2t = NULL;
//double *d_rict_ref = NULL;
double *d_riu = NULL;
//buf
double *d_ifrac = NULL;

//param
int *d_mytid = NULL;

extern "C" void cudainit() {
    //advection_tracer
    CHECK(hipMalloc((void **) &d_at00, dataSize3d));
    CHECK(hipMalloc((void **) &d_atmax, dataSize3d));
    CHECK(hipMalloc((void **) &d_atmin, dataSize3d));

    //tracer
    CHECK(hipMalloc((void **) &d_adv_tt, dataSize3d));
    CHECK(hipMalloc((void **) &d_ori, dataSize3d));
    CHECK(hipMalloc((void **) &d_temp11, dataSize2d));

    //readyc
    CHECK(hipMalloc((void **) &d_hdvk2, dataSize3d));
    CHECK(hipMalloc((void **) &d_hduk2, dataSize3d));

    CHECK(hipMalloc((void **) &d_wp12, dataSize3d));
    CHECK(hipMalloc((void **) &d_wp13, dataSize3d));
    CHECK(hipMalloc((void **) &d_riv1, dataSize3d));
    CHECK(hipMalloc((void **) &d_riv2, dataSize3d));
    CHECK(hipMalloc((void **) &d_diff_back, dataSize2d));
    CHECK(hipMalloc((void **) &d_diff_back_sh, dataSize2d));
    CHECK(hipMalloc((void **) &d_diff_back_nh, dataSize2d));

    CHECK(hipMalloc((void **) &d_uk, dataSize3d));
    CHECK(hipMalloc((void **) &d_vk, dataSize3d));
    CHECK(hipMalloc((void **) &d_u_wface, dataSize3d));
    CHECK(hipMalloc((void **) &d_v_sface, dataSize3d));

    //readyt
    CHECK(hipMalloc((void **) &d_pp, dataSize3d));
    CHECK(hipMalloc((void **) &d_ppa, dataSize3d));
    CHECK(hipMalloc((void **) &d_ppb, dataSize3d));
    CHECK(hipMalloc((void **) &d_ppc, dataSize3d));
    CHECK(hipMalloc((void **) &d_alpha, dataSize3d));
    CHECK(hipMalloc((void **) &d_beta, dataSize3d));

    //bartor
    //operators_div_cu
    CHECK(hipMalloc((void **) &d_div_out, dataSize22d));
    //hmix_del2_hdiffu_del2_cu
    CHECK(hipMalloc((void **) &d_hduk, dataSize22d));
    CHECK(hipMalloc((void **) &d_hdvk, dataSize22d));
    //operators_grad_cu
    CHECK(hipMalloc((void **) &d_gradx, dataSize22d));
    CHECK(hipMalloc((void **) &d_grady, dataSize22d));
    //hmix_del2_hdifft_del2_cu
    CHECK(hipMalloc((void **) &d_hdtk, dataSize22d));

    //bclinc
    //operators_grad_cu
    CHECK(hipMalloc((void **) &d_gradx, dataSize22d));
    CHECK(hipMalloc((void **) &d_grady, dataSize22d));
    //invtriu_device
    CHECK(hipMalloc((void **) &d_a8, dataSize3d));
    CHECK(hipMalloc((void **) &d_b8, dataSize3d));
    CHECK(hipMalloc((void **) &d_c8, dataSize3d));
    CHECK(hipMalloc((void **) &d_d8, dataSize3d));
    CHECK(hipMalloc((void **) &d_e8, dataSize332d));
    CHECK(hipMalloc((void **) &d_f8, dataSize332d));

    CHECK(hipMalloc((void **) &d_a8_1, dataSize3d));
    CHECK(hipMalloc((void **) &d_b8_1, dataSize3d));
    CHECK(hipMalloc((void **) &d_c8_1, dataSize3d));
    CHECK(hipMalloc((void **) &d_d8_1, dataSize3d));
    CHECK(hipMalloc((void **) &d_e8_1, dataSize332d));
    CHECK(hipMalloc((void **) &d_f8_1, dataSize332d));

//	hipMalloc((void **)&d_thisblock, dataSize1block);

    //pconst
//  extern double pconst_mod_mp_dzp_[km];
//  extern double pconst_mod_mp_viv_[max_blocks_clinic][km][jmt][imt];
//	extern double pconst_mod_mp_vit_[max_blocks_clinic][km][jmt][imt];

    CHECK(hipMalloc((void **) &d_ohbu, dataSize2d));
    CHECK(hipMalloc((void **) &d_ebeb, dataSize2d));
    CHECK(hipMalloc((void **) &d_ebea, dataSize2d));
    CHECK(hipMalloc((void **) &d_viv, dataSize3d));
    CHECK(hipMalloc((void **) &d_vit, dataSize3d));
    CHECK(hipMalloc((void **) &d_dzp, dataSize1d));
    CHECK(hipMalloc((void **) &d_dzph, dataSize2d));

    CHECK(hipMalloc((void **) &d_to, dataSize1d));
    CHECK(hipMalloc((void **) &d_so, dataSize1d));
    CHECK(hipMalloc((void **) &d_pconst_dts, dataSizeD));
    CHECK(hipMalloc((void **) &d_c, 9 * dataSize1d));

    CHECK(hipMalloc((void **) &d_adv_tracer, 80 * dataSizeC));
    CHECK(hipMalloc((void **) &d_ist, dataSizeI));
//	hipMalloc((void **)&d_isb, dataSizeI);
    CHECK(hipMalloc((void **) &d_kvt, dataSizeI));

    CHECK(hipMalloc((void **) &d_snlat, dataSize2d));

    CHECK(hipMalloc((void **) &d_ohbt, dataSize2d));
    CHECK(hipMalloc((void **) &d_zkt, dataSize1d));
    CHECK(hipMalloc((void **) &d_epea, dataSize2d));
    CHECK(hipMalloc((void **) &d_epeb, dataSize2d));
    CHECK(hipMalloc((void **) &d_epla, dataSize2d));
    CHECK(hipMalloc((void **) &d_eplb, dataSize2d));
    CHECK(hipMalloc((void **) &d_akmu, dataSize3d));

    CHECK(hipMalloc((void **) &d_ncc, dataSizeI));
    CHECK(hipMalloc((void **) &d_ahv, dataSizeD));
    CHECK(hipMalloc((void **) &d_akt, dataSize4d));
    CHECK(hipMalloc((void **) &d_onbc, dataSizeD));
    CHECK(hipMalloc((void **) &d_oncc, dataSizeD));
    CHECK(hipMalloc((void **) &d_odzp, dataSize1d));
    CHECK(hipMalloc((void **) &d_odzt, dataSize1d));
    CHECK(hipMalloc((void **) &d_gamma, dataSizeD));
    CHECK(hipMalloc((void **) &d_dwndmix, dataSizeD));
    CHECK(hipMalloc((void **) &d_boundary_restore, dataSizeI));

    CHECK(hipMalloc((void **) &d_adv_momentum, 80 * sizeof(char)));

    CHECK(hipMalloc((void **) &d_nss, dataSizeI));

    CHECK(hipMalloc((void **) &d_po, dataSize1d));
    CHECK(hipMalloc((void **) &d_tbice, dataSizeD));
    CHECK(hipMalloc((void **) &d_cp, dataSizeD));

    CHECK(hipMalloc((void **) &d_hbx, dataSize2d));
    CHECK(hipMalloc((void **) &d_hby, dataSize2d));

    CHECK(hipMalloc((void **) &d_zkp, dataSize1112d));
    CHECK(hipMalloc((void **) &d_ak_tide, dataSize3d));
    CHECK(hipMalloc((void **) &d_fz_tide, dataSize3d));
    CHECK(hipMalloc((void **) &d_fztidal, dataSize3d));
    CHECK(hipMalloc((void **) &d_wp3_tidal, dataSize3d));
    CHECK(hipMalloc((void **) &d_akmt, dataSize3d));
    CHECK(hipMalloc((void **) &d_richardson, dataSize3d));

    //grid
    //extern double grid_mp_uarea_[max_blocks_clinic][ny_block][nx_block];
    //extern double grid_mp_tarea_[max_blocks_clinic][ny_block][nx_block];
    CHECK(hipMalloc((void **) &d_au0, dataSize22d));
    CHECK(hipMalloc((void **) &d_aus, dataSize22d));
    CHECK(hipMalloc((void **) &d_auw, dataSize22d));
    CHECK(hipMalloc((void **) &d_ausw, dataSize22d));
    CHECK(hipMalloc((void **) &d_at0, dataSize22d));
    CHECK(hipMalloc((void **) &d_atn, dataSize22d));
    CHECK(hipMalloc((void **) &d_ate, dataSize22d));
    CHECK(hipMalloc((void **) &d_atne, dataSize22d));
    CHECK(hipMalloc((void **) &d_kmtn, dataSize2i));
    CHECK(hipMalloc((void **) &d_kmts, dataSize2i));
    CHECK(hipMalloc((void **) &d_kmte, dataSize2i));
    CHECK(hipMalloc((void **) &d_kmtw, dataSize2i));
    CHECK(hipMalloc((void **) &d_ulat, dataSize22d));
    CHECK(hipMalloc((void **) &d_uarea, dataSize22d));
    CHECK(hipMalloc((void **) &d_tarea, dataSize22d));
    CHECK(hipMalloc((void **) &d_htw, dataSize22d));
    CHECK(hipMalloc((void **) &d_hts, dataSize22d));
    CHECK(hipMalloc((void **) &d_tarea_r, dataSize22d));
    CHECK(hipMalloc((void **) &d_fcor, dataSize22d));
    CHECK(hipMalloc((void **) &d_dxur, dataSize22d));
    CHECK(hipMalloc((void **) &d_dyur, dataSize22d));

    CHECK(hipMalloc((void **) &d_kmu, dataSize2i));
    CHECK(hipMalloc((void **) &d_kmt, dataSize2i));

    CHECK(hipMalloc((void **) &d_area_t, dataSizeD));
    CHECK(hipMalloc((void **) &d_tlat, dataSize22d));
    //hmix_del2
    CHECK(hipMalloc((void **) &d_duc, dataSize22d));
    CHECK(hipMalloc((void **) &d_dum, dataSize22d));
    CHECK(hipMalloc((void **) &d_dun, dataSize22d));
    CHECK(hipMalloc((void **) &d_dus, dataSize22d));
    CHECK(hipMalloc((void **) &d_due, dataSize22d));
    CHECK(hipMalloc((void **) &d_duw, dataSize22d));
    CHECK(hipMalloc((void **) &d_dmc, dataSize22d));
    CHECK(hipMalloc((void **) &d_dmn, dataSize22d));
    CHECK(hipMalloc((void **) &d_dms, dataSize22d));
    CHECK(hipMalloc((void **) &d_dme, dataSize22d));
    CHECK(hipMalloc((void **) &d_dmw, dataSize22d));
    CHECK(hipMalloc((void **) &d_dtn, dataSize22d));
    CHECK(hipMalloc((void **) &d_dts, dataSize22d));
    CHECK(hipMalloc((void **) &d_dte, dataSize22d));
    CHECK(hipMalloc((void **) &d_dtw, dataSize22d));
    CHECK(hipMalloc((void **) &d_ahf, dataSize22d));
    CHECK(hipMalloc((void **) &d_amf, dataSize22d));
    CHECK(hipMalloc((void **) &d_am_factor, dataSize3d)); //jjr
    CHECK(hipMalloc((void **) &d_d2tk, dataSize3d));
    CHECK(hipMalloc((void **) &d_d2uk, dataSize3d));
    CHECK(hipMalloc((void **) &d_d2vk, dataSize3d));

    //dyn
    //extern double(*dyn_mod_mp_u_)[km][jmt][imt];
    //extern double(*dyn_mod_mp_v_)[km][jmt][imt];
    //extern double(*dyn_mod_mp_h0_)[jmt][imt];
    CHECK(hipMalloc((void **) &d_u, dataSize3d));
    CHECK(hipMalloc((void **) &d_v, dataSize3d));
    CHECK(hipMalloc((void **) &d_h0, dataSize2d));
    CHECK(hipMalloc((void **) &d_ub, dataSize2d));
    CHECK(hipMalloc((void **) &d_vb, dataSize2d));
    CHECK(hipMalloc((void **) &d_h0p, dataSize2d));
    CHECK(hipMalloc((void **) &d_h0f, dataSize2d));
    CHECK(hipMalloc((void **) &d_h0bf, dataSize2d));
    CHECK(hipMalloc((void **) &d_ubp, dataSize2d));
    CHECK(hipMalloc((void **) &d_vbp, dataSize2d));
    CHECK(hipMalloc((void **) &d_dlub, dataSize2d));
    CHECK(hipMalloc((void **) &d_dlvb, dataSize2d));

    CHECK(hipMalloc((void **) &d_up, dataSize3d));
    CHECK(hipMalloc((void **) &d_vp, dataSize3d));
    CHECK(hipMalloc((void **) &d_bbcy, dataSize2d));
    CHECK(hipMalloc((void **) &d_dlu, dataSize3d));
    CHECK(hipMalloc((void **) &d_dlv, dataSize3d));
    CHECK(hipMalloc((void **) &d_h0bl, dataSize2d));
    CHECK(hipMalloc((void **) &d_gg, dataSize3d));
    CHECK(hipMalloc((void **) &d_sbcy, dataSize2d));
    CHECK(hipMalloc((void **) &d_utf, dataSize3d));
    CHECK(hipMalloc((void **) &d_vtf, dataSize3d));
    CHECK(hipMalloc((void **) &d_sbcx, dataSize2d));
    CHECK(hipMalloc((void **) &d_bbcx, dataSize2d));

    CHECK(hipMalloc((void **) &d_h0l, dataSize2d));
    CHECK(hipMalloc((void **) &d_utl, dataSize3d));
    CHECK(hipMalloc((void **) &d_vtl, dataSize3d));

    CHECK(hipMalloc((void **) &d_dxu, dataSize22d));
    CHECK(hipMalloc((void **) &d_dyu, dataSize22d));
    CHECK(hipMalloc((void **) &d_hun, dataSize22d));
    CHECK(hipMalloc((void **) &d_hue, dataSize22d));
    CHECK(hipMalloc((void **) &d_uarea_r, dataSize22d));

    CHECK(hipMalloc((void **) &d_ws, dataSize333d));

    //tracer
    //extern double(*tracer_mod_mp_at_)[ntra][km][jmt][imt];
    CHECK(hipMalloc((void **) &d_at, dataSize4d));
    CHECK(hipMalloc((void **) &d_tend, dataSize4d));
    CHECK(hipMalloc((void **) &d_dt_conv, dataSize4d));
    CHECK(hipMalloc((void **) &d_atb, dataSize44d));

    CHECK(hipMalloc((void **) &d_dz, dataSize4d));
    CHECK(hipMalloc((void **) &d_net, dataSize33d));
    CHECK(hipMalloc((void **) &d_dt_diff, dataSize4d));
//    hipMalloc((void **) &d_fw_norm2, dataSizeD);
    CHECK(hipMalloc((void **) &d_penetrate, dataSize3d));
    CHECK(hipMalloc((void **) &d_restore_at, dataSize4d));
    CHECK(hipMalloc((void **) &d_dt_restore, dataSize4d));

    CHECK(hipMalloc((void **) &d_az, dataSize4d));
    CHECK(hipMalloc((void **) &d_ay, dataSize4d));
    CHECK(hipMalloc((void **) &d_ax, dataSize4d));

    CHECK(hipMalloc((void **) &d_pdensity, dataSize332d));

    CHECK(hipMalloc((void **) &d_licomqice, dataSize2d));
    CHECK(hipMalloc((void **) &d_amld, dataSize2d));

    //output
    CHECK(hipMalloc((void **) &d_icmon, 2 * dataSize2d));

    //work
    CHECK(hipMalloc((void **) &d_wka, dataSize3d));
//    CHECK(hipHostMalloc((void **) &h_wka, dataSize3d));
    CHECK(hipMalloc((void **) &d_wka1, dataSize3d));
    CHECK(hipMalloc((void **) &d_wka2, dataSize3d));
    CHECK(hipMalloc((void **) &d_work, dataSize2d));
    CHECK(hipMalloc((void **) &d_work1, dataSize2d));
    CHECK(hipMalloc((void **) &d_work2, dataSize2d));
    CHECK(hipHostMalloc((void **) &h_work, dataSize2d)); 
    CHECK(hipMalloc((void **) &d_wgp, dataSize2d));
    CHECK(hipMalloc((void **) &d_pax, dataSize2d));
    CHECK(hipMalloc((void **) &d_pay, dataSize2d));
    CHECK(hipMalloc((void **) &d_pxb, dataSize2d));
    CHECK(hipMalloc((void **) &d_pyb, dataSize2d));
    CHECK(hipMalloc((void **) &d_whx, dataSize2d));
    CHECK(hipMalloc((void **) &d_why, dataSize2d));
//    hipMalloc((void **) &d_wkk, dataSize11d);

    CHECK(hipMalloc((void **) &d_stf, dataSize2d));
    CHECK(hipMalloc((void **) &d_wkd, dataSize3d));
    CHECK(hipMalloc((void **) &d_wkb, dataSize3d));
    CHECK(hipMalloc((void **) &d_wkc, dataSize3d));
    CHECK(hipMalloc((void **) &d_tf, dataSize3d));

    //forc
    CHECK(hipMalloc((void **) &d_su, dataSize2d));
    CHECK(hipMalloc((void **) &d_sv, dataSize2d));
    CHECK(hipMalloc((void **) &d_psa, dataSize2d));

    CHECK(hipMalloc((void **) &d_buoysol, dataSize2d));
    CHECK(hipMalloc((void **) &d_buoytur, dataSize2d));
    CHECK(hipMalloc((void **) &d_nswv, dataSize2d));


    CHECK(hipMalloc((void **) &d_swv, dataSize2d));
    CHECK(hipMalloc((void **) &d_tsf, dataSize2d));
    CHECK(hipMalloc((void **) &d_ssf, dataSize2d));
    CHECK(hipMalloc((void **) &d_sss, dataSize2d));
    CHECK(hipMalloc((void **) &d_sst, dataSize2d));
    CHECK(hipMalloc((void **) &d_dqdt, dataSize2d));
    CHECK(hipMalloc((void **) &d_seaice, dataSize2d));
    CHECK(hipMalloc((void **) &d_fresh, dataSize2d));
    CHECK(hipMalloc((void **) &d_restore, dataSize4d));

    CHECK(hipMalloc((void **) &d_wave_dis, dataSize2d));
    CHECK(hipMalloc((void **) &d_ustar, dataSize2d));

    //isopyc
    CHECK(hipMalloc((void **) &d_ahisop, dataSize2d));
    CHECK(hipMalloc((void **) &d_k3, dataSize444d));

    //pmix
    CHECK(hipMalloc((void **) &d_pen, dataSize111d));


    CHECK(hipMalloc((void **) &d_rit, dataSize3333d)); //MR. WANG TAO CANCEL THIS CODE 20190724 10:16 A.M.
    CHECK(hipMalloc((void **) &d_ric, dataSize3333d));
    CHECK(hipMalloc((void **) &d_rict, dataSize3333d));
    CHECK(hipMalloc((void **) &d_ricdttms, dataSize3333d));
    CHECK(hipMalloc((void **) &d_ricdt, dataSize3333d));
    CHECK(hipMalloc((void **) &d_rict_replace, dataSize3333d));

    CHECK(hipMalloc((void **) &d_s2u, dataSize3333d));
    CHECK(hipMalloc((void **) &d_s2t, dataSize3333d));
    //CHECK(hipMalloc((void **) &d_rict_ref, dataSize2d));
    CHECK(hipMalloc((void **) &d_riu, dataSize332d));

    //buf
    CHECK(hipMalloc((void **) &d_ifrac, dataSize2d));

    //param
    CHECK(hipMalloc((void **) &d_mytid, dataSizeI));


//
     CHECK(hipMalloc((void **) &d_wk1, dataSize3d));
     CHECK(hipMalloc((void **) &d_wk2, dataSize3d));
     CHECK(hipMalloc((void **) &d_wk3, dataSize3d));
     CHECK(hipMalloc((void **) &d_wp1, dataSize3d));
     CHECK(hipMalloc((void **) &d_wp3, dataSize3d));
     CHECK(hipMalloc((void **) &d_wp7, dataSize3d));
     CHECK(hipMalloc((void **) &d_wp8, dataSize3d));


//readyc.cpp
    CHECK(hipMemcpy(d_utf, dyn_mod_mp_utf_, dataSize3d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_utl, dyn_mod_mp_utl_, dataSize3d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_vtf, dyn_mod_mp_vtf_, dataSize3d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_vtl, dyn_mod_mp_vtl_, dataSize3d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_h0l, dyn_mod_mp_h0l_, dataSize2d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_h0f, dyn_mod_mp_h0f_, dataSize2d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_h0p, dyn_mod_mp_h0p_, dataSize2d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_h0, dyn_mod_mp_h0_, dataSize2d, hipMemcpyHostToDevice));
    //atb copy at steponinit, dlu,dlv is from atb
    //hipMemcpy(d_dlu, dyn_mod_mp_dlu_, dataSize3d, hipMemcpyHostToDevice);
    //hipMemcpy(d_dlv, dyn_mod_mp_dlv_, dataSize3d, hipMemcpyHostToDevice);
    CHECK(hipMemcpy(d_viv, pconst_mod_mp_viv_, dataSize3d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_vit, pconst_mod_mp_vit_, dataSize3d, hipMemcpyHostToDevice));
    //hipMemcpy(d_gg, dyn_mod_mp_gg_, dataSize3d, hipMemcpyHostToDevice);
    CHECK(hipMemcpy(d_psa, forc_mod_mp_psa_, dataSize2d, hipMemcpyHostToDevice));

//readyc.cpp
    CHECK(hipMemcpy(d_up, dyn_mod_mp_up_, dataSize3d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_vp, dyn_mod_mp_vp_, dataSize3d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_u, dyn_mod_mp_u_, dataSize3d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_v, dyn_mod_mp_v_, dataSize3d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_ws, dyn_mod_mp_ws_, dataSize333d, hipMemcpyHostToDevice));
//    hipMemcpy(d_h0bl, dyn_mod_mp_h0bl_, dataSize2d, hipMemcpyHostToDevice);
    CHECK(hipMemcpy(d_h0bf, dyn_mod_mp_h0bf_, dataSize2d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_sbcx, dyn_mod_mp_sbcx_, dataSize2d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_bbcx, dyn_mod_mp_bbcx_, dataSize2d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_sbcy, dyn_mod_mp_sbcy_, dataSize2d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_bbcy, dyn_mod_mp_bbcy_, dataSize2d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_ub, dyn_mod_mp_ub_, dataSize2d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_vb, dyn_mod_mp_vb_, dataSize2d, hipMemcpyHostToDevice));
    hipMemcpy(d_ubp, dyn_mod_mp_ubp_, dataSize2d, hipMemcpyHostToDevice);
    hipMemcpy(d_vbp, dyn_mod_mp_vbp_, dataSize2d, hipMemcpyHostToDevice);

#ifdef  BIHAR
    CHECK(hipMemcpy(d_duc, hmix_del4_mp_duc_, dataSize22d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_dun, hmix_del4_mp_dun_, dataSize22d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_dus, hmix_del4_mp_dus_, dataSize22d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_due, hmix_del4_mp_due_, dataSize22d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_duw, hmix_del4_mp_duw_, dataSize22d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_dum, hmix_del4_mp_dum_, dataSize22d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_dmc, hmix_del4_mp_dmc_, dataSize22d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_dmn, hmix_del4_mp_dmn_, dataSize22d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_dms, hmix_del4_mp_dms_, dataSize22d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_dme, hmix_del4_mp_dme_, dataSize22d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_dmw, hmix_del4_mp_dmw_, dataSize22d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_dte, hmix_del4_mp_dte_, dataSize22d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_dtn, hmix_del4_mp_dtn_, dataSize22d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_dts, hmix_del4_mp_dts_, dataSize22d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_dtw, hmix_del4_mp_dtw_, dataSize22d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_ahf, hmix_del4_mp_ahf_, dataSize22d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_amf, hmix_del4_mp_amf_, dataSize22d, hipMemcpyHostToDevice));
#else
    CHECK(hipMemcpy(d_duc, hmix_del2_mp_duc_, dataSize22d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_dun, hmix_del2_mp_dun_, dataSize22d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_dus, hmix_del2_mp_dus_, dataSize22d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_due, hmix_del2_mp_due_, dataSize22d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_duw, hmix_del2_mp_duw_, dataSize22d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_dum, hmix_del2_mp_dum_, dataSize22d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_dmc, hmix_del2_mp_dmc_, dataSize22d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_dmn, hmix_del2_mp_dmn_, dataSize22d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_dms, hmix_del2_mp_dms_, dataSize22d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_dme, hmix_del2_mp_dme_, dataSize22d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_dmw, hmix_del2_mp_dmw_, dataSize22d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_dte, hmix_del2_mp_dte_, dataSize22d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_dtn, hmix_del2_mp_dtn_, dataSize22d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_dts, hmix_del2_mp_dts_, dataSize22d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_dtw, hmix_del2_mp_dtw_, dataSize22d, hipMemcpyHostToDevice));
#endif

    CHECK(hipMemcpy(d_tarea_r, grid_mp_tarea_r_, dataSize22d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_tarea, grid_mp_tarea_, dataSize22d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_uarea_r, grid_mp_uarea_r_, dataSize22d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_hts, grid_mp_hts_, dataSize22d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_htw, grid_mp_htw_, dataSize22d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_hun, grid_mp_hun_, dataSize22d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_hue, grid_mp_hue_, dataSize22d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_kmu, grid_mp_kmu_, dataSize2i, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_kmt, grid_mp_kmt_, dataSize2i, hipMemcpyHostToDevice));//turb_2 need grid_mp_kmt_
    CHECK(hipMemcpy(d_kmte, grid_mp_kmte_, dataSize2i, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_kmtn, grid_mp_kmtn_, dataSize2i, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_kmts, grid_mp_kmts_, dataSize2i, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_kmtw, grid_mp_kmtw_, dataSize2i, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_dxur, grid_mp_dxur_, dataSize22d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_dyur, grid_mp_dyur_, dataSize22d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_au0, grid_mp_au0_, dataSize22d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_aus, grid_mp_aus_, dataSize22d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_ausw, grid_mp_ausw_, dataSize22d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_auw, grid_mp_auw_, dataSize22d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_at0, grid_mp_at0_, dataSize22d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_atn, grid_mp_atn_, dataSize22d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_ate, grid_mp_ate_, dataSize22d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_atne, grid_mp_atne_, dataSize22d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_tlat, grid_mp_tlat_, dataSize22d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_dxu, grid_mp_dxu_, dataSize22d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_dyu, grid_mp_dyu_, dataSize22d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_fcor, grid_mp_fcor_, dataSize22d, hipMemcpyHostToDevice));

//    hipMemcpy(d_rit, pmix_mod_mp_rit_, dataSize3333d, hipMemcpyHostToDevice);
//    hipMemcpy(d_ric, pmix_mod_mp_ric_, dataSize3333d, hipMemcpyHostToDevice);
//    hipMemcpy(d_rict, pmix_mod_mp_rict_, dataSize3333d, hipMemcpyHostToDevice);
    CHECK(hipMemcpy(d_pen, pmix_mod_mp_pen_, dataSize111d, hipMemcpyHostToDevice));

    CHECK(hipMemcpy(d_ebea, pconst_mod_mp_ebea_, dataSize2d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_ebeb, pconst_mod_mp_ebeb_, dataSize2d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_epea, pconst_mod_mp_epea_, dataSize2d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_epeb, pconst_mod_mp_epeb_, dataSize2d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_epla, pconst_mod_mp_epla_, dataSize2d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_eplb, pconst_mod_mp_eplb_, dataSize2d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_zkt, pconst_mod_mp_zkt_, dataSize1d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_ohbt, pconst_mod_mp_ohbt_, dataSize2d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_snlat, pconst_mod_mp_snlat_, dataSize2d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_odzt, pconst_mod_mp_odzt_, dataSize1d, hipMemcpyHostToDevice));
//    hipMemcpy(d_akmu, pconst_mod_mp_akmu_, dataSize3d, hipMemcpyHostToDevice);
    CHECK(hipMemcpy(d_ohbu, pconst_mod_mp_ohbu_, dataSize2d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_dzp, pconst_mod_mp_dzp_, dataSize1d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_dzph, pconst_mod_mp_dzph_, dataSize2d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_odzp, pconst_mod_mp_odzp_, dataSize1d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_to, pconst_mod_mp_to_, dataSize1d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_so, pconst_mod_mp_so_, dataSize1d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_c, pconst_mod_mp_c_, 9 * dataSize1d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_po, pconst_mod_mp_po_, dataSize1d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_hbx, pconst_mod_mp_hbx_, dataSize2d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_hby, pconst_mod_mp_hby_, dataSize2d, hipMemcpyHostToDevice));
//    hipMemcpy(d_akt, pconst_mod_mp_akt_, dataSize4d, hipMemcpyHostToDevice);
    CHECK(hipMemcpy(d_fz_tide, pconst_mod_mp_fz_tide_, dataSize3d, hipMemcpyHostToDevice));

    CHECK(hipMemcpy(d_restore, forc_mod_mp_restore_, dataSize4d, hipMemcpyHostToDevice));

    CHECK(hipMemcpy(d_pdensity, tracer_mod_mp_pdensity_, dataSize332d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_at, tracer_mod_mp_at_, dataSize4d, hipMemcpyHostToDevice));

//    hipMemcpy(d_ifrac, buf_mod_mp_ifrac_, dataSize2d, hipMemcpyHostToDevice);

//    CHECK(hipMemcpy(d_work, work_mod_mp_work_, dataSize2d, hipMemcpyHostToDevice));
    CHECK(hipMemcpy(d_wka, work_mod_mp_wka_, dataSize3d, hipMemcpyHostToDevice));
    

//        allocate_fortran_();
//move from steponinit
    hipMemcpy(d_ssf, forc_mod_mp_ssf_, dataSize2d, hipMemcpyHostToDevice);
    hipMemcpy(d_sss, forc_mod_mp_sss_, dataSize2d, hipMemcpyHostToDevice);
    hipMemcpy(d_sst, forc_mod_mp_sst_, dataSize2d, hipMemcpyHostToDevice);
    hipMemcpy(d_dqdt, forc_mod_mp_dqdt_, dataSize2d, hipMemcpyHostToDevice);
    hipMemcpy(d_seaice, forc_mod_mp_seaice_, dataSize2d, hipMemcpyHostToDevice);
    hipMemcpy(d_fresh, forc_mod_mp_fresh_, dataSize2d, hipMemcpyHostToDevice);
    hipMemcpy(d_tsf, forc_mod_mp_tsf_, dataSize2d, hipMemcpyHostToDevice);
    hipMemcpy(d_su, forc_mod_mp_su_, dataSize2d, hipMemcpyHostToDevice);
    hipMemcpy(d_sv, forc_mod_mp_sv_, dataSize2d, hipMemcpyHostToDevice);
    hipMemcpy(d_nswv, forc_mod_mp_nswv_, dataSize2d, hipMemcpyHostToDevice);
    hipMemcpy(d_swv, forc_mod_mp_swv_, dataSize2d, hipMemcpyHostToDevice);
    hipMemcpy(d_ifrac, buf_mod_mp_ifrac_, dataSize2d, hipMemcpyHostToDevice);
    hipMemcpy(d_licomqice, tracer_mod_mp_licomqice_, dataSize2d, hipMemcpyHostToDevice);
    return;
}
