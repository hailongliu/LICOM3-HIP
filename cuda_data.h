#ifndef INCLUDE_CUDA_DATA
#define INCLUDE_CUDA_DATA

#define CHECK(call)\
{\
        const hipError_t error = call;\
        if(error != hipSuccess){\
                printf("Error:%s:%d,Code:%d,Reason:%s tid=%d\n",__FILE__,__LINE__,error,hipGetErrorString(error),param_mod_mp_mytid_);\
                exit(1);\
        }\
}
extern size_t dataSizeI;
extern size_t dataSizeD;
extern size_t dataSizeC;
extern size_t dataSizeB;
extern size_t dataSize1d;
extern size_t dataSize11d;
extern size_t dataSize12d;
extern size_t dataSize111d;
extern size_t dataSize1112d;
extern size_t dataSize2d;
extern size_t dataSize22d;
extern size_t dataSize222d;
extern size_t dataSize3d;
extern size_t dataSize33d;
extern size_t dataSize332d;
extern size_t dataSize333d;
extern size_t dataSize3333d;
extern size_t dataSize4d;
extern size_t dataSize44d;
extern size_t dataSize444d;

extern size_t dataSize2i;

//advection_tracer
extern double *d_at00;
extern double *d_atmax;
extern double *d_atmin;

//tracer
extern double *d_adv_tt;
extern double *d_ori;
extern double *d_temp11;

//readyc
extern double *d_wp12;
extern double *d_wp13;
extern double *d_riv1;
extern double *d_riv2;
extern double *d_diff_back;
extern double *d_diff_back_sh;
extern double *d_diff_back_nh;
extern double *d_uk;
extern double *d_vk;
extern double *d_u_wface;
extern double *d_v_sface;
extern double *d_hdvk2;
extern double *d_hduk2;

//readyt
extern double *d_pp;
extern double *d_ppa;
extern double *d_ppb;
extern double *d_ppc;
extern double *d_alpha;
extern double *d_beta;

extern double *d_div_out;
extern double *d_hduk;
extern double *d_hdvk;
extern double *d_gradx;
extern double *d_grady;
extern double *d_hdtk;
//bclinc
extern double *d_a8;
extern double *d_b8;
extern double *d_c8;
extern double *d_d8;
extern double *d_e8;
extern double *d_f8;

extern double *d_a8_1;
extern double *d_b8_1;
extern double *d_c8_1;
extern double *d_d8_1;
extern double *d_e8_1;
extern double *d_f8_1;

//pconst
extern double *d_ohbu;
extern double *d_viv;
extern double *d_vit;
extern double *d_dzp;
extern float *d_f_dzp;
extern double *d_dzph;
extern double *d_to;
extern double *d_so;
extern double *d_pconst_dts;
extern double *d_c;
extern char *d_adv_tracer;

extern int *d_ist;
//extern int *d_isb;
extern int *d_kvt;

extern double *d_ebea;
extern double *d_ebeb;

extern double *d_snlat;

extern double *d_ohbt;
extern double *d_zkt;
extern double *d_epea;
extern double *d_epeb;
extern double *d_epla;
extern double *d_eplb;
extern double *d_akmu;

extern double *d_wk1;
extern double *d_wk2;
extern double *d_wk3;
extern double *d_wp1;
extern double *d_wp3;
extern double *d_wp7;
extern double *d_wp8;

extern double *d_ahf;
extern double *d_amf;
extern double *d_am_factor;
extern double *d_d2tk;
extern double *d_d2uk;
extern double *d_d2vk;


extern int *d_ncc;
extern double *d_ahv;
extern double *d_akt;
extern double *d_onbc;
extern double *d_oncc;
extern double *d_odzp;
extern double *d_odzt;
extern double *d_gamma;
extern double *d_dwndmix;
extern int *d_boundary_restore;

extern char *d_adv_momentum;

extern int *d_nss;

extern double *d_po;
extern double *d_tbice;
extern double *d_cp;

extern double *d_hbx;
extern double *d_hby;

extern double *d_zkp; 
extern double *d_ak_tide; 
extern double *d_fz_tide;
extern double *d_fztidal;
extern double *d_wp3_tidal;
extern double *d_akmt;
extern double *d_richardson;
//grid
extern double *d_au0;
extern double *d_aus;
extern double *d_auw;
extern double *d_ausw;
extern double *d_at0;
extern double *d_atn;
extern double *d_ate;
extern double *d_atne;
extern int *d_kmtn;
extern int *d_kmts;
extern int *d_kmte;
extern int *d_kmtw;
extern double *d_ulat;
extern double *d_uarea;
extern double *d_tarea;
extern double *d_fcor;
extern double *d_dxur;
extern double *d_dyur;
extern double *d_htw;
extern double *d_hts;
extern double *d_tarea_r;
extern double *d_area_t;

extern int *d_kmu;
extern int *d_kmt;

extern double *d_dxu; 
extern double *d_dyu;
extern double *d_hun;
extern double *d_hue;
extern double *d_uarea_r;  
extern double *d_tlat;

//hmix_del2
extern double *d_duc;
extern double *d_dum;
extern double *d_dun;
extern double *d_dus;
extern double *d_due;
extern double *d_duw;
extern double *d_dmc;
extern double *d_dmn;
extern double *d_dms;
extern double *d_dme;
extern double *d_dmw;
extern double *d_dtn;
extern double *d_dts;
extern double *d_dte;
extern double *d_dtw;


//dyn
extern double *d_u;
extern double *d_v;
extern double *d_ub;
extern double *d_vb;
extern double *d_h0;
extern double *d_h0p;
extern double *d_h0f;
extern double *d_h0bf;
extern double *d_ubp;
extern double *d_vbp;
extern double *d_dlub;
extern double *d_dlvb;

extern double *d_up;

//extern double *d_utest; //jjr for test
                             
extern double *d_vp;
extern double *d_bbcy;
extern double *d_dlu;
extern double *d_dlv;
extern double *d_h0bl;
extern double *d_gg;
extern double *d_sbcy;
extern double *d_utf;
extern double *d_vtf;
extern double *d_sbcx;
extern double *d_bbcx;

extern double *d_h0l;
extern double *d_utl;
extern double *d_vtl;

extern double *d_ws;

//tracer
extern double *d_at;
extern double *d_atb;
extern double *d_tend;
extern double *d_dt_conv;

extern double *d_dz;
extern double *d_net;
extern double *d_dt_diff;
//extern double *d_fw_norm2;
extern double *d_penetrate;
extern double *d_restore_at;
extern double *d_dt_restore;

extern double *d_az;
extern double *d_ay; 
extern double *d_ax;  
extern double *d_pdensity;

extern double *d_licomqice;
extern double *d_amld;

//output
extern double *d_icmon;

//work
extern double *d_wka;
extern double *d_wka1;
extern double *d_wka2;
extern double *d_work;
extern double *d_work1;
extern double *d_work2;
extern double *h_wka;
extern double *h_work;
extern double *d_wgp;
extern double *d_pax;
extern double *d_pay;
extern double *d_pxb;
extern double *d_pyb;
extern double *d_whx;
extern double *d_why;
//extern double *d_wkk;

extern double *d_stf;
extern double *d_wkd;
extern double *d_wkb;
extern double *d_wkc;
extern double *d_tf;

//forc
extern double *d_su;
extern double *d_sv;
extern double *d_psa;
extern double *d_swv;
extern double *d_tsf;
extern double *d_ssf;
extern double *d_sss;
extern double *d_seaice;
extern double *d_fresh;
extern double *d_sst;
extern double *d_dqdt; 
extern double *d_restore;

extern double *d_buoysol;
extern double *d_buoytur;
extern double *d_nswv;

extern double *d_wave_dis;      
extern double *d_ustar; 

//isopyc
extern double *d_ahisop;
extern double *d_k3;

//pmix
extern double *d_pen;

extern double *d_rit; 
extern double *d_ric;
extern double *d_rict;
extern double *d_ricdttms;
extern double *d_ricdt;
extern double *d_rict_replace;

extern double *d_s2u; 
extern double *d_s2t;
//extern double *d_rict_ref;
extern double *d_riu;

//buf
extern double *d_ifrac;

//param
extern int *d_mytid;

#endif // !INCLUDE_CUDA_DATA
