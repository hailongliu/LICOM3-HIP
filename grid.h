#ifndef INCLUDE_GRID
#define INCLUDE_GRID

#include "precision_mod.h"
#include "param_mod.h"
#include "blocks.h"
#include "LICOM_Error_mod.h"
#include "pconst_mod.h"

extern double grid_mp_area_t_, grid_mp_volume_t_;
extern double (*grid_mp_uarea_)[ny_block][nx_block];
extern double (*grid_mp_tarea_)[ny_block][nx_block];
extern double (*grid_mp_uarea_r_)[ny_block][nx_block];
extern double (*grid_mp_tarea_r_)[ny_block][nx_block];
extern double (*grid_mp_ulat_)[ny_block][nx_block];
extern double (*grid_mp_ulon_)[ny_block][nx_block];
extern double (*grid_mp_tlat_)[ny_block][nx_block];
extern double (*grid_mp_tlon_)[ny_block][nx_block];
extern double (*grid_mp_fcor_)[ny_block][nx_block];
extern double (*grid_mp_fcort_)[ny_block][nx_block];
extern double (*grid_mp_dxu_)[ny_block][nx_block];
extern double (*grid_mp_dyu_)[ny_block][nx_block];
extern double (*grid_mp_dxur_)[ny_block][nx_block];
extern double (*grid_mp_dyur_)[ny_block][nx_block];
extern double (*grid_mp_dxt_)[ny_block][nx_block];
extern double (*grid_mp_dyt_)[ny_block][nx_block];
extern double (*grid_mp_dxtr_)[ny_block][nx_block];
extern double (*grid_mp_dytr_)[ny_block][nx_block];
extern double (*grid_mp_htw_)[ny_block][nx_block];
extern double (*grid_mp_hts_)[ny_block][nx_block];
extern double (*grid_mp_hue_)[ny_block][nx_block];
extern double (*grid_mp_hun_)[ny_block][nx_block];
extern double (*grid_mp_kmtn_)[ny_block][nx_block];
extern double (*grid_mp_kmts_)[ny_block][nx_block];
extern double (*grid_mp_kmte_)[ny_block][nx_block];
extern double (*grid_mp_kmtw_)[ny_block][nx_block];

extern int (*grid_mp_kmu_)[ny_block][nx_block];
extern int (*grid_mp_kmt_)[ny_block][nx_block];

extern double (*grid_mp_at0_)[ny_block][nx_block];
extern double (*grid_mp_atn_)[ny_block][nx_block];
extern double (*grid_mp_ate_)[ny_block][nx_block];
extern double (*grid_mp_atne_)[ny_block][nx_block];
extern double (*grid_mp_au0_)[ny_block][nx_block];
extern double (*grid_mp_aus_)[ny_block][nx_block];
extern double (*grid_mp_auw_)[ny_block][nx_block];
extern double (*grid_mp_ausw_)[ny_block][nx_block];
//-----------------------------------------------------------------------
//
//  variables which are shared between init_grid1,init_grid2
//
//-----------------------------------------------------------------------
extern char grid_mp_horiz_grid_opt_[char_len];
extern char grid_mp_horiz_grid_opt_[char_len];
extern char grid_mp_vert_grid_opt_[char_len];      
extern char grid_mp_sfc_layer_opt_[char_len];
extern char grid_mp_topography_opt_[char_len];
extern char grid_mp_horiz_grid_file_[char_len];
extern char grid_mp_vert_grid_file_[char_len];
extern char grid_mp_topography_file_[char_len];
extern char grid_mp_region_mask_file_[char_len];
extern char grid_mp_region_info_file_[char_len];
extern char grid_mp_bottom_cell_file_[char_len];
extern char grid_mp_topography_outfile_[char_len];
extern char grid_mp_basin_grid_file_[char_len];

extern void grid_mp_tgrid_to_ugrid_(double[][imt], double[][imt], int *);
extern void grid_mp_ugrid_to_tgrid_(double[][imt], double[][imt], int *, int *);

extern void grid_tgrid_to_ugrid_(double[][imt], double[][imt], int *);
extern void grid_ugrid_to_tgrid_(double[][imt], double[][imt], int *, int *);

#endif // !INCLUDE_GRID

