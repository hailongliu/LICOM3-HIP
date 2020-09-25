/*  by wpf for 419 project, 2018.11 */

#ifndef INCLUDED_PARAM
#define INCLUDED_PARAM //to avoid duplicated include

#include <def-undef.h>
#include "precision_mod.h"
//-------------------------------------------------------------------------------
//
// Author: Yongqiang YU  ( 12 Nov, 2002)
//
//-------------------------------------------------------------------------------
//YU  01/02/2013
#define licom_BlockSizeX   BLCKX // size of block in first  horizontal dimension
#define licom_BlockSizeY   BLCKY // size of block in second horizontal dimension
#define max_blocks_clinic  MXBLCKS
#define max_blocks_tropic  MXBLCKS    //   in each distribution
#define nghost  2       // number of ghost cells around each block
#define nx_block (licom_BlockSizeX + 2*nghost)   //  x,y dir including ghost
#define ny_block (licom_BlockSizeY + 2*nghost)   //  cells
//YU
#define jmt_global NJMT  // Number of the End Grid for Tracer in Latitude.
#define jmm_global jmt_global-1
#define imt_global NIMT    // Number of Grid Points in Longitude
#define km NKM     // Number of Grid Points in Vertical Direction
//
#define num_overlap 2 // Number of overlapping grids for subdomain.

//Nummber of grids in the each subdomain
#define jst 1     // Number of the Strating Grid for Tracer in Latitude.
#define jsm (jst+1) // Number of the Strating Grid for Momentum in Latitude.
//ZHW should Add 1
#define jet ny_block
#define jem (jet-1) // Number of the End Grid for Momentum in Latitude.
#define jmt ny_block
#define imt nx_block
extern int j_loop;             // Loop index of J cycle for the each subdomain.
							   //
#define imm_global (imt_global-1)
#define imm  (imt-1)
#define jmm  (jmt-1)
#define kmp1 (km+1)
#define kmm1 (km-1)
#define ntra 2    // Number of Tracers

extern int param_mod_mp_mytid_ ;
extern int i, j, k, m, n, ierr, mytid;
extern int jj_start, jj_end; // No Overlaped j ,North to South
extern int param_mod_mp_my_task_, param_mod_mp_master_task_;

// double am_hor; // Horizontal laplacian viscosity 
// double ah_hor; // Horizontal laplacian diffusity
// double am_bihar; // Horizontal biharmonic viscosity 
// double ah_bihar; // Horizontal biharmonic diffusity
// extern int num_cpl; // couple frequency 

					//lhl090729
#define s_imt 640
#define s_jmt 320
					//lhl090729
					//YU  01/02/2013

#endif //INCLUDED_PARAM
