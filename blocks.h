/*  by zly 2018.11 */
#ifndef INCLUDE_BLOCKS
#define INCLUDE_BLOCKS

#include "precision_mod.h"
#include "param_mod.h"
#include "LICOM_Error_mod.h"
//Fix block equiv by wpf, 2019.1.18
struct block 
{
	int block_id;
	int local_id;
	int ib, ie, jb, je;
	int iblock, jblock; 
	int i_glob[18], j_glob[18];
};

extern int blocks_mp_nblocks_tot_;
extern struct block (*blocks_mp_all_blocks_);

extern void blocks_mp_create_blocks_(int *, int *, char *, char *);
//extern struct block blocks_mp_get_block_(int *, int *);
extern void blocks_mp_get_block_parameter_(int *, int *, int *, int *, int *, int *, int *, int *, int [], int []);
extern void blocks_mp_destroy_blocks_();

extern struct block* get_block(int, int);

#endif // !INCLUDE_BLOCKS
