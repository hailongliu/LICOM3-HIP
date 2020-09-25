#ifndef INCLUDE_OPERATORS
#define INCLUDE_OPERATORS

#include "precision_mod.h"
#include "param_mod.h"
#include "blocks.h"

extern void operators_mp_div_(int *, double[][nx_block], double[][nx_block], double[][nx_block], struct block *);
extern void operators_mp_grad_(int *, double[][nx_block], double[][nx_block], double[][nx_block], struct block *);
extern void operators_mp_zcurl_(int *, double[][nx_block], double[][nx_block], double[][nx_block], struct block *);

// extern void operators_div_(int *, double[][nx_block], double[][nx_block], double[][nx_block], struct block *);
//extern void operators_div_(int *, double[][nx_block], double[][nx_block], double[][nx_block], int *);
//extern void operators_div_(int *, double[][nx_block], double[][nx_block], double[][nx_block], int *);
extern void operators_grad_(int *, double[][nx_block], double[][nx_block], double[][nx_block], int *);

#endif // !INCLUDE_OPERATORS
