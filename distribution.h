#ifndef INCLUDE_DISTRIBUTION
#define INCLUDE_DISTRIBUTION

#include "param_mod.h"
#include "LICOM_Error_mod.h"
#include "blocks.h"

struct distrb
{
	int nprocs, communicator;
	int *proc, *local_block;
};

#endif // !INCLUDE_DISTRIBUTION
