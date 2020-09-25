#include "param_mod.h"
#include "blocks.h"
#include "pconst_mod.h"
#include "grid.h"
#include "dyn_mod.h"
#include "constant_mod.h"
#include "domain.h"
#include "tracer_mod.h"

extern int msg_mod_mp_mpi_comm_ocn_;
extern void get_blocksinfo_(int*,int*,int*,int*);
extern void output_energy(int*,double*,double*,double*,double*,double*);

void energy() {
	int i, j, k;
	int ierr;

	double ek0, ea0, eb0, et0, es0, volume, volume0, nnp, nnp0,
		ek, ea, eb, et, es, dxx,
		t0, t1, clock0f, nst;
	int nnn, iblock=0;
	int jb,je,ib,ie;

//	struct block *this_block;

	ek = 0.0;
	pconst_mod_mp_month_ = (pconst_mod_mp_iyfm_ - 1) * 12 + pconst_mod_mp_mon0_;

	eb = 0.0;
	nst = 0.0;

	ek = 0.0;
	pconst_mod_mp_month_ = (pconst_mod_mp_iyfm_ - 1) * 12 + pconst_mod_mp_mon0_;

//	for (iblock = 0; iblock < domain_mp_nblocks_clinic_; iblock++) {
//		this_block = get_block(domain_mp_blocks_clinic_[iblock], iblock);
		get_blocksinfo_(&jb,&je,&ib,&ie);
		for (k = 0; k < km; k++) {
			for (j = jb - 1; j < je; j++) {
				for (i = ib - 1; i < ie; i++) {
					ek += pconst_mod_mp_dzp_[k] * grid_mp_uarea_[iblock][j][i] * pconst_mod_mp_viv_[iblock][k][j][i] * (dyn_mod_mp_u_[iblock][k][j][i] * dyn_mod_mp_u_[iblock][k][j][i] + dyn_mod_mp_v_[iblock][k][j][i] * dyn_mod_mp_v_[iblock][k][j][i]);
				}
			}
		}
//	}

	ea = 0.0;
//	for (iblock = 0; iblock < domain_mp_nblocks_clinic_; iblock++) {
//		this_block = get_block(domain_mp_blocks_clinic_[iblock], iblock);

		get_blocksinfo_(&jb,&je,&ib,&ie);
		for (j = jb - 1; j < je; j++) {
			for (i = ib - 1; i < ie; i++) {
				ea += grid_mp_tarea_[iblock][j][i] * pconst_mod_mp_vit_[iblock][0][j][i] * dyn_mod_mp_h0_[iblock][j][i] * dyn_mod_mp_h0_[iblock][j][i];
			}
		}
//	}

	eb = 0.0;
	nst = 0.0;
//	for (iblock = 0; iblock < domain_mp_nblocks_clinic_; iblock++) {
//		this_block = get_block(domain_mp_blocks_clinic_[iblock], iblock);

		get_blocksinfo_(&jb,&je,&ib,&ie);
		for (j = jb - 1; j < je; j++) {
			for (i = ib - 1; i < ie; i++) {
				eb += grid_mp_tarea_[iblock][j][i] * pconst_mod_mp_vit_[iblock][0][j][i] * tracer_mod_mp_at_[iblock][0][0][j][i];
			}
		}
//	}

	et = 0.0;
	es = 0.0;
//	for (iblock = 0; iblock < domain_mp_nblocks_clinic_; iblock++) {
//		this_block = get_block(domain_mp_blocks_clinic_[iblock], iblock);

		get_blocksinfo_(&jb,&je,&ib,&ie);
		for (k = 0; k < km; k++) {
			for (j = jb - 1; j < je; j++) {
				for (i = ib - 1; i < ie; i++) {
					et += pconst_mod_mp_dzp_[k] * grid_mp_tarea_[iblock][j][i] * pconst_mod_mp_vit_[iblock][k][j][i] * tracer_mod_mp_at_[iblock][0][k][j][i];
					es += pconst_mod_mp_dzp_[k] * grid_mp_tarea_[iblock][j][i] * pconst_mod_mp_vit_[iblock][k][j][i] * tracer_mod_mp_at_[iblock][1][k][j][i];
				}
			}
		}
//	}

	mpi1_(&ek, &ek0);
//	mpi_reduce(ek, ek0, 1, MPI_PR, MPI_SUM, 0, msg_mod_mp_mpi_comm_ocn_);
	ek0 = ek0 * 0.5e0 * pconst_mod_mp_d0_;

	mpi1_(&ea, &ea0);
//	mpi_reduce(ea, ea0, 1, MPI_PR, MPI_SUM, 0, msg_mod_mp_mpi_comm_ocn_);
	ea0 = ea0 * 0.5e0 * pconst_mod_mp_d0_ * g;

	mpi1_(&eb, &eb0);
//	mpi_reduce(eb, eb0, 1, MPI_PR, MPI_SUM, 0, msg_mod_mp_mpi_comm_ocn_);
	eb0 = eb0 / grid_mp_area_t_;

	mpi1_(&et, &et0);
	mpi1_(&es, &es0);
//	mpi_reduce(et, et0, 1, MPI_PR, MPI_SUM, 0, msg_mod_mp_mpi_comm_ocn_);
//	mpi_reduce(es, es0, 1, MPI_PR, MPI_SUM, 0, msg_mod_mp_mpi_comm_ocn_);
	et0 = et0 / grid_mp_volume_t_;
	es0 = es0 / grid_mp_volume_t_ * 1.0e+03;
	
	output_energy_(&pconst_mod_mp_month_,&ek0,&ea0,&eb0,&et0,&es0);	
//        WRITE (16,FMT='(I5,I3,6D25.15)') MONTH,IDAY,EK0,EA0,EB0,ET0,ES0,t1-t0
	if(isnan(ek0)) {
	printf("ek0 is NAN\n");
	exit(0);
	}
	return;
}
