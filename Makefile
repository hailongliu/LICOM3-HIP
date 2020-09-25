INC  = -I. -I/public/software/mathlib/netcdf/4.4.1/intel/include -I/opt/hpc/software/mpi/hpcx/v2.4.1/intel-2017.5.239/include
LIB  = -L. -L/public/software/mathlib/netcdf/4.4.1/intel/lib -lnetcdf -lnetcdff -L/opt/rocm/lib -lhip_hcc -lstdc++
#FFLAGS = -O2 -g -traceback -mcmodel=large -convert big_endian -assume byterecl
FFLAGS  = $(CPPDEFS) -O2 -r8 -i4 -free -fp-model precise -convert big_endian -assume byterecl -ftz -traceback
#FFLAGS  = $(CPPDEFS) -O2 -g -r8 -i4 -free -fp-model precise -convert big_endian -assume byterecl -no-vec -ftz -traceback
#FFLAGS = $(CPPDEFS) -O2 -g -r8 -i4 -free -fp-model precise -convert big_endian -assume byterecl -ftz -traceback
#FFLAGS = -O2 -mcmodel=large -shared-intel -convert big_endian -assume byterecl
#FFLAGS = -O2 -mcmodel=large -shared-intel -convert big_endian -assume byterecl -fp-model precise -fp-speculation=safe -no-vec 
FC  = mpif90  -m64 -mcmodel=medium $(FFLAGS)
#OBJS = mk_data_lhl1.o 
#LIB = -L/work1/jjr/soft/netcdf/lib -lnetcdf  -L/soft/mpi/openmpi/1.6.5/intel/lib -lmpi  
#FC  = $(MPIFC)
CC  = mpicc
PRG = run/licom3.exe
#include ./Macros
OBJ = shr_kind_mod.o shr_log_mod.o shr_mpi_mod.o shr_const_mod.o shr_sys_mod.o precision_mod.o param_mod.o msg_mod.o LICOM_Error_mod.o POP_BlocksMod.o \
      POP_CommMod.o POP_SpaceCurveMod.o POP_DistributionMod.o  POP_FieldMod.o  POP_GridHorzMod.o  POP_ReductionsMod.o POP_HaloMod.o \
			control_mod.o broadcast.o blocks.o distribution.o domain.o pconst_mod.o pmix_mod.o constant_mod.o cdf_mod.o gather_scatter.o global_reductions.o \
      grid.o tracer_mod.o forc_mod.o canuto_mod_2002.o licom_drv.o dyn_mod.o work_mod.o boundary.o operators.o output_mod.o diag_mod.o grids.o buf_mod.o hmix_del2.o \
			hmix_del4.o inirun.o isopyc_mod.o smuvh.o barotr.o bclinc.o accumm.o advection.o tracer.o energy.o ssave-cdf-month.o \
			output_1dto4d.o ssave-cdf-instant.o clock0f.o readyt.o icesnow.o convadj.o addps.o nextstep.o vinteg.o utils.o isopyi.o isopyc.o rdriver.o setidealvalue.o \
			upwell.o canuto_ini_2002.o canuto_2002.o readyc.o density.o dens.o thermalexpansion.o yy00.o invtri.o mm00.o k1_3.o k2_3.o k3_123.o isoadv.o isoflux.o jra_daily.o \
			cuda_init.o steponinit.o turb_init.o stepon.o energyMemcpy.o energy_c.o mpi1.o mpi2.o  mpi3.o pop_haloupdate_barotr1.o  pop_haloupdate_bclinc3.o  pop_haloupdate_tracer1.o \
			pop_haloupdate_barotr2.o  pop_haloupdate_readyc.o   pop_haloupdate_tracer2.o pop_haloupdate_bclinc1.o  pop_haloupdate_smts.o \
			pop_haloupdate_bclinc2.o  pop_haloupdate_smuv_3d.o allocate_fortran.o allocate_common.o get_blocksinfo.o  output_energy.o dailyavgini.o
$(PRG):$(OBJ)  
	$(FC) $(INC) -o $@ $(OBJ) $(LIB) 
		      
.SUFFIXES: .F90 .F .f90 .f .c .cpp .o
.c.o:  
	$(CC) $(CC_FLAG) $(INC) -c $*.c -o $*.o  
.cpp.o:  
	hipcc -c -O3 -std=c++11 -fno-gpu-rdc -amdgpu-target=gfx906 $(INC) -c $*.cpp -o $*.o  
.F90.o:				  
	$(FC) $(FFLAGS) $(INC) -c $*.F90 -o $*.o  
clean:  
	rm -f $(OBJ) $(PRG) *.mod
