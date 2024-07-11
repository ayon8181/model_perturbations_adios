FC = mpiifort
# SPECFEM_DIR = ${HOME}/src/mideast_specfem3d_globe
SPECFEM_DIR = ${HOME}/specfem_main/specfem3d_globe
# SPECFEM_DIR = ${HOME}/src/specfem3d_globe_devel
SPECFEM_OBJ_DIR = ${SPECFEM_DIR}/obj
ADIOS1_DIR = /work/07923/orsvuran/shared/adios1/
ADIOS2_DIR = /work/07923/orsvuran/shared/adios2/
ADIOS_LIBS = $(shell ${ADIOS1_DIR}/bin/adios_config -lf)
ADIOS_FLAGS = $(shell ${ADIOS1_DIR}/bin/adios_config -cf)
SPECFEM_INC = -I${SPECFEM_DIR}/setup -I/scratch1/09038/ayon8181/s40rts_deep/OUTPUT_FILES
LIBS =  ${ADIOS_LIBS}  -lstdc++
INCLUDE = ${ADIOS_FLAGS} ${SPECFEM_INC}
FLAGS = -O3 -warn all,noexternal ${INCLUDE} ${LIBS}

SPECFEM_OBJECTS = \
	${SPECFEM_OBJ_DIR}/exit_mpi.shared.o \
	${SPECFEM_OBJ_DIR}/flush_system.shared.o \
	${SPECFEM_OBJ_DIR}/shared_par.shared_module.o \
	${SPECFEM_OBJ_DIR}/parallel.sharedmpi.o	\
	${SPECFEM_OBJ_DIR}/adios_helpers_addons.shared_adios_cc.o \
	${SPECFEM_OBJ_DIR}/adios_helpers.shared_adios.o \
	${SPECFEM_OBJ_DIR}/adios_helpers_definitions.shared_adios.o \
	${SPECFEM_OBJ_DIR}/adios_helpers_readers.shared_adios.o \
	${SPECFEM_OBJ_DIR}/adios_helpers_writers.shared_adios.o \
	${SPECFEM_OBJ_DIR}/adios_manager.shared_adios_module.o \
	${SPECFEM_OBJ_DIR}/param_reader.cc.o \
	${SPECFEM_OBJ_DIR}/read_parameter_file.shared.o \
	${SPECFEM_OBJ_DIR}/read_value_parameters.shared.o


%.o: %.f90
	${FC} -c $< ${FLAGS} -I${SPECFEM_OBJ_DIR}


all: \
	bin/xmodel_add \
	bin/xmodel_add_qmu \
	bin/xmodel_pert

	

#

bin/xmodel_add: src/model_pert_ref.F90
	${FC} -o $@ $^ \
	${SPECFEM_OBJ_DIR}/gll_library.shared.o \
	${SPECFEM_OBJECTS} -I${SPECFEM_OBJ_DIR} ${FLAGS}


bin/xmodel_add_qmu: src/model_pert_ref.F90
	${FC} -DUSE_QMU -o $@ $^ \
	${SPECFEM_OBJ_DIR}/gll_library.shared.o \
	${SPECFEM_OBJECTS} -I${SPECFEM_OBJ_DIR} ${FLAGS}

bin/xmodel_pert: src/model_perturb.f90 
	${FC} -o $@ $^ \
	${SPECFEM_OBJECTS} -I${SPECFEM_OBJ_DIR} ${FLAGS} 

clean:
	rm -rf obj/*.o bin/* *__genmod* *.mod




