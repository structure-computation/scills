#parametres a modifier si necessaire
DIM = 3

DIR_SOURCES_SC = -Ibuild -Ibuild/problem_pb_elast -Isrc -Isrc/DEFINITIONS -Isrc/FORMULATIONS -Isrc/ITERATIONS -Isrc/ITERATIONS/LINEAR -Isrc/ITERATIONS/LOCAL -Isrc/ITERATIONS/ERROR -Isrc/MAILLAGE -Isrc/MATERIAUX -Isrc/MPI -Isrc/OPERATEURS -Isrc/OPERATEURS/INTER -Isrc/OPERATEURS/MACRO -Isrc/OPERATEURS/SST -Isrc/POSTTRAITEMENTS -Isrc/PROBMICRO -Isrc/UTILITAIRES 
DIR_SOURCES_GEOMETRY = -Isrc -Isrc/GEOMETRY -Isrc/COMPUTE -Isrc/UTILS -Isrc/UTILS/hdf -Isrc/UTILS/xdmf -Isrc/UTILS/json_spirit 
#DIR_SOURCES_COMPUTE = -Isrc -Isrc/GEOMETRY -Isrc/COMPUTE -Isrc/UTILS -Isrc/UTILS/hdf -Isrc/UTILS/xdmf  -Isrc/UTILS/json_spirit 

PRG_multi = SC_multi_$(DIM).exe
PRG_create = SC_create_2_cpu_$(DIM).exe
DIR_build_cpu = --comp-dir build/SC_$(DIM)

LOC_MC = metil_comp 
CFLAGS= -LUTIL/metis -lmetis -LUTIL/openmpi/lib -lmpi -lmpi_cxx
DIR_SOURCES_LMT =  -ILMT -ILMT/include -Iusr/include/suitesparse
DIR_SOURCES_CUDA = -Iusr/local/cuda/include -Ihome/ubuntu/driver_toolkit/NVIDIA_GPU_Computing_SDK/C/common/inc 
DIR_SOURCES_MPI = -IUTIL/openmpi -IUTIL/openmpi/include
OPT = -ne -j2 -O3 -ffast-math -fexpensive-optimizations

# all: compact_GEOMETRY 
# all: metil_comp_create_cpu rsync

# all: metil_comp_compute_cpu rsync
all: metil_comp_multi
# all: metil_comp_test
# all: local

# all: compact_FIELD_STRUCTURE 

#codegen local rsync

metil_test :
	metil src/CALCUL/code_metil/test.met

metil_comp_multi :
	$(LOC_MC)  -o  $(PRG_multi) -DCPU  -DDIM=$(DIM) -DCPU  -DTYPE=double -DTYPEREEL=double  -DLDL -Dcrout_alain $(DIR_SOURCES_LMT) $(DIR_SOURCES_SC) $(DIR_SOURCES_GEOMETRY) $(DIR_SOURCES_MPI) $(DIR_build_cpu) $(CFLAGS) $(LIBS) $(OPT)  src/multiscale.cpp

metil_comp_test :
	$(LOC_MC)  -o  $(PRG_multi) -DCPU  -DDIM=$(DIM) -DCPU  -DTYPE=double -DTYPEREEL=double -DLDL -Dcrout_alain $(DIR_SOURCES_LMT) $(DIR_SOURCES_SC) $(DIR_SOURCES_GEOMETRY) $(DIR_SOURCES_MPI) $(DIR_build_cpu) $(CFLAGS) $(LIBS) $(OPT)  src/test.cpp

metil_comp_create_cpu :
	$(LOC_MC)  -o  $(PRG_create) -DDIM=$(DIM) -DCPU  -DTYPE=double -DLDL -DWITH_CHOLMOD -DWITH_UMFPACK $(DIR_SOURCES_LMT) $(DIR_SOURCES_SC) $(DIR_SOURCES_GEOMETRY) $(CFLAGS) $(LIBS) $(OPT)  src/SC_create_2.cpp

# codegen_py:
# 	cd LMT/include/codegen; scons

local:  
	scons -j1 dep_py=1 

# clean:
# 	scons -c
# 	cd LMT/include/codegen; scons -c

rsync : 
# 	rsync -r --exclude '.git' --exclude 'EXEMPLES'  --exclude 'LMT'  --exclude 'UTIL*'  --exclude 'src/COMPUTE'  --exclude 'src/GEOMETRY'  --exclude 'src/SC_create_2.cpp'  --exclude 'src/UTILS'. /home/scproduction/Developpement/SC_code
	rsync SC_multi_*.exe SC_create_* /home/scproduction/Developpement/SC_multi





