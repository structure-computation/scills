#parametres a modifier si necessaire
#dimension du probleme
DIM = 3
#nombre de processeurs pour la compilation
NB_COMP_PROC = 3

PRG_test_load_data = SC_test_load_data_$(DIM).exe
PRG_multi_release = SC_multi_$(DIM).exe
PRG_multi_debug = SC_multi_$(DIM)_debug.exe
DIR_build_release_cpu = --comp-dir build/SC_$(DIM)_release
DIR_build_debug_cpu = --comp-dir build/SC_$(DIM)_debug

#dossiers sources
DIR_SOURCES_SC = -Ibuild -Ibuild/problem_pb_elast -Isrc -Isrc/DEFINITIONS -Isrc/FORMULATIONS -Isrc/ITERATIONS -Isrc/ITERATIONS/LINEAR -Isrc/ITERATIONS/LOCAL -Isrc/ITERATIONS/ERROR -Isrc/MAILLAGE -Isrc/MATERIAUX -Isrc/MPI -Isrc/OPERATEURS -Isrc/OPERATEURS/INTER -Isrc/OPERATEURS/MACRO -Isrc/OPERATEURS/SST -Isrc/POSTTRAITEMENTS -Isrc/PROBMICRO -Isrc/UTILITAIRES 
DIR_SOURCES_GEOMETRY = -Isrc -Isrc/GEOMETRY -Isrc/COMPUTE -Isrc/UTILS -Isrc/UTILS/hdf -Isrc/UTILS/xdmf -Isrc/UTILS/json_spirit 

LOC_MC = metil_comp 
CFLAGS= -LUTIL/metis -lmetis -LUTIL/openmpi/lib -lmpi -lmpi_cxx
DIR_SOURCES_LMT =  -ILMT -ILMT/include -Iusr/include/suitesparse
DIR_SOURCES_CUDA = -Iusr/local/cuda/include -Ihome/ubuntu/driver_toolkit/NVIDIA_GPU_Computing_SDK/C/common/inc 
DIR_SOURCES_MPI = -IUTIL/openmpi -IUTIL/openmpi/include
OPT = -ne -j$(NB_COMP_PROC) -O3 -ffast-math -fexpensive-optimizations
OPT_DBG = -ne -j$(NB_COMP_PROC) -g -g3 -ffast-math -fexpensive-optimizations
GLOB_VAR = -DCPU -DDIM=$(DIM) -DCPU  -DTYPE=double -DTYPEREEL=double  -DLDL -Dcrout_alain


all: clean codegen_py metil_comp_multi

metil_comp_test_load_data :
	$(LOC_MC)  -o  $(PRG_test_load_data) $(GLOB_VAR) $(DIR_SOURCES_LMT) $(DIR_SOURCES_SC) $(DIR_SOURCES_GEOMETRY) $(DIR_SOURCES_MPI) $(DIR_build_release_cpu) $(CFLAGS) $(LIBS) $(OPT)  test/test_load_data.cpp


metil_comp_multi :
	$(LOC_MC)  -o  $(PRG_multi_release) $(GLOB_VAR) $(DIR_SOURCES_LMT) $(DIR_SOURCES_SC) $(DIR_SOURCES_GEOMETRY) $(DIR_SOURCES_MPI) $(DIR_build_release_cpu) $(CFLAGS) $(LIBS) $(OPT)  src/main.cpp

metil_comp_multi_dbg :
	$(LOC_MC)  -o  $(PRG_multi_debug) $(GLOB_VAR) $(DIR_SOURCES_LMT) $(DIR_SOURCES_SC) $(DIR_SOURCES_GEOMETRY) $(DIR_SOURCES_MPI) $(DIR_build_debug_cpu) $(CFLAGS) $(LIBS) $(OPT_DBG)  src/main.cpp

codegen_py:
	cd LMT/include/codegen; scons
	scons -j1 dep_py=1 

clean:
	scons -c
	cd LMT/include/codegen; scons -c
