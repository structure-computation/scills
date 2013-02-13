#parametres a modifier si necessaire
#nombre de processeurs pour la compilation
NB_COMP_PROC = 3
#dimension du probleme
DIM = 3

# emplacement de la libraire MPI
MPI_LIB = /usr/include/openmpi
# type de machine
ARCHITECTURE = CPU


LOC_MC = metil_comp 

#dossiers sources
DIR_SOURCES_SC = -Ibuild -Ibuild/problem_pb_elast -Isrc -Isrc/DEFINITIONS -Isrc/FORMULATIONS -Isrc/ITERATIONS -Isrc/ITERATIONS/LINEAR -Isrc/ITERATIONS/LOCAL -Isrc/ITERATIONS/ERROR -Isrc/MAILLAGE -Isrc/MATERIAUX -Isrc/MPI -Isrc/OPERATEURS -Isrc/OPERATEURS/INTER -Isrc/OPERATEURS/MACRO -Isrc/OPERATEURS/SST -Isrc/POSTTRAITEMENTS -Isrc/PROBMICRO -Isrc/UTILITAIRES 
DIR_SOURCES_GEOMETRY = -Isrc -Isrc/GEOMETRY -Isrc/COMPUTE -Isrc/UTILS -Isrc/UTILS/hdf -Isrc/UTILS/xdmf -Isrc/UTILS/json_spirit 

#dossiers sources pour les client SOCA
DIR_SOURCES_SOCA = -IscillsResultClient -I/home/jbellec/code_dev_scwal/Soca/src -I/home/jbellec/code_dev_scwal/Soca/src/Soca  

#dossiers sources pour scVisu
DIR_SOURCES_SCVISU = -I/home/jbellec/code_dev/scVisu/src

# options communes
GLOB_VAR = -DDIM=$(DIM) -D$(ARCHITECTURE)  -DTYPE=double -DTYPEREEL=double  -DLDL -Dcrout_alain
CFLAGS= -LUTIL/metis -lmetis -L$(MPI_LIB)/lib -lmpi -lmpi_cxx
DIR_SOURCES_LMT =  -ILMT -ILMT/include -Iusr/include/suitesparse
DIR_SOURCES_CUDA = -Iusr/local/cuda/include -Ihome/ubuntu/driver_toolkit/NVIDIA_GPU_Computing_SDK/C/common/inc 

DIR_SOURCES_MPI = -I$(MPI_LIB) -I$(MPI_LIB)/include

# options pour la version release
PRG_multi_release = SC_multi_$(DIM)_$(ARCHITECTURE).exe
DIR_build_release = --comp-dir build/SC_multi_$(DIM)_$(ARCHITECTURE)_release
OPT = -ne -j$(NB_COMP_PROC) -O3 -ffast-math -fexpensive-optimizations -Wno-deprecated

# options pour la version debug
PRG_multi_debug = SC_multi_$(DIM)_$(ARCHITECTURE)_debug.exe
DIR_build_debug = --comp-dir build/SC_multi_$(DIM)_$(ARCHITECTURE)_debug
OPT_DBG = -ne -j$(NB_COMP_PROC) -g -g3 -ffast-math -fexpensive-optimizations


# option pour le test de chargement des donnees
PRG_test_load_data = SC_test_load_data_$(DIM).exe


all: metil_comp_multi

compil_py_files : scons -j1 dep_py=1

update: metil_comp_multi metil_comp_multi_dbg

metil_comp_test_load_data :
	$(LOC_MC)  -o  $(PRG_test_load_data) $(GLOB_VAR) $(DIR_SOURCES_LMT) $(DIR_SOURCES_SC) $(DIR_SOURCES_GEOMETRY) $(DIR_SOURCES_MPI) $(DIR_build_release_cpu) $(CFLAGS) $(LIBS) $(OPT)  test/test_load_data.cpp

metil_comp_multi :

	$(LOC_MC)  -o  $(PRG_multi_release) $(GLOB_VAR) $(DIR_SOURCES_LMT) $(DIR_SOURCES_SC) $(DIR_SOURCES_GEOMETRY) $(DIR_SOURCES_MPI) $(DIR_build_release) $(CFLAGS) $(LIBS) $(OPT)  src/main.cpp

metil_comp_multi_dbg :
	$(LOC_MC)  -o  $(PRG_multi_debug) $(GLOB_VAR) $(DIR_SOURCES_LMT) $(DIR_SOURCES_SC) $(DIR_SOURCES_GEOMETRY) $(DIR_SOURCES_MPI) $(DIR_build_debug) $(CFLAGS) $(LIBS) $(OPT_DBG)  src/main.cpp

codegen_py:
	cd LMT/include/codegen; scons
	scons -j1 dep_py=1 

clean:
	scons -c
	cd LMT/include/codegen; scons -c

codegen_json:
	$(LOC_MC) -o DataUserUpdater.exe $(DIR_SOURCES_LMT) -Ihome/scproduction/code_dev/Metil-test/src $(DIR_SOURCES_SC) $(DIR_SOURCES_GEOMETRY) $(DIR_SOURCES_MPI) --comp-dir build/DataUserUpdater $(CFLAGS) src/SCJSONREADER/run.cpp

scills_result_client:
	$(LOC_MC)  -o scillsResultClient.exe $(GLOB_VAR) $(DIR_SOURCES_LMT) $(DIR_SOURCES_SC) $(DIR_SOURCES_GEOMETRY) $(DIR_SOURCES_MPI) $(CFLAGS) $(LIBS) $(OPT) $(DIR_SOURCES_SCVISU) $(DIR_SOURCES_SOCA) scillsResultClient/main.cpp
