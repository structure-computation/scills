#include "definition_PARAM_AFFICHAGE.h"
#include "definition_PARAM_STRUCTURE.h"
#include "definition_PARAM_MULTI.h"
#include "definition_PARAM_LATIN.h"
#include "definition_PARAM_TEMPS.h"
#include "definition_PARAM_PROPERTY.h"
#include "definition_PARAM_MPI.h"

using namespace LMT;
using namespace std;

inline void allocation_memoire_PARAM(Param &process){
     process.affichage = new AFFICHAGE;
     process.structure = new STRUCTURE;
     process.latin = new LATIN;
     process.multiscale = new MULTI;
     process.temps = new TEMPS;
     process.properties = new PROPERTY;
     process.multi_mpi = new MULTI_MPI;
#ifdef PRINT_ALLOC
     total_allocated[ typeid(AFFICHAGE).name() ] += sizeof(AFFICHAGE);
     total_allocated[ typeid(STRUCTURE).name() ] += sizeof(STRUCTURE);
     total_allocated[ typeid(LATIN).name() ] += sizeof(LATIN);
     total_allocated[ typeid(MULTI).name() ] += sizeof(MULTI);
     total_allocated[ typeid(TEMPS).name() ] += sizeof(TEMPS);
     total_allocated[ typeid(PROPERTY).name() ] += sizeof(PROPERTY);
     total_allocated[ typeid(MULTI_MPI).name() ] += sizeof(MULTI_MPI);
#endif
}

inline void desallocation_memoire_PARAM(Param &process){
   if (process.affichage != NULL) delete process.affichage;
   if (process.structure != NULL) delete process.structure;
   if (process.latin != NULL) delete process.latin;
   if (process.multiscale != NULL) delete process.multiscale;
   if (process.temps != NULL) delete process.temps;
   if (process.properties != NULL) delete process.properties;
   if (process.multi_mpi != NULL) delete process.multi_mpi;
}


