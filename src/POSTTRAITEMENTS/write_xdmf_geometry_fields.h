#include <Metil/BasicVec.h>
#include <Metil/StructCompactor.h>
#include <Metil/Hdf.h>
using namespace Metil;

#include "containers/mat.h"
using namespace LMT;
#include "Process.h"
#include "LatinData.h"
#include "MultiScaleData.h"
#include "TimeData.h"
#include "SavingData.h"

#include "../COMPUTE/DataUser.h"

//ecriture du fichier xdmf a partir de donnees hdf5 pour un maillage contenant des noeuds, elements, des attributs. Utilisable pour les sst ou les interfaces
void write_xdmf_geometry_fields(std::ofstream &f, Process &process, Sc2String name_element, Sc2String name_node, Sc2String generic_grid_name, Sc2String grid_collection_name, bool with_time, BasicVec<Sc2String> &attributs);
//ecriture de 3 fichiers xdmf (elements_0, elements_1, elements_0_skin)  en choisissant les composantes a mettre pour la geometrie
void write_xdmf_file_geometry(Process &process, DataUser &data_user);
//ecriture de 3 fichiers xdmf (elements_0, elements_1, elements_0_skin) en choisissant les composantes a mettre pour les calculs 
void write_xdmf_file_compute(Process &process, DataUser &data_user);
