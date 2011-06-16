#include <Metil/BasicVec.h>
#include <Metil/StructCompactor.h>
#include <Metil/Hdf.h>
using namespace Metil;

#include "containers/mat.h"
using namespace LMT;
#include "definition_PARAM.h"
#include "definition_PARAM_LATIN.h"
#include "definition_PARAM_MULTI.h"
#include "definition_PARAM_TEMPS.h"
#include "definition_PARAM_AFFICHAGE.h"

//ecriture du fichier xdmf a partir de donnees hdf5 pour un maillage contenant des noeuds, elements, des attributs. Utilisable pour les sst ou les interfaces
void write_xdmf_geometry_fields(std::ofstream &f, Param &process, String name_element, String name_node, String generic_grid_name, String grid_collection_name, bool with_time, BasicVec<String> &attributs);
