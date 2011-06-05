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


void write_xdmf_geometry_fields(String output_xdmf, String input_hdf5, String name_geometry, String name_fields, Param &process, bool with_attributes); //Ecriture du fichier xdmf "output_xdmf" avec références aux données du fichier hdf5 "input_hdf5" 

