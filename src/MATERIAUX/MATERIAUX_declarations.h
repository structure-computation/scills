#ifndef MATERIAUX_DECLARATIONS_H
#define MATERIAUX_DECLARATIONS_H

#include "../COMPUTE/DataUser.h"
#include "../COMPUTE/FieldStructureUser.h"

#include "Process.h"
#include "../DEFINITIONS/Sst.h"
#include "../DEFINITIONS/Interface.h"
#include "../DEFINITIONS/SstCarac_InterCarac.h"

#include "../../LMT/include/containers/vec_mt.h"


using namespace LMT;
using namespace Metil;

/** \defgroup Materiaux Lecture et assignation des comportements de sst ou interface
 * 
 * Dans la fonction principale, on appelle successivement :
 * - assignation_materials_property_SST() : fonction principale pour les sst
 *     - read_matprop() : lecture des proprietes materiau et construction du vecteur de proprietes SstCarac à partir du xml
 *     - assignation_material_to_SST : assignation des proprietes aux sst
 * - assignation_materials_property_INTER() : fonction principale pour les interfaces
 *     - read_propinter() : lecture des proprietes d'interfaces particulières à partir du xml
 *     - modif_inter() : application des proprietes aux interfaces
 *
 */
void assignation_materials_property_SST(DataUser &data_user, Vec<SstCarac> &matprops, Vec<Sst> &S, Process &process, FieldStructureUser &field_structure_user);
void assignation_materials_property_INTER(DataUser &data_user, Vec<Interface> &Inter, Process &process, FieldStructureUser &field_structure_user);


#endif //MATERIAUX_DECLARATIONS_H
