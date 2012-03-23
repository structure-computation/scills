#ifndef MATERIAUX_DECLARATIONS_H
#define MATERIAUX_DECLARATIONS_H

#include "../COMPUTE/DataUser.h"
#include "../COMPUTE/FieldStructureUser.h"

#include "../DEFINITIONS/Param.h"
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
 *     - read_matprop() : lecture des propriétés matériau et construction du vecteur de propriétés SstCarac à partir du xml
 *     - assignation_material_to_SST : assignation des propriétés aux sst
 * - assignation_materials_property_INTER() : fonction principale pour les interfaces
 *     - read_propinter() : lecture des propriétés d'interfaces particulières à partir du xml
 *     - modif_inter() : application des propriétés aux interfaces
 *
 */
void assignation_materials_property_SST(DataUser &data_user, Vec<SstCarac> &matprops, Vec<Sst> &S, Param &process, FieldStructureUser &field_structure_user);
void assignation_materials_property_INTER(DataUser &data_user, Vec<Interface> &Inter, Param &process, FieldStructureUser &field_structure_user);


#endif //MATERIAUX_DECLARATIONS_H
