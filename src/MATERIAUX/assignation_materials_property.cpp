// Author: D. Violeau , 24/03/2006

//librairies Hugo
#include "containers/mat.h"
#include "containers/vecpointedvalues.h"
//#include "mesh/mesh.h"
#include "containers/vec_mt.h"

//fichiers de definition des variables 
#include "Param.h"
#include "definition_PARAM_PROPERTY.h"
#include "definition_PARAM_TEMPS.h"
#include "definition_PARAM_COMP_INTER.h"
#include "definition_materials_property.h"
#include "definition_SST_time.h"
#include "definition_INTER_time.h"

#include "utilitaires.h"

#include "read_data_materials_SST_INTER.h"
#include "assignation_materiaux_sst.h"
#include "assignation_comportements_speciaux_inter.h"

#include "definition_fcts.h"
#include "DataUser.h"
#include "FieldStructureUser.h"

using namespace LMT;
using namespace std;

/** \defgroup Materiaux Lecture et assignation des comportements de sst ou interface

Dans la fonction principale, on appelle successivement :
- assignation_materials_property_SST() : fonction principale pour les sst
   - read_material_properties() : lecture des propriétés matériau et construction du vecteur de propriétés SstCarac à partir du xml
   - assignation_material_to_SST : assignation des propriétés aux sst
- assignation_materials_property_INTER() : fonction principale pour les interfaces
   - read_propinter() : lecture des propriétés d'interfaces particulières à partir du xml
   - modif_inter() : application des propriétés aux interfaces
   
*/



void assignation_materials_property_SST(DataUser &data_user, Vec<SstCarac> &matprops, Vec<Sst> &S, Vec<Interface> &Inter,Param &process, FieldStructureUser &field_structure_user){
    if (process.rank == 0) std::cout << "\t Assignation du materiau aux SST" << std::endl;
    // lecture des proprietes materiau des ssts

    BasicVec<BasicVec<TYPEREEL > > mat_prop_temp;
    read_material_properties(matprops,process,data_user, mat_prop_temp); 
    //assignation des proprietes aux SST
    apply_mt(S,process.nb_threads,assignation_material_to_SST(),matprops);        
    //assignation des proprietes aux group_elements (pour GPU)
    field_structure_user.assign_material_id_to_group_elements(data_user);
    field_structure_user.assign_material_properties_to_group_elements(data_user,mat_prop_temp);
}

void assignation_materials_property_INTER(DataUser &data_user, Vec<Interface> &Inter, Vec<Sst> &S, Param &process, FieldStructureUser &field_structure_user){
    if (process.rank == 0) std::cout << "\t Assignation de materiau particulier (Ex : contact) aux Interfaces" << std::endl;
    Vec<InterCarac<DIM,TYPEREEL> > propinter;
    BasicVec<BasicVec<TYPE > > link_prop_temp;
    read_propinter(propinter,data_user, link_prop_temp);
    PRINT("lecture propriete interface ok");

    modif_inter(Inter,propinter,S,process,data_user);
    //idem pour les group_interfaces (pour GPU)
    field_structure_user.assign_link_id_to_group_interfaces(data_user);
    field_structure_user.assign_link_properties_to_group_interfaces(data_user,link_prop_temp);
}
