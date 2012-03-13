// Author: D. Violeau , 24/03/2006

//librairies Hugo
#include "containers/mat.h"
#include "containers/vecpointedvalues.h"
//#include "mesh/mesh.h"
#include "containers/vec_mt.h"

//fichiers de definition des variables 
#include "definition_PARAM.h"
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
   - read_material_properties() : lecture des propri�t�s mat�riau et construction du vecteur de propri�t�s SstCarac � partir du xml
   - assignation_material_to_SST : assignation des propri�t�s aux sst
- assignation_materials_property_INTER() : fonction principale pour les interfaces
   - read_propinter() : lecture des propri�t�s d'interfaces particuli�res � partir du xml
   - modif_inter() : application des propri�t�s aux interfaces
   
*/

//seulement pour le cas d'une formulation endommageable
#ifdef FORMUENDO
#include "create_mesh_meso.h"
#endif

void fake_assignation_materials_property() {
/*    XmlNode n;*/
    Param process;
    DataUser data_user;
    FieldStructureUser field_structure_user;    
//     Vec<Splitted<Sst<DIM,TYPEREEL>, 16> > S3;
    Vec<Sst<DIM,TYPEREEL> > S3;
    Vec<VecPointedValues<Sst<DIM,TYPEREEL> > > SubS3;
    Vec<Interface<DIM,TYPEREEL> > Inter3;
    
    assignation_materials_property_SST(data_user, S3, Inter3,process, field_structure_user);
    assignation_materials_property_INTER(data_user, Inter3, S3, process, field_structure_user);
    
//     assignation_materials_property_SST(n,S3,Inter3,process); 
//     assignation_materials_property_INTER(n, Inter3, S3, process);

}


template <class TV1, class TV2>
void assignation_materials_property_SST(DataUser &data_user, TV1 &S, TV2 &Inter,Param &process, FieldStructureUser &field_structure_user){
       if (process.rank == 0) std::cout << "\t Assignation du materiau aux SST" << std::endl;
        // lecture des proprietes materiau des ssts
        Vec<SstCarac<TV1::template SubType<0>::T::dim,TYPEREEL> > matprop;
        BasicVec<BasicVec<TYPE > > mat_prop_temp;
        read_material_properties(matprop,process,data_user, mat_prop_temp); 
       //assignation des proprietes aux SST
        apply_mt(S,process.nb_threads,assignation_material_to_SST(),matprop);        
        //assignation des proprietes aux group_elements (pour GPU)
        field_structure_user.assign_material_id_to_group_elements(data_user);
        field_structure_user.assign_material_properties_to_group_elements(data_user,mat_prop_temp);

        //pour le mesomodele uniquement
        #ifdef FORMUENDO
        initialisation_mesh_meso(S, Inter, matprop);
        #endif
        
        
}

template <class TV1, class TV2>
void assignation_materials_property_INTER(DataUser &data_user, TV2 &Inter, TV1 &S, Param &process, FieldStructureUser &field_structure_user){
   if (process.rank == 0) std::cout << "\t Assignation de materiau particulier (Ex : contact) aux Interfaces" << std::endl;
   Vec<InterCarac<TV1::template SubType<0>::T::dim,TYPEREEL> > propinter;
   BasicVec<BasicVec<TYPE > > link_prop_temp;
   read_propinter(propinter,data_user, link_prop_temp);
   PRINT("lecture propriete interface ok");
//     for(int i_inter=0; i_inter<propinter.size(); i_inter++){
//         propinter[i_inter].affiche();
//     }
   //assignation des proprietes de liaisons aux interfaces
   modif_inter(Inter,propinter,S,process,data_user);
   //idem pour les group_interfaces (pour GPU)
   field_structure_user.assign_link_id_to_group_interfaces(data_user);
   field_structure_user.assign_link_properties_to_group_interfaces(data_user,link_prop_temp);
   
//    //verification
//     for(int i_inter=0; i_inter<Inter.size(); i_inter++){
//         Inter[i_inter].affiche();
//     }
}
