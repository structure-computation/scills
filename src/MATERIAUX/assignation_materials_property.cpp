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

//seulement pour le cas d'une formulation endommageable
#ifdef FORMUENDO
#include "create_mesh_meso.h"
#endif

void fake_assignation_materials_property() {
    XmlNode n;
    Param process;
    DataUser data_user;
        
    Vec<Sst<DIM,TYPEREEL> > S3;
    Vec<VecPointedValues<Sst<DIM,TYPEREEL> > > SubS3;
    Vec<Interface<DIM,TYPEREEL> > Inter3;
    
    assignation_materials_property_SST(data_user, S3, Inter3,process);
    assignation_materials_property_INTER(data_user, Inter3, S3, process);
    
    assignation_materials_property_SST(n,S3,Inter3,process); 
    assignation_materials_property_INTER(n, Inter3, S3, process);

}


template <class TV1, class TV2>
void assignation_materials_property_SST(DataUser &data_user, TV1 &S, TV2 &Inter,Param &process){
   if (process.rank == 0) std::cout << "\t Assignation du materiau aux SST" << std::endl;
   // lecture des proprietes materiau des ssts
   Vec<SstCarac<TV1::template SubType<0>::T::dim,TYPEREEL> > matprop;
   read_material_properties(matprop,process,data_user); 
   apply_mt(S,process.nb_threads,assignation_material_to_SST(),matprop,process);
   
   //pour le mesomodele uniquement
   #ifdef FORMUENDO
   initialisation_mesh_meso(S, Inter, matprop);
   #endif
}

template <class TV1, class TV2>
void assignation_materials_property_INTER(DataUser &data_user, TV2 &Inter, TV1 &S, Param &process){
   if (process.rank == 0) std::cout << "\t Assignation de materiau particulier (Ex : contact) aux Interfaces" << std::endl;
   Vec<InterCarac<TV1::template SubType<0>::T::dim,TYPEREEL> > propinter;
   read_propinter(propinter,data_user);
//     for(int i_inter=0; i_inter<propinter.size(); i_inter++){
//         propinter[i_inter].affiche();
//     }
   modif_inter(Inter,propinter,S,process);
//    //verification
//     for(int i_inter=0; i_inter<Inter.size(); i_inter++){
//         Inter[i_inter].affiche();
//     }
}




template <class TV1, class TV2>
void assignation_materials_property_SST(const XmlNode &n, TV1 &S, TV2 &Inter,Param &process){
   if (process.rank == 0) std::cout << "\t Assignation du materiau aux SST" << std::endl;
   // lecture des proprietes materiau des ssts
   Vec<SstCarac<TV1::template SubType<0>::T::dim,TYPEREEL> > matprop;
   read_material_properties(matprop,process,n); 
   apply_mt(S,process.nb_threads,assignation_material_to_SST(),matprop,process);
   
   //pour le mesomodele uniquement
   #ifdef FORMUENDO
   initialisation_mesh_meso(S, Inter, matprop);
   #endif
}

template <class TV1, class TV2>
void assignation_materials_property_INTER(const XmlNode &n, TV2 &Inter, TV1 &S, Param &process){
   if (process.rank == 0) std::cout << "\t Assignation de materiau particulier (Ex : contact) aux Interfaces" << std::endl;
   Vec<InterCarac<TV1::template SubType<0>::T::dim,TYPEREEL> > propinter;
   read_propinter(propinter,n);
   modif_inter(Inter,propinter,S,process);
}


