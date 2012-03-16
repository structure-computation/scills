#include "mise_a_jour_quantites.h"
#include <boost/concept_check.hpp>
#include "FieldStructureUser.h"


//Assignation des champs de géométrie d'éléments sur la peau à partir des champs des éléments parents correspondant
struct Projection_data_elements_0_on_skin_sst{
  template<class TE, class TMS, class TM> void operator()(TE &e, const TMS &mskin, const TM &m) const{
    const typename TM::EA *ea = mskin.get_parents_of(e)[0];
    typedef typename TM::TElemList::template SubType<0>::T TypeParent;
    const TypeParent &parent = static_cast<const TypeParent &>( *ea );
    e.numsst_skin = parent.numsst;
    e.num_proc_skin = parent.num_proc;
    e.typmat_skin = parent.typmat;
  }
};

//Assignation des champs de résultat d'éléments sur la peau à partir des champs des éléments parents correspondant
struct Projection_fields_on_skin_sst{
   template<class TE, class TMS, class TM> void operator()(TE &e, const TMS &mskin, const TM &m) const{
      const typename TM::EA *ea = mskin.get_parents_of(e)[0];
      typedef typename TM::TElemList::template SubType<0>::T TypeParent;
      const TypeParent &parent = static_cast<const TypeParent &>( *ea );
      e.sigma_skin = parent.sigma[0];
      e.sigma_mises_skin = parent.sigma_von_mises;
      e.epsilon_skin = parent.epsilon[0];
      e.numsst_skin = parent.numsst;
      e.num_proc_skin = parent.num_proc;
      e.typmat_skin = parent.typmat;

//#warning A modifier pour les comportements orthotrope
      #ifdef FORMUORTHO
        e.sigma_local_skin = parent.sigma_local[0];
      #endif
   }
};

//calcul du champ permettant de réaliser un éclaté de la structure (champ qtrans)
template<class SST>
void calcul_explode_displacements(SST &S){
    typedef Vec<TYPEREEL,DIM> TV;
    TV G=barycenter_constant_rho(*S.mesh.m);
    for(int i=0;i<S.mesh->node_list.size();i++){
        S.mesh->node_list[i].qtrans=S.mesh->node_list[i].pos+G;
    }
}

//Permet de calculer les champs sur la SST et sur sa peau pour un pas de temps donné
template<class TSST>
void calcul_fields_on_sst(TSST &S, Param &process, DataUser &data_user) {
     //assignation des deplacements a partir du deplacement au piquet de temps imic + calcul des champs a partir de ce deplacement
    if(process.nom_calcul=="incr")
        assign_dep_cont_slave(S,S.t[1].q, data_user); 
    else if(process.nom_calcul=="latin"){
        std::cout << "calcul_fields_on_sst non defini pour strategie latin pure " << std::endl; assert(0);
        //assign_dep_cont_slave(S,S.t[1].q);
    }
    else{std::cout << "Type de calcul non reconnu dans save_geometry_sst " << std::endl;assert(0);}
    
    calcul_explode_displacements(S);
    
    S.mesh->update_skin();
/*    apply(S.mesh->skin.elem_list,Projection_elements_0_skin(),S.mesh->skin, *S.mesh.m);*/
    apply(S.mesh->skin.elem_list,Projection_fields_on_skin_sst(),S.mesh->skin,*S.mesh.m);
}

//Conversion des champs d'éléments de peau à la structure de donnée de type BasicVec
struct Extract_fields_on_element_sst_skin{
   template<class TE, class TS> void operator()(TE &e, TS &S ) const{
        int nb_comp=DIM*(DIM+1)/2;
        for(unsigned i_comp=0;i_comp<nb_comp;i_comp++){
            S.sigma_skin[i_comp][e.number]=e.sigma_skin[i_comp];
            S.epsilon_skin[i_comp][e.number]=e.epsilon_skin[i_comp];
            S.sigma_mises_skin[e.number]=e.sigma_mises_skin;
        }
   }
   template<class TE> void operator()(TE &e, FieldStructureGroupElementsUser *field_structure_group_elements ) const{
        int nb_comp=DIM*(DIM+1)/2;
        for(unsigned i_comp=0;i_comp<nb_comp;i_comp++){
            field_structure_group_elements->sigma_skin[i_comp][e.number]=e.sigma_skin[i_comp];
            field_structure_group_elements->epsilon_skin[i_comp][e.number]=e.epsilon_skin[i_comp];
            field_structure_group_elements->sigma_mises_skin[e.number]=e.sigma_mises_skin;
            field_structure_group_elements->num_processor_skin[e.number]=e.num_proc_skin;
            field_structure_group_elements->material_behaviour_skin[e.number]=e.typmat_skin;
        }
   }
};

struct Extract_fields_on_element_sst{
   template<class TE, class TS> void operator()(TE &e, TS &S ) const{
        int nb_comp=DIM*(DIM+1)/2;
        for(unsigned i_comp=0;i_comp<nb_comp;i_comp++){
            S.sigma[i_comp][e.number]=e.sigma[0][i_comp];
            S.epsilon[i_comp][e.number]=e.epsilon[0][i_comp];
            S.sigma_mises[e.number]=e.sigma_von_mises;
        }
   }
   template<class TE> void operator()(TE &e, FieldStructureGroupElementsUser *field_structure_group_elements ) const{
        int nb_comp=DIM*(DIM+1)/2;
        for(unsigned i_comp=0;i_comp<nb_comp;i_comp++){
            field_structure_group_elements->sigma[i_comp][e.number]=e.sigma[0][i_comp];
            field_structure_group_elements->epsilon[i_comp][e.number]=e.epsilon[0][i_comp];
            field_structure_group_elements->sigma_mises[e.number]=e.sigma_von_mises;
            field_structure_group_elements->num_processor[e.number]=e.num_proc;
            field_structure_group_elements->material_behaviour[e.number]=e.typmat;
        }
   }
};

template<class TINTER, class TSST>
void convert_fields_to_field_structure_user(TSST &SubS, TINTER &I, Param &process , DataUser &data_user, FieldStructureUser &field_structure_user, GeometryUser &geometry_user) {
    //sauvegarde des coordonnées des noeuds en entier de la sst
    for(unsigned i_sst=0;i_sst<SubS.size();i_sst++){
        //extraction des champs à partir du resultat de calcul (LMTpp)
        calcul_fields_on_sst(SubS[i_sst],process, data_user);

        int id_sst=SubS[i_sst].id;
        //conversion des deplacements des noeuds de la SST

        for(int d=0;d<DIM;d++){
            for(unsigned i=0;i<field_structure_user.find_group_elements(id_sst)->nb_nodes;i++){
                field_structure_user.find_group_elements(id_sst)->displacement_local_nodes[d][i]=SubS[i_sst].mesh->node_list[i].dep[d];
                field_structure_user.find_group_elements(id_sst)->explode_displacement_local_nodes[d][i]=SubS[i_sst].mesh->node_list[i].qtrans[d];
            }
        }

        //conversion des deplacements des noeuds sur la peau
        for(int d=0;d<DIM;d++){
            for(unsigned i=0;i<field_structure_user.find_group_elements(id_sst)->nb_nodes;i++){
                field_structure_user.find_group_elements(id_sst)->displacement_local_nodes_skin[d][i]=SubS[i_sst].mesh->node_list[i].dep[d];
                field_structure_user.find_group_elements(id_sst)->explode_displacement_local_nodes_skin[d][i]=SubS[i_sst].mesh->node_list[i].qtrans[d];
            }
        }
//         for(int d=0;d<DIM;d++){
//             for(unsigned i=0;i<field_structure_user.find_group_elements(id_sst)->nb_nodes_skin;i++){
//                 field_structure_user.find_group_elements(id_sst)->displacement_local_nodes_skin[d][i]=SubS[i_sst].mesh->skin.node_list[i].dep[d];
//             }
//         }

        //conversion des donnees sur les elements
        //extraction sur la SST
        apply(SubS[i_sst].mesh->elem_list,Extract_fields_on_element_sst(), field_structure_user.find_group_elements(id_sst));
        //extraction sur la peau
        apply(SubS[i_sst].mesh->skin.elem_list,Extract_fields_on_element_sst_skin(), field_structure_user.find_group_elements(id_sst));

//         for(unsigned j=0;j<SubS[i_sst].edge.size();++j) {
//             unsigned q=SubS[i_sst].edge[j].internum;
//             unsigned data=SubS[i_sst].edge[j].datanum;
//             int i_inter=geometry_user.find_group_elements(id_sst)->id_adjacent_group_interfaces[j];
//             int i_side=geometry_user.find_group_elements(id_sst)->side_adjacent_group_interfaces[j];
//             std::cout << i_inter << " " << i_side << std::endl;
//             //if(data==0){
//                 int nbnodeseq=field_structure_user.group_interfaces[i_inter].nb_interfaces;
//                 std::cout << field_structure_user.find_group_interfaces(i_inter)->id << " " << nbnodeseq << " " << I[q].side[data].t[1].F.size() << std::endl;
//                 for(int d=0;d<DIM;d++){
//                     for(int ne=0;ne<nbnodeseq;ne++){
//                         field_structure_user.find_group_interfaces(i_inter)->F[i_side][d][ne]=I[q].side[data].t[1].F[ne*DIM+d];
//                         field_structure_user.find_group_interfaces(i_inter)->W[i_side][d][ne]=I[q].side[data].t[1].W[ne*DIM+d];
//                     }
//                 }       
//             //}
//         }
    }
    
}


struct Extract_connectivities_on_element_sst{
   template<class TE, class TS> void operator()(TE &e, TS &S, int nb_nodes_tot) const{
    for(unsigned i=0;i<e.nb_nodes;i++){
        S.mesh_connectivities[i][e.number]=e.node(i)->number_in_original_mesh()+nb_nodes_tot;
    }
   }
   template<class TE, class TS> void operator()(TE &e, TS &S, int nb_nodes_tot, map<int,int> &correspondance_number_in_original_mesh_local_mesh) const{
      for(unsigned i=0;i<e.nb_nodes;i++){
          S.mesh_connectivities[i][e.number]=correspondance_number_in_original_mesh_local_mesh[e.node(i)->number_in_original_mesh()]+nb_nodes_tot;
      }
      S.num_processor[e.number]=e.num_proc_skin;
      S.num_group[e.number]=e.numsst_skin;
      S.material[e.number]=e.typmat_skin;

   }
};

struct Extract_connectivities_on_element_sst_skin{
   template<class TE, class TS> void operator()(TE &e, TS &S, int nb_nodes_tot) const{
      for(unsigned i=0;i<e.nb_nodes;i++){
          S.mesh_connectivities_skin[i][e.number]=e.node(i)->number_in_original_mesh()+nb_nodes_tot;
      }
      S.num_processor[e.number]=e.num_proc_skin;
      S.num_group[e.number]=e.numsst_skin;
      S.material[e.number]=e.typmat_skin;
   }
};

// struct Extract_connectivities_on_element_inter{
//    template<class TE, class TI> void operator()(TE &e, TI &I, int nb_nodes_tot) const{
//       for(unsigned i=0;i<e.nb_nodes;i++){
//           I.mesh_connectivities[i][e.number]=e.node(i)->number_in_original_mesh()+nb_nodes_tot;
//       }
//       I.number[e.number]=e.num;
//       I.nature[e.number]=e.type;
//    }
// };
// 



// template<class TSST>
// void create_hdf_geometry_data_SST(TSST &S, Param &process, int nb_previous_nodes ) {
//     //sauvegarde des coordonnées des noeuds
//     int nb_nodes=S.mesh->skin.node_list.size();
//     for(int d=0;d<DIM;d++){
//         S.nodes[d].resize(nb_nodes);
//         for(unsigned i=0;i<nb_nodes;i++){
//             S.nodes[d][i]=S.mesh->skin.node_list[i].pos[d];
//         }
//     }
//     
//     //creation map pour connaitre les noeuds locaux
//     map<int,int> correspondance_number_in_original_mesh_local_mesh;
//     for(unsigned i=0;i<nb_nodes;i++)
//           correspondance_number_in_original_mesh_local_mesh[ S.mesh->skin.node_list[i].number_in_original_mesh()]=i;
//     
//     //sauvegarde des connectivites du maillage de peau
//     //definition des tailles
//     S.mesh_connectivities.resize(S.nb_nodes_by_element);
//     for(int ne=0;ne<S.nb_nodes_by_element;ne++)
//         S.mesh_connectivities[ne].resize(S.mesh->skin.elem_list.size());
//     
//     apply(S.mesh->skin.elem_list,Projection_data_elements_0_on_skin_sst(),S.mesh->skin, *S.mesh.m);
//     S.material.resize(S.mesh->skin.elem_list.size());
//     S.num_processor.resize(S.mesh->skin.elem_list.size());
//     S.num_group.resize(S.mesh->skin.elem_list.size());
//     //extraction
//     apply(S.mesh->skin.elem_list,Extract_connectivities_on_element_sst(), S, nb_previous_nodes, correspondance_number_in_original_mesh_local_mesh);
//     
// 
// }

struct Extract_connectivities_on_element_sst_inter{
   template<class TE, class TI, class TV> void operator()(TE &e, TI &I, int nb_nodes_tot, TV &correspondance_ddl_edge_sst) const{
      for(unsigned i=0;i<e.nb_nodes;i++){
          I.mesh_connectivities[i][e.number]=correspondance_ddl_edge_sst[e.node(i)->number_in_original_mesh()*TI::dim]/TI::dim+nb_nodes_tot;
      }
      I.number[e.number]=I.num;
      int type=0;
      if (I.type=="Ext" and I.comp=="depl"){type=0;}
      else if (I.type=="Ext" and I.comp=="effort"){type=1;}
      else if (I.type=="Ext" and ( I.comp=="sym" )){type=2;}
      else if (I.type=="Ext" and ( I.comp=="depl_normal")){type=3;}
      else if (I.type=="Int" and I.comp=="Parfait"){type=4;}
      else if (I.type=="Int" and (I.comp=="Contact" or I.comp=="Contact_jeu" or I.comp=="Contact_jeu_physique" or I.comp=="Contact_ep") ){type=5;}
      else if (I.type=="Int" and I.comp=="Jeu_impose"){type=6;}
      else if (I.type=="Ext" and I.comp=="periodique"){type=7;}
      else {type=8;}
      I.nature[e.number]=type;
   }
};

struct find_nb_elements_with_type{
   template<class TE> void operator()(TE &e, Vec<int,4> &nb_elements_with_type) const{
    if(e.nb_nodes==3){
        nb_elements_with_type[0]+=1;
    }
    else if(e.nb_nodes==4 and e.dim==2){
        nb_elements_with_type[1]+=1;
    }
    else if(e.nb_nodes==4 and e.dim==3){
        nb_elements_with_type[2]+=1;
    }
    else if(e.nb_nodes==8){
        nb_elements_with_type[3]+=1;
    }
    else {std::cerr << "Type d'element non implemente " << std::endl;}
   }
};

template<class TSST>
void test_nb_elements_with_type(TSST &S) {
    Vec<int> ind=find_with_index(S.nb_elements_with_type!=0);
    if(ind.size()!=1){std::cerr << "La procedure de sauvegarde des maillages ne fonctionne qu'avec un type d'element par SST"<< std::endl; assert(0);}
    if(S.nb_elements_with_type[0]!=0){
        S.nb_nodes_by_element_sst=3;
        S.type_elements_sst="Triangle";
        S.nb_nodes_by_element_sst_skin=2;
        S.type_elements_sst_skin="Polygon";
        S.pattern_id=0;
    }
    else if(S.nb_elements_with_type[1]!=0){
        S.nb_nodes_by_element_sst=4;
        S.type_elements_sst="Quadrilateral";
        S.nb_nodes_by_element_sst_skin=2;
        S.type_elements_sst_skin="Polygon";
        S.pattern_id=5;
    }
    else if(S.nb_elements_with_type[2]!=0){
        S.nb_nodes_by_element_sst=4;
        S.type_elements_sst="Tetrahedron";
        S.nb_nodes_by_element_sst_skin=3;
        S.type_elements_sst_skin="Triangle";
        S.pattern_id=1;
    }
    else if(S.nb_elements_with_type[3]!=0){
        S.nb_nodes_by_element_sst=8;
        S.type_elements_sst="Hexahedron";
        S.nb_nodes_by_element_sst_skin=4;
        S.type_elements_sst_skin="Quadrilateral";
        S.pattern_id=6;
    }
}

#include "correspondance_ddl_sst.h"



template<class TSST,class TV2>
void create_hdf_geometry_data_SST_INTER(TSST &S, TV2 &Inter, Param &process, int nb_previous_nodes) {

    //sauvegarde des coordonnées des noeuds en entier de la sst
    int nb_nodes=S.mesh->node_list.size();
    for(int d=0;d<DIM;d++){
        S.nodes[d].resize(nb_nodes);
        for(unsigned i=0;i<nb_nodes;i++){
            S.nodes[d][i]=S.mesh->node_list[i].pos[d];
        }
    }
    
    //assignation d'un champ pour determiner le nombre d'elements et le nombre de type d'elements utilises
    S.nb_elements_with_type[0]=0;
    S.nb_elements_with_type[1]=0;
    S.nb_elements_with_type[2]=0;
    S.nb_elements_with_type[3]=0;
    apply(S.mesh->elem_list,find_nb_elements_with_type(), S.nb_elements_with_type);
    test_nb_elements_with_type(S);
    
    //sauvegarde des connectivites des elements des SST
    S.mesh_connectivities.resize(S.nb_nodes_by_element_sst);
    for(int ne=0;ne<S.nb_nodes_by_element_sst;ne++)
        S.mesh_connectivities[ne].resize(S.mesh->elem_list.size());
    //extraction    
    apply(S.mesh->elem_list,Extract_connectivities_on_element_sst(), S, nb_previous_nodes);
    S.mesh->update_skin();
    //sauvegarde des connectivites des elements de peau des SST
    S.mesh_connectivities_skin.resize(S.nb_nodes_by_element_sst_skin);
    for(int ne=0;ne<S.nb_nodes_by_element_sst_skin;ne++)
        S.mesh_connectivities_skin[ne].resize(S.mesh->skin.elem_list.size());
    
    apply(S.mesh->skin.elem_list,Projection_data_elements_0_on_skin_sst(),S.mesh->skin, *S.mesh.m);
    S.material.resize(S.mesh->skin.elem_list.size());
    S.num_processor.resize(S.mesh->skin.elem_list.size());
    S.num_group.resize(S.mesh->skin.elem_list.size());
    //extraction
    apply(S.mesh->skin.elem_list,Extract_connectivities_on_element_sst_skin(), S, nb_previous_nodes);

    Calc_SST_Correspddl(S,process);

   //sauvegarde des connectivites des elements d'interface
    for(unsigned j=0;j<S.edge.size();++j) {
            unsigned q=S.edge[j].internum;
            unsigned data=S.edge[j].datanum;
            if(data==0){
                //sauvegarde des connectivites en ajoutant un offset du nombre de noeuds des SST precedentes
                Inter[q].mesh_connectivities.resize(S.nb_nodes_by_element_sst_skin);
                for(int ne=0;ne<S.nb_nodes_by_element_sst_skin;ne++)
                    Inter[q].mesh_connectivities[ne].resize(Inter[q].side[data].mesh->elem_list.size());
                //sauvegarde du type d'interface et du numero
                Inter[q].number.resize(Inter[q].side[data].mesh->elem_list.size());
                Inter[q].nature.resize(Inter[q].side[data].mesh->elem_list.size());
                //extraction
                apply(Inter[q].side[data].mesh->elem_list,Extract_connectivities_on_element_sst_inter(),Inter[q], nb_previous_nodes, S.edge[j].repddledge);
            }
    }
    
}

// template<class TINTER>
// void create_hdf_geometry_data_INTER(TINTER &I, Param &process, int nb_previous_nodes ) {
//     //sauvegarde des coordonnées des noeuds
//     int nb_nodes=I.side[0].mesh->node_list.size();
//     for(int d=0;d<DIM;d++){
//         I.nodes[d].resize(nb_nodes);
//         for(unsigned i=0;i<nb_nodes;i++){
//             I.nodes[d][i]=I.side[0].mesh->node_list[i].pos[d];
//         }
//     }
// 
//     //sauvegarde des connectivites du maillage de peau
//     //definition des tailles
//     I.mesh_connectivities.resize(I.nb_nodes_by_element);
//     for(int ne=0;ne<I.nb_nodes_by_element;ne++)
//         I.mesh_connectivities[ne].resize(I.side[0].mesh->elem_list.size());
//     
//     I.nature.resize(I.side[0].mesh->elem_list.size());
//     I.number.resize(I.side[0].mesh->elem_list.size());
//     //extraction
//     apply(I.side[0].mesh->elem_list,Extract_connectivities_on_element_inter(), I, nb_previous_nodes);
//     
// }

// template<class TSST>
// void create_hdf_fields_data(TSST &S, Param &process ) {
//     //sauvegarde des deplacements des noeuds
//     int nb_nodes=S.mesh->skin.node_list.size();
//     for(int d=0;d<DIM;d++){
//         S.dep_nodes[d].resize(nb_nodes);
//         for(unsigned i=0;i<nb_nodes;i++){
//             S.dep_nodes[d][i]=S.mesh->skin.node_list[i].dep[d];
//         }
//     }
//     
//     //sauvegarde des donnees sur les elements
//     //definition des tailles
//     int nb_comp=DIM*(DIM+1)/2;
//     for(unsigned i_comp=0;i_comp<nb_comp;i_comp++){
//         S.sigma[i_comp].resize(S.mesh->skin.elem_list.size());
//         S.epsilon[i_comp].resize(S.mesh->skin.elem_list.size());
//         S.sigma_mises.resize(S.mesh->skin.elem_list.size());
//     }
//     //extraction
//     apply(S.mesh->skin.elem_list,Extract_fields_on_element(), S);
// }


///Ecriture des champs créés lors de la géométrie sur les elements des SST
template<class TSST> void save_data_on_elements_SST(TSST &S, Hdf &hdf_file, String name_list){
    String name_field = name_list + "/num_proc";
    S.num_processor.write_to(hdf_file,name_field);
    name_field = name_list + "/material";
    S.material.write_to(hdf_file,name_field);        
    name_field = name_list + "/num_group";
    S.num_group.write_to(hdf_file,name_field);    
}
///Ecriture des champs créés lors de la géométrie sur les elements des INTERFACES
template<class TSST> void save_data_on_elements_INTER(TSST &S, Hdf &hdf_file, String name_list){
    String name_field = name_list + "/number";
    S.number.write_to(hdf_file,name_field);
    name_field = name_list + "/nature";
    S.nature.write_to(hdf_file,name_field);        
}

/*
///sauvegarde de la geometrie utilise pour l'affichage des champs
template<class TSST>
void save_elements_hdf_sst(TSST &S, Param &process, Hdf &hdf_file, String name_geometry, String name_elements) {
    String name_list ;
    name_list<< name_geometry << "/"<< name_elements <<"/list_" << S.num ;
    for (unsigned i_connect=0;i_connect<S.nb_nodes_by_element;i_connect++) {
        String name_connect;
        name_connect << name_list << "/mesh_c"<<i_connect;
        S.mesh_connectivities[i_connect].write_to( hdf_file, name_connect );
    }
    save_data_on_elements_SST(S, hdf_file, name_list);
    
    String type_elements;
    int pattern_id=0;
#if DIM == 2
    type_elements="Bar";
    pattern_id=0;
#else
    type_elements="Triangle";
    pattern_id=1;
#endif
    hdf_file.add_tag(name_list,"base",type_elements.c_str());
}
*/


///sauvegarde de la geometrie utilise pour l'affichage des champs
template<class TSST, class TI>
void save_elements_hdf_sst_inter(TSST &S, TI &I, Param &process, Hdf &hdf_file, String name_geometry) {
    String name_list ;
    String name_elements;
 
    name_elements="elements_0";
    name_list<< name_geometry << "/"<< name_elements <<"/list_" << S.num ;
    for (unsigned i_connect=0;i_connect<S.nb_nodes_by_element_sst;i_connect++) {
        String name_connect;
        name_connect << name_list << "/local_connectivities_"<<i_connect;
        S.mesh_connectivities[i_connect].write_to( hdf_file, name_connect );
    }    
    hdf_file.add_tag(name_list,"base",S.type_elements_sst.c_str());
    hdf_file.write_tag(name_list,"id",S.id);
    hdf_file.write_tag(name_list,"pattern_base_id",S.id_material);
    
    //sauvegarde des elements_0_skin
    name_elements="elements_0_skin";
    name_list=name_geometry; name_list << "/"<< name_elements <<"/list_" << S.num ;
    for (unsigned i_connect=0;i_connect<S.nb_nodes_by_element_sst_skin;i_connect++) {
        String name_connect;
        name_connect << name_list << "/local_connectivities_"<<i_connect;
        S.mesh_connectivities_skin[i_connect].write_to( hdf_file, name_connect );
    }    
    save_data_on_elements_SST(S, hdf_file, name_list);
    hdf_file.add_tag(name_list,"base",S.type_elements_sst_skin.c_str());
    hdf_file.write_tag(name_list,"id",S.id);
    hdf_file.write_tag(name_list,"pattern_base_id",S.id_material);
    
    //sauvegarde des elements_1
    name_elements="elements_1";
    for(unsigned j=0;j<S.edge.size();++j) {
            unsigned q=S.edge[j].internum;
            unsigned data=S.edge[j].datanum;
            name_list=name_geometry; name_list << "/"<< name_elements <<"/list_" << q ;
            if(data==0){
                for (unsigned i_connect=0;i_connect<S.nb_nodes_by_element_sst_skin;i_connect++) {
                    String name_connect;
                    name_connect << name_list << "/local_connectivities_"<<i_connect;
                    I[q].mesh_connectivities[i_connect].write_to( hdf_file, name_connect );
                }
                save_data_on_elements_INTER(I[q], hdf_file, name_list);
                hdf_file.add_tag(name_list,"base",S.type_elements_sst_skin.c_str());
                hdf_file.add_tag(name_list,"type",I[q].type.c_str());
                hdf_file.write_tag(name_list,"id",I[q].id);
                if(I[q].type=="Ext")
                    hdf_file.write_tag(name_list,"pattern_base_id",I[q].id_bc);
                else hdf_file.write_tag(name_list,"pattern_base_id",I[q].id_link);
            }
    }
}

///sauvegarde de la geometrie utilise pour l'affichage des champs
template<class TI>
void save_elements_hdf_inter(TI &I, Param &process, Hdf &hdf_file, String name_geometry, String name_elements) {
    String name_list ;
    name_list<< name_geometry << "/"<< name_elements <<"/list_" << I.num ;
    for (unsigned i_connect=0;i_connect<I.nb_nodes_by_element;i_connect++) {
        String name_connect;
        name_connect << name_list << "/mesh_c"<<i_connect;
        I.mesh_connectivities[i_connect].write_to( hdf_file, name_connect );
    }
    save_data_on_elements_INTER(I, hdf_file, name_list);
    
    String type_elements;
    int pattern_id=0;
#if DIM == 2
    type_elements="Bar";
    pattern_id=0;
#else
    type_elements="Triangle";
    pattern_id=1;
#endif
    hdf_file.add_tag(name_list,"base",type_elements.c_str());
}



template<class TSST, class TINTER>
void create_hdf_fields_data_SST(TSST &S, TINTER &Inter, Param &process ) {
    //sauvegarde des deplacements des noeuds dans la SST
    int nb_nodes=S.mesh->node_list.size();
    for(int d=0;d<DIM;d++){
        S.dep_nodes[d].resize(nb_nodes);
        for(unsigned i=0;i<nb_nodes;i++){
            S.dep_nodes[d][i]=S.mesh->node_list[i].dep[d];
        }
    }
    //sauvegarde des deplacements des noeuds sur la peau
    nb_nodes=S.mesh->skin.node_list.size();
    for(int d=0;d<DIM;d++){
        S.dep_nodes_skin[d].resize(nb_nodes);
        for(unsigned i=0;i<nb_nodes;i++){
            S.dep_nodes_skin[d][i]=S.mesh->skin.node_list[i].dep[d];
        }
    }
    
    //sauvegarde des donnees sur les elements
    //definition des tailles
    int nb_comp=DIM*(DIM+1)/2;
    for(unsigned i_comp=0;i_comp<nb_comp;i_comp++){
        S.sigma[i_comp].resize(S.mesh->elem_list.size());
        S.epsilon[i_comp].resize(S.mesh->elem_list.size());
        S.sigma_mises.resize(S.mesh->elem_list.size());
        
        S.sigma_skin[i_comp].resize(S.mesh->skin.elem_list.size());
        S.epsilon_skin[i_comp].resize(S.mesh->skin.elem_list.size());
        S.sigma_mises_skin.resize(S.mesh->skin.elem_list.size());
    }
    //extraction sur la SST
    apply(S.mesh->elem_list,Extract_fields_on_element_sst(), S);
    //extraction sur la peau
    apply(S.mesh->skin.elem_list,Extract_fields_on_element_sst_skin(), S);
    
    for(unsigned j=0;j<S.edge.size();++j) {
        unsigned q=S.edge[j].internum;
        unsigned data=S.edge[j].datanum;
        if(data==0){
 /*       Vec<unsigned> &list1=(Inter[q].side[data].ddlcorresp);
    Vec<unsigned> &list2=(Inter.side[1].ddlcorresp);
    Vec<T> Wchap1=Inter.side[0].t[imic].Wpchap[list1];
    Vec<T> Wchap2=Inter.side[1].t[imic].Wpchap[list2];
    Vec<T> Fchap1=Inter.side[0].t[imic].Fchap[list1];
    Vec<T> Fchap2=Inter.side[1].t[imic].Fchap[list2];
    const Vec<T> &Q1=Inter.side[0].t[imic].F[list1];
    const Vec<T> &Q2=Inter.side[1].t[imic].F[list2];
    const Vec<T> &WW1=Inter.side[0].t[imic].Wp[list1];
    const Vec<T> &WW2=Inter.side[1].t[imic].Wp[list2];*/
            int nbnodeseq=Inter[q].side[data].t[1].F.size()/DIM;
            for(int d=0;d<DIM;d++){
                Inter[q].F[d].resize(nbnodeseq);
                Inter[q].W[d].resize(nbnodeseq);
                Inter[q].Fchap[d].resize(nbnodeseq);
                Inter[q].Wchap[d].resize(nbnodeseq);
                for(int ne=0;ne<nbnodeseq;ne++){
                    Inter[q].F[d][ne]=Inter[q].side[data].t[1].F[ne*DIM+d];
                    Inter[q].W[d][ne]=Inter[q].side[data].t[1].W[ne*DIM+d];
                    Inter[q].Fchap[d][ne]=Inter[q].side[data].t[1].Fchap[ne*DIM+d];
                    Inter[q].Wchap[d][ne]=Inter[q].side[data].t[1].Wchap[ne*DIM+d];
                }
            }       
        }
    }
    
    
}


///sauvegarde de la geometrie utilise pour l'affichage des champs
template<class TSST>
void save_fields_hdf_SST(TSST &S, Param &process, Hdf &hdf , String name_group_fields) {
    String name_fields; name_fields << name_group_fields <<"/pt_"<< process.temps->pt_cur;
    String name_sigma, name_epsilon ;
    name_sigma<< name_fields << "/sigma/list_" << S.num ;
    name_epsilon<< name_fields << "/epsilon/list_" << S.num ;
#if DIM==2
    BasicVec<String> tensor_comp= BasicVec<String>("/xx","/yy","/xy");
#else
    BasicVec<String> tensor_comp= BasicVec<String>("/xx","/yy","/zz","/xy","/xz","/yz");        
#endif
    for(unsigned i_comp=0;i_comp<tensor_comp.size();i_comp++){
        String name_sigma_field, name_epsilon_field;
        name_sigma_field=name_sigma+tensor_comp[i_comp];
        name_epsilon_field=name_epsilon+tensor_comp[i_comp];
        S.sigma[i_comp].write_to(hdf,name_sigma_field.c_str());
        S.epsilon[i_comp].write_to(hdf,name_epsilon_field.c_str());
    }
    String name_sigma_mises;
    name_sigma_mises <<name_fields<<"/sigma_von_mises/list_" << S.num ;
    S.sigma_mises.write_to(hdf,name_sigma_mises.c_str());
}

///sauvegarde de la geometrie utilise pour l'affichage des champs
template<class TSST, class TINTER>
void save_fields_hdf_SST_INTER(TSST &S, TINTER &I, Param &process, Hdf &hdf , String name_group_fields) {
    String name_fields; name_fields << name_group_fields <<"/pt_"<< process.temps->pt_cur;
#if DIM==2
    BasicVec<String> tensor_comp= BasicVec<String>("/xx","/yy","/xy");
#else
    BasicVec<String> tensor_comp= BasicVec<String>("/xx","/yy","/zz","/xy","/xz","/yz");        
#endif
    //sauvegarde des champs des SST
    String name_sigma, name_epsilon ;
    name_sigma<< name_fields << "/elements_0/sigma/list_" << S.num ;
    name_epsilon<< name_fields << "/elements_0/epsilon/list_" << S.num ;
    for(unsigned i_comp=0;i_comp<tensor_comp.size();i_comp++){
        String name_sigma_field, name_epsilon_field;
        name_sigma_field=name_sigma+tensor_comp[i_comp];
        name_epsilon_field=name_epsilon+tensor_comp[i_comp];
        S.sigma[i_comp].write_to(hdf,name_sigma_field.c_str());
        S.epsilon[i_comp].write_to(hdf,name_epsilon_field.c_str());
    }
    String name_sigma_mises;
    name_sigma_mises <<name_fields<<"/elements_0/sigma_von_mises/list_" << S.num ;
    S.sigma_mises.write_to(hdf,name_sigma_mises.c_str());
     //sauvegarde des champs de peau
    String name_sigma_skin, name_epsilon_skin ;
    name_sigma_skin<< name_fields << "/elements_0_skin/sigma_skin/list_" << S.num ;
    name_epsilon_skin<< name_fields << "/elements_0_skin/epsilon_skin/list_" << S.num ;
    for(unsigned i_comp=0;i_comp<tensor_comp.size();i_comp++){
        String name_sigma_field, name_epsilon_field;
        name_sigma_field=name_sigma_skin+tensor_comp[i_comp];
        name_epsilon_field=name_epsilon_skin+tensor_comp[i_comp];
        S.sigma_skin[i_comp].write_to(hdf,name_sigma_field.c_str());
        S.epsilon_skin[i_comp].write_to(hdf,name_epsilon_field.c_str());
    }
    String name_sigma_mises_skin;
    name_sigma_mises_skin <<name_fields<<"/elements_0_skin/sigma_von_mises_skin/list_" << S.num ;
    S.sigma_mises_skin.write_to(hdf,name_sigma_mises_skin.c_str());
      
    //sauvegarde des champs sur les interfaces
    BasicVec<String> name_direction("x","y","z");
    BasicVec<String> list_name_field("F","W","Fchap","Wchap");
    for(unsigned j=0;j<S.edge.size();++j) {
            unsigned q=S.edge[j].internum;
            unsigned data=S.edge[j].datanum;
            String name_F; name_F<< name_fields << "/elements_1/F/list_" << q ;
            String name_W; name_W<< name_fields << "/elements_1/W/list_" << q ;
            String name_Fchap; name_Fchap<< name_fields << "/elements_1/Fchap/list_" << q ;
            String name_Wchap; name_Wchap<< name_fields << "/elements_1/Wchap/list_" << q ;
            if(data==0){
                for (unsigned d=0;d<DIM;d++) {
                    String name_Fdim=name_F+"/"+name_direction[d]; 
                    I[q].F[d].write_to(hdf,name_Fdim);
                    String name_Wdim=name_W+"/"+name_direction[d]; 
                    I[q].W[d].write_to(hdf,name_Wdim);
                    String name_Fchapdim=name_Fchap+"/"+name_direction[d]; 
                    I[q].Fchap[d].write_to(hdf,name_Fchapdim);
                    String name_Wchapdim=name_Wchap+"/"+name_direction[d]; 
                    I[q].Wchap[d].write_to(hdf,name_Wchapdim);
                }
            }
    }
   
}

#include "utils_2.h"

template<class TSST>
void write_hdf_geometry_SST(TSST &SubS, Param &process ) {
    
    //chaque processeur calcul stocke les noeuds de ces sst
    BasicVec<int> nb_previous_nodes;
    nb_previous_nodes.push_back(0);
    for(unsigned i=0;i<SubS.size();i++) {
        SubS[i].mesh->update_skin();
        nb_previous_nodes.push_back(SubS[i].mesh->skin.node_list.size()+nb_previous_nodes[i]);
    }
    //creation des donnees hdf et sauvegarde d'un fichier pour chaque processeur
    for(unsigned i=0;i<SubS.size();i++) create_hdf_geometry_data_SST(SubS[i],process,nb_previous_nodes[i]);
    String name_hdf ; name_hdf << process.affichage->name_hdf <<"_"<< process.rank<<".h5";
    if(FileExists(name_hdf.c_str())){ String command = "rm -rf "+name_hdf; int syst_rm=system(command.c_str());}
    Hdf hdf_file( name_hdf.c_str() );
    //ecriture des noeuds (concatenation en attendant la possibilite d'ecrire à la suite en hdf)
    BasicVec<BasicVec<TYPE>,DIM> nodes;
    BasicVec<String> name_direction("x","y","z");
    for(unsigned d=0;d<DIM;d++) {
        nodes[d].resize(nb_previous_nodes[SubS.size()]);
        for(unsigned i=0;i<SubS.size();i++) 
            for(unsigned j=0;j<SubS[i].nodes[d].size();j++)
                nodes[d][j+nb_previous_nodes[i]]=SubS[i].nodes[d][j];
        String name_dim;  name_dim << process.affichage->name_geometry << "/local_nodes_0/" << name_direction[d];
        nodes[d].write_to(hdf_file,name_dim);
    }
    //ecriture des elements et des champs aux elements créés lors de la géométrie
    for(unsigned i=0;i<SubS.size();i++) {
        SubS[i].mesh->update_skin();
        save_elements_hdf_sst(SubS[i], process, hdf_file, process.affichage->name_geometry, "elements_0_skin");
    }
}

template<class TINTER>
void write_hdf_geometry_INTER(TINTER &SubI, Param &process ) {
    
    //chaque processeur calcul stocke les noeuds de ces interfaces
    BasicVec<int> nb_previous_nodes;
    nb_previous_nodes.push_back(0);
    for(unsigned i=0;i<SubI.size();i++) {
        nb_previous_nodes.push_back(SubI[i].side[0].mesh->node_list.size()+nb_previous_nodes[i]);
    }
    //creation des donnees hdf et sauvegarde d'un fichier pour chaque processeur
    for(unsigned i=0;i<SubI.size();i++) create_hdf_geometry_data_INTER(SubI[i],process,nb_previous_nodes[i]);
    String name_hdf ; name_hdf << process.affichage->name_hdf <<"_"<< process.rank<<".h5";
    Hdf hdf_file( name_hdf.c_str() );
    //ecriture des noeuds (concatenation en attendant la possibilite d'ecrire à la suite en hdf)
    BasicVec<BasicVec<TYPE>,DIM> nodes;
    BasicVec<String> name_direction("x","y","z");
    for(unsigned d=0;d<DIM;d++) {
        nodes[d].resize(nb_previous_nodes[SubI.size()]);
        for(unsigned i=0;i<SubI.size();i++) 
            for(unsigned j=0;j<SubI[i].nodes[d].size();j++)
                nodes[d][j+nb_previous_nodes[i]]=SubI[i].nodes[d][j];
        String name_dim;  name_dim << process.affichage->name_geometry << "/local_nodes_1/" << name_direction[d];
        nodes[d].write_to(hdf_file,name_dim);
    }
    //ecriture des elements et des champs aux elements créés lors de la géométrie
    for(unsigned i=0;i<SubI.size();i++) {
        save_elements_hdf_inter(SubI[i], process, hdf_file, process.affichage->name_geometry, "elements_1");
    }
}



template<class TINTER, class TSST>
void write_hdf_geometry_SST_INTER(TSST &SubS, TINTER &I, Param &process , GeometryUser &geometry_user) {
    
    
    //chaque processeur calcul stocke les noeuds de ces sst
    BasicVec<int> nb_previous_nodes;
    nb_previous_nodes.push_back(0);
    for(unsigned i=0;i<SubS.size();i++) {
        nb_previous_nodes.push_back(SubS[i].mesh->node_list.size()+nb_previous_nodes[i]);
    }
    //creation des donnees hdf et sauvegarde d'un fichier pour chaque processeur
    for(unsigned i=0;i<SubS.size();i++) {
        create_hdf_geometry_data_SST_INTER(SubS[i],I, process,nb_previous_nodes[i]);
    }
    String name_hdf ; name_hdf << process.affichage->name_hdf <<"_"<< process.rank<<".h5";
    if(FileExists(name_hdf.c_str())){ String command = "rm -rf "+name_hdf; int syst_rm=system(command.c_str());}
    Hdf hdf_file( name_hdf.c_str() );
    //ecriture des noeuds (concatenation en attendant la possibilite d'ecrire à la suite en hdf)
    BasicVec<BasicVec<TYPE>,DIM> nodes;
    BasicVec<String> name_direction("x","y","z");
    for(unsigned d=0;d<DIM;d++) {
        nodes[d].resize(nb_previous_nodes[SubS.size()]);
        for(unsigned i=0;i<SubS.size();i++) 
            for(unsigned j=0;j<SubS[i].nodes[d].size();j++)
                nodes[d][j+nb_previous_nodes[i]]=SubS[i].nodes[d][j];
        String name_dim;  name_dim << process.affichage->name_geometry << "/nodes/" << name_direction[d];
        nodes[d].write_to(hdf_file,name_dim);
    }
    //ecriture des elements et des champs aux elements créés lors de la géométrie
    
    for(unsigned i=0;i<SubS.size();i++) {
        save_elements_hdf_sst_inter(SubS[i], I, process, hdf_file, process.affichage->name_geometry);
    }
    
    


}

template<class TSST>
void write_hdf_fields_SST(TSST &SubS, Param &process ) {
    //chaque processeur calcul stocke les noeuds de ces sst
    BasicVec<int> nb_previous_nodes;
    nb_previous_nodes.push_back(0);
    for(unsigned i=0;i<SubS.size();i++) {
        SubS[i].mesh->update_skin();
        nb_previous_nodes.push_back(SubS[i].mesh->skin.node_list.size()+nb_previous_nodes[i]);
    }
    //ouverture d'un fichier pour chaque processeur
    String name_hdf ; name_hdf << process.affichage->name_hdf <<"_"<< process.rank<<".h5";
    Hdf hdf_file( name_hdf.c_str() );
    //calcul des champs sur le maillage a partir de la solution et écriture des champs hdf
    for(unsigned i=0;i<SubS.size();i++){
        calcul_fields_on_sst(SubS[i],process);
        create_hdf_fields_data(SubS[i],process);
    }
    
    //concatenation des noeuds et ecriture dans le hdf
    BasicVec<BasicVec<TYPE>, DIM> dep_nodes;
    String name_displacements; name_displacements<< process.affichage->name_fields <<"/pt_"<<process.temps->pt_cur <<"/displacements";
    BasicVec<String> displacements_coor= BasicVec<String>("/x","/y","/z");
    for(unsigned d=0;d<DIM;d++) {
        dep_nodes[d].resize(nb_previous_nodes[SubS.size()]);
        for(unsigned i=0;i<SubS.size();i++) 
            for(unsigned j=0;j<SubS[i].nodes[d].size();j++)
                dep_nodes[d][j+nb_previous_nodes[i]]=SubS[i].dep_nodes[d][j];
        String name_displacement_coor; name_displacement_coor=name_displacements+displacements_coor[d];
        dep_nodes[d].write_to(hdf_file,name_displacement_coor.c_str());
    }

    //ecriture des champs par elements dans le hdf
    for(unsigned i=0;i<SubS.size();i++) {
        save_fields_hdf_SST(SubS[i],process, hdf_file , process.affichage->name_fields);
    }
    String name_fields ;
    name_fields<< process.affichage->name_fields <<"/pt_"<< process.temps->pt_cur ;
    int i_step=process.temps->step_cur;
    int i_pt=process.temps->time_step[i_step].pt_cur;
    TYPE val_time=process.temps->time_step[i_step].t_ini+(i_pt+1)*process.temps->time_step[i_step].dt ;
    hdf_file.write_tag(name_fields,"time",val_time);

}

template<class TSST, class TINTER>
void write_hdf_fields_SST_INTER(TSST &SubS, TINTER &Inter,Param &process , DataUser &data_user) {
    //chaque processeur calcul stocke les noeuds de ces sst
    BasicVec<int> nb_previous_nodes;
    nb_previous_nodes.push_back(0);
    for(unsigned i=0;i<SubS.size();i++) {
        //SubS[i].mesh->update_skin();
        nb_previous_nodes.push_back(SubS[i].mesh->node_list.size()+nb_previous_nodes[i]);
    }
    //ouverture d'un fichier pour chaque processeur
    String name_hdf ; name_hdf << process.affichage->name_hdf <<"_"<< process.rank<<".h5";
    Hdf hdf_file( name_hdf.c_str() );
    //calcul des champs sur le maillage a partir de la solution et écriture des champs hdf
    for(unsigned i=0;i<SubS.size();i++){
        calcul_fields_on_sst(SubS[i],process, data_user);
        create_hdf_fields_data_SST(SubS[i],Inter, process);
    }
    
    //concatenation des noeuds et ecriture dans le hdf
    BasicVec<BasicVec<TYPE>, DIM> dep_nodes;
    String name_displacements; name_displacements<< process.affichage->name_fields <<"/pt_"<<process.temps->pt_cur <<"/elements_0/displacements";
    BasicVec<String> displacements_coor= BasicVec<String>("/x","/y","/z");
    for(unsigned d=0;d<DIM;d++) {
        dep_nodes[d].resize(nb_previous_nodes[SubS.size()]);
        for(unsigned i=0;i<SubS.size();i++) 
            for(unsigned j=0;j<SubS[i].mesh->node_list.size();j++)
                dep_nodes[d][j+nb_previous_nodes[i]]=SubS[i].dep_nodes[d][j];
        String name_displacement_coor; name_displacement_coor=name_displacements+displacements_coor[d];
        dep_nodes[d].write_to(hdf_file,name_displacement_coor.c_str());
    }

    BasicVec<int> nb_previous_nodes_skin;
    nb_previous_nodes_skin.push_back(0);
    for(unsigned i=0;i<SubS.size();i++) {
        //SubS[i].mesh->update_skin();
        nb_previous_nodes_skin.push_back(SubS[i].mesh->skin.node_list.size()+nb_previous_nodes_skin[i]);
    }
    //PRINT(nb_previous_nodes_skin);
    BasicVec<BasicVec<TYPE>, DIM> dep_nodes_skin;
    String name_displacements_skin; name_displacements_skin<< process.affichage->name_fields <<"/pt_"<<process.temps->pt_cur <<"/elements_0_skin/displacements_skin";
    for(unsigned d=0;d<DIM;d++) {
        dep_nodes_skin[d].resize(nb_previous_nodes_skin[SubS.size()]);
        for(unsigned i=0;i<SubS.size();i++) 
            for(unsigned j=0;j<SubS[i].mesh->skin.node_list.size();j++)
                dep_nodes_skin[d][j+nb_previous_nodes_skin[i]]=SubS[i].dep_nodes_skin[d][j];
        String name_displacement_coor; name_displacement_coor=name_displacements_skin+displacements_coor[d];
        dep_nodes_skin[d].write_to(hdf_file,name_displacement_coor.c_str());
    }
    
    hdf_file.add_tag(name_displacements,"type","Node");
    hdf_file.add_tag(name_displacements_skin,"type","Node");
    hdf_file.write_tag(name_displacements,"nb_comp",DIM);
    hdf_file.write_tag(name_displacements_skin,"nb_comp",DIM);

    //ecriture des champs par elements dans le hdf
    for(unsigned i=0;i<SubS.size();i++) {
        save_fields_hdf_SST_INTER(SubS[i] , Inter, process, hdf_file , process.affichage->name_fields);
    }
    
    
    String name_fields ;
    name_fields<< process.affichage->name_fields <<"/pt_"<< process.temps->pt_cur ;
    int i_step=process.temps->step_cur;
    int i_pt=process.temps->time_step[i_step].pt_cur;
    TYPE val_time=process.temps->time_step[i_step].t_ini+(i_pt+1)*process.temps->time_step[i_step].dt ;
    hdf_file.write_tag(name_fields,"time",val_time);
    
    int nb_comp_tensor=DIM*(DIM+1)/2, nb_comp_vector=DIM;
    BasicVec< BasicVec<String>, 3 > list_name_field_elements;
    BasicVec<BasicVec<int>, 3 > list_nb_comp;
    list_name_field_elements[0]=BasicVec<String>("sigma","epsilon","sigma_von_mises");
    list_nb_comp[0]=BasicVec<int,3>(nb_comp_tensor,nb_comp_tensor,1);
    list_name_field_elements[1]=BasicVec<String>("sigma_skin","epsilon_skin","sigma_von_mises_skin");
    list_nb_comp[1]=BasicVec<int,3>(nb_comp_tensor,nb_comp_tensor,1);
    list_name_field_elements[2]=BasicVec<String>("F","W","Fchap","Wchap");
    list_nb_comp[2]=BasicVec<int,3>(nb_comp_vector,nb_comp_vector,nb_comp_vector);
    BasicVec<String,3> group_type("elements_0","elements_0_skin","elements_1");
    for(int i=0;i<3;i++){
        for(int j=0;j<list_name_field_elements[i].size();j++){
            String name_field=name_fields+"/"+group_type[i]+"/"+list_name_field_elements[i][j];
            //PRINT(name_field);
            hdf_file.add_tag(name_field,"type","Elem");
            hdf_file.write_tag(name_field,"nb_comp",list_nb_comp[i][j]);
            
        }
    }
 
}
