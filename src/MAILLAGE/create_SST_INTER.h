#include<fstream>
#include <sstream>
#include "mesh/read_avs.h"
#include "mesh/read_geof.h"
#include "mesh/write_avs.h"
#include "mesh/read_mshadv.h"
#include "containers/algo.h"
#include <map>
#include "find_entity.h"

#include "GeometryUser.h"
#include "DataUser.h"
using namespace Metil;
// #include "definition_PARAM_AFFICHAGE.h"
// #include "affichage_mesh_SST.h"

//#include "mesh/ordering.h"
//#include "containers/evaluate_nb_cycles.h"

#include "mpi.h"

extern "C" {
    // #include "metis.h"
//     void METIS_PartGraphRecursive(int *, long long int *, long long int *, long long int *, long long int *, int *, int *, int *, int *, int *, long long int *);
    void METIS_PartGraphRecursive(int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *);
}

using namespace LMT;
using namespace std;

/** \defgroup Maillage_geometrie_sst Géométrie et maillage des Sous-structures
\ingroup Maillage_geometrie
*/
/** \defgroup Maillage_geometrie_inter Géométrie et maillage des Interfaces
\ingroup Maillage_geometrie
*/


/** \ingroup Maillage_geometrie_sst
\brief Création du nombre de Sst et affectation du numero du materiau par l'intermédiaire du data_user (fichier json)
 
La structure data_user contient une liste de group_elements correspondant à chaque sst et une id_material. Cet id_material correspond à l'identificateur des matériaux donnés dans le fichier json.
A cette etape, on spécifie la taille du vecteur de Sous-structures par l'intermédiaire geometry_user.nb_group_elements et on assigne le numéro du matériau pour chaque sous-structure.
*/
template<class TV1>
void create_SST_typmat(DataUser &data_user, GeometryUser &geometry_user,TV1 &S,Param &process) {
    
    //initialisation de la taille des sst et de leur id == num   
    //S.resize(geometry_user.nb_group_elements);
    
    for(unsigned i=0;i<S.size();i++) {
        S[i].num = geometry_user.group_elements[i].id;
        S[i].id = geometry_user.group_elements[i].id;
    }
    
    //assignation des numero de materiaux aux sst
    for(int i_data_group=0; i_data_group<data_user.group_elements.size(); i_data_group++){
        int id_mat = data_user.find_behaviour_materials(data_user.group_elements[i_data_group].id_material)->id;
        int index_mat = data_user.find_index_behaviour_materials(data_user.group_elements[i_data_group].id_material);
        find_sst(S,data_user.group_elements[i_data_group].id)->id_material = id_mat ;
        find_sst(S,data_user.group_elements[i_data_group].id)->typmat = index_mat ; 
    }
}


/**\ingroup Maillage_geometrie_sst
\brief Lecture et assignation d'un maillage par sst
 
Le maillage de chaque sst est lu à partir des données chargée dans geometry_user. Ces données viennent de SC_create_2
 
*/
///Fonction generique
template<class TE, class TM, class TR>
void add_new_elem(TE &e, TM &m, TR &rep_nodes) {}
///Fonctions specialisees
template<class TNB,class TN,class TD,unsigned NET, class TM, class TR>
void add_new_elem(Element<Tetra,TNB,TN,TD,NET> &e, TM &m, TR &rep_nodes) {
    m.add_element(Tetra(),DefaultBehavior(),rep_nodes.ptr() );
}
template<class TNB,class TN,class TD,unsigned NET, class TM, class TR>
void add_new_elem(Element<Wedge,TNB,TN,TD,NET> &e, TM &m, TR &rep_nodes) {
    m.add_element(Wedge(),DefaultBehavior(),rep_nodes.ptr() );
}
template<class TNB,class TN,class TD,unsigned NET, class TM, class TR>
void add_new_elem(Element<Hexa,TNB,TN,TD,NET> &e, TM &m, TR &rep_nodes) {
    m.add_element(Hexa(),DefaultBehavior(),rep_nodes.ptr() );
}
template<class TNB,class TN,class TD,unsigned NET, class TM, class TR>
void add_new_elem(Element<Triangle,TNB,TN,TD,NET> &e, TM &m, TR&rep_nodes) {
    m.add_element(Triangle(),DefaultBehavior(),rep_nodes.ptr() );
}
template<class TNB,class TN,class TD,unsigned NET, class TM, class TR>
void add_new_elem(Element<Quad,TNB,TN,TD,NET> &e, TM &m, TR &rep_nodes) {
    m.add_element(Quad(),DefaultBehavior(),rep_nodes.ptr() );
}
template<class TNB,class TN,class TD,unsigned NET, class TM, class TR>
void add_new_elem(Element<Bar,TNB,TN,TD,NET> &e, TM &m, TR &rep_nodes) {
    m.add_element(Bar(),DefaultBehavior(),rep_nodes.ptr() );
}

//Structure utilisee pour realiser une operation sur les elements du maillage de peau d'une sous-structure : 
//      copie des connectivites des elements de peau du maillage local d'une sous-structure (LMT::mesh) vers l'objet geometry_user (GPU)
struct ConvertMeshConnectivitiesSkin{
   template<class TE> void operator()(TE &e, GeometryUser &geometry_user, int id_sst) const{
      for(unsigned i=0;i<e.nb_nodes;i++){
          geometry_user.find_group_elements(id_sst)->local_connectivities_skin[i][e.number]=e.node(i)->number_in_original_mesh();
      }
   }
};

//conversion du maillage de peau d'une sous-structure (LMT::mesh) aux group_elements
template<class TS>
void convert_mesh_skin_to_geometry_user(TS &S, GeometryUser &geometry_user){
    int id_sst=S.id;
    int nb_nodes_by_element_skin=(geometry_user.patterns.find_type(geometry_user.find_group_elements(id_sst)->pattern_base_id)).nb_nodes_by_sides;
    geometry_user.find_group_elements(id_sst)->local_connectivities_skin.resize(nb_nodes_by_element_skin);
    for(int ne=0;ne<nb_nodes_by_element_skin;ne++){
        geometry_user.find_group_elements(id_sst)->local_connectivities_skin[ne].resize(S.mesh->skin.elem_list.size());
    }
    apply(S.mesh->skin.elem_list,ConvertMeshConnectivitiesSkin(),geometry_user, id_sst);
   
}


template<class TV1>
void create_maillage_SST(DataUser &data_user, GeometryUser &geometry_user, TV1 &S, Param &process) {
    for(unsigned i=0;i<S.size();++i) {
        std::string namein = data_user.find_group_elements(S[i].num)->name;
        S[i].mesh.name=namein;
        S[i].mesh.load(geometry_user, S[i].id);
        S[i].mesh.load();
        S[i].mesh->update_skin();
        geometry_user.find_group_elements(S[i].id)->nb_elements_skin=S[i].mesh->skin.elem_list.size();
        geometry_user.find_group_elements(S[i].id)->nb_nodes_skin=S[i].mesh->skin.node_list.size();
        convert_mesh_skin_to_geometry_user(S[i], geometry_user);
        S[i].mesh.unload();

        //creation de la boite englobant le maillage de la Sst (utile pour la suite)
//         if (process.rank == 0) S[i].box=create_box_mesh(*S[i].mesh.m);
//         if (process.size > 1)
//             for(unsigned j=0 ;j<S[i].box.size() ;j++ ){
//                 MPI_Bcast(S[i].box[j].ptr(),S[i].box[j].size(),MPI_DOUBLE, 0,MPI_COMM_WORLD);
//             }
//         if (process.rank == 0) S[i].mesh.unload();
    }
}


/** \ingroup Maillage_geometrie_inter
\brief Création des interfaces comprises entre les Sst. Par défaut elles sont considérées comme parfaites.
 
Le maillage de l'interface est obtenu par la lecture de geometry_user issue de SC_create_2
 
*/

template<class TM>
void read_mesh_interface_geometry_user(TM &mesh, GeometryUser &geometry_user, int num_inter) throw(std::runtime_error) {
    //TM mesh;
    typedef typename TM::Tpos T;
    typedef typename TM::Pvec Pvec;
    typedef typename TM::TNode TNode;
    typedef typename TM::EA EA;
  
    // obtaining nbnode, nbelem
    unsigned nbnode = geometry_user.find_group_interfaces(num_inter)->map_global_nodes.size();
    unsigned nbelem = geometry_user.find_group_interfaces(num_inter)->nb_interfaces;
  
    //ajout des noeuds au maillage
    map<int,TNode *> map_num_node;
    Vec<TYPE,DIM> vec;
    for(int i_node=0; i_node<nbnode; i_node++){
        for(unsigned d=0; d<DIM; d++){
            vec[d] = geometry_user.find_group_interfaces(num_inter)->local_nodes[d][i_node];
        }
        map_num_node[i_node] = mesh.add_node(vec);
    }
    
    //ajout des elements
    switch (geometry_user.find_group_interfaces(num_inter)->interface_base_id){
        //for Bar
        case 0 :{
            int nb_node_elem = 2;
            Vec<TNode *> vn;
            vn.resize(nb_node_elem);
            for(int i_elem=0; i_elem<nbelem; i_elem++) {
                for(int i_node=0; i_node<nb_node_elem; i_node++) {
                    vn[i_node] = map_num_node[geometry_user.find_group_interfaces(num_inter)->local_connectivities[i_node][i_elem]];
                }
                typename TM::EA *ne = reinterpret_cast<typename TM::EA *>(mesh.add_element(Bar(),DefaultBehavior(),&vn[0]));
                ne->group = num_inter;
            }
            break;
        }
        //for Bar_3
        case 1 :{
            int nb_node_elem = 3;
            Vec<TNode *> vn;
            vn.resize(nb_node_elem);
            for(int i_elem=0; i_elem<nbelem; i_elem++) {
                for(int i_node=0; i_node<nb_node_elem; i_node++) {
                    vn[i_node] = map_num_node[geometry_user.find_group_interfaces(num_inter)->local_connectivities[i_node][i_elem]];
                }
                typename TM::EA *ne = reinterpret_cast<typename TM::EA *>(mesh.add_element(Bar_3(),DefaultBehavior(),&vn[0]));
                ne->group = num_inter;
            }
            break;
        }
        //for Triangle
        case 2 :{
            int nb_node_elem = 3;
            Vec<TNode *> vn;
            vn.resize(nb_node_elem);
            for(int i_elem=0; i_elem<nbelem; i_elem++) {
                for(int i_node=0; i_node<nb_node_elem; i_node++) {
                    vn[i_node] = map_num_node[geometry_user.find_group_interfaces(num_inter)->local_connectivities[i_node][i_elem]];
                }
                typename TM::EA *ne = reinterpret_cast<typename TM::EA *>(mesh.add_element(Triangle(),DefaultBehavior(),&vn[0]));
                ne->group = num_inter;
            }
            break;
        }
        //for Triangle_6
        case 3 :{
            int nb_node_elem = 6;
            Vec<TNode *> vn;
            vn.resize(nb_node_elem);
            for(int i_elem=0; i_elem<nbelem; i_elem++) {
                for(int i_node=0; i_node<nb_node_elem; i_node++) {
                    vn[i_node] = map_num_node[geometry_user.find_group_interfaces(num_inter)->local_connectivities[i_node][i_elem]];
                }
                typename TM::EA *ne = reinterpret_cast<typename TM::EA *>(mesh.add_element(Triangle_6(),DefaultBehavior(),&vn[0]));
                ne->group = num_inter;
            }
            break;
        }
        default :{
            std::cerr << "type de pattern non implemente" << std::endl; assert(0);                    
        }
    }
}

template<class TV1, class TV2>
void create_perfect_interfaces(DataUser &data_user, GeometryUser &geometry_user, TV1 &S, TV2 &Inter, Param process) {  
    // initialisation de la taille du vecteur d'interfaces
    int nb_inter = 0;
    BasicVec< int > rep_id_inter;
    for(int i_group=0; i_group<geometry_user.group_interfaces.size(); i_group++){
        if(geometry_user.group_interfaces[i_group].type == 2){
            nb_inter += 1;
            rep_id_inter.push_back(geometry_user.group_interfaces[i_group].id);
        }
    }
    
    // assignation des id d'interfaces et des reférence vers un comportement
    Inter.resize(nb_inter);
    for(int i_inter=0; i_inter<nb_inter; i_inter++){
        Inter[i_inter].num = rep_id_inter[i_inter];
        Inter[i_inter].id = rep_id_inter[i_inter];
        int id_link = data_user.find_group_interfaces(Inter[i_inter].id)->id_link;
        int index_link = data_user.find_index_behaviour_links(id_link);
        Inter[i_inter].id_link = id_link ; 
        
        //ajout des numeros des Sst voisines et cotes correspondants
        typename TV1::template SubType<0>::T::Edge edge;
        edge.internum=i_inter;

        BasicVec < int > id_sst;
        BasicVec < int > index_sst;
        Inter[i_inter].side.resize(2);
        id_sst.resize(2);
        index_sst.resize(2);
        for(int i_group=0; i_group<2; i_group++){
            id_sst[i_group] = geometry_user.find_group_interfaces(Inter[i_inter].num)->group_elements_id[i_group];
            index_sst[i_group] = find_index_sst(S, id_sst[i_group]);
        }
        edge.datanum=0;
        find_sst(S, id_sst[0])->edge.push_back(edge);
        find_sst(S, id_sst[0])->vois.push_back(index_sst[1]);
        edge.datanum=1;
        find_sst(S, id_sst[1])->edge.push_back(edge);
        find_sst(S, id_sst[1])->vois.push_back(index_sst[0]);
        
        Inter[i_inter].vois=Vec<unsigned>(index_sst[0], S[index_sst[0]].edge.size()-1, index_sst[1], S[index_sst[1]].edge.size()-1);
        Inter[i_inter].side[0].vois=Vec<unsigned>(index_sst[0], S[index_sst[0]].edge.size()-1);
        Inter[i_inter].side[1].vois=Vec<unsigned>(index_sst[1], S[index_sst[1]].edge.size()-1);
        
        //modification des proprietes de l'interface
        Inter[i_inter].type="Int";
        Inter[i_inter].comp="Parfait";
        
        Inter[i_inter].num=geometry_user.find_group_interfaces(Inter[i_inter].id)->nb_interfaces; // pourquoi ???
    }
    
    
    // construction du maillage des interfaces
//     for(int i_inter=0; i_inter<nb_inter; i_inter++){
//         typename TV2::template SubType<0>::T::TMESH meshnew;
//         read_mesh_interface_geometry_user(meshnew, geometry_user, Inter[i_inter].id); 
//         
//         Inter[i_inter].side[0].mesh=new typename TV2::template SubType<0>::T::TMESH;
//         Inter[i_inter].side[0].mesh->append(meshnew);
//         Inter[i_inter].side[0].mesh->sub_mesh(LMT::Number<1>()).elem_list.change_hash_size( *Inter[i_inter].side[0].mesh,1);
//         Inter[i_inter].side[0].mesh->sub_mesh(LMT::Number<2>()).elem_list.change_hash_size( *Inter[i_inter].side[0].mesh,1);
//         Inter[i_inter].side[0].mesh->sub_mesh(LMT::Number<0>()).elem_list.change_hash_size( *Inter[i_inter].side[0].mesh,1);
//         Inter[i_inter].side[0].mesh->elem_list.change_hash_size( *Inter[i_inter].side[0].mesh,1);
//         Inter[i_inter].side[1].mesh=Inter[i_inter].side[0].mesh;
//         
//         Inter[i_inter].num=meshnew.elem_list.size();
//     }
}

/** \ingroup Maillage_geometrie_inter
\brief Création des interfaces contenues dans une condition aux limites donnée
 
On recherche les Ssts ayant une intersection avec chacune des conditions aux limites (intersection de boites). On extrait alors du maillage de peau de chaque Ssts sélectionnée le maillage contenu dans une condition aux limites donnée. 
*/
template<class TV1, class TV2, class TV5>
void create_interfaces_CL(DataUser &data_user, GeometryUser &geometry_user, TV1 &S, TV2 &Inter, TV5 &CL, Param &process) {
    typedef Vec<TYPEREEL, TV1::template SubType<0>::T::dim> Pvec;

    // initialisation de la taille du vecteur d'interfaces
    int nb_inter = 0;
    BasicVec< int > rep_id_inter;
    for(int i_group=0; i_group<geometry_user.group_interfaces.size(); i_group++){
        if(geometry_user.group_interfaces[i_group].type == 0 and geometry_user.group_interfaces[i_group].is_splited == 0 and geometry_user.group_interfaces[i_group].edge_id != -1 and geometry_user.group_interfaces[i_group].edge_id != -2){
            nb_inter += 1;
            rep_id_inter.push_back(geometry_user.group_interfaces[i_group].id);
        }
    }
    int nb_inter_actuel =  Inter.size();
    int nb_inter_total =  nb_inter_actuel + nb_inter;
    Inter.resize(nb_inter_total);
    
    PRINT(nb_inter_actuel);
    PRINT(nb_inter);
    PRINT(Inter.size());
    // assignation des id d'interfaces et des reférence vers un comportement
    for(int i_inter=0; i_inter<nb_inter; i_inter++){

        int num_inter = nb_inter_actuel + i_inter;
        Inter[num_inter].num = rep_id_inter[i_inter];
        Inter[num_inter].id = rep_id_inter[i_inter]; 
        Inter[num_inter].type="Ext";
        Inter[num_inter].edge_id = geometry_user.find_group_interfaces(Inter[num_inter].id)->edge_id; 
        int id_bc = data_user.find_group_edges( Inter[num_inter].edge_id )->id_CL ;
        int index_bc = data_user.find_index_behaviour_bc(id_bc);
        Inter[num_inter].id_bc = id_bc;
        Inter[num_inter].refCL = index_bc;
        Inter[num_inter].comp = data_user.behaviour_bc[index_bc].type;
        
        //ajout des numeros des Sst voisines et cotes correspondants
        typename TV1::template SubType<0>::T::Edge edge;
        edge.internum=num_inter; 
        BasicVec < int > id_sst;
        BasicVec < int > index_sst;
        Inter[num_inter].side.resize(1);
        id_sst.resize(1);
        index_sst.resize(1);
        for(int i_group=0; i_group<1; i_group++){
            id_sst[i_group] = geometry_user.find_group_interfaces(Inter[num_inter].num)->group_elements_id[i_group];
            index_sst[i_group] = find_index_sst(S, id_sst[i_group]);
        }
        edge.datanum=0;
        find_sst(S, id_sst[0])->edge.push_back(edge);
        find_sst(S, id_sst[0])->vois.push_back(-1);
        Inter[num_inter].vois=Vec<unsigned>(index_sst[0], S[index_sst[0]].edge.size()-1);
        Inter[num_inter].side[0].vois=Vec<unsigned>(index_sst[0], S[index_sst[0]].edge.size()-1);
        Inter[num_inter].num=geometry_user.find_group_interfaces(Inter[num_inter].id)->nb_interfaces; // pourquoi ???
    }
    
    // construction du maillage des interfaces
//     for(int i_inter=0; i_inter<nb_inter; i_inter++){
//         int num_inter = nb_inter_actuel + i_inter;
//         typename TV2::template SubType<0>::T::TMESH meshnew;
//         read_mesh_interface_geometry_user(meshnew, geometry_user, Inter[num_inter].id); 
//         
//         Inter[num_inter].side[0].mesh=new typename TV2::template SubType<0>::T::TMESH;
//         Inter[num_inter].side[0].mesh->append(meshnew);
//         Inter[num_inter].side[0].mesh->sub_mesh(LMT::Number<1>()).elem_list.change_hash_size( *Inter[num_inter].side[0].mesh,1);
//         Inter[num_inter].side[0].mesh->sub_mesh(LMT::Number<2>()).elem_list.change_hash_size( *Inter[num_inter].side[0].mesh,1);
//         Inter[num_inter].side[0].mesh->sub_mesh(LMT::Number<0>()).elem_list.change_hash_size( *Inter[num_inter].side[0].mesh,1);
//         Inter[num_inter].side[0].mesh->elem_list.change_hash_size( *Inter[num_inter].side[0].mesh,1);
//         
//         Inter[num_inter].num=meshnew.elem_list.size();
//     }
}

struct make_interface_inter {
    template <class TV2>
    void operator()(TV2 &SubI, GeometryUser &geometry_user) const {
        if (SubI.comp=="Parfait"){
            if (SubI.side[0].mesh == NULL){             
                SubI.side[0].mesh=new typename TV2::TMESH;
                read_mesh_interface_geometry_user(*SubI.side[0].mesh, geometry_user, SubI.id); 
                SubI.side[0].mesh->sub_mesh(LMT::Number<1>()).elem_list.change_hash_size( *SubI.side[0].mesh,1);
                SubI.side[0].mesh->sub_mesh(LMT::Number<2>()).elem_list.change_hash_size( *SubI.side[0].mesh,1);
                SubI.side[0].mesh->sub_mesh(LMT::Number<0>()).elem_list.change_hash_size( *SubI.side[0].mesh,1);
                SubI.side[0].mesh->elem_list.change_hash_size( *SubI.side[0].mesh,1);
                SubI.side[1].mesh=SubI.side[0].mesh;
            }
        }
    }
};

struct make_interface_CL {
    template <class TV2>
    void operator()(TV2 &SubI, GeometryUser &geometry_user) const {
        if (SubI.type=="Ext"){
            if (SubI.side[0].mesh==NULL){
                SubI.side[0].mesh=new typename TV2::TMESH;
                read_mesh_interface_geometry_user(*SubI.side[0].mesh, geometry_user, SubI.id); 
                SubI.side[0].mesh->sub_mesh(LMT::Number<1>()).elem_list.change_hash_size( *SubI.side[0].mesh,1);
                SubI.side[0].mesh->sub_mesh(LMT::Number<2>()).elem_list.change_hash_size( *SubI.side[0].mesh,1);
                SubI.side[0].mesh->sub_mesh(LMT::Number<0>()).elem_list.change_hash_size( *SubI.side[0].mesh,1);
                SubI.side[0].mesh->elem_list.change_hash_size( *SubI.side[0].mesh,1);
            }
        }
    }
};

struct mesh_unload {
    template <class TV1>
    void operator()(TV1 &S) const {
        S.mesh.unload();
    }
};

//******************************************************************************************************
// Programme permettant la construction des SST et Interfaces aveclecture de maillage issu de castem
//******************************************************************************************************
/** \ingroup Maillage_geometrie
\brief Création des Sous-structures et des Interfaces (maillages).
 
On construit ici le vecteur des sous-structures et des interfaces à partir des informations du fichier xml.
 
- Assignation du type de matériau aux sous-structures create_SST_typmat().
- Lecture du maillage des sous-structures create_maillage_SST().
- creation des interfaces parfaites create_perfect_interfaces().
- creation des interfaces avec jeu physique create_gap_interfaces().
- creation des interfaces de conditions aux limites. create_interfaces_CL().
*/
                     
#ifndef INFO_TIME
#define INFO_TIME
#endif 
#include "containers/evaluate_nb_cycles.h"

template<class TV1, class TV2, class TV3, class TV4,class TV5>
void create_SST_INTER(DataUser &data_user, GeometryUser &geometry_user, TV1 &S,TV2 &Inter, TV5 &CL, Param &process, TV3 &Stot,TV3 &SubS,TV4 &SubI) {
#ifdef INFO_TIME
    if (process.size>1) MPI_Barrier(MPI_COMM_WORLD);
    TicToc tic1;
    if (process.rank==0) tic1.start();
    
    
#endif
    if (process.rank == 0)
        std::cout << "Creation de la geometrie des SST" <<std::endl;
    create_SST_typmat(data_user, geometry_user,S,process);      // to be TEST

#ifdef INFO_TIME
    if (process.size>1) MPI_Barrier(MPI_COMM_WORLD);
    if (process.rank==0) std::cout << "create_SST_typmat : " <<std::endl;
    if (process.rank==0) tic1.stop();;
    if (process.rank==0) tic1.start();
#endif
    
    
    if (process.rank == 0)
        std::cout << "\t Lecture maillages des SST" <<std::endl;
    create_maillage_SST(data_user, geometry_user,S,process);    // to be TEST
#ifdef INFO_TIME
    if (process.size>1) MPI_Barrier(MPI_COMM_WORLD);
    if (process.rank==0) std::cout << "Creation maillage : " ;
    if (process.rank==0) tic1.stop();;
    if (process.rank==0) tic1.start();
#endif


    if (process.rank == 0)
        std::cout << "Creation des Interfaces" <<std::endl;
    if (process.rank == 0)
        std::cout << "\tParfaites" << std::endl;
    create_perfect_interfaces(data_user, geometry_user, S, Inter, process);     // to be TEST
#ifdef INFO_TIME
    if (process.size>1) MPI_Barrier(MPI_COMM_WORLD);
    if (process.rank==0) std::cout << "Creation interfaces parfaites : " ;
    if (process.rank==0) tic1.stop();;
    if (process.rank==0) tic1.start();
#endif

    //create_other_interfaces(structure,S,Inter);
    if (process.rank == 0)
        std::cout << "\tConditions aux limites" << std::endl;
    create_interfaces_CL(data_user, geometry_user, S,Inter,CL, process);        // to be TEST
#ifdef INFO_TIME
    if (process.size>1) MPI_Barrier(MPI_COMM_WORLD);
    if (process.rank==0) std::cout << "Creation interfaces conditions limites : " ;
    if (process.rank==0) tic1.stop();
    if (process.rank==0) tic1.start();
#endif
    
    //make subs et subi
    //assignation_mpi...
    mpi_repartition(S, Inter,process,Stot,SubS,SubI);                            // to be TEST
#ifdef INFO_TIME
    if (process.size>1) MPI_Barrier(MPI_COMM_WORLD);
    if (process.rank==0) std::cout << "MPI Repartition : " ;
    if (process.rank==0) tic1.stop();
    if (process.rank==0) tic1.start();
#endif

    //create mesh of interface perfect (only necessary)
//     if (process.size > 1) {
        apply(SubI,make_interface_inter(), geometry_user );                              // to be TEST
        apply(SubI,make_interface_CL(), geometry_user);                                  // to be TEST
        apply(S,mesh_unload());                                                 
//     }

    if (process.size>1) MPI_Barrier(MPI_COMM_WORLD);
    if (process.rank==0) std::cout << "\t Assignation des numeros aux interfaces " << std::endl;
    //assignation du numero de l'interface
    for(unsigned i=0;i<Inter.size();i++){
        Inter[i].num=i;
    }

};

