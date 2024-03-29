#include "create_SST_INTER.h"
#include <fstream>
#include <sstream>
#include "../../LMT/include/mesh/read_avs.h"
#include "../../LMT/include/mesh/read_geof.h"
#include "../../LMT/include/mesh/write_avs.h"
#include "../../LMT/include/mesh/read_mshadv.h"
#include "../../LMT/include/containers/algo.h"
#include <map>
#include "../../LMT/include/containers/vecpointedvalues.h"


#include "../DEFINITIONS/GeneralData.h"
#include "../MPI/assignation_mpi.h"
#include "../DEFINITIONS/SavingData.h"

void create_SST_typmat(DataUser &data_user, GeometryUser &geometry_user,Vec<Sst> &S,Process &process) {
    
    ///initialisation de la taille des sst et de leur id == num   
    //S.resize(geometry_user.nb_group_elements);
    
    for(unsigned i=0;i<S.size();i++) {
        S[i].num = geometry_user.group_elements[i].id;
        S[i].id = geometry_user.group_elements[i].id;
    }
    
    /// assignation des numero de materiaux aux sst
    for(int i_data_group=0; i_data_group<data_user.group_elements.size(); i_data_group++){
        int id_mat = data_user.find_materials_index(data_user.group_elements[i_data_group].id_material)->id;
        int index_mat = data_user.find_materials_index(data_user.group_elements[i_data_group].id_material);
        Sst::find_sst(S,data_user.group_elements[i_data_group].id)->id_material = id_mat ;
        Sst::find_sst(S,data_user.group_elements[i_data_group].id)->typmat = index_mat ; 
    }
}



void convert_mesh_skin_to_geometry_user(Sst &S, GeometryUser &geometry_user){
    int id_sst=S.id;
    int nb_nodes_by_element_skin=(geometry_user.patterns.find_type(geometry_user.find_group_elements(id_sst)->pattern_base_id)).nb_nodes_by_sides;
    geometry_user.find_group_elements(id_sst)->local_connectivities_skin.resize(nb_nodes_by_element_skin);
    for(int ne=0;ne<nb_nodes_by_element_skin;ne++){
        geometry_user.find_group_elements(id_sst)->local_connectivities_skin[ne].resize(S.mesh->skin.elem_list.size());
    }
    apply(S.mesh->skin.elem_list,ConvertMeshConnectivitiesSkin(),geometry_user, id_sst);
    
}


void create_maillage_SST(DataUser &data_user, GeometryUser &geometry_user, Vec<Sst> &S, Process &process) {
    for(unsigned i=0;i<S.size();++i) {
        Sc2String namein = data_user.find_group_elements_pointer(S[i].num)->name;
        S[i].mesh.name=namein;
        S[i].mesh.load(geometry_user, S[i].id);
        S[i].mesh.load();
        S[i].mesh->update_skin();
        geometry_user.find_group_elements_pointer(S[i].id)->nb_elements_skin=S[i].mesh->skin.elem_list.size();
        geometry_user.find_group_elements_pointer(S[i].id)->nb_nodes_skin=S[i].mesh->skin.node_list.size();
        convert_mesh_skin_to_geometry_user(S[i], geometry_user);
        S[i].mesh.unload();
    }
}


void read_mesh_interface_geometry_user(Interface::TMESH &mesh, GeometryUser &geometry_user, int num_inter) throw(std::runtime_error) {
    //TM mesh;
    typedef Interface::TMESH::Tpos T;
    typedef Interface::TMESH::Pvec Pvec;
    typedef Interface::TMESH::TNode TNode;
    typedef Interface::TMESH::EA EA;
    
    // obtaining nbnode, nbelem
    unsigned nbnode = geometry_user.find_group_interfaces(num_inter)->map_global_nodes.size();
    unsigned nbelem = geometry_user.find_group_interfaces(num_inter)->nb_interfaces;
    
    //ajout des noeuds au maillage
    map<int,TNode *> map_num_node;
    Vec<TYPEREEL,DIM> vec;
    for(int i_node=0; i_node<nbnode; i_node++){
        for(unsigned d=0; d<DIM; d++){
            vec[d] = geometry_user.find_group_interfaces(num_inter)->local_nodes[d][i_node];
        }
        map_num_node[i_node] = mesh.add_node(vec);
    }
    
    ///ajout des elements
    switch (geometry_user.find_group_interfaces(num_inter)->interface_base_id){
        ///for Bar
        case 0 :{
            int nb_node_elem = 2;
            Vec<TNode *> vn;
            vn.resize(nb_node_elem);
            for(int i_elem=0; i_elem<nbelem; i_elem++) {
                for(int i_node=0; i_node<nb_node_elem; i_node++) {
                    vn[i_node] = map_num_node[geometry_user.find_group_interfaces(num_inter)->local_connectivities[i_node][i_elem]];
                }
                EA *ne = reinterpret_cast<EA *>(mesh.add_element(Bar(),DefaultBehavior(),&vn[0]));
                ne->group = num_inter;
            }
            break;
        }
        ///for Bar_3
        case 1 :{
            int nb_node_elem = 3;
            Vec<TNode *> vn;
            vn.resize(nb_node_elem);
            for(int i_elem=0; i_elem<nbelem; i_elem++) {
                for(int i_node=0; i_node<nb_node_elem; i_node++) {
                    vn[i_node] = map_num_node[geometry_user.find_group_interfaces(num_inter)->local_connectivities[i_node][i_elem]];
                }
                EA *ne = reinterpret_cast<EA *>(mesh.add_element(Bar_3(),DefaultBehavior(),&vn[0]));
                ne->group = num_inter;
            }
            break;
        }
        ///for Triangle
        case 2 :{
            int nb_node_elem = 3;
            Vec<TNode *> vn;
            vn.resize(nb_node_elem);
            for(int i_elem=0; i_elem<nbelem; i_elem++) {
                for(int i_node=0; i_node<nb_node_elem; i_node++) {
                    vn[i_node] = map_num_node[geometry_user.find_group_interfaces(num_inter)->local_connectivities[i_node][i_elem]];
                }
                EA *ne = reinterpret_cast<EA *>(mesh.add_element(Triangle(),DefaultBehavior(),&vn[0]));
                ne->group = num_inter;
            }
            break;
        }
        ///for Triangle_6
        case 3 :{
            int nb_node_elem = 6;
            Vec<TNode *> vn;
            vn.resize(nb_node_elem);
            for(int i_elem=0; i_elem<nbelem; i_elem++) {
                for(int i_node=0; i_node<nb_node_elem; i_node++) {
                    vn[i_node] = map_num_node[geometry_user.find_group_interfaces(num_inter)->local_connectivities[i_node][i_elem]];
                }
                EA *ne = reinterpret_cast<EA *>(mesh.add_element(Triangle_6(),DefaultBehavior(),&vn[0]));
                ne->group = num_inter;
            }
            break;
        }
        default :{
            std::cerr << "type de pattern non implemente" << std::endl; assert(0);                    
        }
    }
}


void create_perfect_interfaces(DataUser &data_user, GeometryUser &geometry_user, Vec<Sst> &S, Vec<Interface> &Inter, Process &process) {  
    // initialisation de la taille du vecteur d'interfaces
    int nb_inter = 0;
    BasicVec< int > rep_id_inter;
    for(int i_group=0; i_group<geometry_user.group_interfaces.size(); i_group++){
        if(geometry_user.group_interfaces[i_group].type == 2){
            nb_inter += 1;
            rep_id_inter.push_back(geometry_user.group_interfaces[i_group].id);
        }
    }
    
    // assignation des id d'interfaces et des reference vers un comportement
    Inter.resize(nb_inter);
    for(int i_inter=0; i_inter<nb_inter; i_inter++){
        Inter[i_inter].num = rep_id_inter[i_inter];
        Inter[i_inter].id = rep_id_inter[i_inter];
        int id_link = data_user.find_interfaces_pointer(Inter[i_inter].id)->link_id;
        int index_link = data_user.find_links_index(id_link);
        Inter[i_inter].id_link = id_link ; 
        
        //ajout des numeros des Sst voisines et cotes correspondants
        Sst::Edge edge;
        edge.internum=i_inter;
        
        BasicVec < int > id_sst;
        BasicVec < int > index_sst;
        Inter[i_inter].side.resize(2);
        id_sst.resize(2);
        index_sst.resize(2);
        for(int i_group=0; i_group<2; i_group++){
            id_sst[i_group] = geometry_user.find_group_interfaces(Inter[i_inter].num)->group_elements_id[i_group];
            index_sst[i_group] = data_user.find_pieces_index(id_sst[i_group]);
        }
        edge.datanum=0;
        S[index_sst[0]].edge.push_back(edge);
        S[index_sst[0]].vois.push_back(index_sst[1]);
        edge.datanum=1;
        S[index_sst[1]].edge.push_back(edge);
        S[index_sst[1]].vois.push_back(index_sst[0]);
        
        Inter[i_inter].vois=Vec<unsigned>(index_sst[0], S[index_sst[0]].edge.size()-1, index_sst[1], S[index_sst[1]].edge.size()-1);
        Inter[i_inter].side[0].vois=Vec<unsigned>(index_sst[0], S[index_sst[0]].edge.size()-1);
        Inter[i_inter].side[1].vois=Vec<unsigned>(index_sst[1], S[index_sst[1]].edge.size()-1);
        
        //modification des proprietes de l'interface
        Inter[i_inter].type="Int";
        Inter[i_inter].comp="Parfait";
    }
}


void create_interfaces_CL(DataUser &data_user, GeometryUser &geometry_user, Vec<Sst> &S, Vec<Interface> &Inter, Vec<Boundary> &CL, Process &process) {
    // initialisation de la taille du vecteur d'interfaces
    int nb_inter = 0;
    BasicVec< int > rep_id_inter;
    for(int i_group=0; i_group<geometry_user.group_interfaces.size(); i_group++){
        if(geometry_user.group_interfaces[i_group].type == 0 and geometry_user.group_interfaces[i_group].is_splited == 0 and geometry_user.group_interfaces[i_group].edge_id != -1 and geometry_user.group_interfaces[i_group].edge_id != -2){
            nb_inter += 1;
            rep_id_inter.push_back(geometry_user.group_interfaces[i_group].id);
        }
    }
    const int nb_inter_actuel =  Inter.size();
    const int nb_inter_total =  nb_inter_actuel + nb_inter;
    Inter.resize(nb_inter_total);
    // assignation des id d'interfaces et des reference vers un comportement
    for(int i_inter=0; i_inter<nb_inter; i_inter++){
        
        int num_inter = nb_inter_actuel + i_inter;
        Inter[num_inter].num = rep_id_inter[i_inter];
        Inter[num_inter].id = rep_id_inter[i_inter]; 
        Inter[num_inter].type="Ext";
        Inter[num_inter].edge_id = geometry_user.find_group_interfaces(Inter[num_inter].id)->edge_id; 
        int id_bc = data_user.find_edges_pointer( Inter[num_inter].edge_id )->boundary_condition_id ;
        int index_bc = data_user.find_boundary_conditions_index(id_bc);
        Inter[num_inter].id_bc = id_bc;
        Inter[num_inter].refCL = index_bc;
        Inter[num_inter].comp = data_user.boundary_conditions_vec[index_bc].condition_type;
        
        //ajout des numeros des Sst voisines et cotes correspondants
        Sst::Edge edge;
        edge.internum=num_inter; 
        BasicVec < int > id_sst;
        BasicVec < int > index_sst;
        Inter[num_inter].side.resize(1);
        id_sst.resize(1);
        index_sst.resize(1);
        for(int i_group=0; i_group<1; i_group++){
            id_sst[i_group] = geometry_user.find_group_interfaces(Inter[num_inter].num)->group_elements_id[i_group];
            index_sst[i_group] = data_user.find_pieces_index(id_sst[i_group]);
        }
        edge.datanum=0;
        S[index_sst[0]].edge.push_back(edge);
        S[index_sst[0]].vois.push_back(-1);
        Inter[num_inter].vois=Vec<unsigned>(index_sst[0], S[index_sst[0]].edge.size()-1);
        Inter[num_inter].side[0].vois=Vec<unsigned>(index_sst[0], S[index_sst[0]].edge.size()-1);
    }
}


void make_interface_inter::operator()(Interface &SubI, GeometryUser &geometry_user) const {
    if (SubI.comp=="Parfait"){
        if (SubI.side[0].mesh == NULL){             
            SubI.side[0].mesh=new Interface::TMESH;
            read_mesh_interface_geometry_user(*SubI.side[0].mesh, geometry_user, SubI.id); 
            SubI.side[0].mesh->sub_mesh(LMT::Number<1>()).elem_list.change_hash_size( *SubI.side[0].mesh,1);
            SubI.side[0].mesh->sub_mesh(LMT::Number<2>()).elem_list.change_hash_size( *SubI.side[0].mesh,1);
            SubI.side[0].mesh->sub_mesh(LMT::Number<0>()).elem_list.change_hash_size( *SubI.side[0].mesh,1);
            SubI.side[0].mesh->elem_list.change_hash_size( *SubI.side[0].mesh,1);
            SubI.side[1].mesh=SubI.side[0].mesh;
        }
    }
}


void make_interface_CL::operator()(Interface &SubI, GeometryUser &geometry_user) const {
    if (SubI.type=="Ext"){
        if (SubI.side[0].mesh==NULL){
            SubI.side[0].mesh=new Interface::TMESH;
            read_mesh_interface_geometry_user(*SubI.side[0].mesh, geometry_user, SubI.id); 
            SubI.side[0].mesh->sub_mesh(LMT::Number<1>()).elem_list.change_hash_size( *SubI.side[0].mesh,1);
            SubI.side[0].mesh->sub_mesh(LMT::Number<2>()).elem_list.change_hash_size( *SubI.side[0].mesh,1);
            SubI.side[0].mesh->sub_mesh(LMT::Number<0>()).elem_list.change_hash_size( *SubI.side[0].mesh,1);
            SubI.side[0].mesh->elem_list.change_hash_size( *SubI.side[0].mesh,1);
        }
    }
}


void mesh_unload::operator()(Sst &S) const {
    S.mesh.unload();
}


void create_SST_INTER(DataUser                          &data_user, 
                      GeometryUser                      &geometry_user, 
                      Vec<Sst>                          &S,
                      Vec<Interface>                    &Inter, 
                      Vec<Boundary>                     &CL, 
                      Process                           &process, 
                      Vec<VecPointedValues<Sst> >       &Stot,
                      Vec<VecPointedValues<Sst> >       &SubS,
                      Vec<VecPointedValues<Interface> > &SubI) {
    /// Timer pour le chronometrage
    #ifdef INFO_TIME
    process.parallelisation->synchronisation();
    TicTac tic1;
    if (process.parallelisation->is_master_cpu()) tic1.start();
    #endif
    
    /// Assignation des id de materiaux des SST
    process.print(" - Assignation des id de materiaux des SST : ... ",true);
    create_SST_typmat(data_user, geometry_user,S,process);      // to be TEST
    #ifdef INFO_TIME
    process.print_duration(tic1);
    #endif
    
    /// Creation des maillages des SST
    process.print(" - Creation des maillages des SST : ... ",true);
    create_maillage_SST(data_user, geometry_user,S,process);    // to be TEST
    #ifdef INFO_TIME
    process.print_duration(tic1);
    #endif
    
    /// Creation des maillages des interfaces ...
    process.print(" - Creation des Interfaces : ");
    /// ... interieures ...
    process.print("     - Interieures parfaites : ... ",true);
    create_perfect_interfaces(data_user, geometry_user, S, Inter, process);     // to be TEST
    #ifdef INFO_TIME
    process.print_duration(tic1);
    #endif
    //create_other_interfaces(structure,S,Inter);
    /// ... et exterieures
    process.print("     - Exterieures avec CL : ... ",true);
    create_interfaces_CL(data_user, geometry_user, S,Inter,CL, process);        // to be TEST
    #ifdef INFO_TIME
    process.print_duration(tic1);
    #endif
    
    /// Repartition MPI
    //make subs et subi
    //assignation_mpi...
    //process.parallelisation->repartition_sst.free();  // un bug memoire se produit sinon
    process.print(" - Repartition MPI : ...",true);
    mpi_repartition(S, Inter,process,Stot,SubS,SubI, geometry_user);                            // to be TEST
    #ifdef INFO_TIME
    process.print_duration(tic1);
    #endif
    
    apply(SubI,make_interface_inter(), geometry_user );                              // to be TEST
    apply(SubI,make_interface_CL(), geometry_user);                                  // to be TEST
    apply(S,mesh_unload());
        
    process.parallelisation->synchronisation();
    process.print(" - Assignation des numeros aux interfaces ",true);
    ///assignation du numero de l'interface
    for(unsigned i=0;i<Inter.size();i++){
        Inter[i].num=i;
    }
    ///copie des "id" et "side" des interfaces adjacentes aux groupes d'elements
    for(unsigned i_sst=0;i_sst<S.size();i_sst++){
        int id_sst=S[i_sst].id;
        geometry_user.find_group_elements(id_sst)->id_adjacent_group_interfaces.resize(S[i_sst].edge.size());
        geometry_user.find_group_elements(id_sst)->side_adjacent_group_interfaces.resize(S[i_sst].edge.size());
        for(unsigned i_side=0;i_side<S[i_sst].edge.size();i_side++) {
            int id_interface=S[i_sst].edge[i_side].internum;
            int side_interface=S[i_sst].edge[i_side].datanum;
            geometry_user.find_group_elements(id_sst)->id_adjacent_group_interfaces[i_side]=Inter[id_interface].id;
            geometry_user.find_group_elements(id_sst)->side_adjacent_group_interfaces[i_side]=side_interface;
        }
    }
        
};
