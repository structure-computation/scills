//
// C++ Interface: GEOMETRY
//
// Description:
//
//
// Author: Jeremie Bellec <j.bellec@structure-computation.com>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include <iostream>
#include <fstream>
#include <cassert>

#include <Metil/BasicVec.h>
#include <Metil/StructCompactor.h>
#include <Metil/Hdf.h>

#include "json_spirit.h"
#include "GeometryUser.h"
#include "DataUser.h"
#include "utils_2.h"

using namespace json_spirit;
using namespace Metil;


//Methode génériques **************************************************************************************************************************************
GroupElementsUser* GeometryUser::find_group_elements(int id_) {
    for (int i_group=0; i_group<group_elements.size(); i_group++) {
        if (group_elements[i_group].id == id_) {
            return &group_elements[i_group];
            break;
        }
    }
}

int GeometryUser::find_index_group_elements(int id_) {
    for (int i_group=0; i_group<group_elements.size(); i_group++) {
        if (group_elements[i_group].id == id_) {
            return i_group;
            break;
        }
    }
}

//recherche d'un group_interfaces particulier avec son id---------------
GroupInterfacesUser* GeometryUser::find_group_interfaces(int id_) {
    for (int i_group=0; i_group<group_interfaces.size(); i_group++) {
        if (group_interfaces[i_group].id == id_) {
            return &group_interfaces[i_group];
            break;
        }
    }
}

int GeometryUser::find_index_group_interfaces(int id_) {
    for (int i_group=0; i_group<group_interfaces.size(); i_group++) {
        if (group_interfaces[i_group].id == id_) {
            return i_group;
            break;
        }
    }
}



//Methode d'initialisation*******************************************************
GeometryUser::GeometryUser() {
    dim = DIM;
}

GeometryUser::GeometryUser(MeshUser &mesh_) {
    std::cout << "**********************************************************************************************************************" << std::endl;
    std::cout << "** create geometry_user  *********************************************************************************************" << std::endl;
    std::cout << "**********************************************************************************************************************" << std::endl;
    dim = DIM;

    std::cout << "** create nodes          *********************************************************************************************" << std::endl;
    nodes = mesh_.elements_pos_nodes;

    std::cout << "** create group_elements *********************************************************************************************" << std::endl;
    initialize_group_elements_from_MeshUser(mesh_);

    std::cout << "** create group_interfaces *******************************************************************************************" << std::endl;
    initialize_group_interfaces_from_MeshUser(mesh_);

    std::cout << "** écriture du json       ********************************************************************************************" << std::endl;
    write_json( mesh_ );

    std::cout << "** écriture du hdf5       ********************************************************************************************" << std::endl;
    write_visu_hdf5( mesh_ );
/*
    std::cout << "** écriture du xdmf       ********************************************************************************************" << std::endl;
    String name_geometry_hdf5 = "/Level_0/Geometry";
    String name_fields_hdf5 = "/Level_0/Structure_Fields";
    write_xdmf(mesh_.name_visu_xdmf, mesh_.name_visu_hdf, name_geometry_hdf5,0, name_fields_hdf5, BasicVec<String>("none"));*/

    std::cout << "** écriture du xdmf       ********************************************************************************************" << std::endl;
    String name_geometry_hdf5 = "/Level_0/Geometry";
    write_xdmf(mesh_.name_visu_xdmf, mesh_.name_visu_hdf, name_geometry_hdf5,0);

}


void GeometryUser::initialize_GPU() {
    std::cout << "** initialisation pour GPU     ***************************************************************************************" << std::endl;
    for(int i_group=0; i_group < group_elements.size(); i_group++){
        group_elements[i_group].initialize_GPU(patterns);
    }
    for(int i_group=0; i_group < group_interfaces.size(); i_group++){
        group_interfaces[i_group].initialize_GPU(patterns);
    }
    std::cout << "** fin initialisation pour GPU ***************************************************************************************" << std::endl;
}


void GeometryUser::initialize_group_elements_from_MeshUser(MeshUser &mesh) {
    for (int i_elem=0; i_elem<mesh.nb_elements; i_elem++) {
        //test pour l'ajout d'un group_elements
        bool add_group_elements = true;
        for (int i_group=0; i_group<group_elements.size(); i_group++) {
            if (group_elements[i_group].add_entity_element(mesh.list_elements[i_elem])) {
                add_group_elements = false;
                break;
            }
        }
        if (add_group_elements) {
            GroupElementsUser current_group(group_elements.size(), mesh.list_elements[i_elem]);
            group_elements.push_back(current_group);
            nb_group_elements = group_elements.size();
        }
    }
}


void GeometryUser::initialize_group_interfaces_from_MeshUser(MeshUser &mesh) {
    for (int i_inter=0; i_inter<mesh.nb_interfaces; i_inter++) {
        mesh.list_interfaces[i_inter].calcul_interface_connectivity(mesh.list_elements, mesh.patterns);
        //test pour l'ajout d'un group_interfaces
        bool add_group_interfaces = true;
        for (int i_group=0; i_group<group_interfaces.size(); i_group++) {
            if (group_interfaces[i_group].add_entity_interface( mesh.list_interfaces[i_inter], mesh.list_elements, group_elements, mesh.patterns )) {
                add_group_interfaces = false;
                break;
            }
        }
        if (add_group_interfaces) {
            GroupInterfacesUser current_group(group_interfaces.size(), mesh.list_interfaces[i_inter], mesh.list_elements, group_elements, mesh.patterns);
            group_interfaces.push_back(current_group);
            nb_group_interfaces = group_interfaces.size();
        }
    }
}


//Methode de niveau sup**********************************************************

bool GeometryUser::do_respect_geometry(int i_group, int num_edge, DataUser::Geometry &geom){
    BasicVec< BasicVec<TYPE,DIM> > vertex_point;
    vertex_point.resize(interfaces.find_type(group_interfaces[i_group].interface_id).nb_vertex_nodes);
    for(int i_node=0; i_node<vertex_point.size(); i_node++){
        for(int d=0; d<DIM; d++){
            vertex_point[i_node][d] = nodes[d][group_interfaces[i_group].connectivities[i_node][num_edge]];
        }
    }
    
    bool S_respect_geometry = true;
    if(geom.type=="all"){
        S_respect_geometry = true;
        return S_respect_geometry;
    }else if(geom.type=="is_in"){
        if(geom.nature=="box"){
            for(int num_point=0; num_point<vertex_point.size(); num_point++){
                if(!pt_in_box(vertex_point[num_point],geom.points)){
                    S_respect_geometry = false;
                    break;
                }
            }
            return S_respect_geometry;
        }else if(geom.nature=="cylinder"){
            for(int num_point=0; num_point<vertex_point.size(); num_point++){
                if(!pt_in_cylinder(vertex_point[num_point],geom.points,geom.radius)){
                    S_respect_geometry = false;
                    break;
                }
            }
            return S_respect_geometry;
        }else if(geom.nature=="sphere"){
            for(int num_point=0; num_point<vertex_point.size(); num_point++){
                if(!pt_in_sphere(vertex_point[num_point],geom.points[0],geom.radius)){
                    S_respect_geometry = false;
                    break;
                }
            }
            return S_respect_geometry;
        }
    }else if(geom.type=="is_on"){
        if(geom.nature=="plan"){
            for(int num_point=0; num_point<vertex_point.size(); num_point++){
                if(!pt_on_plan(vertex_point[num_point],geom.points[0],geom.pdirection)){
                    S_respect_geometry = false;
                    break;
                }
            }
            return S_respect_geometry;
        }else if(geom.nature=="disc"){
            for(int num_point=0; num_point<vertex_point.size(); num_point++){
                if(!pt_on_disc(vertex_point[num_point],geom.points[0],geom.pdirection,geom.radius)){
                    S_respect_geometry = false;
                    break;
                }
            }
            return S_respect_geometry;
        }else if(geom.nature=="cylinder"){
            for(int num_point=0; num_point<vertex_point.size(); num_point++){
                if(!pt_on_cylinder(vertex_point[num_point],geom.points,geom.radius)){
                    S_respect_geometry = false;
                    break;
                }
            }
            return S_respect_geometry;
        }else if(geom.nature=="sphere"){
            for(int num_point=0; num_point<vertex_point.size(); num_point++){
                if(!pt_on_sphere(vertex_point[num_point],geom.points[0],geom.radius)){
                    S_respect_geometry = false;
                    break;
                }
            }
            return S_respect_geometry;
        }else if(geom.nature=="equation"){
            for(int num_point=0; num_point<vertex_point.size(); num_point++){
               if(!pt_match_equation(vertex_point[num_point],geom.equation)){
                   S_respect_geometry = false;
                   break;
               }
            }
            return S_respect_geometry;
//             cout << "equation non implementé" << endl;
        }
    }
};



void GeometryUser::split_group_edges_within_geometry(DataUser &data_user) {
    PRINT(nb_group_interfaces);
    for(int i_group=0; i_group<nb_group_interfaces; i_group++){
        if(group_interfaces[i_group].type == 0){
            BasicVec< int > edge_assigned;
            edge_assigned.resize(group_interfaces[i_group].nb_interfaces,-2);
            BasicVec< bool > find_group_edge;
            find_group_edge.resize(data_user.group_edges.size(),false);
            for(int i_data_group=0; i_data_group<data_user.group_edges.size(); i_data_group++){
                for(int i_edge=0; i_edge<group_interfaces[i_group].nb_interfaces; i_edge++){
                    if(edge_assigned[i_edge]==-2 and do_respect_geometry(i_group, i_edge, data_user.group_edges[i_data_group].geom)){
                        edge_assigned[i_edge] = data_user.group_edges[i_data_group].id;
                        find_group_edge[i_data_group] = true;
                    }
                }
            }
            PRINT(data_user.group_edges.size());
            PRINT(data_user.group_edges[0].id);
            //PRINT(edge_assigned);
            
            for(int i_data_group=0; i_data_group<data_user.group_edges.size(); i_data_group++){
                if(find_group_edge[i_data_group]){
                    PRINT(i_data_group);
                    int id = group_interfaces.size();
                    GroupInterfacesUser group_interface_temp(group_interfaces[i_group], edge_assigned, data_user.group_edges[i_data_group].id, id);
                    group_interfaces.push_back(group_interface_temp);
                }
            }
        }
    }
    nb_group_interfaces = group_interfaces.size();
    PRINT(nb_group_interfaces);
}


//Methode de niveau courrant ****************************************************
//ecriture du fichier json pour l'interface utilisateur javascript---------------
void GeometryUser::write_json(MeshUser &mesh_user) {
    //ecriture du nom du fichier json
    std::string path=mesh_user.name_directory+"/MESH/mesh.txt";
    std::ofstream os( path.c_str() );
    //ecriture du debut du fichier------------------
    Object mesh;
    mesh.push_back(Pair( "model_directory", mesh_user.name_directory ));
    mesh.push_back(Pair( "mesh_directory", mesh_user.mesh_directory ));
    mesh.push_back(Pair( "mesh_name", mesh_user.name_mesh_user ));
    mesh.push_back(Pair( "extension", mesh_user.extension ));

    mesh.push_back(Pair( "nb_sst", mesh_user.nb_elements ));
    mesh.push_back(Pair( "nb_inter", mesh_user.nb_interfaces ));
    mesh.push_back(Pair( "nb_groups_elem", nb_group_elements ));
    mesh.push_back(Pair( "nb_groups_inter", nb_group_interfaces ));

    Object mesh_o;
    mesh_o.push_back( Pair( "mesh", mesh ) );

    //ecriture des groupes d'elements---------------
    Array group_elem_a;
    //Object group_elem_a;
    for (int i_group=0; i_group< nb_group_elements; i_group++) {
        Object group;
        group.push_back(Pair( "id", group_elements[i_group].id ));
        group.push_back(Pair( "origine", "from_bulkdata" ));
        group.push_back(Pair( "identificateur", group_elements[i_group].id ));
        std::ostringstream name_group;
        name_group << "piece " << group_elements[i_group].id;
        group.push_back(Pair( "name", name_group.str() ));
        group.push_back(Pair( "assigned", -1 ));
        group.push_back(Pair( "group", -1 ));

        std::ostringstream num;
        num << i_group;
        //group_elem_a.push_back(Pair(num.str(),group));
        group_elem_a.push_back(group);
    }
    Object group_elem_o;
    group_elem_o.push_back( Pair( "groups_elem", group_elem_a ) );

    //ecriture des groupes d'interfaces------------
    Array group_inter_a;
    //Object group_inter_a;
    for (int i_group=0; i_group< nb_group_interfaces; i_group++) {
        if (group_interfaces[i_group].type==2) {
            Object group;
            group.push_back(Pair( "id", group_interfaces[i_group].id ));
            group.push_back(Pair( "origine", "from_bulkdata" ));
            group.push_back(Pair( "type", "between_group_elem" ));
            std::ostringstream name_group;
            name_group << "interface " << group_interfaces[i_group].id;
            group.push_back(Pair( "name", name_group.str() ));
            group.push_back(Pair( "assigned", -1 ));
            group.push_back(Pair( "group", -1 ));
            std::ostringstream adj;
            adj << group_interfaces[i_group].group_elements_id[0] << " " << group_interfaces[i_group].group_elements_id[1];
            group.push_back(Pair( "adj_num_group", adj.str() ));
            std::ostringstream num;
            num << i_group;
            //group_inter_a.push_back(Pair(num.str(),group));
            group_inter_a.push_back(group);
        }
    }
    Object group_inter_o;
    group_inter_o.push_back( Pair( "groups_inter", group_inter_a ) );

    //ecriture des proprietes_interfaces par defaut---------------
    Array prop_inter_a;
    //Object prop_inter_a;
    for (unsigned i=0;i< 1;i++) {
        Object prop;
        prop.push_back(Pair( "id", 0 ));
        prop.push_back(Pair( "type", "parfait" ));
        prop.push_back(Pair( "coef_frottement", 0. ));
        std::ostringstream num;
        num << i;
        //prop_inter_a.push_back(Pair(num.str(),prop));
        prop_inter_a.push_back(prop);
    }
    Object prop_inter_o;
    prop_inter_o.push_back( Pair( "proprietes_interfaces", prop_inter_a ) );


    //ecriture des CL par defaut-------------------------------
    Array cl_a;
    //Object cl_a;
    for (unsigned i=0;i< 1;i++) {
        Object cl;
        cl.push_back(Pair( "id", 0 ));
        cl.push_back(Pair( "type", "effort" ));
        cl.push_back(Pair( "fcts_spatiales", "0 0" ));
        cl.push_back(Pair( "fcts_temporelles", "0 0" ));
        std::ostringstream num;
        num << i;
        //cl_a.push_back(Pair(num.str(),cl));
        cl_a.push_back(cl);
    }
    Object cl_o;
    cl_o.push_back( Pair( "CL", cl_a ) );

    //regroupement des donnees pour ecriture du fichier de sortie
    Array output;
    output.push_back(mesh_o);
    output.push_back(group_elem_o);
    output.push_back(group_inter_o);
    output.push_back(prop_inter_o);
    output.push_back(cl_o);

    write_formatted( output, os );

    os.close();

}

//ecriture du fichier hdf5 -----------------------------------------------------------------------

void GeometryUser::write_visu_hdf5(MeshUser &mesh_user) {

    //ecriture des champs dans le fichier hdf5
//     std::string name_hdf =mesh_user.name_directory + "/MESH/visu_geometry.h5";

    if (FileExists(mesh_user.name_visu_hdf.c_str())) {
        String command = "rm -rf "; command << mesh_user.name_visu_hdf;
        int syst_rm=system(command.c_str());
    }

    Hdf hdf( mesh_user.name_visu_hdf.c_str() );
    String name;
    int num_level=0;
    name << "/Level_" << num_level << "/Geometry";

    //ecriture des noeuds
    nodes[0].write_to( hdf, name + "/nodes/x" );
    nodes[1].write_to( hdf, name + "/nodes/y" );
#if DIM==3
    nodes[2].write_to( hdf, name + "/nodes/z" );
#endif

    //ecriture des groupes d'elements
    String name_group_0;
    name_group_0 << name << "/elements_0";
    for (unsigned ng=0;ng<group_elements.size();ng++) {
        String name_list ;
        name_list<< name_group_0 << "/list_" << group_elements[ng].id ;
        for (unsigned nb_connect=0;nb_connect<group_elements[ng].connectivities.size();nb_connect++) {
            String name_connect;
            name_connect << name_list << "/c"<<nb_connect;
            group_elements[ng].connectivities[nb_connect].write_to( hdf, name_connect );
        }
        hdf.add_tag(name_list,"piece",name_list.c_str());
        //TODO utiliser plutôt le flag reperant le type de pattern
        String type_elements;       
//         String type_elements=(patterns.find_type(group_elements[ng].pattern_id)).name.c_str();       

#if DIM == 2
        type_elements="Triangle";
#else
        type_elements="Tetra";
#endif
        hdf.add_tag(name_list,"base",type_elements.c_str());
        hdf.write_tag(name_list,"pattern_id",group_elements[ng].pattern_id);
        hdf.write_tag(name_list,"id",group_elements[ng].id);
        hdf.write_tag(name_list,"nb_elements",group_elements[ng].nb_elements);
        //sauvegarde des sides pour chaque element (repere du groupe et du numéro local dans le groupe de l'interface voisine)
        for(unsigned nside=0;nside<group_elements[ng].interface_group_id.size();nside++){
            String name_side ;
            name_side<< name_group_0 << "/list_" << group_elements[ng].id << "/sides/side_" << nside ;
            String name_side_group_id;
            name_side_group_id << name_side << "/interface_group_id";
            group_elements[ng].interface_group_id[nside].write_to( hdf, name_side_group_id );
            String name_side_num_in_group;
            name_side_num_in_group << name_side << "/num_in_group" ;
            group_elements[ng].interface_num_in_group[nside].write_to( hdf, name_side_num_in_group );
        }
        //sauvegarde des id des groupes d'interfaces voisins
        BasicVec<String> name_type("/id_group_interfaces_edges","/id_group_interfaces_interiors","/id_group_interfaces_links");
        for(unsigned nside_g=0;nside_g<group_elements[ng].group_interfaces_id.size();nside_g++){
            if(group_elements[ng].group_interfaces_id[nside_g].size() != 0){
                String name_side_g ;
                name_side_g<< name_group_0 << "/list_" << group_elements[ng].id << name_type[nside_g] ;
                group_elements[ng].group_interfaces_id[nside_g].write_to( hdf, name_side_g );
            }
        }
    }


    //ecriture des groupes d'interfaces
    String name_group_1;
    name_group_1 << name << "/elements_1";

    
    PRINT(group_interfaces.size());
    for (unsigned ng=0;ng<group_interfaces.size();ng++) {
        String name_list ;
        name_list<< name_group_1 << "/list_" << group_interfaces[ng].id ;
        
        for (unsigned nb_connect=0;nb_connect<group_interfaces[ng].connectivities.size();nb_connect++) {
            String name_connect;
            name_connect << name_list << "/c"<<nb_connect;
            group_interfaces[ng].connectivities[nb_connect].write_to( hdf, name_connect );
        }
        if (group_interfaces[ng].type==0) {//bord
            //PRINT(name_list);
            String name_piece;
            name_piece << name_group_1 << name << "/elements_0/list_" << group_interfaces[ng].group_elements_id[0];
            hdf.add_tag(name_list,"edge_of_0",name_piece.c_str());
            hdf.add_tag(name_list,"edge_of_1","");
            hdf.write_tag(name_list,"id_edge_of_0",group_interfaces[ng].group_elements_id[0]);
            hdf.write_tag(name_list,"type",group_interfaces[ng].type);
        }
        else if (group_interfaces[ng].type==1 or group_interfaces[ng].type==2) {//interieur
            String name_piece_0, name_piece_1;
            name_piece_0 << name_group_1 << name << "/elements_0/list_" << group_interfaces[ng].group_elements_id[0];
            name_piece_1 << name_group_1 << name << "/elements_0/list_" << group_interfaces[ng].group_elements_id[1];
            hdf.add_tag(name_list,"edge_of_0",name_piece_0.c_str());
            hdf.add_tag(name_list,"edge_of_1",name_piece_1.c_str());
            hdf.write_tag(name_list,"id_edge_of_0",group_interfaces[ng].group_elements_id[0]);
            hdf.write_tag(name_list,"id_edge_of_1",group_interfaces[ng].group_elements_id[1]);            
            hdf.write_tag(name_list,"type",group_interfaces[ng].type);
        }
        //ecriture des id des groupes d'elements voisins ainsi que de l'id du pattern des groupes adjacents
        for(unsigned nside_g=0;nside_g<group_interfaces[ng].group_elements_id.size();nside_g++){
            String name_side ;
            name_side<< name_group_1 << "/list_" << group_interfaces[ng].id << "/sides/side_" << nside_g ;
            String name_side_elements_id;
            name_side_elements_id << name_side << "/element_num_in_group";
            group_interfaces[ng].element_num_in_group[nside_g].write_to( hdf, name_side_elements_id );
            String name_side_num_side;
            name_side_num_side << name_side << "/element_num_side";
            group_interfaces[ng].element_num_side[nside_g].write_to( hdf, name_side_num_side );
            String name_id="id_group_elements_";
            name_id << nside_g;
            hdf.add_tag(name_side,name_id.c_str(),group_interfaces[ng].group_elements_id[nside_g]);
            String name_pattern="id_pattern_group_elements_";
            name_pattern << nside_g;            
            hdf.add_tag(name_side,name_pattern.c_str(),group_interfaces[ng].patterns_id[nside_g]);
        } 
        String type_elements;
#if DIM == 2
        type_elements="Bar";
#else
        type_elements="Triangle";
#endif
        hdf.add_tag(name_list,"base",type_elements.c_str());
        hdf.write_tag(name_list,"interface_id",patterns.find_type(group_interfaces[ng].patterns_id[0]).interface_id);
        hdf.write_tag(name_list,"id",group_interfaces[ng].id);
        hdf.write_tag(name_list,"nb_interfaces",group_interfaces[ng].nb_interfaces);
    }
}



void GeometryUser::read_hdf5(String name_file_hdf5) {
    if (FileExists(name_file_hdf5.c_str())==0) {
        std::cerr << "Le fichier hdf5 "<< name_file_hdf5.c_str() << " n'existe pas " << std::endl;
        assert(0);
    }
    Hdf hdf( name_file_hdf5.c_str() );

    dim=DIM;
    
    String name;
    num_level=0;
    name << "/Level_" << num_level << "/Geometry";
    //lecture des noeuds
    nodes.resize(DIM);
    nodes[0].read_from( hdf, name + "/nodes/x" );
    nodes[1].read_from( hdf, name + "/nodes/y" );
#if DIM==3
    nodes[2].read_from( hdf, name + "/nodes/z" );
#endif
    
    //lecture des groupes d'elements
    PRINT("lecture des group_elements");
    String name_group_0; name_group_0 << name << "/elements_0";
    BasicVec<String> list_groups;
    list_groups=hdf.list_dir( name_group_0 );
    group_elements.resize(list_groups.size());
    nb_group_elements=group_elements.size();
    for (unsigned ng=0;ng<group_elements.size();ng++) { 
        String name_list ;
        name_list<< name_group_0 << "/"<< list_groups[ng];
        hdf.read_tag(name_list,"id",group_elements[ng].id,1);
        hdf.read_tag(name_list,"pattern_id",group_elements[ng].pattern_id,1);
        hdf.read_tag(name_list,"nb_elements",group_elements[ng].nb_elements,1);
        
        int nb_nodes = patterns.find_type(group_elements[ng].pattern_id).nb_nodes;
        int nb_sides = patterns.find_type(group_elements[ng].pattern_id).nb_sides;
        
        group_elements[ng].connectivities.resize(nb_nodes);
        
        for (unsigned i_connect=0;i_connect<group_elements[ng].connectivities.size();i_connect++) {
            String name_connect;
            name_connect << name_list << "/c"<<i_connect;
            group_elements[ng].connectivities[i_connect].read_from( hdf, name_connect );
        }

        //lecture des sides pour chaque element (repere du groupe et du numéro local dans le groupe de l'interface voisine)
        BasicVec<String> list_side;
        list_side=hdf.list_dir( name_list + "/sides");
        group_elements[ng].interface_group_id.resize(list_side.size());
        group_elements[ng].interface_num_in_group.resize(list_side.size());
        for(unsigned nside=0;nside<group_elements[ng].interface_group_id.size();nside++){
            String name_side_group_id;
            name_side_group_id << name_list << "/sides/" <<list_side[nside] << "/interface_group_id";
            group_elements[ng].interface_group_id[nside].read_from( hdf, name_side_group_id );
            String name_side_num_in_group;
            name_side_num_in_group << name_list << "/sides/" <<list_side[nside] << "/num_in_group" ;
            group_elements[ng].interface_num_in_group[nside].read_from( hdf, name_side_num_in_group );
        }
        //lecture des id des groupes d'interfaces voisins
        BasicVec<String> name_type("/id_group_interfaces_edges","/id_group_interfaces_interiors","/id_group_interfaces_links");
        group_elements[ng].group_interfaces_id.resize(3);
        for(unsigned nside_g=0;nside_g<group_elements[ng].group_interfaces_id.size();nside_g++){
            String name_side_g ;
            name_side_g<< name_group_0 << "/list_" << group_elements[ng].id << name_type[nside_g] ;
           if ( hdf.dataset_exist(name_side_g) ){
                group_elements[ng].group_interfaces_id[nside_g].read_from( hdf, name_side_g );
           }
        }
    }

    PRINT("lecture des group_interfaces");
    //lecture des groupes d'interfaces
    String name_group_1;
    name_group_1 << name << "/elements_1";
    list_groups=hdf.list_dir( name_group_1 );
    group_interfaces.resize(list_groups.size());
    nb_group_interfaces=group_interfaces.size();
    for (unsigned ng=0;ng<group_interfaces.size();ng++) {
        String name_list ;
        name_list<< name_group_1 << "/"<< list_groups[ng];
        hdf.read_tag(name_list,"id",group_interfaces[ng].id,1);
        hdf.read_tag(name_list,"type",group_interfaces[ng].type,1);
        hdf.read_tag(name_list,"interface_id",group_interfaces[ng].interface_id,1);
        hdf.read_tag(name_list,"nb_interfaces",group_interfaces[ng].nb_interfaces,1);
        
        int nb_nodes = interfaces.find_type(group_interfaces[ng].interface_id).nb_nodes;
        int nb_nodes_eq = interfaces.find_type(group_interfaces[ng].interface_id).nb_nodes_eq;
        
        group_interfaces[ng].connectivities.resize(nb_nodes);
        
        for (unsigned nb_connect=0;nb_connect<group_interfaces[ng].connectivities.size();nb_connect++) {
            String name_connect;
            name_connect << name_list << "/c"<<nb_connect;
            group_interfaces[ng].connectivities[nb_connect].read_from( hdf, name_connect );
        }

        if (group_interfaces[ng].type==0) {//bord
            group_interfaces[ng].group_elements_id.resize(1);
            hdf.read_tag(name_list,"id_edge_of_0",group_interfaces[ng].group_elements_id[0]);
        }
        else if (group_interfaces[ng].type==1 or group_interfaces[ng].type==2) {//interieur
            group_interfaces[ng].group_elements_id.resize(2);
            hdf.read_tag(name_list,"id_edge_of_0",group_interfaces[ng].group_elements_id[0]);
            hdf.read_tag(name_list,"id_edge_of_1",group_interfaces[ng].group_elements_id[1]);
        }

        //lecture des id des groupes d'elements voisins ainsi que de l'id du pattern des groupes adjacents
        //lecture des sides pour chaque element (repere du groupe et du numéro local dans le groupe de l'interface voisine)
        
        group_interfaces[ng].element_num_in_group.resize(group_interfaces[ng].group_elements_id.size());
        group_interfaces[ng].element_num_side.resize(group_interfaces[ng].group_elements_id.size());
        group_interfaces[ng].patterns_id.resize(group_interfaces[ng].group_elements_id.size());
        nb_group_elements=group_interfaces[ng].group_elements_id.size();
        for(unsigned nside_g=0;nside_g<group_interfaces[ng].group_elements_id.size();nside_g++){
            String name_side ;
            name_side<< name_group_1 << "/list_" << group_interfaces[ng].id << "/sides/side_" << nside_g ;
            String name_side_elements_id;
            name_side_elements_id << name_side << "/element_num_in_group";
            group_interfaces[ng].element_num_in_group[nside_g].read_from( hdf, name_side_elements_id );
            String name_side_num_side;
            name_side_num_side << name_side << "/element_num_side";            
            group_interfaces[ng].element_num_side[nside_g].read_from( hdf, name_side_num_side );
//             String name_pattern="id_pattern_group_elements_";
//             name_pattern << nside_g;            
//             hdf.read_tag(name_side,name_pattern.c_str(),group_interfaces[ng].patterns_id[nside_g]);
        }   
    }
    
    
}

#include "write_xdmf.h"

/** * Ecriture des DataItems : coordonnées des noeuds, connectivités classées par groupes, attributs
* Ecriture des groupes pour les éléemnts parents : Geometrie_0 : pieces_i
* Ecriture des groupes pour les éléments enfants en séparant les groupes de bord (Edges) et les liaisons entre pièces (Liaisons) : Geometrie_1 : Edges/Liaisons : Edge_i/Interface_i)
L'entier skin (0 ou 1) permet de sortir les groupes d'élements parents et enfants (skin=0) ou uniquement les éléments enfants (skin=1) **/
void GeometryUser::write_xdmf(String output_xdmf, String input_hdf5, String name_geometry, int skin){
    std::ofstream f(output_xdmf.c_str());
    f << "<Xdmf Version=\"2.0\">" << std::endl;
    f << "  <Domain>" << std::endl;
    
    String  name_nodes = name_geometry + "/nodes";
    Hdf hdf(input_hdf5.c_str());
    int nb_nodes;
    String data_node_x = name_nodes + "/x";
    hdf.read_size(data_node_x,nb_nodes);
    
    //ecriture des dataitems noeuds + attributs
    write_nodes_datasets(f,input_hdf5, name_nodes,nb_nodes);
    
    //ecriture des dataitems list d'elements connectivites + champs par element
    for (unsigned el=skin;el<2;el++) {
        int nb_list;
        String name_group_elements; name_group_elements<<name_geometry<<"/elements_"<<el;
        hdf.read_group_size(name_group_elements,nb_list);
        write_groups_datasets(f,input_hdf5, name_group_elements,nb_list,el);
    }
    
    
    //ecriture de la geometrie_0 (elements parents)
    BasicVec<String> attributs=("none");
    String name_fields="none";
    if (skin==0) {
        f<< "           <Grid Name=\"Geometry_0\" GridType=\"Tree\">" << std::endl;
        int nb_list;
        String name_elements_0;
        name_elements_0<<name_geometry<<"/elements_0";
        hdf.read_group_size(name_elements_0,nb_list);
        String generic_grid_name="piece";
        BasicVec<int> list;
        list.resize(nb_list);
        for(unsigned i=0;i<nb_list;i++) list[i]=i;
        write_grids(f,input_hdf5, name_elements_0,list,nb_nodes,generic_grid_name,name_fields,attributs,0,0,0);
        f<<"            </Grid>"<< std::endl;
    }
       
    //ecriture de la geometrie_1 (elements enfants)
    f<< "           <Grid Name=\"Geometry_1\" GridType=\"Tree\">" << std::endl;
    int nb_list;
    String name_elements_1;
    name_elements_1<<name_geometry<<"/elements_1";
    hdf.read_group_size(name_elements_1,nb_list);

     
    int nb_links=0, nb_edges=0;
    BasicVec<BasicVec<int>,3 > lists ;
    for (unsigned i=0; i< nb_list; i++) {
        String name_list;
        name_list << name_elements_1<<  "/list_" <<  i;
        int type=-1;
        hdf.read_tag(name_list,"type",type);
        if (type==0) {//bord
            lists[0].push_back(i);
        }
        else if (type==2) {//liaison entre 2 pièces
            lists[2].push_back(i);
        }
    }
 
    for (unsigned tl=0;tl<lists.size();tl++) {
        String generic_grid_name;
        if (tl==0) {//bord
            f<<"                    <Grid Name=\" Edges \" GridType=\"Tree\">" << std::endl;
            generic_grid_name="Edge";
            write_grids(f,input_hdf5, name_elements_1,lists[tl],nb_nodes,generic_grid_name,name_fields,attributs,0,0,1);
        }
        else if (tl==1) {
            continue;
        }
        else if (tl==2) {//liaison entre 2 pièces
            f<<"                    <Grid Name=\" Liaisons \" GridType=\"Tree\" >" << std::endl;
            generic_grid_name="Interface";
            write_grids(f,input_hdf5, name_elements_1,lists[tl],nb_nodes,generic_grid_name,name_fields,attributs,0,0,1);
        }
        
        f<<"                 </Grid>"<< std::endl;
    }
    f<<"            </Grid>"<<std::endl;
    
    f<<"    </Domain>"<< std::endl;
    f<<"</Xdmf>"<<std::endl;

    f.close();
}

