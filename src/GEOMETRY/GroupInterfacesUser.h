//
// C++ Interface: Définition des group_elements et group_interfaces
//
// Description: 
//
//
// Author: Jeremie Bellec <j.bellec@structure-computation.com>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//


#include <Metil/BasicVec.h>
#include <Metil/StructCompactor.h>

#include "Patterns.h"
// #include <boost/concept_check.hpp>
//#include "containers/mat.h"

using namespace Metil;
//using namespace LMT;

#ifndef GROUP_INTERFACES_USER_H
#define GROUP_INTERFACES_USER_H
              
//group_interfaces---------------------------------------------------------------------------------------------------------------------------
struct GroupInterfacesUser{ 
    // attributs du groupe------------------------------------------
    int id;                                                          //identité du group
    int interface_id;                                                //id du type d'interface
    int edge_id;                                                     //id de l'edge géometrique définit par l'utilisateur pour créer ce nouveau group
    int type;                                                        //type du group 0 : bord | 1 : interieur a un group_elements  |  2 : entre deux group_elements
    int nb_group_elements;
    int nb_interfaces;
    int nb_interfaces_for_GPU;
    int from_group_id;                                              // de quel group vient il (si c'est un split de group)
    
    BasicVec< int > group_elements_id;                            //identité des groupes d'elements adjacents (max 2)
    BasicVec< int > patterns_id;                                  //identité patterns pour chaque groupes d'elements adjacents (max 2)
    BasicVec< BasicVec< int > > element_num_in_group;             //numéros locals des elements adjacents à chaque interface dans leur group_elements (max 2, nb_interfaces) !meme ordre que id_group_elements
    BasicVec< BasicVec< int > > element_num_side;                 //numéros des sides des elements adjacents à chaque interface dans leur group_elements (max 2, nb_interfaces)
    
    BasicVec< BasicVec< int > > connectivities;                   //liste des connectivité des interfaces (nb_node_by_side_pattern, nb_interfaces)
    BasicVec< BasicVec< int > > correspondance_between_nodes;      //correspondance des noeuds du side 2 par rapport au side 1 (nb_node_by_side_pattern, nb_interfaces)
    
    BasicVec< BasicVec< TYPE > > Ne;                                    //operateur de passage entre W (depl sur interface) et U(déplacement sur le bord) (nb_node_side,nb_interfaces)
    BasicVec< TYPE > Me;                                                //mesure de chaque interface, utilise pour faire l'integration FW sur un element (nb_interfaces) 
    TYPE measure;               //mesure du groupe d'interface (surface ou longueur)
//     

    // méthodes du groupe------------------------------------------
    template<class TB,class TP>
    void apply_bs( TB &res, TP ) const {
        res.set_type( "GroupInterfacesCompactClass" );
        APPLY_WN( res, id    );
        APPLY_WN( res, nb_group_elements );
        APPLY_WN( res, nb_interfaces  );
        APPLY_WN( res, group_elements_id    );
        APPLY_WN( res, patterns_id    );
        APPLY_WN( res, element_num_in_group    );
        APPLY_WN( res, element_num_side    );
        APPLY_WN( res, connectivities    );  
        APPLY_WN( res, correspondance_between_nodes    );
        APPLY_WN( res, Ne    );
        APPLY_WN( res, Me    );
        APPLY_WN( res, measure    );
    }
    //création à partir du maillage
    GroupInterfacesUser(int id_, EntityInterfaceUser &entity_interface, BasicVec< EntityElementUser > &list_elements, BasicVec<GroupElementsUser> &group_elements, Patterns &patterns){
        //std::cout << "ajout d'un group_interfaces n° " << id_ << std::endl;
        id = id_;
        from_group_id = -1;
        edge_id = -2;
        nb_interfaces = 0;
        nb_group_elements = entity_interface.nb_sides;
        group_elements_id.resize(nb_group_elements);
        patterns_id.resize(nb_group_elements);
        element_num_in_group.resize(nb_group_elements);
        element_num_side.resize(nb_group_elements);
        connectivities.resize(entity_interface.correspondance_between_nodes.size());
        correspondance_between_nodes.resize(entity_interface.correspondance_between_nodes.size());
        
        for(int i_group_elements=0; i_group_elements<nb_group_elements; i_group_elements++){
            //list_elements[entity_interface.elements_global_num[i_group_elements]].affiche();
            group_elements_id[i_group_elements] = list_elements[entity_interface.elements_global_num[i_group_elements]].group_id; 
            patterns_id[i_group_elements] = list_elements[entity_interface.elements_global_num[i_group_elements]].pattern_id; 
        }
        type = entity_interface.mesh_type;
        if(group_elements_id.size() == 2 && group_elements_id[0] != group_elements_id[1]){
            type = 2;
        }
        
        //ajout de l'id du group_interfaces a la liste des group_elements adjacents
        for(int i_group_elements=0; i_group_elements<nb_group_elements; i_group_elements++){
            group_elements[group_elements_id[i_group_elements]].group_interfaces_id[type].push_back(id) ;
        }
        //ajout de l'interface
        
        add_entity_interface(entity_interface, list_elements, group_elements, patterns);    
    }
    //création à partir du split d'un autre group
    GroupInterfacesUser( GroupInterfacesUser &base_group_interfaces, BasicVec< int > edge_assigned, int edge_id_, int id_){
        //std::cout << "ajout d'un group_interfaces n° " << id_ << std::endl;
        id = id_;
        from_group_id = base_group_interfaces.id;
        interface_id = base_group_interfaces.interface_id;
        nb_group_elements = base_group_interfaces.nb_group_elements;
        type = base_group_interfaces.type;
        edge_id = edge_id_;
        group_elements_id = base_group_interfaces.group_elements_id;
        patterns_id = base_group_interfaces.patterns_id;
        
        nb_interfaces = 0;
        BasicVec< int > copy_interface_from_base;
        for(int i_inter=0; i_inter<base_group_interfaces.nb_interfaces; i_inter++){
            if(edge_assigned[i_inter]==edge_id_){
                nb_interfaces+=1;
                copy_interface_from_base.push_back(i_inter);
            }
        }
        
        //initialisation des taille
        group_elements_id.resize(base_group_interfaces.group_elements_id.size());
        patterns_id.resize(base_group_interfaces.patterns_id.size());
        element_num_in_group.resize(base_group_interfaces.element_num_in_group.size());
        element_num_side.resize(base_group_interfaces.element_num_side.size());
        connectivities.resize(base_group_interfaces.connectivities.size());
        correspondance_between_nodes.resize(base_group_interfaces.correspondance_between_nodes.size());
        
        for(int i_connect=0; i_connect<connectivities.size(); i_connect++){
            connectivities[i_connect].resize(nb_interfaces);
        }
        for(int i=0; i<correspondance_between_nodes.size(); i++){
            correspondance_between_nodes[i].resize(nb_interfaces);
        }
        for(int i=0; i<element_num_in_group.size(); i++){
            element_num_in_group[i].resize(nb_interfaces);
            element_num_side[i].resize(nb_interfaces);
        }
        
        //assignation des interface a partir du groupe de base_group_interfaces
        for(int i_connect=0; i_connect<connectivities.size(); i_connect++){
            for(int i_inter=0; i_inter<nb_interfaces; i_inter++){
                connectivities[i_connect][i_inter]=base_group_interfaces.connectivities[i_connect][copy_interface_from_base[i_inter]];
            }
        }
        for(int i=0; i<correspondance_between_nodes.size(); i++){
            for(int i_inter=0; i_inter<nb_interfaces; i_inter++){
                correspondance_between_nodes[i][i_inter]=base_group_interfaces.correspondance_between_nodes[i][copy_interface_from_base[i_inter]];
            }  
        }
        for(int i=0; i<element_num_in_group.size(); i++){
            for(int i_inter=0; i_inter<nb_interfaces; i_inter++){
                element_num_in_group[i][i_inter]=base_group_interfaces.element_num_in_group[i][copy_interface_from_base[i_inter]];
                element_num_side[i][i_inter]=base_group_interfaces.element_num_side[i][copy_interface_from_base[i_inter]];
            } 
        }      
    }
    GroupInterfacesUser(){ 
        edge_id = -2;
    } 
    // ajout d'une entity_interface au group
    bool add_entity_interface(EntityInterfaceUser &entity_interface, BasicVec< EntityElementUser > &list_elements, BasicVec<GroupElementsUser> &group_elements, Patterns &patterns){
        bool entity_interface_added = false;
        if(entity_interface.valid_group_elements( list_elements, group_elements_id, patterns )){
            
            for(int i_connect=0; i_connect<correspondance_between_nodes.size(); i_connect++){
                correspondance_between_nodes[i_connect].push_back(entity_interface.correspondance_between_nodes[i_connect]);
                connectivities[i_connect].push_back(entity_interface.interface_connectivity[i_connect]);
            }
            for(int i_group_elements=0; i_group_elements<group_elements_id.size(); i_group_elements++){
                int id_group_elements = group_elements_id[i_group_elements];
                int num_side = entity_interface.elements_num_side[i_group_elements];
                int num_in_group = entity_interface.elements_num_in_group[i_group_elements];
                element_num_side[i_group_elements].push_back(num_side);
                element_num_in_group[i_group_elements].push_back(num_in_group);
                group_elements[id_group_elements].interface_group_id[num_side][num_in_group] = id;
                group_elements[id_group_elements].interface_num_in_group[num_side][num_in_group] = nb_interfaces;
            }
            entity_interface.group_id = id;
            entity_interface.num_in_group = nb_interfaces;  
            nb_interfaces += 1;
            //entity_interface.affiche();
            entity_interface_added = true;
        }
        return entity_interface_added;
    }
    // resize des liste pour GPU
    void initialize_GPU(Patterns &patterns){
        int threadsPerBlock = 48;
        int n_group_blocks = std::ceil( nb_interfaces / threadsPerBlock );  // obtenir l'entier superieur
        nb_interfaces_for_GPU = (n_group_blocks + 1) * threadsPerBlock;
        
        for(int i_connect=0; i_connect<connectivities.size(); i_connect++){
            connectivities[i_connect].resize(nb_interfaces_for_GPU);
        }
        for(int i=0; i<correspondance_between_nodes.size(); i++){
            correspondance_between_nodes[i].resize(nb_interfaces_for_GPU);
        }
        for(int i=0; i<element_num_in_group.size(); i++){
            element_num_in_group[i].resize(nb_interfaces_for_GPU);
            element_num_side[i].resize(nb_interfaces_for_GPU);
        }

        int nb_nodes_by_sides = patterns.find_type(patterns_id[0]).nb_nodes_by_sides;
        Ne.resize(nb_nodes_by_sides);
        for(int i_node=0; i_node < nb_nodes_by_sides; i_node++){
            Ne[i_node].resize(nb_interfaces_for_GPU,1);
        }
        Me.resize(nb_interfaces_for_GPU,1);
    }
    // affichage pour vérification
    void affiche(){
        PRINT("group_interfaces-----------------------------------------------------------");
        PRINT(id);
        PRINT(type);
        PRINT(nb_interfaces);
        PRINT(group_elements_id);
        PRINT(element_num_in_group);
        PRINT(element_num_side);
        PRINT(connectivities);
    }
};

#endif //GROUP_INTERFACES_USER_H


