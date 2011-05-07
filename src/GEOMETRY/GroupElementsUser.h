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

#include <cmath>
#include <Metil/BasicVec.h>
#include <Metil/StructCompactor.h>

#include "Patterns.h"
// #include <boost/concept_check.hpp>

using namespace Metil;

#ifndef GROUP_ELEMENTS_USER_H
#define GROUP_ELEMENTS_USER_H


//group elements---------------------------------------------------------------------------------------------------------------------------
struct GroupElementsUser{
    //flag du group
    BasicVec< int > flags;
    BasicVec< String > flags_names;
  
    // attributs du groupe------------------------------------------
    int id;
    int nb_elements;
    int nb_elements_for_GPU;
    int pattern_id;
    
    // attributs des elements du groupe
    BasicVec< BasicVec< int > > connectivities;                        //liste des connectivité des elements (nb_node_by_pattern, nb_elements)  
    BasicVec< BasicVec< int > > interface_group_id;                    //numero local du groupe des interfaces adjactentes à chaque motif rangés dans l'ordre des sides (nb_side_by_pattern, nb_elements)
    BasicVec< BasicVec< int > > interface_num_in_group;                //numero local des interfaces adjactentes à chaque motif dans leur groupe rangés dans l'ordre des sides (nb_side_by_pattern, nb_elements)
    
    // attributs des groupes d'interfaces adjacents au groupe
    BasicVec< BasicVec< int > > group_interfaces_id;                   //identité des groupes d'interfaces adjactents taille 3 * nb_group_interfaces adjacents en fonction du type : 0 : bord | 1 : interieur au group_elements  |  2 : entre deux group_elements
    
    // operateur geometrique des elements 
    BasicVec< BasicVec< BasicVec< BasicVec< TYPE > > > > side_N;       // operateur de passage entre W et U(nb_node_side,nb_node_eq_side,nb_side,nb_elements)  
    BasicVec< BasicVec< BasicVec< TYPE > > > side_M;       // operateur pour l'intégration sur un side entre W et F (nb_node_eq_side,nb_side,nb_elements)  
    
    // méthodes du groupe------------------------------------------
    
    //Méthode permettant de spécifier les champs à inclure dans la structure compactée dont le nom est donné par set_type. Cette méthode est appelée lors de la génération des classes compactées.
    template<class TB,class TP>
    void apply_bs( TB &res, TP ) const {
        res.set_type( "GroupElementsCompactClass" );
        APPLY_WN( res, id    );
//         APPLY_WN( res, flags );
//         APPLY_WN( res, flags_names );
        APPLY_WN( res, nb_elements  );
        APPLY_WN( res, pattern_id  );
        APPLY_WN( res, connectivities    );
        APPLY_WN( res, interface_group_id    );
        APPLY_WN( res, interface_num_in_group    );
        APPLY_WN( res, group_interfaces_id    );
        APPLY_WN( res, side_N    );
        APPLY_WN( res, side_M    );
    }
    
    //Constructeur
    GroupElementsUser(int id_, EntityElementUser &entity_element){
        //std::cout << "ajout d'un group_elements n° " << id_ << std::endl;
        id = id_;
        nb_elements = 0;
        flags = entity_element.flags;
        flags_names = entity_element.flags_names;
        for( int i_flag=0; i_flag<flags_names.size(); i_flag++){
            if(flags_names[i_flag] == "pattern_id"){
                pattern_id = flags[i_flag];
            }
        }
        group_interfaces_id.resize(3);
        connectivities.resize(entity_element.element_connectivity.size());
        interface_group_id.resize(entity_element.nb_sides);
        interface_num_in_group.resize(entity_element.nb_sides);
        add_entity_element(entity_element);
    }
    //Constructeur simplifié
    GroupElementsUser(){
        group_interfaces_id.resize(3);
    }
    // Ajout d'un entity_element au group
    bool add_entity_element(EntityElementUser &entity_element){
        bool entity_element_added = false;
        if(entity_element.valid_flags(flags, flags_names)){
            for(int i_connect=0; i_connect<connectivities.size(); i_connect++){
                connectivities[i_connect].push_back(entity_element.element_connectivity[i_connect]);
            }
            for(int i_side=0; i_side<entity_element.nb_sides; i_side++){
                interface_group_id[i_side].push_back(-1);
                interface_num_in_group[i_side].push_back(-1);
            }
            entity_element.group_id = id;
            entity_element.num_in_group = nb_elements;  
            nb_elements += 1;
            entity_element_added = true;
        }
        return entity_element_added;
    }
    // Définition des tailles des listes pour utilisation sous GPU
    void initialize_GPU(Patterns &patterns){
        int threadsPerBlock = 48;
        int n_group_blocks = std::ceil( nb_elements / threadsPerBlock );  // obtenir l'entier superieur
        nb_elements_for_GPU = (n_group_blocks + 1) * threadsPerBlock;
        
        for(int i_connect=0; i_connect<connectivities.size(); i_connect++){
            connectivities[i_connect].resize(nb_elements_for_GPU);
        }
        for(int i_side=0; i_side<interface_group_id.size(); i_side++){
            interface_group_id[i_side].resize(nb_elements_for_GPU);
            interface_num_in_group[i_side].resize(nb_elements_for_GPU);
        }
        //(nb_node_side,nb_node_eq_side,nb_side,num_motif) 
        int nb_sides = patterns.find_type(pattern_id).nb_sides;
        int nb_nodes_by_sides = patterns.find_type(pattern_id).nb_nodes_by_sides;
        int nb_nodes_eq_by_sides = patterns.find_type(pattern_id).nb_nodes_eq_by_sides;
        
        side_N.resize(nb_nodes_by_sides);
        for(int i_node=0; i_node < nb_nodes_by_sides; i_node++){
            side_N[i_node].resize(nb_nodes_eq_by_sides);
            for(int i_node_eq=0; i_node_eq < nb_nodes_eq_by_sides; i_node_eq++){
                side_N[i_node][i_node_eq].resize(nb_sides);
                for(int i_side=0; i_side < nb_sides; i_side++){
                    side_N[i_node][i_node_eq][i_side].resize(nb_elements_for_GPU,1);
                }
            }
        }
        
        side_M.resize(nb_nodes_eq_by_sides);
        for(int i_node_eq=0; i_node_eq < nb_nodes_eq_by_sides; i_node_eq++){
            side_M[i_node_eq].resize(nb_sides);
            for(int i_side=0; i_side < nb_sides; i_side++){
                side_M[i_node_eq][i_side].resize(nb_elements_for_GPU,1);
            }
        }        
        
    }
    // Affichage pour vérification
    void affiche(){
        PRINT("group_elements-----------------------------------------------------------");
        PRINT(id);
        PRINT(pattern_id);
        PRINT(nb_elements);
        PRINT(interface_group_id);
        PRINT(interface_num_in_group);
        PRINT(connectivities);
    }
};


#endif //GROUP_ELEMENTS_USER_H


