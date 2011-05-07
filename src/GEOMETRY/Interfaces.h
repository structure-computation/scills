//
// Description: fichier généré par GEOM_Interfaces.met 
//
//
// Author: Jeremie Bellec <j.bellec@structure-computation.com>, (C) 2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef INTERFACES_H
#define INTERFACES_H

#include <Metil/BasicVec.h>

using namespace Metil;

struct Interfaces{
    //atributs **********************************************************************************************
    struct Types{
        std::string name ;
        int id ;
        int dim ;
        int nb_nodes ;
        int nb_nodes_eq ; 
        int nb_vertex_nodes ; 
        
        template<class TB,class TP>
        void apply_bs( TB &res, TP ) const {
            res.set_type( "Types" );
            APPLY_WN( res, id );
            APPLY_WN( res, dim );
            APPLY_WN( res, nb_nodes );
            APPLY_WN( res, nb_nodes_eq );
            APPLY_WN( res, nb_vertex_nodes );
        }
    };
    BasicVec< Types > types;
    
    
    //méthodes **********************************************************************************************
    Types find_type(int id_){
        for(int num_type=0; num_type<4; num_type++){
            if(types[num_type].id == id_){
                return types[num_type];
                break;
            }
        }
    }
    
    int find_index_type(int id_){
        for(int num_type=0; num_type<4; num_type++){
            if(types[num_type].id == id_){
                return num_type ;
                break;
            }
        }
    }
    
    template<class TB,class TP>
    void apply_bs( TB &res, TP ) const {
        res.set_type( "Interfaces" );
        APPLY_WN( res, types );
    }
          
    Interfaces(){
        types.resize(4);

        //
        // type 0 = Bar --------------------------------------
        types[0].name = "Bar";
        types[0].id = 0;
        types[0].dim = 2;
        types[0].nb_nodes = 2;
        types[0].nb_nodes_eq = 2;
        types[0].nb_vertex_nodes = 2;

        //
        // type 0 = Bar_3 --------------------------------------
        types[0].name = "Bar_3";
        types[0].id = 1;
        types[0].dim = 2;
        types[0].nb_nodes = 3;
        types[0].nb_nodes_eq = 2;
        types[0].nb_vertex_nodes = 2;

        //
        // type 0 = Triangle --------------------------------------
        types[0].name = "Triangle";
        types[0].id = 2;
        types[0].dim = 3;
        types[0].nb_nodes = 3;
        types[0].nb_nodes_eq = 4;
        types[0].nb_vertex_nodes = 3;

        //
        // type 0 = Triangle_6 --------------------------------------
        types[0].name = "Triangle_6";
        types[0].id = 3;
        types[0].dim = 3;
        types[0].nb_nodes = 6;
        types[0].nb_nodes_eq = 4;
        types[0].nb_vertex_nodes = 3;

    }
};
#endif //INTERFACES_H
