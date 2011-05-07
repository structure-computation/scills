//
// C++ Interface: DataUser
//
// Description: 
//
//
// Author: Jeremie Bellec <j.bellec@structure-computation.com>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef DATA_USER_H
#define DATA_USER_H


#include <iostream>
#include <fstream>

#include <Metil/Level1/CompilationEnvironment.h>
#include <Metil/DynamicLibrary.h>
#include <Metil/BasicVec.h>

#include "Patterns.h"
#include "Properties.h"
#include "json_spirit.h"

using namespace json_spirit;  
using namespace Metil;

class DataUser{
    //******************************************************************************************
    //Attributs
    //******************************************************************************************
    public :
    DataUser(std::string model_path_, std::string calcul_path_, const char *id_calcul_){
        id_calcul << id_calcul_;
        model_path = model_path_;
        calcul_path = calcul_path_;
        PRINT(id_calcul);
    }
    DataUser(){
    }
    //nom et chemins de fichier--------------------------------
    std::string metil_comp_path;                //chemin vers metil_comp pour la compilation à la volée
    std::string model_path; 		        //chemin d'acces au repertoire model
    std::string name_mesh_user; 		//nom du fichier mesh utilisateur
    std::string extension; 		        //extention pour le nom du fichier de mesh
    std::string name_directory; 		//nom du repertoire pour sauvegarder le motif
    std::string mesh_directory; 		//nom du repertoire pour sauvegarder les meshs  
    std::string calcul_path;                    //chemin d'acces au repertoire model
    String id_calcul;                           // id du calcul en cours
    Patterns patterns;
    Properties properties;  
    
    // attribut des groupes ----------------------------------------------------------------
    // DATA pour les groupes d'elements--------------------
    struct GroupElements{
        int id;
        std::string name;
        int num_in_mesh_file;
        std::string origine; //origine du groupe d'element
        int id_material;
        
        template<class TB,class TP>
        void apply_bs( TB &res, TP ) const {
            res.set_type( "GroupElements" );
            //pour SC_create-------------------------
            //info du groupe
            APPLY_WN( res, id    );
//             APPLY_WN( res, name    );
            APPLY_WN( res, num_in_mesh_file    );
            APPLY_WN( res, id_material    );
        }
    };
    BasicVec<GroupElements> group_elements;
    
    // DATA pour les groupes d'interfaces------------------
    struct GroupInterfaces{
        int id;
        std::string name;
        int num_in_mesh_file;
        BasicVec<GroupInterfaces *,2> adj_group_elements;
        BasicVec<int> adj_num_group;
        std::string origine; //origine du groupe d'interface
        std::string type;
        int id_link;
        float assigned;
        
        GroupInterfaces(){
            adj_num_group.resize(2);
        }
        
        template<class TB,class TP>
        void apply_bs( TB &res, TP ) const {
            res.set_type( "GeometryEdges" );
            //pour SC_create-------------------------
            //info du groupe
            APPLY_WN( res, id    );
//             APPLY_WN( res, name    );
            APPLY_WN( res, num_in_mesh_file    );
            APPLY_WN( res, adj_num_group    );
            APPLY_WN( res, id_link    );
        }
    };
    BasicVec<GroupInterfaces> group_interfaces;
    
    // DATA pour les groupes de bords --------------------
    struct Geometry{
        std::string origine;
        std::string type;
        std::string nature;
        double radius;
        BasicVec<BasicVec<double, DIM> > points;
        BasicVec<double, DIM> pdirection;
        std::string equation;
        
        Geometry(){
            points.resize(DIM);
            pdirection.resize(DIM);
        }
        
        template<class TB,class TP>
        void apply_bs( TB &res, TP ) const {
            res.set_type( "GeometryEdges" );
            //pour SC_create-------------------------
            //info du groupe
//             APPLY_WN( res, origine    );
//             APPLY_WN( res, type    );
//             APPLY_WN( res, nature    );
            APPLY_WN( res, radius    );
            APPLY_WN( res, points    );
            APPLY_WN( res, pdirection    );
//             APPLY_WN( res, equation    );
        }
        
    };
    struct GroupEdges{
        int id;
        std::string name;
        int num_in_mesh_file;
        BasicVec<GroupElements *> adj_group_elements;
        BasicVec<int> adj_num_group;
        Geometry geom;
        int id_CL;
        float assigned;
        
        GroupEdges(){
            adj_num_group.resize(1);
        }
        
        template<class TB,class TP>
        void apply_bs( TB &res, TP ) const {
            res.set_type( "GroupEdges" );
            //pour SC_create-------------------------
            //info du groupe
            APPLY_WN( res, id    );
//             APPLY_WN( res, name    );
            APPLY_WN( res, num_in_mesh_file    );
            APPLY_WN( res, adj_num_group    );
            APPLY_WN( res, id_CL    );
            APPLY_WN( res, assigned    );
            APPLY_WN( res, geom    );
        }
    };
    BasicVec<GroupEdges> group_edges;
    
    
    // DATA pour les comportement materiaux --------------------
    struct BehaviourMaterial {
        BehaviourMaterial(){
            mat_prop.resize(properties.mat_prop_name.size(),"0");
            mat_prop_name = properties.mat_prop_name;
        }
	std::string name;
	std::string familly;
	std::string name_select;
	std::string description;
	std::string type;
	std::string comp;
        std::string resolution;
	int    type_num;
	int    id;
    
        Properties properties;
        BasicVec< std::string > mat_prop;
        BasicVec< std::string > mat_prop_name;
        
        template<class TB,class TP>
        void apply_bs( TB &res, TP ) const {
            res.set_type( "BehaviourMaterial" );
            //pour SC_create-------------------------
            //info du groupe
            APPLY_WN( res, id    );
//             APPLY_WN( res, name    );
//             APPLY_WN( res, familly    );
//             APPLY_WN( res, type    );
//             APPLY_WN( res, comp    );
            APPLY_WN( res, type_num    );
//             APPLY_WN( res, resolution    );
        }
	
    };
    BasicVec<BehaviourMaterial > behaviour_materials;
    BasicVec< BasicVec< int > > num_materials_id_group_elements;
    
    
    // DATA pour les comportement liaison --------------------
    struct BehaviourLink {
        BehaviourLink(){
            link_prop.resize(properties.link_prop_name.size(),"0");
            link_prop_name = properties.link_prop_name;
        }
        int id;
        std::string  name;
        std::string  familly;
        int  company_id;
        int  reference;
        int  id_select;
        std::string  name_select;
        std::string  type;
        std::string  comp_complexe;
        int  type_num;

        Properties properties;
        BasicVec< std::string > link_prop;
        BasicVec< std::string > link_prop_name;
        
        template<class TB,class TP>
        void apply_bs( TB &res, TP ) const {
            res.set_type( "BehaviourLink" );
            //pour SC_create-------------------------
            //info du groupe
            APPLY_WN( res, id    );
//             APPLY_WN( res, name    );
//             APPLY_WN( res, familly    );
            APPLY_WN( res, id_select    );
//             APPLY_WN( res, type    );
//             APPLY_WN( res, comp_complexe    );
            APPLY_WN( res, type_num    );
        }
        
    };
    BasicVec<BehaviourLink > behaviour_links;
    BasicVec< BasicVec< int > > num_links_id_group_interfaces;
    
    
    // DATA pour les condition limite --------------------
    struct StepBc{
        StepBc(){
            CL_step_prop.resize(properties.BC_step_prop_name.size(),"0");
            CL_step_prop_name = properties.BC_step_prop_name;
        }
        
        Properties properties;
        BasicVec< std::string > CL_step_prop;
        BasicVec< std::string > CL_step_prop_name;   
    };
    struct BehavBc {
        int id;
        std::string type;
        BasicVec<StepBc> step;
        
        template<class TB,class TP>
        void apply_bs( TB &res, TP ) const {
            res.set_type( "BehaviourBc" );
            //pour SC_create-------------------------
            //info du groupe
            APPLY_WN( res, id    );
//             APPLY_WN( res, type    );
        }
    }; 
    BasicVec<BehavBc > behaviour_bc;
   
    
    // DATA pour les condition limite en volume-----------
    struct StepBcVolume{
        StepBcVolume(){
            CLv_step_prop.resize(properties.BCv_step_prop_name.size(),"0");
            CLv_step_prop_name = properties.BCv_step_prop_name;  
        }
        
        Properties properties;
        BasicVec< std::string > CLv_step_prop;
        BasicVec< std::string > CLv_step_prop_name;   
        
    };
    struct BehavBcVolume{
        std::string name;
        BasicVec<StepBcVolume> step;
        std::string type;
        int select;
        int ref;
        
        template<class TB,class TP>
        void apply_bs( TB &res, TP ) const {
            res.set_type( "BehaviourBcVolume" );
            //pour SC_create-------------------------
            //info du groupe
//             APPLY_WN( res, name    );
//             APPLY_WN( res, type    );
            APPLY_WN( res, select    );
            APPLY_WN( res, ref    );
        }
    };
    BasicVec<BehavBcVolume> behaviour_bc_volume;
    
    
    // DATA pour les options de calcul-----------------------
    struct OptionsCalcul {
        template<class TB,class TP>
        void apply_bs( TB &res, TP ) const {
            res.set_type( "OptionsCalcul" );
            //pour SC_create-------------------------
            //info du groupe
            APPLY_WN( res, dissipation    );
            APPLY_WN( res, nb_option    );
            APPLY_WN( res, nb_level    );
//             APPLY_WN( res, Temp_statique   );
//             APPLY_WN( res, mode    );
//             APPLY_WN( res, resolution_2D    );
            APPLY_WN( res, LATIN_nb_iter_max    );
            APPLY_WN( res, LATIN_crit_error    );
            APPLY_WN( res, PREC_erreur    );
            APPLY_WN( res, multiechelle    );
            APPLY_WN( res, LATIN_current_iter    );
            APPLY_WN( res, save_depl    );
        }
        OptionsCalcul(){
            Multiresolution_on = false;
            Multiresolution_nb_cyle = 1;
        }
        std::string mode;
        std::string resolution_2D;
        int nb_option;
        std::string Temp_statique;
        int LATIN_nb_iter_max;
        TYPE LATIN_crit_error;
        int PREC_erreur;
        int nb_level;
        int multiechelle;
        int save_depl;        
        int Multiresolution_on;
        int Multiresolution_nb_cyle;
        int dissipation;
        int LATIN_current_iter; 
    };
    OptionsCalcul options;
    
    
    
    
    template<class TB,class TP>
    void apply_bs( TB &res, TP ) const {
        res.set_type( "DataCompactClass" );
        //pour SC_create-------------------------
        //info du groupe
        APPLY_WN( res, options    );
        APPLY_WN( res, group_elements    );
        APPLY_WN( res, group_interfaces    );
        APPLY_WN( res, group_edges    );
        APPLY_WN( res, behaviour_materials    );
        APPLY_WN( res, behaviour_links    );
        APPLY_WN( res, behaviour_bc    );
        APPLY_WN( res, behaviour_bc_volume    );
    }
 
    //******************************************************************************************
    //Méthodes
    //******************************************************************************************
    //lecture du json pour les differente structure de données------------------------------------------------------------
    void read_json_groups_elements(const Object& gr);	                        //lecture des groupes d'elements
    void read_json_groups_interfaces(const Object& gr);	                        //lecture des groupes d'interface
    void read_json_groups_edges(const Object& gr);	                        //lecture des groupes de bord
    
    void read_json_behaviour_materials(const Object& gr);	                //lecture des données matériaux
    void read_json_behaviour_interfaces(const Object& gr);                      //lecture des données liaisons
    
    void read_step_bc_volume(const Object& gr, BasicVec<StepBcVolume> &step);	//lecture des step BCvolumes
    void read_json_behaviour_bc_volume(const Object& gr);		        //lecture des BCvolumes
    
    void read_step_bc(const Object& gr, BasicVec<StepBc> &step);		//lecture des step BC
    void read_json_behaviour_bc(const Object& gr);		                //lecture des BC
    
    void read_json_calcul(std::string file_calcul); 
    
};



#endif //DATA_USER_H