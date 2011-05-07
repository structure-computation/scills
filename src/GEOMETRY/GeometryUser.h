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


#include <Metil/BasicVec.h>
#include <Metil/StructCompactor.h>

#include "Patterns.h"
#include "Interfaces.h"
#include "EntityElementUser.h"
#include "EntityInterfaceUser.h"
#include "GroupElementsUser.h"
#include "GroupInterfacesUser.h"
#include "MeshUser.h"
#include "DataUser.h"

using namespace Metil;

#ifndef GEOMETRYUSER_H
#define GEOMETRYUSER_H

class GeometryUser{
    public:
      //******************************************************************************************
      //Attributs
      //******************************************************************************************
      //Attributs génériques--------------------------------------------------------------------------------------------------------------------------
      int num_level;
      int dim;
      Patterns patterns;
      Interfaces interfaces;
      //Attributs level sup---------------------------------------------------------------------------------------------------------------------------
      
      //Attributs level courant-----------------------------------------------------------------------------------------------------------------------
      
      BasicVec< BasicVec< TYPE > > nodes; //coordonnées des noeuds (dim , nb_nodes)
      
      //Attributs des elements
      int nb_group_elements;
      BasicVec<GroupElementsUser> group_elements;
          
      //Attributs des interfaces
      int nb_group_interfaces;
      BasicVec< GroupInterfacesUser > group_interfaces;
       
      //***********************************************************************************
      //Methode pour générer la class compacté GEOMETRY.h
      //***********************************************************************************   
      template<class TB,class TP>
      void apply_bs( TB &res, TP ) const {
          res.set_type( "GeometryCompactClass" );
          APPLY_WN( res, num_level );
          APPLY_WN( res, dim );
          APPLY_WN( res, patterns );
          APPLY_WN( res, nodes );
          APPLY_WN( res, group_elements );
          APPLY_WN( res, group_interfaces );
          
      }
      //***********************************************************************************
      //Methode pour SC_create
      //*********************************************************************************** 
      //Methode générique  ---------------------------------------------------------------------------------------------------------------------------
      
      GroupElementsUser* find_group_elements(int id_);              // recherche d'un group_elements particulier avec son id
      int find_index_group_elements(int id_);                   // recherche de l'index d'un group_elements particulier avec son id
      
      GroupInterfacesUser* find_group_interfaces(int id_);          // recherche d'un group_interfaces particulier avec son id
      int find_index_group_interfaces(int id_);                 // recherche de l'index d'un group_interfaces particulier avec son id
      

      //Methode d'initialisation----------------------------------------------------------------------------------------------------------------------
      GeometryUser();
      GeometryUser(MeshUser &mesh_);
      void initialize_group_elements_from_MeshUser(MeshUser &mesh_);                    // initialisation à partir des elements du maillage
      void initialize_group_interfaces_from_MeshUser(MeshUser &mesh_);                  // initialisation à partir des interfaces du maillage
      void write_json(MeshUser &mesh_user);                                             // ecriture du fichier json pour l'interface
      void write_visu_hdf5(MeshUser &mesh_user);                                        // ecriture du fichier hdf5 à partir de la géométrie
      void write_xdmf(String output_xdmf, String input_hdf5, String name_geometry, int skin); //Ecriture du fichier xdmf "output_xdmf" avec références aux données du fichier hdf5 "input_hdf5" 
      void read_hdf5(String name_file_hdf5);                                            //lecture des données du fichier hdf5 et assignation des champs de la classe GEOMETRY_USER
      
      void initialize_GPU();                                                            // initialisation pour calcul sur GPU

      //Methode de niveau sup-------------------------------------------------------------------------------------------------------------------------
      bool do_respect_geometry(int i_group, int num_edge, DataUser::Geometry &geom);    // repere si un element d'interface est dans une geometry donnée
      void split_group_edges_within_geometry(DataUser &data_user);                      // séparation d'un group_interfaces de type 0 suivant un critère géometrique au moment du calcul

      //Methode de niveau courrant -------------------------------------------------------------------------------------------------------------------
      
      
      
       
};

#endif //GEOMETRY_H


