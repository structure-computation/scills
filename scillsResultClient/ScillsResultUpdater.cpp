#include <Soca/Model/TypedArray.h>
#include <QtCore/QTemporaryFile>
#include <QtCore/QDataStream>
#include <containers/vec.h>
#include "ScillsResultUpdater.h"

#include <string.h>


#include "../COMPUTE/DataUser.h"
#include "../COMPUTE/FieldStructureUser.h"
#include "../GEOMETRY/GeometryUser.h"


typedef LMT::Vec<double,3> Pvec;

struct AutoRm {
    AutoRm( QString f ) : f( f ) {}
    ~AutoRm() { QFile::remove( f ); }
    QString f;
};

bool ScillsResultUpdater::run( MP mp ) {
    // add_message( mp, ET_Info, "Test info msg" );
    
      double id_model = mp[ "id_model" ];
      if ( id_model == -1 ) {
          add_message( mp, ET_Error, "Model is not available" );
          //return false;
      }
    
      double id_calcul = mp[ "id_calcul" ];
      if ( id_calcul <= -2 ) {
          add_message( mp, ET_Error, "Data are not available" );
      }else if( id_calcul == -1 ){
          add_message( mp, ET_Info, "You choose to visualize the base mesh" );
//           using namespace Metil;
//           using namespace LMT;
          
          
          // version 1------------------------------------------------------------------
          Sc2String str_id_model;
          str_id_model << id_model;
          Sc2String str_id_calcul;
          str_id_calcul << id_calcul;
          
          GeometryUser        geometry_user;                                    /// structure de stockage des informations du fichier HDF5
          geometry_user.initialisation(str_id_model, str_id_calcul);
          geometry_user.read_hdf5(false,true,"test");                           /// true si on lit les info micro, true si on lit toutes les infos
           
          MP oec = mp[ "_children[ 2 ]" ];            // output Edge collection
          //qDebug() << oec[ "_edge_profile" ]; 
    
          for (int i_group=0; i_group<geometry_user.group_interfaces.size(); i_group++) {
              if(geometry_user.group_interfaces[i_group].type == 0){
                  // MP edge_ = oec[ "_edge_profile" ];
                  MP edge = MP::new_obj( "ScillsEdgeItem" );
                  Sc2String str_name_edge;
                  str_name_edge << "edge_" << geometry_user.group_interfaces[i_group].id;
                  // edge[ "_ico" ] = edge_[ "_ico" ];
                  // edge[ "_name" ] = str_name_edge.c_str();
                  
                  edge[ "_ico" ] = oec[ "_edge_profile._ico" ];
                  edge[ "_name" ] = str_name_edge.c_str();
                  edge[ "_viewable" ] = true;
                  edge[ "_children" ] = MP::new_lst();
                  edge[ "_output" ] = MP::new_lst();
                  edge[ "_name_class" ] = "";
                  edge[ "_allow_vmod" ] = true;
                  //edge[ "_mesh" ] = oec[ "_edge_profile._mesh" ];
                  //edge[ "visualization" ] = oec[ "_edge_profile.visualization" ];
                  
                  edge[ "_mesh" ] = MP::new_obj( "Mesh" );

                  edge[ "_mesh.visualization" ] = MP::new_obj();
                  edge[ "_mesh.visualization.display_style" ] = MP::new_obj( "Choice" );
                  edge[ "_mesh.visualization.display_style.num" ] = 1;
                  edge[ "_mesh.visualization.display_style.lst" ] = MP::new_lst();
                  edge[ "_mesh.visualization.display_style.lst" ] << "Points";
                  edge[ "_mesh.visualization.display_style.lst" ] << "Wireframe";
                  edge[ "_mesh.visualization.display_style.lst" ] << "Surface";
                  edge[ "_mesh.visualization.display_style.lst" ] << "Surface with Edges";
                  edge[ "_mesh.visualization.point_edition" ] = false;

                  //qDebug() << edge ;

                  //ajout du maillage de peau
                  MP om = edge[ "_mesh" ];
                  om[ "points" ] = MP::new_lst( "Lst_Point" );
                  om[ "_elements" ] = MP::new_lst();
                  om[ "_selected_points" ] = MP::new_lst();
                  om[ "_pelected_points" ] = MP::new_lst();
                  om[ "_selected_elements" ] = MP::new_lst();
                  om[ "_pelected_elements" ] = MP::new_lst();

                  //
                  for(int i=0; i<geometry_user.group_interfaces[i_group].local_nodes[0].size(); i++){
                      MP pos = MP::new_lst( "Vec_3" );
                      for(int d=0; d<geometry_user.group_interfaces[i_group].local_nodes.size(); d++){
                          pos << geometry_user.group_interfaces[i_group].local_nodes[d][i];
                      }
                      //PRINT(geometry_user.group_interfaces[i_group].local_nodes[0][i]);
                      MP pnt = MP::new_obj( "Point" );
                      pnt[ "pos" ] = pos;

                      om[ "points" ] << pnt;
                  }

                  //
                  TypedArray<int> *connectivity = new TypedArray<int>;
                  for(int i=0; i<geometry_user.group_interfaces[i_group].local_connectivities[0].size(); i++){

                      for(int d=0; d<geometry_user.group_interfaces[i_group].local_connectivities.size(); d++){
                          connectivity->_data << geometry_user.group_interfaces[i_group].local_connectivities[d][i];
                      }

                  }
                  connectivity->_size.resize( 2 );
                  connectivity->_size[ 0 ] = 3;
                  connectivity->_size[ 1 ] = connectivity->_data.size() / 3;

                  MP triangles = MP::new_obj( "Element_TriangleList" );
                  triangles[ "indices" ] = connectivity;
                  om[ "_elements" ] << triangles;
                  
                  
                  //                  qDebug() <<  oec[ "_edge_profile" ];
                  //                  qDebug() <<  edge ;
                  
                  
                  oec[ "_children" ] << edge;
                  
              }
          }
          mp.flush();
          
          
      }else{
          add_message( mp, ET_Info, "You choose to visualize a result" );
      }
      add_message( mp, ET_Info, "ScillsResult just finish" );
      
}


