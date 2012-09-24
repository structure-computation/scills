#include <Soca/Model/TypedArray.h>
#include <QtCore/QTemporaryFile>
#include <QtCore/QDataStream>
#include <containers/vec.h>
#include "ScillsResultUpdater.h"



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
          
          MP om = mp[ "_mesh" ];
          om[ "points" ].clear();
          om[ "_elements" ].clear();
          
          MP oec = mp[ "_children[ 0 ].ScillsEdgeSetItem" ];            // output Edge collection
          PRINT(oec.ok());
          //if(oec.ok()){
          //PRINT(geometry_user.mesh_nodes.size());
          //PRINT(geometry_user.mesh_nodes[0].size());
          
              for (int i_group=0; i_group<geometry_user.group_interfaces.size(); i_group++) {
                  if(geometry_user.group_interfaces[i_group].type == 0){
                      MP edge = MP::new_obj( "ScillsEdgeItem" );
                      oec[ "_children[ 0 ]" ] << edge;
                    
                    
                      PRINT(geometry_user.group_interfaces[i_group].local_nodes.size());
                      PRINT(geometry_user.group_interfaces[i_group].local_nodes[0].size());
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

                      break;
                  }
              }
          //}
          
          // version 2------------------------------------------------------------------
          
          
//           Ptr<DisplayScene> scene = NEW( DisplayScene );
//           scene->load_hdf( "/share/sc2/Developpement/MODEL/model_217/MESH/visu_geometry.h5" );
//           
//           
//           MP opc = mp[ "_children[ 0 ].ScillsPartSetItem" ];            // output Part collection
//           MP oic = mp[ "_children[ 0 ].ScillsInterfaceSetItem" ];       // output Interface collection
//           MP oec = mp[ "_children[ 0 ].ScillsEdgeSetItem" ];            // output Edge collection
//           
//           MP om = mp[ "_mesh" ];
//           om[ "points" ].clear();
//           om[ "_elements" ].clear();
//           
//           for(int i=0; i<scene->group_data.size(); i++){
//               if(scene->group_data[i].base == "Triangle"){
//                   add_message( mp, ET_Info, "add group Triangle" );
//                   MP part = MP::new_obj( "ScillsPartItem" );
//                   opc["_children"] << part; 
//                   //opc.flush();
//               }
//           }
          //opc.flush();
          mp.flush();
          
          
      }else{
          add_message( mp, ET_Info, "You choose to visualize a result" );
      }
      add_message( mp, ET_Info, "ScillsResult just finish" );
      
}



