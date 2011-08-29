#include "write_xdmf_geometry_fields.h"
#include <iostream>
#include <fstream>
#include <cassert>
using namespace Metil;

#include "write_xdmf.h"

void write_xdmf_geometry_fields(std::ofstream &f, Param &process, String name_element, String name_node, String generic_grid_name, String grid_collection_name, bool with_time, BasicVec<String> &attributs){
    String name_hdf5=process.affichage->name_hdf;
    String name_geometry=process.affichage->name_geometry;
    String name_fields=process.affichage->name_fields;
    
    
    String name_nodes = name_geometry + "/" + name_node;
    String name_elements = name_geometry + "/" + name_element;
    
    
    int i_proc_ini=1;
    if(process.size==1) i_proc_ini=0;        
    
    //lecture du fichier hdf5 créé par chaque processeur
    for(unsigned i_proc=i_proc_ini; i_proc < process.size ;i_proc++){
        String input_hdf5; input_hdf5 << name_hdf5 <<"_"<<i_proc<<".h5";
        Hdf hdf(input_hdf5.c_str());    

        PRINT(input_hdf5);
        int nb_nodes;
        String data_node_x = name_nodes + "/x";
        hdf.read_size(data_node_x,nb_nodes);
    
        //ecriture des dataitems noeuds
        //write_nodes_datasets(f,input_hdf5, name_nodes,nb_nodes,i_proc);
    
        //ecriture des dataitems list d'elements connectivites
        PRINT("dataset");
        write_groups_datasets_2(f,input_hdf5, name_elements);
        
        //ecriture des dataitems attributs (noeuds et elements) pour chaque piquet de temps ou avant calcul
        if(with_time==1){
            process.temps->pt_cur=0;
            for(unsigned i_step=0;i_step<process.temps->nb_step;i_step++){
                for(unsigned i_pt = 0 ; i_pt < process.temps->time_step[i_step].nb_time_step; i_pt++){
                    process.temps->pt_cur+=1;
                    PRINT("nodes");
                    write_nodes_attributs_datasets(f,input_hdf5, nb_nodes,name_fields, attributs,process.temps->pt_cur, i_proc);
                    PRINT("groups");
                    write_groups_attributs_datasets(f,input_hdf5, name_elements, name_fields, attributs, process.temps->pt_cur);
                }
            }
        }
        else{
            write_groups_attributs_datasets(f,input_hdf5, name_elements,name_fields, attributs, 0);
        }
    }
    
    //ecriture des grids
    //grilles temporelles avec champs associes ou grilles simples
    if(with_time==1){
        f<< "               <Grid Name=\""<< grid_collection_name.c_str() <<"\" GridType=\"Collection\" CollectionType=\"Temporal\" >" << std::endl;
        process.temps->pt_cur=0;
        for(unsigned i_step=0;i_step<process.temps->nb_step;i_step++){
            for(unsigned i_pt = 0 ; i_pt < process.temps->time_step[i_step].nb_time_step; i_pt++){
                process.temps->pt_cur+=1;
                TYPE val_time=process.temps->time_step[i_step].t_ini+(i_pt+1)*process.temps->time_step[i_step].dt ;
                f<<"                    <Grid Name=\"Time "<<process.temps->pt_cur<<"\"  GridType=\"Collection\" >" << endl;
                f<<"                            <Time Value=\""<<val_time<<"\" />" << endl;
                for(unsigned i_proc=i_proc_ini; i_proc < process.size ;i_proc++){
                    write_grids_2(f,name_hdf5, name_elements,name_nodes, name_fields, generic_grid_name, attributs, process.temps->pt_cur, i_proc);
                }
            f<<"                        </Grid>"<< std::endl;
            }
        }
        f<<"                </Grid>"<< std::endl;        
    }
    else{
        f<< "               <Grid Name=\""<< grid_collection_name.c_str() <<"\" GridType=\"Tree\" >" << std::endl;
        for(unsigned i_proc=i_proc_ini; i_proc < process.size ;i_proc++){
            write_grids_2(f,name_hdf5, name_elements,name_nodes,name_fields, generic_grid_name,attributs, 0, i_proc);
        }    
        f<<"                </Grid>"<< std::endl;
    }           
}

void write_nodes(std::ofstream &f, Param &process, String name_node){
    String name_hdf5=process.affichage->name_hdf;
    String name_nodes = process.affichage->name_geometry + "/" + name_node;
    
    int i_proc_ini=1;
    if(process.size==1) i_proc_ini=0;        
    std::cout << process.size << endl;
    
    //lecture du fichier hdf5 créé par chaque processeur
    for(unsigned i_proc=i_proc_ini; i_proc < process.size ;i_proc++){
        String input_hdf5; input_hdf5 << name_hdf5 <<"_"<<i_proc<<".h5";
        Hdf hdf(input_hdf5.c_str());    

        int nb_nodes;
        String data_node_x = name_nodes + "/x";
        hdf.read_size(data_node_x,nb_nodes);
    
        //ecriture des dataitems noeuds
        write_nodes_datasets(f,input_hdf5, name_nodes,nb_nodes,i_proc);
    }
}

void write_xdmf_file_geometry(Param &process, DataUser &data_user){
/*        process.affichage->name_xdmf_geometry << data_user.name_directory.c_str() << "/calcul_" << data_user.id_calcul.c_str() << "/results/geometry.xdmf";*/
//         process.affichage->name_geometry = "/Level_0/Geometry";
//         process.affichage->name_fields = "/Level_0/Fields";
        String name_xdmf_file; 
        name_xdmf_file << data_user.name_directory.c_str() << "/calcul_" << data_user.id_calcul.c_str() << "/results/geometry";
        BasicVec<String> xdmf_file_type=BasicVec<String>("_elements_0.xdmf","_elements_0_skin.xdmf","_elements_1.xdmf");
        for(unsigned i=0;i<xdmf_file_type.size();i++){
            process.affichage->name_xdmf_geometry =name_xdmf_file+xdmf_file_type[i];
            std::ofstream f(process.affichage->name_xdmf_geometry.c_str());
            f << "<Xdmf Version=\"2.0\">" << std::endl;
            f << "  <Domain>" << std::endl;
            
            //ecriture des noeuds (les mêmes pour tous les maillages)
            write_nodes(f, process, "local_nodes");
        
            if(xdmf_file_type[i]=="_elements_0.xdmf"){
                //ecriture du maillage des sst
                BasicVec<String> attributs("none");
                write_xdmf_geometry_fields(f, process,"elements_0", "local_nodes", "piece_volume" ,"Geometry_0",0, attributs);
            }
            else if(xdmf_file_type[i]=="_elements_0_skin.xdmf"){
                //ecriture du maillage de peau des sst
                BasicVec<String> attributs=BasicVec<String>("num_proc", "material", "num_group");
                write_xdmf_geometry_fields(f, process,"elements_0_skin", "local_nodes", "piece_peau" ,"Geometry_0_skin",0, attributs);
            }
            else if(xdmf_file_type[i]=="_elements_1.xdmf"){
                //ecriture du maillage des interfaces
                BasicVec<String> attributs=BasicVec<String>("number", "nature");
                write_xdmf_geometry_fields(f, process,"elements_1", "local_nodes", "interface" ,"Geometry_1",0, attributs);
            }
            f<<"    </Domain>"<< std::endl;
            f<<"</Xdmf>"<<std::endl;
            f.close();     
        }
}

void write_xdmf_file_compute(Param &process, DataUser &data_user){
    //Ecriture du fichier contenant les donnees sur les sst
        String name_xdmf_file; 
        name_xdmf_file << data_user.name_directory.c_str() << "/calcul_" << data_user.id_calcul.c_str() << "/results/geometry_fields";
        BasicVec<String> xdmf_file_type=BasicVec<String>("_elements_0.xdmf","_elements_0_skin.xdmf","_elements_1.xdmf");
        for(unsigned i=0;i<xdmf_file_type.size();i++){
            process.affichage->name_xdmf_fields =name_xdmf_file+xdmf_file_type[i];
            std::ofstream f(process.affichage->name_xdmf_fields.c_str());
            f << "<Xdmf Version=\"2.0\">" << std::endl;
            f << "  <Domain>" << std::endl;
            //ecriture des noeuds (les mêmes pour tous les maillages)
            write_nodes(f, process, "local_nodes");
            if(xdmf_file_type[i]=="_elements_0.xdmf"){
                //ecriture des donnees de calcul au maillage des sst
                BasicVec<String> attributs("displacements","sigma","epsilon","sigma_von_mises");
                write_xdmf_geometry_fields(f,process, "elements_0", "local_nodes", "piece_volume" ,"Geometry_0",1, attributs);
            }
            else if(xdmf_file_type[i]=="_elements_0_skin.xdmf"){
                //ecriture des donnees de calcul au maillage de peau des sst
                BasicVec<String> attributs=BasicVec<String>("sigma_skin","epsilon_skin","sigma_von_mises_skin");
                write_xdmf_geometry_fields(f,process, "elements_0_skin", "local_nodes", "piece_peau" ,"Geometry_0_skin",1, attributs);
            }
            else if(xdmf_file_type[i]=="_elements_1.xdmf"){
                //ecriture des donnees de calcul au maillage des interfaces
                BasicVec<String> attributs=BasicVec<String>("number", "nature","F","W", "Fchap", "Wchap");
                write_xdmf_geometry_fields(f, process,"elements_1", "local_nodes", "interface" ,"Geometry_1",1, attributs);
            }
            f<<"    </Domain>"<< std::endl;
            f<<"</Xdmf>"<<std::endl;
            f.close();     
        }
}
