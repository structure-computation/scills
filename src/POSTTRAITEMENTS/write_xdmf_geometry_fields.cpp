#include "write_xdmf_geometry_fields.h"
#include <fstream>
#include <cassert>
using namespace Metil;

#include "write_xdmf.h"

void write_xdmf_geometry_fields(std::ofstream &f, Process &process, Sc2String name_element, Sc2String name_node, Sc2String generic_grid_name, Sc2String grid_collection_name, bool with_time, BasicVec<Sc2String> &attributs){
    Sc2String name_hdf5=process.affichage->name_hdf;
    Sc2String name_geometry=process.affichage->name_geometry;
    Sc2String name_fields=process.affichage->name_fields;
    
    
    Sc2String name_nodes = name_geometry + "/" + name_node;
    Sc2String name_elements = name_geometry + "/" + name_element;
    
    
    int i_proc_ini=1;
    if(not process.parallelisation->is_multi_cpu()) i_proc_ini=0;        
    
    //lecture du fichier hdf5 cree par chaque processeur
    for(unsigned i_proc=i_proc_ini; i_proc < process.parallelisation->size ;i_proc++){
        Sc2String input_hdf5; input_hdf5 << name_hdf5 <<"_"<<i_proc<<".h5";
        Hdf hdf(input_hdf5.c_str());    

//         PRINT(input_hdf5);
        int nb_nodes;
        Sc2String data_node_x = name_nodes + "/x";
        hdf.read_size(data_node_x.c_str(),nb_nodes);
    
        //ecriture des dataitems noeuds
        //write_nodes_datasets(f,input_hdf5, name_nodes,nb_nodes,i_proc);
    
        //ecriture des dataitems list d'elements connectivites
        //PRINT("dataset");
        write_groups_datasets_2(f,input_hdf5, name_elements);
        
        //ecriture des dataitems attributs (noeuds et elements) pour chaque piquet de temps ou avant calcul
        if(with_time==1){
            process.temps->pt_cur=0;
            for(unsigned i_step=0;i_step<process.temps->nb_step;i_step++){
                for(unsigned i_pt = 0 ; i_pt < process.temps->time_step[i_step].nb_time_step; i_pt++){
                    process.temps->pt_cur+=1;
                    //PRINT("nodes");
                    write_nodes_attributs_datasets(f,input_hdf5, nb_nodes,name_fields, attributs,process.temps->pt_cur, i_proc);
                    //PRINT("groups");
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
                f<<"                    <Grid Name=\"Time "<<process.temps->pt_cur<<"\"  GridType=\"Collection\" >" << std::endl;
                f<<"                            <Time Value=\""<<val_time<<"\" />" << std::endl;
                for(unsigned i_proc=i_proc_ini; i_proc < process.parallelisation->size ;i_proc++){
                    write_grids_2(f,name_hdf5, name_elements,name_nodes, name_fields, generic_grid_name, attributs, process.temps->pt_cur, i_proc);
                }
            f<<"                        </Grid>"<< std::endl;
            }
        }
        f<<"                </Grid>"<< std::endl;        
    }
    else{
        f<< "               <Grid Name=\""<< grid_collection_name.c_str() <<"\" GridType=\"Tree\" >" << std::endl;
        for(unsigned i_proc=i_proc_ini; i_proc < process.parallelisation->size ;i_proc++){
            write_grids_2(f,name_hdf5, name_elements,name_nodes,name_fields, generic_grid_name,attributs, 0, i_proc);
        }    
        f<<"                </Grid>"<< std::endl;
    }           
}

void write_nodes(std::ofstream &f, Process &process, Sc2String name_node){
    Sc2String name_hdf5=process.affichage->name_hdf;
    Sc2String name_nodes = process.affichage->name_geometry + "/" + name_node;
    
    int i_proc_ini=1;
    if(not process.parallelisation->is_multi_cpu()) i_proc_ini=0;
    
    //lecture du fichier hdf5 cree par chaque processeur
    for(unsigned i_proc=i_proc_ini; i_proc < process.parallelisation->size ;i_proc++){
        Sc2String input_hdf5; input_hdf5 << name_hdf5 <<"_"<<i_proc<<".h5";
        Hdf hdf(input_hdf5.c_str());
        std::cout << "Test write_nodes 2 ************************************************************************************" << std::endl;
        
        int nb_nodes;
        Sc2String data_node_x = name_nodes + "/x";
        std::cout << "DATA_NODE_X" << data_node_x << std::endl << std::endl;
        hdf.read_size(data_node_x.c_str(),nb_nodes);
        std::cout << "Test write_nodes 3 ************************************************************************************" << std::endl;
    
        //ecriture des dataitems noeuds
        write_nodes_datasets(f,input_hdf5, name_nodes,nb_nodes,i_proc);
        std::cout << "Test write_nodes 4 ************************************************************************************" << std::endl;
    }
}

void write_xdmf_file_geometry(Process &process, DataUser &data_user){
/*        process.affichage->name_xdmf_geometry << data_user.name_directory.c_str() << "/calcul_" << data_user.id_calcul.c_str() << "/results/geometry.xdmf";*/
//         process.affichage->name_geometry = "/Level_0/Geometry";
//         process.affichage->name_fields = "/Level_0/Fields";
        Sc2String name_xdmf_file; 
        name_xdmf_file << data_user.name_directory.c_str() << "/calcul_" << data_user.id_calcul.c_str() << "/results/geometry";
        BasicVec<Sc2String> xdmf_file_type=BasicVec<Sc2String>("_elements_0.xdmf","_elements_0_skin.xdmf","_elements_1.xdmf");
        for(unsigned i=0;i<xdmf_file_type.size();i++){
            process.affichage->name_xdmf_geometry =name_xdmf_file+xdmf_file_type[i];
            std::ofstream f(process.affichage->name_xdmf_geometry.c_str());
            f << "<Xdmf Version=\"2.0\">" << std::endl;
            f << "  <Domain>" << std::endl;
            
            //ecriture des noeuds (les mêmes pour tous les maillages)
            write_nodes(f, process, "local_nodes");
            
            if(xdmf_file_type[i]=="_elements_0.xdmf"){
                //ecriture du maillage des sst
                BasicVec<Sc2String> attributs("none");
                std::cout << "Test 3a************************************************************************************************************************************************************************************************" << std::endl;
                write_xdmf_geometry_fields(f, process,"elements_0", "local_nodes", "piece_volume" ,"Geometry_0",0, attributs);
                std::cout << "Test 4a************************************************************************************************************************************************************************************************" << std::endl;
                
            }
            else if(xdmf_file_type[i]=="_elements_0_skin.xdmf"){
                //ecriture du maillage de peau des sst
                BasicVec<Sc2String> attributs=BasicVec<Sc2String>("num_proc", "material", "num_group");
                std::cout << "Test 3b************************************************************************************************************************************************************************************************" << std::endl;
                write_xdmf_geometry_fields(f, process,"elements_0_skin", "local_nodes", "piece_peau" ,"Geometry_0_skin",0, attributs);
                std::cout << "Test 4b************************************************************************************************************************************************************************************************" << std::endl;
                
            }
            else if(xdmf_file_type[i]=="_elements_1.xdmf"){
                //ecriture du maillage des interfaces
                BasicVec<Sc2String> attributs=BasicVec<Sc2String>("number", "nature");
                std::cout << "Test 3c************************************************************************************************************************************************************************************************" << std::endl;
                write_xdmf_geometry_fields(f, process,"elements_1", "local_nodes", "interface" ,"Geometry_1",0, attributs);
                std::cout << "Test 4c************************************************************************************************************************************************************************************************" << std::endl;
                
            }
            f<<"    </Domain>"<< std::endl;
            f<<"</Xdmf>"<<std::endl;
            f.close();
            std::cout << "Test 5************************************************************************************************************************************************************************************************" << std::endl;
            
        }
}

void write_xdmf_file_compute(Process &process, DataUser &data_user){
    //Ecriture du fichier contenant les donnees sur les sst
        Sc2String name_xdmf_file; 
        name_xdmf_file << data_user.name_directory.c_str() << "/calcul_" << data_user.id_calcul.c_str() << "/results/geometry_fields";
        BasicVec<Sc2String> xdmf_file_type=BasicVec<Sc2String>("_elements_0.xdmf","_elements_0_skin.xdmf","_elements_1.xdmf");
        for(unsigned i=0;i<xdmf_file_type.size();i++){
            process.affichage->name_xdmf_fields =name_xdmf_file+xdmf_file_type[i];
            std::ofstream f(process.affichage->name_xdmf_fields.c_str());
            f << "<Xdmf Version=\"2.0\">" << std::endl;
            f << "  <Domain>" << std::endl;
            //ecriture des noeuds (les mêmes pour tous les maillages)
            write_nodes(f, process, "local_nodes");
            if(xdmf_file_type[i]=="_elements_0.xdmf"){
                //ecriture des donnees de calcul au maillage des sst
                BasicVec<Sc2String> attributs("displacements","sigma","epsilon","sigma_von_mises");
                write_xdmf_geometry_fields(f,process, "elements_0", "local_nodes", "piece_volume" ,"Geometry_0",1, attributs);
            }
            else if(xdmf_file_type[i]=="_elements_0_skin.xdmf"){
                //ecriture des donnees de calcul au maillage de peau des sst
                BasicVec<Sc2String> attributs=BasicVec<Sc2String>("sigma_skin","epsilon_skin","sigma_von_mises_skin");
                write_xdmf_geometry_fields(f,process, "elements_0_skin", "local_nodes", "piece_peau" ,"Geometry_0_skin",1, attributs);
            }
            else if(xdmf_file_type[i]=="_elements_1.xdmf"){
                //ecriture des donnees de calcul au maillage des interfaces
                BasicVec<Sc2String> attributs=BasicVec<Sc2String>("number", "nature","F","W", "Fchap", "Wchap");
                write_xdmf_geometry_fields(f, process,"elements_1", "local_nodes", "interface" ,"Geometry_1",1, attributs);
            }
            f<<"    </Domain>"<< std::endl;
            f<<"</Xdmf>"<<std::endl;
            f.close();     
        }
}
