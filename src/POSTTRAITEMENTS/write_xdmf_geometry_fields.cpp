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
    
        //ecriture des dataitems list d'elements connectivites
        write_groups_datasets_2(f,input_hdf5, name_elements);
        
        //ecriture des dataitems attributs (noeuds et elements) pour chaque piquet de temps ou avant calcul
        if(with_time==1){
            process.temps->pt_cur=0;
            for(unsigned i_step=0;i_step<process.temps->nb_step;i_step++){
                for(unsigned i_pt = 0 ; i_pt < process.temps->time_step[i_step].nb_time_step; i_pt++){
                    process.temps->pt_cur+=1;
                    write_nodes_attributs_datasets(f,input_hdf5, nb_nodes,name_fields, attributs,process.temps->pt_cur, i_proc);
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
