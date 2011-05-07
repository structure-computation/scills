//
// C++ Implementation: write_in_xdmf
//
// Description:
//
//
// Author: David Violeau <dvioleau@structure-computation.com>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//


#include <Metil/Hdf.h>
#include <iostream>
#include <fstream>
#include <cassert>

using namespace Metil;
using namespace std;

void write_nodes_datasets(std::ofstream &f,String input_hdf5, String name_nodes,int nb_nodes) {

    //ecriture des coordonnées des noeuds
    f <<"           <DataItem Name=\"X\" Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" Dimensions=\" "<< nb_nodes << " 1 \">" << input_hdf5.c_str() <<":"<< name_nodes.c_str() << "/x </DataItem>" << endl;
    f <<"           <DataItem Name=\"Y\" Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" Dimensions=\" "<< nb_nodes << " 1 \">" << input_hdf5.c_str() <<":"<< name_nodes.c_str() << "/y </DataItem>" << endl;
#if DIM==3
    f <<"           <DataItem Name=\"Z\" Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" Dimensions=\" "<< nb_nodes << " 1 \">" << input_hdf5.c_str() <<":"<< name_nodes.c_str() << "/z </DataItem>" << endl;
#endif

}

void write_nodes_attributs_datasets(std::ofstream &f,String input_hdf5, String name_nodes,int nb_nodes,String name_group_fields, BasicVec<String>  attributs, int time) {
    //ecriture des attributs aux noeuds
    for (unsigned attr=0;attr<attributs.size();attr++) {
        if (attributs[attr]=="displacements") {
            f <<"           <DataItem Name=\"dx_"<<time<<"\" Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" Dimensions=\" "<< nb_nodes << " 1 \">" << input_hdf5.c_str() <<":"<< name_group_fields.c_str() <<"/displacements/pt_"<<time<<"/x </DataItem>" << endl;
            f <<"           <DataItem Name=\"dy_"<<time<<"\" Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" Dimensions=\" "<< nb_nodes << " 1 \">" << input_hdf5.c_str() <<":"<< name_group_fields.c_str() <<"/displacements/pt_"<<time<<"/y </DataItem>" << endl;
#if DIM==3
            f <<"           <DataItem Name=\"dz_"<<time<<"\" Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" Dimensions=\" "<< nb_nodes << " 1 \">" << input_hdf5.c_str() <<":"<< name_group_fields.c_str() <<"/displacements/pt_"<<time<<"/z </DataItem>" << endl;
#endif
        }
    }
}

void write_groups_datasets(std::ofstream &f,String input_hdf5, String name_group_elements,int nb_list, int type_element) {
    Hdf hdf(input_hdf5.c_str());

    for (unsigned i=0; i< nb_list; i++) {
        String name_list;
        name_list << name_group_elements <<"/list_" <<  i;
        String name_list_c0=name_list+"/c0";
        int nb_elems;
        hdf.read_size(name_list_c0,nb_elems);
        int nb_nodes_elem;
        int pattern_id;
        if(type_element==0){
            hdf.read_tag(name_list,"pattern_id",pattern_id,1);
            if (pattern_id==0) nb_nodes_elem=3;             // Triangle
            else if (pattern_id==1) nb_nodes_elem=3;        // Triangle_6, on ne prend que les 3 premiers noeuds
            else if (pattern_id==2) nb_nodes_elem=4;        // Tetra
            else if (pattern_id==3) nb_nodes_elem=4;        // Tetra_10,   on ne prend que les 4 premiers noeuds
            else {
                std::cerr << "Type d'element non implementé - voir le tag \"base\" dans le fichier hdf5" ;
                assert(0);
            }
        }else if(type_element==1){ 
            hdf.read_tag(name_list,"interface_id",pattern_id,1);
            if (pattern_id==0) nb_nodes_elem=2;             // Bar
            else if (pattern_id==1) nb_nodes_elem=3;        // Bar_3
            else if (pattern_id==2) nb_nodes_elem=3;        // Triangle
            else if (pattern_id==3) nb_nodes_elem=3;        // Triangle_6, on ne prend que les 3 premiers noeuds
            else {
                std::cerr << "Type d'element non implementé - voir le tag \"base\" dans le fichier hdf5" ;
                assert(0);
            }
        }
        for (unsigned j=0;j<nb_nodes_elem ; j++) {
            String data_item;
            data_item << name_list << "/c" <<j ;
            f <<"           <DataItem Name=\"" << data_item.c_str() << "\" Format=\"HDF\" NumberType=\"Int\" Dimensions=\" "<< nb_elems << " 1\">" <<  input_hdf5.c_str() <<":"<< data_item.c_str() <<" </DataItem>" << endl;
        }
    }
}

void write_groups_attributs_datasets(std::ofstream &f,String input_hdf5, String name_group_elements,int nb_list,String name_group_fields, BasicVec<String>  attributs, int time) {
    Hdf hdf(input_hdf5.c_str());

    for (unsigned i=0; i< nb_list; i++) {
        String name_list;
        name_list << name_group_elements <<"/list_" <<  i;
        String name_list_c0=name_list+"/c0";
        int nb_elems;
        hdf.read_size(name_list_c0,nb_elems);

        //Ecriture des attributs en fonction de leur nom
        for (unsigned attr=0;attr<attributs.size();attr++) {
            if (attributs[attr]=="sigma") {
                f <<"           <DataItem Name=\"" << name_group_elements.c_str() << "_sxx_"<<time<<"\" Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" Dimensions=\" "<< nb_elems << " \">" << input_hdf5.c_str() <<":"<< name_group_fields.c_str() <<"/sigma/pt_"<<time<<"/xx </DataItem>" << endl;
                f <<"           <DataItem Name=\"" << name_group_elements.c_str() << "_syy_"<<time<<"\" Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" Dimensions=\" "<< nb_elems << " \">" << input_hdf5.c_str() <<":"<< name_group_fields.c_str() <<"/sigma/pt_"<<time<<"/yy </DataItem>" << endl;
                f <<"           <DataItem Name=\"" << name_group_elements.c_str() << "_szz_"<<time<<"\" Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" Dimensions=\" "<< nb_elems << " \">" << input_hdf5.c_str() <<":"<< name_group_fields.c_str() <<"/sigma/pt_"<<time<<"/zz </DataItem>" << endl;
                f <<"           <DataItem Name=\"" << name_group_elements.c_str() << "_sxy_"<<time<<"\" Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" Dimensions=\" "<< nb_elems << " \">" << input_hdf5.c_str() <<":"<< name_group_fields.c_str() <<"/sigma/pt_"<<time<<"/xy </DataItem>" << endl;
                f <<"           <DataItem Name=\"" << name_group_elements.c_str() << "_sxz_"<<time<<"\" Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" Dimensions=\" "<< nb_elems << " \">" << input_hdf5.c_str() <<":"<< name_group_fields.c_str() <<"/sigma/pt_"<<time<<"/xz </DataItem>" << endl;
                f <<"           <DataItem Name=\"" << name_group_elements.c_str() << "_syz_"<<time<<"\" Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" Dimensions=\" "<< nb_elems << " \">" << input_hdf5.c_str() <<":"<< name_group_fields.c_str() <<"/sigma/pt_"<<time<<"/yz </DataItem>" << endl;
            }
            else if (attributs[attr]=="sigma_von_mises") {
                f <<"           <DataItem Name=\"" << name_group_elements.c_str() << "_smises_"<<time<<"\" Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" Dimensions=\" "<< nb_elems << " \">" << input_hdf5.c_str() <<":"<< name_group_fields.c_str() <<"/sigma_von_mises/pt_"<<time<<" </DataItem>" << endl;
            }
            else if (attributs[attr]=="epsilon") {
                f <<"           <DataItem Name=\"" << name_group_elements.c_str() << "_exx_"<<time<<"\" Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" Dimensions=\" "<< nb_elems << " \">" << input_hdf5.c_str() <<":"<< name_group_fields.c_str() <<"/epsilon/pt_"<<time<<"/xx </DataItem>" << endl;
                f <<"           <DataItem Name=\"" << name_group_elements.c_str() << "_eyy_"<<time<<"\" Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" Dimensions=\" "<< nb_elems << " \">" << input_hdf5.c_str() <<":"<< name_group_fields.c_str() <<"/epsilon/pt_"<<time<<"/yy </DataItem>" << endl;
                f <<"           <DataItem Name=\"" << name_group_elements.c_str() << "_ezz_"<<time<<"\" Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" Dimensions=\" "<< nb_elems << " \">" << input_hdf5.c_str() <<":"<< name_group_fields.c_str() <<"/epsilon/pt_"<<time<<"/zz </DataItem>" << endl;
                f <<"           <DataItem Name=\"" << name_group_elements.c_str() << "_exy_"<<time<<"\" Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" Dimensions=\" "<< nb_elems << " \">" << input_hdf5.c_str() <<":"<< name_group_fields.c_str() <<"/epsilon/pt_"<<time<<"/xy </DataItem>" << endl;
                f <<"           <DataItem Name=\"" << name_group_elements.c_str() << "_exz_"<<time<<"\" Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" Dimensions=\" "<< nb_elems << " \">" << input_hdf5.c_str() <<":"<< name_group_fields.c_str() <<"/epsilon/pt_"<<time<<"/xz </DataItem>" << endl;
                f <<"           <DataItem Name=\"" << name_group_elements.c_str() << "_eyz_"<<time<<"\" Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" Dimensions=\" "<< nb_elems << " \">" << input_hdf5.c_str() <<":"<< name_group_fields.c_str() <<"/epsilon/pt_"<<time<<"/yz </DataItem>" << endl;
            }
        }
    }
}

void write_grids(std::ofstream &f,String input_hdf5, String name_group_elements, BasicVec<int> list, int nb_nodes, String generic_grid_name, String name_fields, BasicVec<String>  attributs, int nt, TYPE val_time, int type_element) {
    Hdf hdf(input_hdf5.c_str());

    for (unsigned i=0; i< list.size(); i++) {
        String name_list;
        name_list << name_group_elements <<  "/list_" <<  list[i];
        String name_list_c0=name_list+"/c0";
        int nb_elems;
        hdf.read_size(name_list_c0,nb_elems);
        int nb_nodes_elem;
        String JOIN, name_element_xdmf;
        
        if(type_element==0){
            int pattern_id;
            hdf.read_tag(name_list,"pattern_id",pattern_id);
            if (pattern_id==0 or pattern_id==1) {
                nb_nodes_elem=3;
                JOIN="JOIN($0 , $1 , $2  )";
                name_element_xdmf="Triangle";
            }
            else if (pattern_id==2 or pattern_id==3) {
                nb_nodes_elem=4;
                JOIN="JOIN($0 , $1 , $2  , $3)";
                name_element_xdmf="Tetrahedron";
            }
            else {
                std::cerr << "Type d'element enfant non implementé" ;
                assert(0);
            }
        }else if(type_element==1){  
            int interface_id;
            hdf.read_tag(name_list,"interface_id",interface_id);
            if (interface_id==0 or interface_id==1) {
                nb_nodes_elem=2;
                JOIN="JOIN($0 , $1 )";
                name_element_xdmf="Bar";
            }
            else if (interface_id==2 or interface_id==3) {
                nb_nodes_elem=3;
                JOIN="JOIN($0 , $1 , $2  )";
                name_element_xdmf="Triangle";
            }
            else {
                std::cerr << "Type d'element enfant non implementé" ;
                assert(0);
            }
        }
        
        f<<"                    <Grid Name=\""<< generic_grid_name.c_str() << "_" << i <<"\">" << endl;
        f<<"                            <Time Value=\""<<val_time<<"\" />" << endl;
        f<<"                            <Topology Type=\""<<name_element_xdmf.c_str()<<"\" NumberOfElements=\" "<< nb_elems << " \" >" << endl;
        f<<"                                    <DataItem Dimensions=\" "<< nb_elems << "  " << nb_nodes_elem <<" \" ItemType=\"Function\" Function=\""<<JOIN.c_str()<<"\"> " <<endl;
        for (unsigned j=0;j<nb_nodes_elem ; j++) {
            String data_item;
            data_item << name_list << "/c" <<j ;
            f<<"                                            <DataItem Reference=\"XML\"> /Xdmf/Domain/DataItem[@Name=\""<< data_item.c_str() <<"\"] </DataItem> " << endl;
        }
        f<<"                                     </DataItem> " <<endl;
        f<<"                            </Topology>" <<endl;
#if DIM==2
        f<<"                            <Geometry Type=\"XY_Z\">"<<endl;
#else
        f<<"                            <Geometry Type=\"X_Y_Z\">"<<endl;
#endif
        f<<"                                    <DataItem Reference=\"XML\"> /Xdmf/Domain/DataItem[@Name=\"X\"] </DataItem> " << endl;
        f<<"                                    <DataItem Reference=\"XML\"> /Xdmf/Domain/DataItem[@Name=\"Y\"] </DataItem> " << endl;
#if DIM==3
        f<<"                                    <DataItem Reference=\"XML\"> /Xdmf/Domain/DataItem[@Name=\"Z\"] </DataItem> " << endl;
#endif
        f<<"                            </Geometry>" <<endl;
        for (unsigned attr=0;attr<attributs.size();attr++) {
            if (attributs[attr]=="displacements") {
                f<<"                            <Attribute Name=\"displacements\"  AttributeType=\"Vector\" Center=\"Node\">" << endl;
                f<<"                                    <DataItem Dimensions=\" "<< nb_nodes << " "<< DIM <<" \" ItemType=\"Function\" Function=\""<<JOIN.c_str()<<"\"> " <<endl;
                f<<"                                            <DataItem Reference=\"XML\"> /Xdmf/Domain/DataItem[@Name=\"dx_"<<nt<<"\"] </DataItem> " << endl;
                f<<"                                            <DataItem Reference=\"XML\"> /Xdmf/Domain/DataItem[@Name=\"dy_"<<nt<<"\"] </DataItem> " << endl;
#if DIM==3
                f<<"                                            <DataItem Reference=\"XML\"> /Xdmf/Domain/DataItem[@Name=\"dz_"<<nt<<"\"] </DataItem> " << endl;
#endif
                f<<"                                     </DataItem> " <<endl;
                f<<"                            </Attribute>"<<endl;
            }
            else if(attributs[attr]=="sigma") {
#if DIM==2
                int nb_comp_tensor=3;
                String JOIN_tensor="JOIN($0 , $1 , $2 )";
                BasicVec<String> tensor_comp_sigma=BasicVec<String>("sxx","syy","sxy");
#else
                int nb_comp_tensor=6;
                String JOIN_tensor="JOIN($0 , $1 , $2, $3 , $4 , $5 )";
                BasicVec<String> tensor_comp_sigma=BasicVec<String>("sxx","syy","szz","sxy","sxz","syz");
#endif
                f<<"                            <Attribute AttributeType=\"Tensor6\" Name=\"sigma\" Center=\"Cell\">" << endl;
                f<<"                                    <DataItem Dimensions=\" "<< nb_elems << " " << nb_comp_tensor << " \" ItemType=\"Function\" Function=\""<<JOIN_tensor.c_str()<<"\"> " <<endl;
                for(unsigned nc=0;nc<nb_comp_tensor;nc++)
                    f<<"                                    <DataItem Reference=\"XML\"> /Xdmf/Domain/DataItem[@Name=\""<< name_group_elements.c_str() << "_"<<tensor_comp_sigma[nc].c_str()<<"_"<<nt<<"\"] </DataItem> " << endl;
                f<<"                                    </DataItem> " <<endl;
                f<<"                            </Attribute>"<<endl;
            }
            else if(attributs[attr]=="epsilon") {
#if DIM==2
                int nb_comp_tensor=3;
                String JOIN_tensor="JOIN($0 , $1 , $2 )";
                BasicVec<String> tensor_comp_sigma=BasicVec<String>("exx","eyy","exy");
#else
                int nb_comp_tensor=6;
                String JOIN_tensor="JOIN($0 , $1 , $2, $3 , $4 , $5 )";
                BasicVec<String> tensor_comp_sigma=BasicVec<String>("exx","eyy","ezz","exy","exz","eyz");
#endif
                f<<"                            <Attribute AttributeType=\"Tensor6\" Name=\"epsilon\" Center=\"Cell\">" << endl;
                f<<"                                    <DataItem Dimensions=\" "<< nb_elems << " " << nb_comp_tensor << " \" ItemType=\"Function\" Function=\""<<JOIN_tensor.c_str()<<"\"> " <<endl;
                for(unsigned nc=0;nc<nb_comp_tensor;nc++)
                    f<<"                                    <DataItem Reference=\"XML\"> /Xdmf/Domain/DataItem[@Name=\""<< name_group_elements.c_str() << "_"<<tensor_comp_sigma[nc].c_str()<<"_"<<nt<<"\"] </DataItem> " << endl;
                f<<"                                    </DataItem> " <<endl;
                f<<"                            </Attribute>"<<endl;
            }
            else if(attributs[attr]=="sigma_von_mises") {
                f<<"                            <Attribute AttributeType=\"Scalar\" Name=\"sigma_von_mises\" Center=\"Cell\">" << endl;
                f<<"                                    <DataItem Dimensions=\" "<< nb_elems << " 1 \" Reference=\"XML\"> /Xdmf/Domain/DataItem[@Name=\""<< name_group_elements.c_str() << "_smises_"<< nt<<"\"] </DataItem> " << endl;
                f<<"                            </Attribute>"<<endl;               
            }
        }
        f<<"                    </Grid>"<< endl;
    }
}

