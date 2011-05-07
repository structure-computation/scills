#include "DataUser.h"
#include <string>
#include <Metil/BasicVec.h>

#include <Metil/Level1/CompilationEnvironment.h>
#include <Metil/DynamicLibrary.h>
#include <Metil/MathBasicVec.h>
#include <Metil/System.h>
#include <Metil/Math.h>
#include <Metil/Md5.h>
#include <set>
#include <map>

using namespace Metil;
using namespace json_spirit;  


//*************************************************************************************************************
// fonction de lecture du fichier de calcul
//*************************************************************************************************************
//lecture des groupes d'elements dans un fichier json--------------------------------------
void DataUser::read_json_groups_elements(const Object& gr){
    for(Object::size_type i=0;i != gr.size() ;i++){
        const Pair& pair1 = gr[i];        
        const Value&  value1 = pair1.value_;
        Object obj=value1.get_obj();
        for( Object::size_type j = 0; j != obj.size(); ++j )
        {
            const Pair& pair = obj[j];
            const std::string& name  = pair.name_;
            const Value&  value = pair.value_;

            if( name == "id" )
            {
                group_elements[i].id= value.get_int();
            }
            else if( name == "origine" )
            {
                group_elements[i].origine = value.get_str() ;
            }
            else if( name == "identificateur" )
            {
                group_elements[i].num_in_mesh_file = value.get_int() ;
            }
            else if( name == "name" )
            {
                group_elements[i].name = value.get_str();
            }
            else if( name == "group" )
            {
                //TODO
            }
            else if( name == "assigned" )
            {
                //TODO
            }
            else if(name=="id_material")
            {
                group_elements[i].id_material= value.get_int();
            }
            else
            {
                assert( "Donnees groups_elem non implementee" );
            }
        }
    }

}


///lecture des groupes d'interfaces dans un fichier json
void DataUser::read_json_groups_interfaces( const Object& gr ){
    for(Object::size_type i=0;i != gr.size() ;i++){
        const Pair& pair1 = gr[i];        
        const Value&  value1 = pair1.value_;
        Object obj=value1.get_obj();
        for( Object::size_type j = 0; j != obj.size(); ++j )
        {
            const Pair& pair = obj[j];
            const std::string& name  = pair.name_;
            const Value&  value = pair.value_;

            if( name == "id" )
            {
                group_interfaces[i].id= value.get_int();
            }
            else if( name == "origine" )
            {
                group_interfaces[i].origine= value.get_str();
            }
            else if( name == "name" )
            {
                group_interfaces[i].name= value.get_str();
            }
            else if( name == "type" )
            {
                group_interfaces[i].type= value.get_str();
            }
            else if( name == "adj_num_group" )
            {
                std::string adj=value.get_str();
                std::istringstream s(adj);
                BasicVec<int,2> adjnum;
                for(unsigned k=0;k< 2 ;k++){
                    s>>adjnum[k];
                }
                group_interfaces[i].adj_num_group=adjnum;
            }
            else if( name == "id_link" )
            {
                group_interfaces[i].id_link= value.get_int();
            }
            else if( name == "group" )
            {
                //TODO
            }
            else if( name == "assigned" )
            {
                group_interfaces[i].assigned= value.get_real();
            }
            else
            {
                assert( "Donnee groups_inter non implementee" );
            }
        }
    }
}

///lecture des groupes de bord dans un fichier json
void DataUser::read_json_groups_edges( const Object& gr){
    for(Object::size_type i=0;i != gr.size() ;i++){
        group_edges[i].geom.points.resize(2);
        for(int k=0; k<2; k++){
            group_edges[i].geom.points[k].resize(DIM);
        }
        const Pair& pair1 = gr[i];        
        const Value&  value1 = pair1.value_;
        Object obj=value1.get_obj();
        for( Object::size_type j = 0; j != obj.size(); ++j )
        {

            const Pair& pair = obj[j];
            const std::string& name  = pair.name_;
            const Value&  value = pair.value_;
            if( name == "id" )
            {
                group_edges[i].id= value.get_int();
            }
            else if( name == "origine" )
            {
                group_edges[i].geom.origine= value.get_str();
            }
            else if( name == "name" )
            {
                group_edges[i].name= value.get_str();
            }
            else if( name == "type" )
            {
                group_edges[i].geom.type= value.get_str();
            }
            else if( name == "id_CL" )
            {
                group_edges[i].id_CL= value.get_int();
            }
            else if( name == "assigned" )
            {
                group_edges[i].assigned= value.get_real();
            }
            else if(name == "pdirection_x" )
            {
                std::string temp=value.get_str();
                std::istringstream s(temp);
                s >> group_edges[i].geom.pdirection[0];
            }
            else if(name == "pdirection_y" )
            {
                std::string temp=value.get_str();
                std::istringstream s(temp);
                s >> group_edges[i].geom.pdirection[1];
            }
#if DIM==3
            else if(name == "pdirection_z" )
            {
                std::string temp=value.get_str();
                std::istringstream s(temp);
                s >> group_edges[i].geom.pdirection[2];
            }
#endif
            else if( name == "geometry" )
            {
                group_edges[i].geom.nature=value.get_str();
                group_edges[i].geom.points.resize(2);
                for(int k=0; k<2; k++){
                    group_edges[i].geom.points[k].resize(DIM);
                }
            }
            else if(name == "radius" )
            {
                std::string temp=value.get_str();
                std::istringstream is(temp);
                is >>group_edges[i].geom.radius;
//                 group_edges[i].geom.radius=value.get_real();
            }
            else if(name == "equation" )
            {
                group_edges[i].geom.equation=value.get_str();
            }
            else if(name == "point_1_x" )
            {
                std::string temp=value.get_str();
                std::istringstream s(temp);
                s >> group_edges[i].geom.points[0][0];
            }
            else if(name == "point_1_y" )
            {
                std::string temp=value.get_str();
                std::istringstream s(temp);
                s >> group_edges[i].geom.points[0][1];
            }
#if DIM==3
            else if(name == "point_1_z" )
            {
                std::string temp=value.get_str();
                std::istringstream s(temp);
                s >> group_edges[i].geom.points[0][2];                
            }
#endif
            else if(name == "point_2_x" )
            {
                std::string temp=value.get_str();
                std::istringstream s(temp);
                s >> group_edges[i].geom.points[1][0];
            }
            else if(name == "point_2_y" )
            {
                std::string temp=value.get_str();
                std::istringstream s(temp);
                s >> group_edges[i].geom.points[1][1];
            }
#if DIM==3
            else if(name == "point_2_z" )
            {
                std::string temp=value.get_str();
                std::istringstream s(temp);
                s >> group_edges[i].geom.points[1][2];
            }
#endif
            else
            {
                std::cout << "Donnee groups_edge non implementee : " << name << std::endl;
            }
        }
    }
    //Ajout du groupe en fin des groupes d'edge (comportement generique pour tous les edges non assignees
    group_edges[gr.size()].id=-1;
    group_edges[gr.size()].geom.type="all";
    group_edges[gr.size()].id_CL=-1; 
    

}

///lecture des proprietes materiaux donnees dans un fichier json
void DataUser::read_json_behaviour_materials(const Object& gr){

    for(Object::size_type i=0;i != gr.size() ;i++){
        //Object obj=pr[j].get_obj();;
    
        const Pair& pair1 = gr[i];        
        const Value&  value1 = pair1.value_;
        Object obj=value1.get_obj();  
                
        for( Object::size_type j = 0; j != obj.size(); ++j ){
            const Pair& pair = obj[j];
            const std::string& name  = pair.name_;
            const Value&  value = pair.value_;
            if( name == "id" )
            {
                behaviour_materials[i].id= value.get_int();
            }
            else if( name == "mtype" )
            {
                behaviour_materials[i].type= value.get_str();
            }
            else if( name == "type_num" )
            {
                behaviour_materials[i].type_num= value.get_int();
            }
            else if( name == "comp" )
            {
                behaviour_materials[i].comp= value.get_str();
            }
            else if( name == "resolution" )
            {
                behaviour_materials[i].resolution= value.get_str();
            }
            else if( name == "elastic_modulus" )
            {
                //behaviour_materials[i].elastic_modulus=value.get_str();
                behaviour_materials[i].mat_prop[0]=value.get_str();
            }
            else if( name == "poisson_ratio" )
            {
                //behaviour_materials[i].poisson_ratio=value.get_str();
                behaviour_materials[i].mat_prop[1]=value.get_str();
            }
            else if( name == "alpha" )
            {
                //behaviour_materials[i].alpha=value.get_str();
                behaviour_materials[i].mat_prop[2]=value.get_str();
            }
            else if( name == "rho" )
            {
                //behaviour_materials[i].rho=value.get_str();
                behaviour_materials[i].mat_prop[3]=value.get_str();
            }
            else if( name == "viscosite" )
            {
                //behaviour_materials[i].viscosite=value.get_str();
                behaviour_materials[i].mat_prop[4]=value.get_str();
            }
            else if( name == "dir_1_x" )
            {
                //behaviour_materials[i].dir_1_x = value.get_str();
                behaviour_materials[i].mat_prop[5]=value.get_str();
            }
            else if( name == "dir_1_y" )
            {
                //behaviour_materials[i].dir_1_y  = value.get_str();
                behaviour_materials[i].mat_prop[6]=value.get_str();
            }
#if DIM==3
            else if( name == "dir_1_z" )
            {
                //behaviour_materials[i].dir_1_z  = value.get_str();
                behaviour_materials[i].mat_prop[7]=value.get_str();
            }
#endif
            else if( name == "dir_2_x" )
            {
                //behaviour_materials[i].dir_2_x = value.get_str();
                behaviour_materials[i].mat_prop[8]=value.get_str();
            }
            else if( name == "dir_2_y" )
            {
                //behaviour_materials[i].dir_2_y = value.get_str();
                behaviour_materials[i].mat_prop[9]=value.get_str();
            }
#if DIM==3
            else if( name == "dir_2_z" )
            {
                //behaviour_materials[i].dir_2_z = value.get_str();
                behaviour_materials[i].mat_prop[10]=value.get_str();
            }
#endif
            else if( name == "dir_3_x" )
            {
                //behaviour_materials[i].dir_3_x = value.get_str();
                behaviour_materials[i].mat_prop[11]=value.get_str();
            }
            else if( name == "dir_3_y" )
            {
                //behaviour_materials[i].dir_3_y = value.get_str();
                behaviour_materials[i].mat_prop[12]=value.get_str();
            }
#if DIM==3
            else if( name == "dir_3_z" )
            {
                //behaviour_materials[i].dir_3_z = value.get_str();
                behaviour_materials[i].mat_prop[13]=value.get_str();
            }
#endif
            else if( name == "E1" )
            {
                //behaviour_materials[i].E1=value.get_str();
                behaviour_materials[i].mat_prop[14]=value.get_str();
            }
            else if( name == "E2" )
            {
                //behaviour_materials[i].E2=value.get_str();
                behaviour_materials[i].mat_prop[15]=value.get_str();
            }
            else if( name == "E3" )
            {
                //behaviour_materials[i].E3=value.get_str();
                behaviour_materials[i].mat_prop[16]=value.get_str();
            }
            else if( name == "G12" )
            {
                //behaviour_materials[i].G12=value.get_str();
                behaviour_materials[i].mat_prop[17]=value.get_str();
            }
            else if( name == "G23" )
            {
                //behaviour_materials[i].G23=value.get_str();
                behaviour_materials[i].mat_prop[18]=value.get_str();
            }
            else if( name == "G13" )
            {
                //behaviour_materials[i].G13=value.get_str();
                behaviour_materials[i].mat_prop[19]=value.get_str();
            }
            else if( name == "nu12" )
            {
                //behaviour_materials[i].nu12=value.get_str();
                behaviour_materials[i].mat_prop[20]=value.get_str();
            }
            else if( name == "nu23" )
            {
                //behaviour_materials[i].nu23=value.get_str();
                behaviour_materials[i].mat_prop[21]=value.get_str();
            }
            else if( name == "nu13" )
            {
                //behaviour_materials[i].nu13=value.get_str();
                behaviour_materials[i].mat_prop[22]=value.get_str();
            }
            else if( name == "alpha_1" )
            {
                //behaviour_materials[i].alpha_1=value.get_str();
                behaviour_materials[i].mat_prop[23]=value.get_str();
            }
            else if( name == "alpha_2" )
            {
                //behaviour_materials[i].alpha_2=value.get_str();
                behaviour_materials[i].mat_prop[24]=value.get_str();
            }
            else if( name == "alpha_3" )
            {
                //behaviour_materials[i].alpha_3=value.get_str();
                behaviour_materials[i].mat_prop[25]=value.get_str();
            }
            else if( name == "Yo" )
            {
                //behaviour_materials[i].Yo = value.get_str();
                behaviour_materials[i].mat_prop[26]=value.get_str();
            }
            else if( name == "Ysp" )
            {
                //behaviour_materials[i].Ysp = value.get_str();
                behaviour_materials[i].mat_prop[27]=value.get_str();
            }
            else if( name == "Yop" )
            {
                //behaviour_materials[i].Yop = value.get_str();
                behaviour_materials[i].mat_prop[28]=value.get_str();
            }
            else if( name == "Yc" )
            {
                //behaviour_materials[i].Yc = value.get_str();
                behaviour_materials[i].mat_prop[29]=value.get_str();
            }
            else if( name == "Ycp" )
            {
                //behaviour_materials[i].Ycp = value.get_str();
                behaviour_materials[i].mat_prop[30]=value.get_str();
            }
            else if( name == "b" )
            {
                //behaviour_materials[i].b = value.get_str();
                behaviour_materials[i].mat_prop[31]=value.get_str();
            }
            else
            {
                std::cerr << "Donnee materials non implementee"<< std::endl;
            }
        }
    }
}

///Lecture des proprietes des interfaces donnees dans un fichier json
void DataUser::read_json_behaviour_interfaces(const Object& gr){
   
    for(Object::size_type i=0;i != gr.size() ;i++){
        //Object obj=pr[i].get_obj();;
        const Pair& pair1 = gr[i];        
        const Value&  value1 = pair1.value_;
        Object obj=value1.get_obj();
        for( Object::size_type j = 0; j != obj.size(); ++j )
        {
            const Pair& pair = obj[j];
            const std::string& name  = pair.name_;
            const Value&  value = pair.value_;

            if( name == "id" )
            {
                behaviour_links[i].id= value.get_int();
            }
            else if( name == "id_select" )
            {
                behaviour_links[i].id_select= value.get_int();
            }
            else if( name == "type_num" )
            {
                behaviour_links[i].type_num= value.get_int();
            }
            else if( name == "type" )
            {
                behaviour_links[i].type= value.get_str();
            }
            else if( name == "comp_complexe" )
            {
                behaviour_links[i].comp_complexe= value.get_str();
            }
            else if( name == "coef_frottement" )
            {
                //behaviour_links[i].coef_frottement= value.get_str();
                behaviour_links[i].link_prop[0]=value.get_str();
            }
            else if( name == "Ep" )
            {
                //behaviour_links[i].Ep= value.get_str();
                behaviour_links[i].link_prop[1]=value.get_str();
            }
            else if( name == "jeux" )
            {
                //behaviour_links[i].jeux= value.get_str();
                behaviour_links[i].link_prop[2]=value.get_str();
            }
            else if( name == "R" )
            {
                //behaviour_links[i].R= value.get_str();
                behaviour_links[i].link_prop[3]=value.get_str();
            }
            else if( name == "Lp" )
            {
                //behaviour_links[i].Lp= value.get_str();
                behaviour_links[i].link_prop[4]=value.get_str();
            }
            else if( name == "Dp" )
            {
                //behaviour_links[i].Dp= value.get_str();
                behaviour_links[i].link_prop[5]=value.get_str();
            }
            else if( name == "p" )
            {
                //behaviour_links[i].p= value.get_str();
                behaviour_links[i].link_prop[6]=value.get_str();
            }
            else if( name == "Lr" )
            {
                //behaviour_links[i].Lr= value.get_str();
                behaviour_links[i].link_prop[7]=value.get_str();
            }
            else
            {
                std::cout << "Champ proprietes_interface non implementee" << std::endl;
            }
        }
    }  
    //Ajout du comportement parfait en fin des groupes d'interface (comportement generique pour toutes les interfaces non assignees
    behaviour_links[gr.size()].type="Parfait";
    behaviour_links[gr.size()].id=-1; 
}

void DataUser::read_step_bc_volume(const Object& gr, BasicVec<StepBcVolume> &step){
    for( Object::size_type k = 0; k != gr.size(); ++k )
    {
        const Pair& pair1 = gr[k];
        //const std::string& name1  = pair1.name_;
        const Value&  value1 = pair1.value_;
        Object obj=value1.get_obj();
        for( Object::size_type l = 0; l != obj.size(); ++l )
        {
            const Pair& pair = obj[l];
            const std::string& name  = pair.name_;
            const Value&  value = pair.value_;
            
            if( name == "pdirection_x" )
            {
                step[k].CLv_step_prop[0]=value.get_str();
            }
            else if( name == "pdirection_y" )
            {
                step[k].CLv_step_prop[1]=value.get_str();
            }
            else if( name == "pdirection_z" )
            {
                step[k].CLv_step_prop[2]=value.get_str();
            }
            else if( name == "point_1_x" )
            {
                step[k].CLv_step_prop[3]=value.get_str();
            }
            else if( name == "point_1_y" )
            {
                step[k].CLv_step_prop[4]=value.get_str();
            }            
            else if( name == "point_1_z" )
            {
                step[k].CLv_step_prop[5]=value.get_str();
            }
            else if( name == "gravity" )
            {
                step[k].CLv_step_prop[6]=value.get_str();
            }
            else if( name == "wrotation" )
            {
                step[k].CLv_step_prop[7]=value.get_str();
            }
        }
    }
}

///lecture des CL volumiques dans un fichier json
void DataUser::read_json_behaviour_bc_volume(const Object& gr){
    for(Object::size_type i=0;i != gr.size() ;i++){
        const Pair& pair1 = gr[i];        
        const Value&  value1 = pair1.value_;
        Object obj=value1.get_obj();
        for( Object::size_type j = 0; j != obj.size(); ++j )
        {
            const Pair& pair = obj[j];
            const std::string& name  = pair.name_;
            const Value&  value = pair.value_;

            if( name == "name" )
            {
                behaviour_bc_volume[i].name= value.get_str();
            }
            else if( name == "type" )
            {
                behaviour_bc_volume[i].type = value.get_str() ;
            }
            else if( name == "select" )
            {
                std::string temp=value.get_str() ;
                std::istringstream s(temp);
                s >> behaviour_bc_volume[i].select;
            }
            else if( name == "ref" )
            {
                behaviour_bc_volume[i].ref = value.get_int();
            }
            else if( name == "step" )
            {
                Object obj2=value.get_obj();
                behaviour_bc_volume[i].step.resize(obj2.size());
                read_step_bc_volume(obj2,behaviour_bc_volume[i].step);
            }
            else
            {
                assert( "Donnee behaviour_bc_volume non implementee" );
            }
        }
    }

}


void DataUser::read_step_bc(const Object& gr, BasicVec<StepBc> &step){
    for( Object::size_type k = 0; k != gr.size(); ++k )
    {
        const Pair& pair1 = gr[k];
        //const std::string& name1  = pair1.name_;
        const Value&  value1 = pair1.value_;
        Object obj=value1.get_obj();
        for( Object::size_type l = 0; l != obj.size(); ++l )
        {
            const Pair& pair = obj[l];
            const std::string& name  = pair.name_;
            const Value&  value = pair.value_;
            if( name == "fct_spatiale_x" )
            {
                step[k].CL_step_prop[0]=value.get_str();
            }
            else if( name == "fct_spatiale_y" )
            {
                step[k].CL_step_prop[1]=value.get_str();
            }
            else if( name == "fct_spatiale_z" )
            {
                step[k].CL_step_prop[2]=value.get_str();
            }
            else if( name == "fct_temporelle_x" )
            {
                step[k].CL_step_prop[3]=value.get_str();
            }
            else if( name == "fct_temporelle_y" )
            {
                step[k].CL_step_prop[4]=value.get_str();
            }
            else if( name == "fct_temporelle_z" )
            {
                step[k].CL_step_prop[5]=value.get_str();
            }
        }
    }
}
///lecture des CL donnees dans un fichier json
void DataUser::read_json_behaviour_bc(const Object& gr){
    for(Object::size_type i=0;i != gr.size() ;i++){
        //Object obj=pr[i].get_obj();;
    const Pair& pair1 = gr[i];        
        const Value&  value1 = pair1.value_;
    Object obj=value1.get_obj();

        for( Object::size_type j = 0; j != obj.size(); ++j )
        {
            const Pair& pair = obj[j];
            const std::string& name  = pair.name_;
            const Value&  value = pair.value_;

            if( name == "id" )
            {
                behaviour_bc[i].id= value.get_int();
            }
            else if( name == "type" )
            {
                behaviour_bc[i].type= value.get_str();
            }
            else if( name == "step" )
            {
                Object obj2=value.get_obj();
                behaviour_bc[i].step.resize(obj2.size());
                read_step_bc(obj2,behaviour_bc[i].step);
            }
            else
            {
                std::cout << "Champ CL non implementee : " << name << std::endl;
            }        
        }
    }   
    //Ajout d'une CL en effort nul pour les interfaces non assignees sur les cotes
    behaviour_bc[gr.size()].type="effort";
    behaviour_bc[gr.size()].id=-1;
    
    StepBc step_0;
    step_0.CL_step_prop[0] = "0";
    step_0.CL_step_prop[1] = "0";
    step_0.CL_step_prop[2] = "0";
    step_0.CL_step_prop[3] = "0";
    step_0.CL_step_prop[4] = "0";
    step_0.CL_step_prop[5] = "0";

    behaviour_bc[gr.size()].step.resize(1);
    behaviour_bc[gr.size()].step[0] = step_0;
}


void DataUser::read_json_calcul(std::string file_calcul){
    std::ifstream is( file_calcul.c_str() );
    Value value_i;

    read( is, value_i );
    const Object& input = value_i.get_obj();
    
    for( Object::size_type i = 0; i != input.size(); ++i )
    {
        const Pair& pair_groups = input[i];
        const std::string& name_groups  = pair_groups.name_;
        const Value& value_groups= pair_groups.value_;
        
        if(name_groups=="mesh"){//lecture des donnees de maillage
            const Object obj=value_groups.get_obj();
            for( Object::size_type j = 0; j != obj.size(); ++j ){
                const Pair& pair = obj[j];
                const std::string& name  = pair.name_;
                const Value&  value = pair.value_;
                if( name == "model_directory" ){name_directory=value.get_str();}
                else if( name == "mesh_directory" ){mesh_directory=value.get_str();}
                else if( name == "mesh_name" ){name_mesh_user= value.get_str();}
                else if( name == "extension" ){extension= value.get_str();}
                else if( name == "nb_sst" ){std::cout << "nb_sst : " <<  value.get_int() << std::endl;}
                else if( name == "nb_inter" ){std::cout << "nb_inter : " <<  value.get_int() << std::endl;}
                else if( name == "nb_groups_elem" ){std::cout << "nb_groups_elem : " <<  value.get_int() << std::endl;}
                else if( name == "nb_groups_inter" ){std::cout << "nb_groups_inter : " <<  value.get_int() << std::endl;}
                else{assert( false );}
            }        
        }
        
        if(name_groups=="groups_elem"){//lecture des groupes d'elements
            std::cout << "in groups_elem " <<  std::endl;
            const Object& gr = value_groups.get_obj();
            group_elements.resize(gr.size());
            read_json_groups_elements(gr);
        }
        
        if(name_groups=="groups_inter"){//lecture des groupes d'interfaces
            std::cout << "in groups_inter " <<  std::endl;
            const Object& gr = value_groups.get_obj();
            group_interfaces.resize(gr.size());
            read_json_groups_interfaces(gr);
        }
        
        if(name_groups=="groups_edge"){//lecture des groupes d'interfaces
            std::cout << "in groups_edge " <<  std::endl;
            const Object& gr = value_groups.get_obj();
            group_edges.resize(gr.size()+1);
            read_json_groups_edges( gr);
        }
        
        if(name_groups=="links"){//lecture des proprietes d'interfaces et creation de behaviour_links
            std::cout << "in links " <<  std::endl;
            const Object& gr = value_groups.get_obj();
            behaviour_links.resize(gr.size()+1);
            read_json_behaviour_interfaces(gr);
        }
        
        if(name_groups=="materials"){//lecture des caracteristiques materiaux et creation de behaviour_materials
            std::cout << "in materials " <<  std::endl;
            const Object& mat = value_groups.get_obj();
            behaviour_materials.resize(mat.size());
            read_json_behaviour_materials(mat);
        }
        
        if(name_groups=="CL"){//lecture des proprietes d'interfaces et creation de behaviour_links
            std::cout << "in CL " <<  std::endl;
            const Object& cl = value_groups.get_obj();
            behaviour_bc.resize(cl.size()+1);
            read_json_behaviour_bc(cl);
        }
        if(name_groups=="CLvolume"){//lecture des proprietes d'interfaces et creation de behaviour_links
            std::cout << "in CLvolume " <<  std::endl;
            const Object& clvol = value_groups.get_obj();
            behaviour_bc_volume.resize(clvol.size());
            read_json_behaviour_bc_volume(clvol);
        }
        if(name_groups=="options"){//lecture des donnees de maillage
            const Object obj=value_groups.get_obj();
            for( Object::size_type j = 0; j != obj.size(); ++j ){
                const Pair& pair = obj[j];
                const std::string& name  = pair.name_;
                const Value&  value = pair.value_;
                if( name == "dissipation" ){std::string temp=value.get_str();std::istringstream s(temp); s>> options.dissipation;}
                else if( name == "nb_options" ){options.nb_option=value.get_int();}
                else if( name == "PREC_nb_niveaux" ){std::string temp=value.get_str();std::istringstream s(temp); s>>options.nb_level;}
                else if( name == "Temp_statique" ){options.Temp_statique=value.get_str();}
                else if( name == "LATIN_nb_iter" ){std::string temp=value.get_str();std::istringstream s(temp); s>>options.LATIN_nb_iter_max;}
                else if( name == "LATIN_conv" ){std::string temp=value.get_str();std::istringstream s(temp); s>>options.LATIN_crit_error;}
                else if( name == "2D_resolution" ){std::string temp=value.get_str();std::istringstream s(temp); s>>options.resolution_2D;}
                else if( name == "Multiresolution_on" ){options.Multiresolution_on=value.get_int();}
                else if( name == "Multiresolution_nb_cyle" ){options.Multiresolution_nb_cyle=value.get_int();}    
            }        
        }
    }
}





