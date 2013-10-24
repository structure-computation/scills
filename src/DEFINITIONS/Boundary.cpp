#include "Boundary.h"
//#include "../UTILITAIRES/utilitaires.h"


ParameterGroup Boundary::CL_parameters;
MainParameter Boundary::x("x",&(Boundary::CL_parameters));
MainParameter Boundary::y("y",&(Boundary::CL_parameters));
#if DIM == 3
MainParameter Boundary::z("z",&(Boundary::CL_parameters));
#endif

Boundary::Boundary():
#if DIM == 2
fcts_spatiales(Vec<UserParameter>(UserParameter("F1"),UserParameter("F2"))),
#elif DIM == 3
fcts_spatiales(Vec<UserParameter>(UserParameter("F1"),UserParameter("F2"),UserParameter("F3"))),
#endif
vitesse(Vec<UserParameter>(UserParameter("R1"),UserParameter("R2"),UserParameter("R3"))),
rotation(Vec<UserParameter>(UserParameter("M1"),UserParameter("M2"),UserParameter("M3")))
{
    CL_parameters.addParameter(&(fcts_spatiales[0]));
    CL_parameters.addParameter(&(fcts_spatiales[1]));
    CL_parameters.addParameter(&(vitesse[0]));
    CL_parameters.addParameter(&(vitesse[1]));
    CL_parameters.addParameter(&(vitesse[2]));
    CL_parameters.addParameter(&(rotation[0]));
    CL_parameters.addParameter(&(rotation[1]));  
    CL_parameters.addParameter(&(rotation[2]));
#if DIM == 3
    CL_parameters.addParameter(&(fcts_spatiales[2])); 
#endif
    ///initialisation de la base
    for( unsigned d=0;d<3;d++ ) {
        for( unsigned dd=0;dd<3;dd++ ) {
            if(d==dd) Base[d][dd] = 1.;
            else Base[d][dd] = 0.; 
        }
    }
    
    
}

Boundary::~Boundary(){}

void Boundary::read_data_user(int index,DataUser &data_user){
    const DataUser::Json_boundary_conditions &json_boundary_condition = data_user.boundary_conditions_vec[index];
    id = json_boundary_condition.id_in_calcul;
    comp = json_boundary_condition.condition_type;
    fcts_spatiales[0] = json_boundary_condition.spatial_function_x;
    if(comp.find("normal")<comp.size()){/// Pas besoin des autres si c'est selon la normale
        return;
    }
    fcts_spatiales[1] = json_boundary_condition.spatial_function_y;
#if DIM == 3
    fcts_spatiales[2] = json_boundary_condition.spatial_function_z;
#endif
    
    if(comp == "cinetic_torseur"){/// Pas besoin des autres si c'est selon la normale
        PRINT("read data user : cinetic_torseur");
        Centre[0] = json_boundary_condition.point_1_x;
        Centre[1] = json_boundary_condition.point_1_y;
        Centre[2] = json_boundary_condition.point_1_z; 
        
        Base[0][0] = json_boundary_condition.dir_1_x;
        Base[0][1] = json_boundary_condition.dir_1_y;
        Base[0][2] = json_boundary_condition.dir_1_z;
        
        Base[1][0] = json_boundary_condition.dir_2_x;
        Base[1][1] = json_boundary_condition.dir_2_y;
        Base[1][2] = json_boundary_condition.dir_2_z;
        
        Base[2][0] = json_boundary_condition.dir_3_x;
        Base[2][1] = json_boundary_condition.dir_3_y;
        Base[2][2] = json_boundary_condition.dir_3_z;
        
        vitesse[0] = json_boundary_condition.R_0;
        vitesse[1] = json_boundary_condition.R_1;
        vitesse[2] = json_boundary_condition.R_2;
        rotation[0] = json_boundary_condition.M_0;
        rotation[1] = json_boundary_condition.M_1;
        rotation[2] = json_boundary_condition.M_2;
        
        imp_vitesse[0] = json_boundary_condition.imp_R_0;
        imp_vitesse[1] = json_boundary_condition.imp_R_1;
        imp_vitesse[2] = json_boundary_condition.imp_R_2;
        imp_rotation[0] = json_boundary_condition.imp_M_0;
        imp_rotation[1] = json_boundary_condition.imp_M_1;
        imp_rotation[2] = json_boundary_condition.imp_M_2;    
    }
    
    
}

void Boundary::prepareParameters(){
    CL_parameters.prepareParameters();
}

void Boundary::updateParameters(){
    /// Remettre ici les fonctions dans les fichiers prelocalstage
    std::cout << "update boudaries ok" << std::endl;
}

void Boundary::affiche(){
    std::cout << std::endl;
    std::cout << "*************************************************************" << std::endl;
    std::cout << "********************** Debug Boundary : *********************" << std::endl;
    std::cout << "*************************************************************" << std::endl;
    std::cout << "id : " << id << std::endl;
    std::cout << "comportement : " << comp << std::endl;
    std::cout << "valeur courante :" << std::endl;
    if(comp.find("normal")<comp.size()){
        std::cout << "    fn :" << fcts_spatiales[0].str_expr << " = " << fcts_spatiales[0];
    }else{
        std::cout << "    fx :" << fcts_spatiales[0].str_expr << " = " << fcts_spatiales[1];
        std::cout << " ;     fy :" << fcts_spatiales[1].str_expr << " = " << fcts_spatiales[2];
        if(DIM == 3)
            std::cout << " ;     fz :" << fcts_spatiales[2].str_expr << " = " << fcts_spatiales[2];
    }
    std::cout << std::endl;
    std::cout << "*************************************************************" << std::endl;
    std::cout << "*************************************************************" << std::endl;
    std::cout << std::endl;
}