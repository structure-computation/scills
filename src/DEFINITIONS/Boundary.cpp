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
fcts_spatiales(Vec<UserParameter>(UserParameter("F1"),UserParameter("F2")))
#elif DIM == 3
fcts_spatiales(Vec<UserParameter>(UserParameter("F1"),UserParameter("F2"),UserParameter("F3")))
#endif
{
    CL_parameters.addParameter(&(fcts_spatiales[0]));
    CL_parameters.addParameter(&(fcts_spatiales[1]));
#if DIM == 3
    CL_parameters.addParameter(&(fcts_spatiales[2]));
#endif
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
}

void Boundary::prepareParameters(){
    CL_parameters.prepareParameters();
}

void Boundary::updateParameters(){
    /// Remettre ici les fonctions dans les fichiers prelocalstage
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