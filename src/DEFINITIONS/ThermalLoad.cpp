#include "ThermalLoad.h"

ThermalLoad::ThermalLoad():
parameters(),
x("x"),
deltaT(UserParameter("deltaT")),
y("y"),
z("z")
{
    parameters.addParameter(&x);
    parameters.addParameter(&y);
    parameters.addParameter(&z);
    parameters.addParameter(&deltaT);
}

void ThermalLoad::read_data_user(const DataUser &data_user){
    const DataUser::Json_thermal_load &json_thermal_load = data_user.thermal_load;
    deltaT = json_thermal_load.function;
}

void ThermalLoad::prepareParameters(){
    parameters.prepareParameters();
}

void ThermalLoad::updateParameters(){
    /// Remettre ici le code de maj des CL
}

void ThermalLoad::affiche() const{
    std::cout << "*************************************************************" << std::endl;
    std::cout << "******************* Debug ThermalLoad : *********************" << std::endl;
    std::cout << "*************************************************************" << std::endl;

    std::cout << "    delatT : " << deltaT.str_expr << " = " << deltaT           << std::endl;
    
    std::cout << "*************************************************************" << std::endl;
    std::cout << "*************************************************************" << std::endl;
}