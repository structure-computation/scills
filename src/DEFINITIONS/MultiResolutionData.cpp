#include "MultiResolutionData.h"

MultiResolutionData::MultiResolutionData():
parameters(),
m("m")
{
    parameters.addParameter(&m);
    
}

void MultiResolutionData::free()
{
}

void MultiResolutionData::read_data_user(const Metil::DataUser &data_user){
    const DataUser::Json_multiresolution_parameters &multiresolution_parameters = data_user.multiresolution_parameters;
    /// Recuperation du type
    type = multiresolution_parameters.multiresolution_type;
    if(type == "off"){
        nb_calculs = 1;
    }else if(type == "sequential"){
        PRINT(multiresolution_parameters.multiresolution_type);
        PRINT(multiresolution_parameters.resolution_number);
        nb_calculs = multiresolution_parameters.resolution_number;
        const int nb_parameters = multiresolution_parameters.collection_vec.size();
        parameters_data.resize(nb_parameters);
        for(unsigned i = 0; i < nb_parameters; i++){
            parameters_data[i].name = multiresolution_parameters.collection_vec[i].name;
            Sc2String type = multiresolution_parameters.collection_vec[i].type;
            if(type == "function") {
                parameters_data[i].type = ParameterData::Function;
                parameters_data[i].values << multiresolution_parameters.collection_vec[i].function;
            }
            else if(type == "list") {
                parameters_data[i].type = ParameterData::List;
                parameters_data[i].values = multiresolution_parameters.collection_vec[i].values;
            }
            parameters_data[i].user_parameter = new UserParameter(multiresolution_parameters.collection_vec[i].name,&parameters);
        }
    }else{
        std::cerr << "Mauvais type de multi-resolution : " << type << std::endl;
        assert(0);
    }
}

void MultiResolutionData::prepareParameters(){
    for(int i = 0; i < parameters_data.size(); i++) {
        switch(parameters_data[i].type) {
            case ParameterData::Function:
                parameters_data[i].user_parameter->setExpression(parameters_data[i].values[0]);
                break;
            case ParameterData::List:
                parameters_data[i].user_parameter->setExpression(parameters_data[i].values[calcul_cur]);
                break;
        }
    }
    parameters.prepareParameters();
}

void MultiResolutionData::updateParameters(){
    if(type == "off"){
        /// il n'y a rien a faire
    }else if(type == "sequential"){
        //affiche();
        prepareParameters();            /// Necessaire pour la gestion des parametres de type "List"
        parameters.updateParameters();  /// Une seule fonction, pour tous les calculs
    } else{
        std::cerr << "Mauvais type de calcul parametrique : '" << type << "'" << std::endl;
        assert(0);
    }
    //* DEBUG : affichage des parametres multi-resolution apres mise a jour
    std::cout << "Mise a jour des parametres de multi-resolution *******************************" << std::endl;
    std::cout << m.symbol << " = " << m.value << std::endl;
    for(int i = 0; i < parameters.user_parameters.size(); i++){
        std::cout << parameters.user_parameters[i]->symbol << " = " << parameters.user_parameters[i]->value << std::endl;
    }
    std::cout << "******************************************************************************" << std::endl;
    //*/
}

void MultiResolutionData::init(){
    calcul_cur = 0;
    m.value = calcul_cur;
}

void MultiResolutionData::next(){
    calcul_cur++;
    m.value = calcul_cur;
}

bool MultiResolutionData::has_next(){
    return (calcul_cur < nb_calculs);
}

void MultiResolutionData::affiche(){
    std::cout << "*************************************************************" << std::endl;
    std::cout << "**************** Debug MultiResolutionData : ****************" << std::endl;
    std::cout << "*************************************************************" << std::endl;
    std::cout << "Type : " << type << std::endl;
    std::cout << "nb_calculs : " << nb_calculs << std::endl;
    std::cout << "calcul_cur : " << calcul_cur << std::endl;
    std::cout << "m : " << std::endl; m.affiche();
    std::cout << "parameters : " << std::endl; parameters.affiche();
    std::cout << "*************************************************************" << std::endl;
    std::cout << "*************************************************************" << std::endl;
}
