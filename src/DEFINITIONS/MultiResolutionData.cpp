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
    if(type == "Off"){
        nb_calculs = 1;
    }else if(type == "plan d'experience"){
        std::cerr << "Type de multi-resolution : " << type << " a implementer" << std::endl;
        assert(0);
    }else if(type == "function"){
        nb_calculs = multiresolution_parameters.resolution_number;
        for(unsigned i = 0; i < multiresolution_parameters.collection_vec.size(); i++){
            /// Creation de la fonction parametrique (pour le moment...)
            Sc2String function,minvalue,maxvalue,nbvalues;
            minvalue << multiresolution_parameters.collection_vec[i].min_value;
            maxvalue << multiresolution_parameters.collection_vec[i].max_value;
            nbvalues << multiresolution_parameters.collection_vec[i].nb_value;
            if(multiresolution_parameters.collection_vec[i].nb_value>1){
                function = minvalue + "+(" + maxvalue + "-" + minvalue + ")*m/(" + nbvalues + "-1)";
            }else{
                function << multiresolution_parameters.collection_vec[i].nominal_value;
            }
            /// Creation du parametre
            UserParameter *PM = new UserParameter(multiresolution_parameters.collection_vec[i].name,&parameters);
            PM->setExpression(function);
        }
    }else{
        std::cerr << "Mauvais type de multi-resolution : " << type << std::endl;
        assert(0);
    }
}

void MultiResolutionData::prepareParameters(){
    parameters.prepareParameters();
}

void MultiResolutionData::updateParameters(){
    if(type == "Off"){
        /// il n'y a rien a faire
    } else if(type == "plan d'experience"){
        //TODO
        //int tmp = nb_calculs;
        //Codegen::Ex::MapExNum values = user_parameters.getValuesMap();
    } else if(type == "liste"){
        parameters.updateParameters();  /// Tous les parametres sont mis a l'expression indicee calcul_cur
    }else if(type == "function"){
        parameters.updateParameters();  /// Une seule fonction, pour tous les calculs
    } else{
        std::cerr << "Mauvais type de calcul parametrique : '" << type << "'" << std::endl;
        assert(0);
    }
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
