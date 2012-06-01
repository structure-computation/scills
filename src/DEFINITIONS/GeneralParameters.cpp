#include "GeneralParameters.h"
#include "../UTILS/utils_2.h"

GeneralParameters::GeneralParameters(){
    extension=".avs";
}

void GeneralParameters::read_data_user(DataUser &data_user) {
    repertoire_des_maillages = data_user.model_path + "/MESH/";
    nom_fichier_qualification_materiaux = "";
    nom_des_maillages = "";
    nb_maillages = 0;
    extension = "";
    jeu_physique = 0;   
}

void GeneralParameters::display_all_data(){
    std::cout << std::endl;
    std::cout << "*************************************************************" << std::endl;
    std::cout << "***************** Debug GeneralParameters : *****************" << std::endl;
    std::cout << "*************************************************************" << std::endl;
    debug("Sc2String nom_fichier_qualification_materiaux",nom_fichier_qualification_materiaux);
    debug("Sc2String repertoire_des_maillages           ",repertoire_des_maillages);
    debug("unsigned nb_maillages                        ",nb_maillages);
    debug("Sc2String nom_des_maillages                  ",nom_des_maillages);
    debug("Sc2String extension                          ",extension);
    debug("unsigned jeu_physique                        ",jeu_physique);
    debug("Sc2String nom_maillages_jeu                  ",nom_maillages_jeu);
    debug("Vec< Vec<int> > inter_jeu                    ",inter_jeu);
    debug("TYPEREEL volumetot",volumetot);
    std::cout << "*************************************************************" << std::endl;
    std::cout << "*************************************************************" << std::endl;
}