#include "DataUser.h"

using namespace Metil;
using namespace LMT;
using namespace std;

inline void read_data_structure(Param &process, DataUser &data_user) {
    process.structure->repertoire_des_maillages = data_user.model_path + "/MESH/";
    process.structure->nom_fichier_qualification_materiaux = "";
    process.structure->nom_des_maillages = "";
    process.structure->nb_maillages = 0;
    process.structure->extension = "";
    process.structure->jeu_physique = 0;   
};
