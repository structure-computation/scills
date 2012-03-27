#ifndef MULTISCALEPARAMETERS_H
#define MULTISCALEPARAMETERS_H

#include "../COMPUTE/DataUser.h"

#include "../../LMT/include/containers/vec.h"
using namespace LMT;


/** Parametres multiechelles
*/
struct MultiScaleParameters{
    //attributs==============================================================================================
    int multiechelle;           /// multiechelle = 1, monoechelle 0
    unsigned type_base_macro;   /// type de fct macro : 1 : Resultantes, 2: Moments, 3: partie lineaire
    unsigned sizeM;             /// nombre d'inconnues macro total
    Vec<unsigned> nb_fcts_macro_par_inter;  ///< nombre de fonctions macro par interface dans le cas particulier où le nombre est différent pour certaines interfaces
    Vec<unsigned> inter_correspondantes;    ///< interfaces dont le nombre de fcts macro est different du nombre par defaut
    bool nbmacro_identique;     ///< booleen indiquant si le nombre de fonctions de base macro est identique par interface
    bool opti_multi;            ///< optimisation du calcul en stoppant le macro a partir d'une erreur  donnée
    double erreur_macro;        ///< niveau d'erreur sur la solution macro à partir de laquelle on passe en monoéchelle
    
    //methodes===============================================================================================
    MultiScaleParameters();
    void read_data_user(DataUser &data_user);
    void display_all_data();
};

#endif //MULTISCALEPARAMETERS_H

