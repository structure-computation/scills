#ifndef TIMEPARAMETERS_H
#define TIMEPARAMETERS_H

#include "../COMPUTE/DataUser.h"
#include "../UTILS/Sc2String.h"
#include "../../LMT/include/containers/vec.h"
using namespace LMT;

/** \ingroup Paramètres_calcul
Parametres temporels
*/

struct STEP{
    STEP(){
        dt=1;
        t_ini=0;
        t_fin=1;
        nb_time_step=1;
        pt_cur=0;
    }
    
    double dt;          ///< pas de temps dans le step
    double t_ini;       ///< piquet de temps initial du step
    double t_fin;       ///< piquet de temps final du step
    int nb_time_step;   ///< nombre de pas de temps dans le step
    int pt_cur;         ///< pas de temps courant dans le step
};

struct TimeParameters{
    
    //attributs==============================================================================================
    
    int pt;                 ///< piquet de temps de calcul (0 ou 1)
    int pt_cur;             ///< piquet de temps (ou intervalle de temps) courant global
    int step_cur;           ///< step courant
    double current_time;    ///< valeur courante du piquet de temps de calcul
    double dt;              ///< pas de temps
    double theta;           ///< paramètre de la theta_methode
    unsigned nbpastemps;    ///< nb de pas de temps total
    Sc2String type_de_calcul;   ///< type de calcul, choix entre : stat (statique (1 pas de temps)), qstat (plusieurs pas de temps)
    Vec<STEP> time_step;        ///< steps de calcul
    int nb_step;                ///< nombre de step de calcul
    
    
    //methodes===============================================================================================
    
    TimeParameters();
    void read_data_user(DataUser &data_user);
    void display_all_data();
};

#endif // TIMEPARAMETERS_H

