#ifndef PARAM_TEMPS_H
#define PARAM_TEMPS_H

#include "../UTILS/Sc2String.h"
using namespace LMT;

/** \ingroup Paramètres_calcul
Parametres temporels
*/

struct STEP{
    STEP(){
        dt=1;           ///< pas de temps dans le step
        t_ini=0;        ///< piquet de temps initial du step
        t_fin=1;        ///< piquet de temps final du step
        nb_time_step=1; ///< nombre de pas de temps dans le step
        pt_cur=0;       ///< pas de temps courant dans le step
    }
    double dt, t_ini, t_fin;
    int nb_time_step;
    int pt_cur;
};

struct TEMPS{
    TEMPS(){
        pt=0;
        pt_cur=0;
        dt=1;
        theta = 1;
        nbpastemps=1;
        type_de_calcul="stat";
        nb_step=1;
        current_time=0;
    }
    void read_data_user(DataUser &data_user){
        if(data_user.options.Temp_statique == "statique"){
            type_de_calcul= "stat";
        }else if(data_user.options.Temp_statique == "quasistatique"){
            type_de_calcul= "Qstat";
        }
        
        if (type_de_calcul=="stat") {
            nbpastemps=1;
            dt=0;
            nb_step = 1;
            time_step.resize(nb_step);
            nbpastemps=nb_step;
            for(unsigned i_step=0;i_step<nb_step;i_step++){
                time_step[i_step].dt=1;
                time_step[i_step].t_ini=0;
                time_step[i_step].t_fin=0;
                time_step[i_step].nb_time_step=1;
            }
            
        }else if(type_de_calcul=="Qstat") {
            nb_step = data_user.time_step.size();
            time_step.resize(nb_step);
            nbpastemps= 0;
            for(unsigned i_step=0;i_step<nb_step;i_step++){
                time_step[i_step].dt = data_user.time_step[i_step].dt;
                time_step[i_step].t_ini = data_user.time_step[i_step].ti;
                time_step[i_step].t_fin = data_user.time_step[i_step].tf;
                time_step[i_step].nb_time_step = data_user.time_step[i_step].nb_time_step;
                nbpastemps += time_step[i_step].nb_time_step;
            }
        }
    }
    int pt;                 ///< piquet de temps de calcul (0 ou 1)
    int pt_cur;             ///< piquet de temps (ou intervalle de temps) courant global
    int step_cur;           ///< step courant
    double current_time;    ///< valeur courante du piquet de temps de calcul
    double dt;              ///< pas de temps
    double theta;           ///< paramètre de la theta_methode
    unsigned nbpastemps;    ///< nb de pas de temps total
    Sc2String type_de_calcul;     ///< type de calcul, choix entre : stat (statique (1 pas de temps)), qstat (plusieurs pas de temps)
    Vec<STEP> time_step;            ///< steps de calcul
    int nb_step;                    ///< nombre de step de calcul
};

#endif // PARAM_TEMPS_H

