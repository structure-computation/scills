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
    
    TYPEREEL dt;        /// pas de temps dans le step
    TYPEREEL t_ini;     /// piquet de temps initial du step
    TYPEREEL t_fin;     /// piquet de temps final du step
    int nb_time_step;   /// nombre de pas de temps dans le step
    int pt_cur;         /// pas de temps courant dans le step
};

struct TimeParameters{
    
    //attributs==============================================================================================
    
    Sc2String type_de_calcul;   /// type de calcul, choix entre : stat (statique (1 pas de temps)), qstat (plusieurs pas de temps)
    unsigned nbpastemps;        /// nb de pas de temps total
    Vec<STEP> time_step;        /// steps de calcul
    int nb_step;                /// nombre de step de calcul
    int step_cur;               /// step courant
    int pt;                     /// piquet de temps de calcul (0 ou 1)
    int pt_cur;                 /// piquet de temps (ou intervalle de temps) courant global
    TYPEREEL current_time;      /// valeur courante du piquet de temps de calcul
    TYPEREEL dt;                /// duree du pas de temps courant
    TYPEREEL theta;             /// paramètre de la theta_methode (A REVOIR)
    bool __step_changed;
    
    
    //methodes===============================================================================================
    
    TimeParameters();                           /// Constructeur
    void read_data_user(DataUser &data_user);   /// Lecture du DataUser
    void display_all_data();                    /// Affichage pour le debug
    
    void init();            /// (Re)Initialiser les variables au 1er pas de temps
    void next();            /// Passer au pas de temps (step si necessaire) suivant
    bool has_next();        /// Indique s'il reste d'autres pas de temps a calculer
    bool step_changed();    /// Indique si on a change de step (ne peut etre appele qu'une seule fois!!! raz du booleen de stockage dans la fonction)
};

#endif // TIMEPARAMETERS_H

