#ifndef TIMEDATA_H
#define TIMEDATA_H

#include "main_typedef.h"
#include "../PARAMETERS/Parameters.h"

/** \ingroup Paramètres_calcul
Parametres temporels
*/

struct STEP{
    STEP(){
        id = -1;
        dt=1;
        t_ini=0;
        t_fin=1;
        nb_time_step=1;
        pt_cur=0;
    }
    int id;
    Scalar dt;        /// pas de temps dans le step
    Scalar t_ini;     /// piquet de temps initial du step
    Scalar t_fin;     /// piquet de temps final du step
    int nb_time_step;   /// nombre de pas de temps dans le step
    int pt_cur;         /// pas de temps courant dans le step
};

struct TimeData{
    
    //attributs==============================================================================================
    
    Sc2String type_de_calcul;   /// type de calcul, choix entre : stat (statique (1 pas de temps)), qstat (plusieurs pas de temps)
    unsigned nbpastemps;        /// nb de pas de temps total
    Vec<STEP> time_step;        /// steps de calcul
    int nb_step;                /// nombre de step de calcul
    int step_old;               /// step courant
    int step_cur;               /// step courant
    int pt;                     /// piquet de temps de calcul (0 ou 1)
    int pt_cur;                 /// piquet de temps (ou intervalle de temps) courant global
    Scalar dt;                  /// duree du pas de temps courant
    Scalar theta;               /// paramètre de la theta_methode (A REVOIR)
    MainParameter t_cur;        /// Parametre representant le temps
    ParameterGroup parameters;  /// structure contenant les parametres temporels
    LMT::Vec<LMT::Vec<Sc2String> > expressions;   /// Expression de tous les parametres a chaque step de temps
    
    
    //methodes===============================================================================================
    
    TimeData();                                     /// Constructeur
    ~TimeData();                                    /// Destructeur
    void read_data_user(const DataUser &data_user); /// Lecture du DataUser
    void prepareParameters();                       /// Preparation des parametres pour leur future evaluation (assignation des expressions)
    void updateParameters();                        /// Mis a jour des expressions (si changement de step) et valeurs des parametres
    void affiche();                                 /// Affichage pour le debug
    
    void init();                /// (Re)Initialiser les variables au 1er pas de temps
    void next();                /// Passer au pas de temps (step si necessaire) suivant
    bool has_next();            /// Indique s'il reste d'autres pas de temps a calculer
    bool step_changed();        /// Indique si on a change de step (ne peut etre appele qu'une seule fois!!! raz du booleen de stockage dans la fonction)
};

#endif // TIMEDATA_H

